# ==============================================================================
# Sex-matched subgroup selection - BARIA cohort
# ==============================================================================
# Goal: select a subgroup of male and female BARIA participants that is
# balanced on age and medication use, so that downstream comparisons between
# sexes are not confounded by these variables.
#
# Matching variables: age (v0_age), lipid-lowering drug use (v0_lld),
# metformin use (v0_metformine)
#
# Approach: optimal matching (MatchIt::matchit) using sex as the "treatment"
# variable, matching each male to multiple females (see ratio below).
# ==============================================================================

# --- Libraries ----------------------------------------------------------------
library(tidyverse)
library(MatchIt)
library(ggplot2)
library(tableone)

# --- Load data ------------------------------------------------------------------
serum       <- readRDS("data/serum_ids.RDS")            # IDs of participants with serum samples available
rna         <- readRDS("data/idsrnaseq.RDS")             # IDs of participants with RNAseq data available
baria_all   <- readRDS("data/baria_metadata.RDS")         # full BARIA metadata (all participants)
bariadiet   <- readRDS("data/bariatot.RDS")               # participants with processed dietary data

# Some participants have dietary questionnaires that were not yet processed
# into `bariadiet` at the time of matching, so their raw dietary lists are
# read in separately and combined by ID below.
otherdiet <- readxl::read_xlsx("data/extradieetlijsten.xlsx") # at the time of matching, I did not have the processed dietary data of these participants yet.
otherdiet$ID <- str_c("BARIA_", otherdiet$studienummer)

# --- Build eligible cohort ------------------------------------------------------
# Step 1: keep participants that have (raw or processed) dietary data AND a
# serum sample available.
baria <- baria_all %>% filter(ID %in% c(bariadiet$ID, otherdiet$ID)) %>%
        filter(ID %in% serum)
dim(baria)

# Step 2: require a follow-up (visit 4) BMI measurement.
baria <- baria %>% dplyr::filter(!is.na(v4_bmi))
dim(baria)
summary(baria$sex)

# Step 3: exclude participants with recent antibiotic use or diarrhoea at
# baseline (both can disturb gut/immune measures), and re-apply the serum
# filter for safety.
baria <- baria %>% filter(v0_antibiotics != "yes" & v0_diarrhoea != "yes") %>%
        filter(ID %in% serum)
dim(baria)
summary(baria$sex)          # sex distribution of the final eligible cohort
summary(baria$ID %in% rna)  # how many eligible participants also have RNAseq data

# --- Matching: males vs. females --------------------------------------------
# All eligible males are used as the pool for matching (previously a random
# subsample of 50 was used instead - see commented-out `sample()` call).
malepart <- baria %>% filter(sex == "male") %>% pull(ID)
malesubgroup <- malepart # sample(malepart, 50)

# Restrict to the matching variables and drop anyone with missing values on
# any of them (MatchIt cannot match on incomplete cases).
matchData <- baria %>% select(ID, v0_age, sex, v0_dmmed, v0_lld, v0_metformine) %>%
        filter(ID %in% malesubgroup | sex == "female") %>% # all females vs 50 randomly selected males
        drop_na()

# Optimal matching of females to males on age, LLD use and metformin use.
# ratio = 2.12 -> up to ~2 females matched per male (approximates matching
# ALL eligible females against the male pool); max.controls caps this at 10
# females per male.
matchRes <- matchit(sex ~ v0_age + v0_lld + v0_metformine, data = matchData, method="optimal",
                    ratio = 2.12, max.controls = 10) # 2.12 because we get the right number out
matchedCases <- match.data(matchRes)

# --- Add descriptive variables for Table 1 --------------------------------------
# These variables were not part of the matching itself, but are added
# afterwards to check balance / describe the matched sample.
extradata <- baria %>% select(ID, v0_dm, v0_diabp, v0_sysbp, v0_bmi, v0_leuko, v4_leuko, v0_smoking,
                        v4_bmi, v5_bmi, v0_dm)
matchedCases <- matchedCases %>% left_join(extradata, by="ID")

# --- Table 1: check balance between matched males and females -------------------
tableMatched <- CreateTableOne(vars = c('v0_age','sex', 'v0_bmi', 'v0_sysbp', 'v0_diabp', 'v0_smoking',
                                                'v0_dm','v0_lld','v0_dmmed', 'v0_metformine',
                                                'v0_leuko', 'v4_leuko', 'v4_bmi', 'v5_bmi'
                                                ),
                            data = matchedCases,
                            factorVars = c('sex','v0_dm','v0_smoking', 'v0_lld', 'v0_dmmed',
                                           'v0_metformine'),
                            strata = 'sex')
tab <- print(tableMatched,missing = T, noSpaces = TRUE)

# --- Save outputs ---------------------------------------------------------------
write.csv(tab, "results/matching/matchingtable.csv", quote = F, row.names = T)
matchedCases$ID
write.csv(matchedCases$ID, "results/matching/matched_ids.csv", row.names = F)
