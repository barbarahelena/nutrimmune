# Matching of men and women in BARIA cohort
# Vars to match: age, medication use
# Medication: LLD, metformine, antihypertensiva

# Library
library(tidyverse)
library(MatchIt)
library(ggplot2)
library(tableone)

# Open data
serum <- readRDS("data/serum_ids.RDS")
rna <- readRDS("data/idsrnaseq.RDS")
bariadiet <- readRDS("data/bariatot.RDS")
baria_all <- readRDS("data/baria_metadata.RDS")
otherdiet <- readxl::read_xlsx("data/extradieetlijsten.xlsx")
otherdiet$ID <- str_c("BARIA_", otherdiet$studienummer)
baria <- baria_all %>% filter(ID %in% c(bariadiet$ID, otherdiet$ID)) %>% 
        filter(ID %in% serum)
dim(baria)
# summary(serum %in% baria$ID)
# notdiet <- serum[!serum %in% baria$ID]

baria <- baria %>% dplyr::filter(!is.na(v4_bmi))
dim(baria)
summary(baria$sex)
baria <- baria %>% filter(v0_antibiotics != "yes" & v0_diarrhoea != "yes") %>% 
        filter(ID %in% serum)
dim(baria)
summary(baria$sex)
summary(baria$ID %in% rna)

## Matching: 50 males vs 50 females
malepart <- baria %>% filter(sex == "male") %>% pull(ID)
malesubgroup <- malepart # sample(malepart, 50)
matchData <- baria %>% select(ID, v0_age, sex, v0_dmmed, v0_lld, v0_metformine) %>% 
        filter(ID %in% malesubgroup | sex == "female") %>% # all females vs 50 randomly selected males
        drop_na()
matchRes <- matchit(sex ~ v0_age + v0_lld + v0_metformine, data = matchData, method="optimal", 
                    ratio = 2.12, max.controls = 10)
matchedCases <- match.data(matchRes)
extradata <- baria %>% select(ID, v0_dm, v0_diabp, v0_sysbp, v0_bmi, v0_leuko, v4_leuko, v0_smoking,
                        v4_bmi, v5_bmi, v0_dm) 
matchedCases <- matchedCases %>% left_join(extradata, by="ID")
tableMatched <- CreateTableOne(vars = c('v0_age','sex', 'v0_bmi', 'v0_sysbp', 'v0_diabp', 'v0_smoking',
                                                'v0_dm','v0_lld','v0_dmmed', 'v0_metformine', 
                                                'v0_leuko', 'v4_leuko', 'v4_bmi', 'v5_bmi'
                                                ), 
                            data = matchedCases, 
                            factorVars = c('sex','v0_dm','v0_smoking', 'v0_lld', 'v0_dmmed', 
                                           'v0_metformine'), 
                            strata = 'sex')
tab <- print(tableMatched,missing = T, noSpaces = TRUE)
write.csv(tab, "results/matching/matchingtable.csv", quote = F, row.names = T)
matchedCases$ID
write.csv(matchedCases$ID, "results/matching/matched_ids.csv", row.names = F)
