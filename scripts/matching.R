# Matching of men and women in BARIA cohort
# Vars to match: age, medication use
# Medication: LLD, metformine, antihypertensiva

# Library
library(tidyverse)
library(MatchIt)
library(ggplot2)
library(tableone)

# Open data
baria <- readRDS("data/bariatot.RDS")
baria <- baria %>% dplyr::filter(!is.na(v4_bmi) & !is.na(TotalCal))
dim(baria)
summary(baria$sex)
baria <- baria %>% filter(v0_antibiotics != "yes" & v0_diarrhoea != "yes")
dim(baria)
summary(baria$sex)

## Matching
malepart <- baria %>% filter(sex == "male") %>% pull(ID)
malesubgroup <- sample(malepart, 50)
matchData <- baria %>% select(ID, v0_age, sex, v0_dmmed, v0_lld, v0_metformine) %>% 
        filter(ID %in% malesubgroup | sex == "female") %>% # all females vs 50 randomly selected males
        drop_na()
matchRes <- matchit(sex ~ v0_age + v0_dmmed + v0_lld + v0_metformine, data = matchData, method="optimal")
matchedCases <- match.data(matchRes)
extradata <- baria %>% select(ID, v0_dm, v0_diabp, v0_sysbp, v0_bmi, v0_smoking) 
matchedCases <- matchedCases %>% left_join(extradata, by="ID")
tableMatched <- CreateTableOne(vars = c('v0_age','sex','v0_sysbp', 'v0_diabp', 'v0_smoking',
                                                'v0_dm','v0_lld','v0_dmmed', 'v0_metformine'), 
                            data = matchedCases, 
                            factorVars = c('sex','v0_dm','v0_smoking', 'v0_lld', 'v0_dmmed', 
                                           'v0_metformine'), 
                            strata = 'sex')
tab <- print(tableMatched,missing = T)
write.csv(tab, "results/matching/matchingtable.csv")
