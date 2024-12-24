## Table 1
## Barbara Verhaar

library(tableone)
library(tidyverse)

## Data 
baria <- readRDS("data/bariatot.RDS")

table1 <- baria %>% filter(!is.na(v0_age)) %>% 
    dplyr::select(v0_age, sex, v0_dm, v0_hypertension, v0_hyperchol, v0_myocardInf,
                  v0_smoking, v0_alcohol, v0_alcoholunits, v0_bmi, v0_sysbp, v0_diabp,
                  v0_fatperc, v0_crp, v0_leuko, bmi_v1v4, bmi_v1v5,
                  bmiperc_v1v4, bmiperc_v1v5) %>% 
    CreateTableOne(data=., test = FALSE) %>% 
    print(nonnormal=c("v0_crp", "v0_alcoholunits", "v0_fatperc"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)
write.csv(table1, "results/table1/table1.csv")

## BMI changes after surgery
bmiv0v4 <- baria %>% select(ID, v0_bmi, v2_bmi, v3_bmi, v4_bmi, v5_bmi) %>% 
    pivot_longer(., cols = c(v0_bmi, v2_bmi, v3_bmi, v4_bmi, v5_bmi), names_to = c("timepoint", "var"), 
                 names_sep = "_", values_to = "bmi") %>% 
    mutate(timepoint = case_when(timepoint == "v0" ~ "baseline", 
                                 timepoint == "v2" ~ "6 weeks",
                                 timepoint == "v3" ~ "6 months",
                                 timepoint == "v4" ~ "1 year postop",
                                 timepoint == "v5" ~ "2 years postop"),
           timepoint = fct_relevel(timepoint, "baseline", after = 0L),
           timepoint = fct_relevel(timepoint, "6 weeks", after = 1L),
           timepoint = fct_relevel(timepoint, "6 months", after = 2L),
           )
ggplot(data = bmiv0v4, aes(x = timepoint, y = bmi, fill = timepoint)) +
    geom_violin() +
    geom_boxplot(fill = "white", width = 0.25, outlier.shape = NA) +
    theme_Publication() +
    scale_fill_bmj(guide = "none") +
    labs(x = "", y = "BMI (kg/m2)", title = "Change in BMI")

bmiv0v4 <- baria %>% select(ID, v0_bmi, v2_bmi, v3_bmi, v4_bmi, v5_bmi) %>% 
    mutate(across(c(v0_bmi, v2_bmi, v3_bmi, v4_bmi, v5_bmi), ~ (.x / v0_bmi) *100)) %>% 
    pivot_longer(., cols = c(v0_bmi, v2_bmi, v3_bmi, v4_bmi, v5_bmi), names_to = c("timepoint", "var"), 
                 names_sep = "_", values_to = "bmi") %>% 
    mutate(timepoint = case_when(timepoint == "v0" ~ "baseline", 
                                 timepoint == "v2" ~ "6 weeks",
                                 timepoint == "v3" ~ "6 months",
                                 timepoint == "v4" ~ "1 year postop",
                                 timepoint == "v5" ~ "2 years postop"),
           timepoint = fct_relevel(timepoint, "baseline", after = 0L),
           timepoint = fct_relevel(timepoint, "6 weeks", after = 1L),
           timepoint = fct_relevel(timepoint, "6 months", after = 2L),
    )
ggplot(data = bmiv0v4, aes(x = timepoint, y = bmi, fill = timepoint)) +
    geom_violin() +
    geom_boxplot(fill = "white", width = 0.25, outlier.shape = NA) +
    theme_Publication() +
    scale_fill_bmj(guide = "none") +
    labs(x = "", y = "BMI % of baseline", title = "Percentage change in BMI")


## CRP changes after surgery
crp <- baria %>% select(ID, v0_crp, v2_crp, v3_crp, v4_crp, v5_crp) %>% 
    pivot_longer(., cols = c(v0_crp, v2_crp, v3_crp, v4_crp, v5_crp), names_to = c("timepoint", "var"), 
                 names_sep = "_", values_to = "crp") %>% 
    mutate(timepoint = case_when(timepoint == "v0" ~ "baseline", 
                                 timepoint == "v2" ~ "6 weeks",
                                 timepoint == "v3" ~ "6 months",
                                 timepoint == "v4" ~ "1 year postop",
                                 timepoint == "v5" ~ "2 years postop"),
           timepoint = fct_relevel(timepoint, "baseline", after = 0L),
           timepoint = fct_relevel(timepoint, "6 weeks", after = 1L),
           timepoint = fct_relevel(timepoint, "6 months", after = 2L),
    ) %>% filter(crp < 50 & timepoint != "2 years postop")
ggplot(data = crp, aes(x = timepoint, y = log10(crp+0.1), fill = timepoint)) +
    geom_violin() +
    geom_boxplot(fill = "white", width = 0.25, outlier.shape = NA) +
    theme_Publication() +
    scale_fill_bmj(guide = "none") +
    labs(x = "", y = "log10(CRP)", title = "Change in CRP")

leuko <- baria %>% select(ID, v0_leuko, v2_leuko, v3_leuko, v4_leuko, v5_leuko) %>% 
    pivot_longer(., cols = c(v0_leuko, v2_leuko, v3_leuko, v4_leuko, v5_leuko), 
                 names_to = c("timepoint", "var"), 
                 names_sep = "_", values_to = "leuko") %>% 
    mutate(timepoint = case_when(timepoint == "v0" ~ "baseline", 
                                 timepoint == "v2" ~ "6 weeks",
                                 timepoint == "v3" ~ "6 months",
                                 timepoint == "v4" ~ "1 year postop",
                                 timepoint == "v5" ~ "2 years postop"),
           timepoint = fct_relevel(timepoint, "baseline", after = 0L),
           timepoint = fct_relevel(timepoint, "6 weeks", after = 1L),
           timepoint = fct_relevel(timepoint, "6 months", after = 2L),
    ) %>% filter(leuko < 40)
ggplot(data = leuko, aes(x = timepoint, y = leuko, fill = timepoint)) +
    geom_violin() +
    geom_boxplot(fill = "white", width = 0.25, outlier.shape = NA) +
    theme_Publication() +
    scale_fill_bmj(guide = "none") +
    labs(x = "", y = "Leukocytes", title = "Change in leukocytes")


homair <- baria %>% select(ID, homaIR_v1, homaIR_v4, homaIR_v5) %>% 
    pivot_longer(., cols = c(homaIR_v1, homaIR_v4, homaIR_v5), 
                 names_to = c("var", "timepoint"), 
                 names_sep = "_", values_to = "homair") %>% 
    mutate(timepoint = case_when(timepoint == "v1" ~ "baseline", 
                                 timepoint == "v4" ~ "1 year postop",
                                 timepoint == "v5" ~ "2 years postop"),
           timepoint = fct_relevel(timepoint, "baseline", after = 0L)
    )
ggplot(data = homair, aes(x = timepoint, y = log10(homair), fill = timepoint)) +
    geom_violin() +
    geom_boxplot(fill = "white", width = 0.25, outlier.shape = NA) +
    theme_Publication() +
    scale_fill_bmj(guide = "none") +
    labs(x = "", y = "log10(HOMA-IR)", title = "Change in HOMA-IR")
ggsave("results/table1/homair.pdf", width = 4, height = 5)
