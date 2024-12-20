## Data cleaning
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggpubr)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

# Data
datevars <- c("v1_date", "v2_date", "v3_date", "v4_date",
              "v5_date", "v6_date", "v7_date")
clin <- readRDS("data/BARIA_clinical_20241209.RDS")
clin1 <- clin%>% 
    rename(ID = Subject_ID) %>% 
    mutate(ID = str_c("BARIA_", ID)) %>% 
    select(ID, v0_date = `Participant Creation Date`, 
           v1_date = v1_date, v2_date = V2_date, v3_date = V3_date, 
           v4_date = V4_date, v5_date = V5_date, v6_date = V6_date, v7_date = V7_date,
           # V0 data
           v0_age = Age, sex, v0_dm = dm, v0_hypertension = hypertension,
           v0_hyperchol = dyslipidemia_v0, v0_myocardInf = myocardial_infarction_v0,
           v0_ischCVA = CNS_infarction_v0, v0_smoking = smoking, 
           v0_packyears = packyears, v0_alcohol = alcohol, v0_alcoholunits = alcoholunits,
           v0_bmi = bmi, v0_sysbp = sysbp, v0_diabp = diabp, 
           v0_fatperc = tbf_percent, v0_leuko = leukocytes, v0_crp = crp, 
           v0_tnfa = tnfa, v0_ifng = ifng, v0_il1b = il1b, v0_il1ra = il1ra, 
           v0_il6 = il6, v0_il8 = il8, v0_il10 = il10, 
           v0_leptin = leptin, v0_adiponectin = adiponectin,
           v0_glp1 = glp1, v0_pyy = PYY,
           # V4 data
           v4_bmi = V4_bmi, v4_sysbp = V4_sysbp, v4_diabp = V4_diabp,
           v4_leuko = V4_leukocytes, v4_crp = V4_crp, v4_hba1c = V4_hba1c,
           v4_tnfa = V4_tnfa, v4_ifng = V4_ifng, v4_il1b = V4_il1b, v4_il1ra = V4_il1ra, 
           v4_il6 = V4_il6, v4_il8 = V4_il8, v4_il10 = V4_il10, 
           v4_leptin = V4_leptin, v4_adiponectin = V4_adiponectin,
           v4_glp1 = V4_glp1, v4_pyy = V4_PYY,
           # V5 data
           v5_tnfa = V5_tnfa, v5_ifng = V5_ifng, v5_il1b = V5_il1b, v5_il1ra = V5_il1ra, 
           v5_il6 = V5_il6, v5_il8 = V5_il8, v5_il10 = V5_il10, 
           v5_leptin = V5_leptin, v5_adiponectin = V5_adiponectin,
           v5_glp1 = V5_glp1, v5_crp = V5_crp, v5_bmi = V5_bmi
           ) %>% 
    mutate(across(all_of(datevars), 
                  ~ {
                      # Replace invalid dates like "2999-01-01" with NA
                      .x <- ifelse(.x == "01-01-2999", NA, .x)
                      # Convert to Date using dmy after replacing invalid dates
                      dmy(.x)
                  }, 
                  .names = "cv_{.col}"),
           v0_date = as_date(dmy_hms(v0_date)),
           across(all_of(datevars), ~ifelse(.x > ymd("2025-01-01"), NA, .x)),
           v1v3_diff = difftime(cv_v3_date, cv_v1_date, units = "days"),
           v1v2_diff = difftime(cv_v2_date, cv_v1_date, units = "days"),
           v1v4_diff = difftime(cv_v4_date, cv_v1_date, units = "days"),
           v1v5_diff = difftime(cv_v5_date, cv_v1_date, units = "days"),
           crp_v1v4 = v4_crp - v0_crp,
           crp_v1v5 = v5_crp - v0_crp,
           bmi_v1v4 = v4_bmi - v0_bmi,
           bmi_v1v5 = v5_bmi - v0_bmi
           )

clinv4 <- clin %>% select(contains("V4"), contains("v4"))
names(clinv4)
clinv1 <- clin %>% select(contains("V2"), contains("v1"))
names(clinv1)
clin[1:5,1:5]
Amelia::missmap(clin1)

baria <- rio::import("data/baria_dieet.xlsx") 
    select(ID = naam, Date = datum, Meal = moment,
           TotalCal = `energie (kcal)`, Protein = `eiwit totaal (g)`,
           Carbs = `koolhydraten totaal (g)`, Fat = `vet totaal (g)`,
           SatFat = `vetzuren verzadigd (g)`, UnsatFat = `vetzuren onverzadigd (g)`,
           TransFat = `trans vetzuren (g)`, LinolicAcid = `linolzuur (g)`,
           Cholesterol = `cholesterol (mg)`, Alcohol = `alcohol (g)`,
           Fibers = `voedingsvezels (g)`, Water = `water (ml)`, MonoDiSacch = `tot mono disach (g)`,
           PolySacch = `polysacchariden (g)`, PolyunsatFat = `vetzuren meerv onverzadigd (g)`) %>% 
    mutate(
        ID = str_c("BARIA_", ID),
        DateTime = ymd_hms(Date),
        Time = format(DateTime, "%H:%M:%S"),
        Date = as_date(DateTime),
        Year = format(Date, "%Y"),
        Water = case_when(Water > 25000 ~ Water / 1000, .default = Water)
    )
            
diabar <- rio::import("data/diabar_dieet.xlsx") %>% 
    select(ID = naam, Date = datum, Meal = moment,
           TotalCal = `energie (kcal)`, Protein = `eiwit totaal (g)`,
           Carbs = `koolhydraten totaal (g)`, Fat = `vet totaal (g)`,
           SatFat = `vetzuren verzadigd (g)`, UnsatFat = `vetzuren onverzadigd (g)`,
           TransFat = `trans vetzuren (g)`, LinolicAcid = `linolzuur (g)`,
           Cholesterol = `cholesterol (mg)`, Alcohol = `alcohol (g)`,
           Fibers = `voedingsvezels (g)`, Water = `water (ml)`, MonoDiSacch = `tot mono disach (g)`,
           PolySacch = `polysacchariden (g)`, PolyunsatFat = `vetzuren meerv onverzadigd (g)`) %>% 
    mutate(
        ID = str_c("DIABAR_", ID),
        DateTime = ymd_hms(Date),
        Time = format(DateTime, "%H:%M:%S"),
        Date = as_date(DateTime),
        Year = format(Date, "%Y")
    )


bariaday <- baria %>% 
    group_by(ID, Year, Date) %>% 
    summarise(across(TotalCal:PolyunsatFat, sum)) %>% 
    group_by(ID, Year) %>% 
    summarise(across(TotalCal:PolyunsatFat, mean))

bariayear <- baria %>%
    group_by(ID, Year, Date) %>%
    summarise(across(TotalCal:PolyunsatFat, ~ sum(.x, na.rm = TRUE), .names = "daily_sum_{col}"), .groups = "drop") %>%
    group_by(ID, Year) %>%
    mutate(across(starts_with("daily_sum"), ~ slider::slide_dbl(.x, mean, .before = 2, .complete = FALSE), 
                  .names = "rolling_avg_{col}"))
head(bariayear)

bariabase <- bariayear %>% group_by(ID) %>% 
    filter(Year == min(Year)) %>%
    slice_tail(n = 1) %>% 
    ungroup(.) %>% 
    select(-starts_with("daily_sum_")) %>% 
    rename_with(~gsub("^rolling_avg_daily_sum_", "", .), starts_with("rolling_avg_"))

bariatot <- left_join(bariabase, clin1)    

crp <- bariatot %>% filter(v0_crp < 100 & v4_crp < 100) %>% filter(TotalCal < 5000)
bmi <- bariatot %>% filter(TotalCal < 5000)

plist <- c()
for(a in names(crp)[4:18]) {
    crp$var <- crp[[a]]
    print(gghistogram(crp$var, title = a, fill = ggsci::pal_bmj()(1)))
    pl <- ggplot(data = crp, aes(x = var, y = crp_v1v4)) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "crp change at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])

plist <- c()
for(a in names(crp)[4:18]) {
    crp$var <- crp[[a]]
    pl <- ggplot(data = crp, aes(x = var, y = v4_crp)) +
        geom_point(alpha = 0.7) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi change at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])

plist <- c()
for(a in names(bmi)[4:18]) {
    bmi$var <- bmi[[a]]
    pl <- ggplot(data = bmi, aes(x = var, y = bmi_v1v4)) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi change at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])

plist <- c()
for(a in names(bmi)[4:18]) {
    bmi$var <- bmi[[a]]
    pl <- ggplot(data = bmi, aes(x = var, y = bmi_v1v5)) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi change at 2 years",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
