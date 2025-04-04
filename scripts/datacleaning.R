## Data cleaning
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)

# Data
datevars <- c("v1_date", "v2_date", "v3_date", "v4_date",
              "v5_date", "v6_date", "v7_date")
clin <- readRDS("data/BARIA_clinical_20241209.RDS")
clin1 <- clin %>% 
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
           v0_glp1 = glp1, v0_pyy = PYY, v0_ldl = ldlchol, v0_hdl = hdlchol, 
           v0_totchol = totchol,
           # V1 data?
           # V2 data
           v2_bmi = V2_bmi, v2_crp = V2_crp, v2_leuko = V2_leukocytes,
           # V3 data
           v3_bmi = V3_bmi, v3_crp = V3_crp, v3_leuko = V3_leukocytes,
           v3_totchol = V3_totchol, v3_ldl = vis3_ldl, v3_totchol = V3_totchol,
           # V4 data
           v4_bmi = V4_bmi, v4_sysbp = V4_sysbp, v4_diabp = V4_diabp,
           v4_leuko = V4_leukocytes, v4_crp = V4_crp, v4_hba1c = V4_hba1c,
           v4_tnfa = V4_tnfa, v4_ifng = V4_ifng, v4_il1b = V4_il1b, v4_il1ra = V4_il1ra, 
           v4_il6 = V4_il6, v4_il8 = V4_il8, v4_il10 = V4_il10, 
           v4_leptin = V4_leptin, v4_adiponectin = V4_adiponectin,
           v4_glp1 = V4_glp1, v4_pyy = V4_PYY, v4_ldl = V4_ldlchol, 
           v4_totchol = V4_totchol, v4_hdlchol = V4_hdlchol,
           # V5 data
           v5_tnfa = V5_tnfa, v5_ifng = V5_ifng, v5_il1b = V5_il1b, v5_il1ra = V5_il1ra, 
           v5_il6 = V5_il6, v5_il8 = V5_il8, v5_il10 = V5_il10, v5_leuko = V5_leukocytes,
           v5_leptin = V5_leptin, v5_adiponectin = V5_adiponectin,
           v5_glp1 = V5_glp1, v5_crp = V5_crp, v5_bmi = V5_bmi, 
           v5_ldlchol = V5_ldlchol, v5_hdlchol = V5_hdlchol, v5_totchol = V5_totchol,
           contains("insulin"),
           contains("gluc")
           ) %>% 
    mutate(across(all_of(datevars), 
                  ~ {.x <- case_when(.x == "01-01-2999" | .x == "01-01-2995" ~ NA, 
                                     .default = .x)
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
           bmi_v1v5 = v5_bmi - v0_bmi,
           v2_bmi = case_when(v2_bmi < 10 ~ NA, .default = v2_bmi),
           bmiperc_v1v4 = (bmi_v1v4/v0_bmi) *100,
           bmiperc_v1v5 = (bmi_v1v5/v0_bmi) *100,
           across(contains("insulin"), ~ .x / 6),
           homaIR_v1 = (min0insulin * min0gluc) / 22.5,
           homaIR_v4 = (V4_min0insulin * V4_min0gluc) / 22.5,
           homaIR_v5 = (V5_min0insulin * V5_min0gluc) / 22.5,
           homaIR_v1v4 = homaIR_v4 - homaIR_v1,
           homaIR_v1v5 = homaIR_v5 - homaIR_v1,
           sex = case_when(sex == 2 ~ "female", sex == 1 ~ "male"),
           v0_dm = case_when(v0_dm == 1 ~ "yes", v0_dm == 2 ~ "no"),
           v0_hypertension = case_when(v0_hypertension == 1 ~ "yes", v0_hypertension == 2 ~ "no"),
           v0_hyperchol = case_when(v0_hyperchol == 1 ~ "yes", v0_hyperchol == 2 ~ "no"),
           v0_myocardInf = case_when(v0_myocardInf == 1 ~ "yes", v0_myocardInf == 2 ~ "no"),
           v0_smoking = case_when(v0_smoking == 1 ~ "current smoking", v0_smoking == 2 ~ "former smoking", 
                                  v0_smoking == 3 ~ "never"),
           v0_alcohol = case_when(v0_alcohol == 1 ~ "yes", v0_alcohol == 2 ~ "no"),
           across(where(is.character), as.factor)
           )
saveRDS(clin1, "data/baria_metadata.RDS")
write.csv(clin1, "data/baria_metadata.csv")

clinv4 <- clin %>% select(contains("V4"), contains("v4"))
names(clinv4)
clinv1 <- clin %>% select(contains("V2"), contains("v1"))
names(clinv1)
clin[1:5,1:5]
Amelia::missmap(clin1)

baria <- rio::import("data/baria_dieet.xlsx") %>% 
    select(ID = naam, Date = datum, Meal = moment,
           TotalCal = `energie (kcal)`, 
           Protein = `eiwit totaal (g)`,
           Carbs = `koolhydraten totaal (g)`, 
           Fat = `vet totaal (g)`,
           SatFat = `vetzuren verzadigd (g)`, 
           UnsatFat = `vetzuren onverzadigd (g)`,
           TransFat = `trans vetzuren (g)`, 
           LinolicAcid = `linolzuur (g)`,
           Cholesterol = `cholesterol (mg)`, 
           Alcohol = `alcohol (g)`,
           Fibers = `voedingsvezels (g)`, 
           Water = `water (ml)`, 
           MonoDiSacch = `tot mono disach (g)`,
           PolySacch = `polysacchariden (g)`, 
           PolyunsatFat = `vetzuren meerv onverzadigd (g)`
           ) %>% 
    mutate(
        ID = str_c("BARIA_", ID),
        DateTime = ymd_hms(Date),
        Time = format(DateTime, "%H:%M:%S"),
        Date = as_date(DateTime),
        Year = format(Date, "%Y"),
        Water = case_when(Water > 25000 ~ Water / 1000, .default = Water), # correct unit
    )

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
        ID = str_c("BARIA_", ID),
        DateTime = ymd_hms(Date),
        Time = format(DateTime, "%H:%M:%S"),
        Date = as_date(DateTime),
        Year = format(Date, "%Y")
    )

diabaryear <- diabar %>%
    group_by(ID, Year, Date) %>%
    summarise(across(TotalCal:PolyunsatFat, ~ sum(.x, na.rm = TRUE), .names = "daily_sum_{col}"), .groups = "drop") %>%
    group_by(ID, Year) %>%
    mutate(across(starts_with("daily_sum"), ~ slider::slide_dbl(.x, mean, .before = 2, .complete = FALSE), 
                  .names = "rolling_avg_{col}"))
head(diabaryear)

diabarbase <- diabaryear %>% group_by(ID) %>% 
    filter(Year == min(Year)) %>%
    slice_tail(n = 1) %>% 
    ungroup(.) %>% 
    select(-starts_with("daily_sum_")) %>% 
    rename_with(~gsub("^rolling_avg_daily_sum_", "", .), starts_with("rolling_avg_"))

diettot <- rbind(bariabase, diabarbase)
bariatot <- left_join(diettot, clin1) %>% filter(TotalCal < 5000)  %>% 
    mutate(
        baselineYear = format(dmy(v1_date), "%Y"),
        diffDate = time_length(interval(start = Date, end = dmy(v1_date)), unit = "days"),
        TotalCalBin = case_when(TotalCal < median(TotalCal, na.rm = TRUE) ~ 
                                    str_c("<", round(median(TotalCal, na.rm = TRUE), 0), " kcal"),
                                TotalCal >= median(TotalCal, na.rm = TRUE) ~ 
                                    str_c(">=", round(median(TotalCal, na.rm = TRUE), 0), " kcal")),
        ProteinBin = case_when(Protein < median(Protein, na.rm = TRUE) ~ 
                                   str_c("<", round(median(Protein, na.rm = TRUE), 0), " g"),
                               Protein >= median(Protein, na.rm = TRUE) ~ 
                                   str_c(">=", round(median(Protein, na.rm = TRUE), 0), " g")),
        FibersBin = case_when(Fibers < median(Fibers, na.rm = TRUE) ~ 
                                  str_c("<", round(median(Fibers, na.rm = TRUE),0), " g"),
                              Fibers >= median(Fibers, na.rm = TRUE) ~ 
                                  str_c(">=", round(median(Fibers, na.rm = TRUE), 0), " g")),
        CarbsBin = case_when(Carbs < median(Carbs, na.rm = TRUE) ~ 
                                 str_c("<", round(median(Carbs, na.rm = TRUE), 0), " g"),
                             Carbs >= median(Carbs, na.rm = TRUE) ~ 
                                 str_c(">=", round(median(Carbs, na.rm = TRUE),0), " g")),
        FatBin = case_when(Carbs < median(Fat, na.rm = TRUE) ~ 
                               str_c("<", round(median(Fat, na.rm = TRUE), 0), " g"),
                           Carbs >= median(Fat, na.rm = TRUE) ~ 
                               str_c(">=", round(median(Fat, na.rm = TRUE), 0), " g")),
        SatFatBin = case_when(SatFat < median(SatFat, na.rm = TRUE) ~ 
                                  str_c("<", round(median(SatFat, na.rm = TRUE),0), " g"),
                              SatFat >= median(SatFat, na.rm = TRUE) ~ 
                                  str_c(">=", round(median(SatFat, na.rm = TRUE),0), " g")),
        UnsatFatBin = case_when(UnsatFat < median(UnsatFat, na.rm = TRUE) ~ 
                                    str_c("<", round(median(UnsatFat, na.rm = TRUE),0), " g"),
                                UnsatFat >= median(UnsatFat, na.rm = TRUE) ~ 
                                    str_c(">=", round(median(UnsatFat, na.rm = TRUE),0), " g")),
        TransFatBin = case_when(TransFat < median(TransFat, na.rm = TRUE) ~ 
                                    str_c("<", round(median(TransFat, na.rm = TRUE),0), " g"),
                                TransFat >= median(TransFat, na.rm = TRUE) ~ 
                                    str_c(">=", round(median(TransFat, na.rm = TRUE),0), " g")),
        LinolicAcidBin = case_when(LinolicAcid < median(LinolicAcid, na.rm = TRUE) ~ 
                                       str_c("<", round(median(LinolicAcid, na.rm = TRUE),0), " g"),
                                   LinolicAcid >= median(LinolicAcid, na.rm = TRUE) ~ 
                                       str_c(">=", round(median(LinolicAcid, na.rm = TRUE),0), " g")),
        CholesterolBin = case_when(Cholesterol < median(Cholesterol, na.rm = TRUE) ~ 
                                       str_c("<", round(median(Cholesterol, na.rm = TRUE),0), " mg"),
                                   Cholesterol >= median(Cholesterol, na.rm = TRUE) ~ 
                                       str_c(">=", round(median(Cholesterol, na.rm = TRUE),0), " mg")),
        AlcoholBin = case_when(Alcohol == 0 ~ "0 g", Alcohol > 0 ~ ">0 g"),
        MonoDiSacchBin = case_when(MonoDiSacch < median(MonoDiSacch, na.rm = TRUE) ~ 
                                       str_c("<", round(median(MonoDiSacch, na.rm = TRUE),0), " g"),
                                   MonoDiSacch >= median(MonoDiSacch, na.rm = TRUE) ~ 
                                       str_c(">=", round(median(MonoDiSacch, na.rm = TRUE),0), " g")),
        PolySacchBin = case_when(PolySacch < median(PolySacch, na.rm = TRUE) ~ 
                                     str_c("<", round(median(PolySacch, na.rm = TRUE),0), " g"),
                                 PolySacch >= median(PolySacch, na.rm = TRUE) ~ 
                                     str_c(">=", round(median(PolySacch, na.rm = TRUE),0), " g")),
        PolyunsatFatBin = case_when(PolyunsatFat < median(PolyunsatFat, na.rm = TRUE) ~ 
                                        str_c("<", round(median(PolyunsatFat, na.rm = TRUE),0), " g"),
                                    PolyunsatFat >= median(PolyunsatFat, na.rm = TRUE) ~ 
                                        str_c(">=", round(median(PolyunsatFat, na.rm = TRUE),0), " g")),
        across(contains("Bin"), as.factor)
    ) %>% filter(diffDate >= 0 & diffDate < 365)

saveRDS(bariatot, "data/bariatot.RDS")

## Gene list
gene_list <- read.csv("data/ICP_list.csv", sep = ';')
gene_list$ICP_symbol <- trimws(gene_list$ICP_symbol)
head(gene_list)
saveRDS(gene_list, "data/gene_list.RDS")
