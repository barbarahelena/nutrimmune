## Data cleaning
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(dplyr)
library(stringr)
library(lubridate)
library(ggpubr)

# Data
datevars <- c("v1_date", "v2_date", "v3_date", "v4_date",
              "v5_date", "v6_date", "v7_date")
lipidlowering <- c(
  "atorvastatine", "lipitor",
  "rosuvastatine", "crestor",
  "simvastatine", "zocor",
  "pravastatine", "pravachol",
  "fluvastatine", "lescol",
  "ezetimib", "ezetrol",
  "alirocumab", "praluent",
  "evolocumab", "repatha",
  "inclisiran", "leqvio",
  "bempedoïnezuur", "nexletol", "nustendi",
  "colestyramine", "questran",
  "fenofibraat", "lipanthyl", "tricor",
  "gemfibrozil", "lopid",
  "omega-3-vetzuren", "omacor", "epanova",
  "nicotinezuur", "niaspan"
)
lld_pattern <- str_c("\\b(", str_c(tolower(lipidlowering), collapse = "|"), ")\\b")
dmmed <- c(
  # Metformine
  "metformine", "glucophage",
  
  # Sulfonylureumderivaten
  "gliclazide", "diamicron",
  "glibenclamide", "daonil",
  "glimepiride", "amaryl",
  "glipizide", "minodiab",
  
  # DPP-4-remmers
  "sitagliptine", "januvia",
  "vildagliptine", "galvus",
  "saxagliptine", "onglyza",
  "linagliptine", "trajenta",
  "alogliptine", "vipidia",
  
  # GLP-1-receptoragonisten
  "liraglutide", "victoza",
  "semaglutide", "ozempic", "rybelsus",
  "exenatide", "byetta", "bydureon",
  "dulaglutide", "trulicity",
  
  # SGLT2-remmers
  "dapagliflozine", "forxiga",
  "empagliflozine", "jardiance",
  "canagliflozine", "invokana",
  "ertugliflozine", "steeglatro",
  
  # Insuline (kort, middellang en langwerkend)
  "insuline", "lantus", "levemir", "novorapid", "apidra", "toujeo", "tresiba", "humalog", "novomix", "fiasp", "actrapid"
)
dmmed_pattern <- str_c("\\b(", str_c(tolower(dmmed), collapse = "|"), ")\\b")
metforminemed <- c("metformine", "glucophage", "glucient", "metformax")
metforminemed_pattern <- str_c("\\b(", str_c(tolower(metforminemed), collapse = "|"), ")\\b")

clin <- readRDS("data/BARIA_clinical_20241209.RDS")
clin1 <- clin %>% 
    dplyr::rename(ID = Subject_ID) %>% 
    dplyr::mutate(ID = str_c("BARIA_", ID)) %>% 
    dplyr::select(ID, v0_date = `Participant Creation Date`, 
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
           v0_totchol = totchol, v0_antibiotics = scre_ab, v0_diarrhoea = scre_chrondiarree,
           v0_probiotics = scre_probiotic, v0_medication = medication_v0_freetext,
           # V1 data
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
    mutate(
           across(all_of(datevars), ~ case_when(
                .x == "01-01-2999" | .x == "01-01-2995" ~ NA_character_,
                TRUE ~ .x
            ) %>% dmy(), .names = "cv_{.col}"),
           v0_date = as_date(dmy_hms(v0_date)),
           across(all_of(datevars), ~ifelse(.x > ymd("2025-01-01"), NA, .x)),
           v0_antibiotics = case_when(v0_antibiotics == "1" ~ "yes", 
                                        v0_antibiotics == "2" ~ "no"),
           v0_diarrhoea = case_when(v0_diarrhoea == "1" ~ "yes", 
                                    v0_diarrhoea == "2" ~ "no"),
           v0_probiotics = case_when(v0_probiotics == "1" ~ "yes", 
                                        v0_probiotics == "2" ~ "no"),
           v1v3_diff = difftime(cv_v3_date, cv_v1_date, units = "days"),
           v1v2_diff = difftime(cv_v2_date, cv_v1_date, units = "days"),
           v1v4_diff = difftime(cv_v4_date, cv_v1_date, units = "days"),
           v1v5_diff = difftime(cv_v5_date, cv_v1_date, units = "days"),
           crp_v1v4 = v4_crp - v0_crp,
           crp_v1v5 = v5_crp - v0_crp,
           bmi_v1v4 = v4_bmi - v0_bmi,
           bmi_v1v5 = v5_bmi - v0_bmi,
           v2_bmi = case_when(
                v2_bmi < 10 ~ NA_real_,
                .default = v2_bmi
            ),
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
           v0_lld = case_when(str_detect(tolower(v0_medication), lld_pattern) ~ "yes", .default = "no"),
           v0_dmmed = case_when(str_detect(tolower(v0_medication), dmmed_pattern) ~ "yes", .default = "no"),
           v0_metformine = case_when(str_detect(tolower(v0_medication), metforminemed_pattern) ~ "yes", .default = "no"),
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

baria <- rio::import("data/baria_ids_nutrimmune_b.xlsx") %>% 
    select(ID = naam, Date = datum, Visit = visite, 
           Meal = moment,
           TotalCal_kcal = `energie (kcal)`, 
           Protein_g = `eiwit totaal (g)`,
           Carbs_g = `koolhydraten totaal (g)`, 
           Fat_g = `vet totaal (g)`,
           SatFat_g = `vetzuren verzadigd (g)`, 
           UnsatFat_g = `vetzuren onverzadigd (g)`,
           TransFat_g = `trans vetzuren (g)`, 
           LinolicAcid_g = `linolzuur (g)`,
           Cholesterol_mg = `cholesterol (mg)`, 
           Alcohol_g = `alcohol (g)`,
           Fibers_g = `voedingsvezels (g)`, 
           Water_ml = `water (ml)`, 
           MonoDiSacch_g = `tot mono disach (g)`,
           PolySacch_g = `polysacchariden (g)`, 
           PolyunsatFat_g = `vetzuren meerv onverzadigd (g)`) %>% 
    mutate(
        ID = str_c("BARIA_", ID),
        Visit = str_remove_all(Visit, "\\–"),
        Visit = as.factor(Visit),
        DateTime = ymd_hms(Date),
        Time = format(DateTime, "%H:%M:%S"),
        Date = as_date(DateTime),
        Year = format(Date, "%Y")
    ) %>% 
    filter(!is.na(Date) & !is.na(Visit))
names(baria)
dim(baria)
baria$Visit[1:50]

bariayear <- baria %>%
    group_by(ID, Visit, Date) %>%
    summarise(across(TotalCal_kcal:PolyunsatFat_g, 
        ~ sum(.x, na.rm = TRUE), .names = "daily_sum_{col}"), .groups = "drop") %>%
    ungroup() %>% 
    group_by(ID, Visit) %>%
    mutate(across(starts_with("daily_sum"), ~ mean(.x, na.rm = TRUE), .names = "mean_{col}")) |> 
    mutate(n_days = n())
head(bariayear)
dim(bariayear)

# baria <- rio::import("data/baria_dieet.xlsx") %>% 
#     select(ID = naam, Date = datum, Meal = moment,
#            TotalCal = `energie (kcal)`, 
#            Protein = `eiwit totaal (g)`,
#            Carbs = `koolhydraten totaal (g)`, 
#            Fat = `vet totaal (g)`,
#            SatFat = `vetzuren verzadigd (g)`, 
#            UnsatFat = `vetzuren onverzadigd (g)`,
#            TransFat = `trans vetzuren (g)`, 
#            LinolicAcid = `linolzuur (g)`,
#            Cholesterol = `cholesterol (mg)`, 
#            Alcohol = `alcohol (g)`,
#            Fibers = `voedingsvezels (g)`, 
#            Water = `water (ml)`, 
#            MonoDiSacch = `tot mono disach (g)`,
#            PolySacch = `polysacchariden (g)`, 
#            PolyunsatFat = `vetzuren meerv onverzadigd (g)`
#            ) %>% 
#     mutate(
#         ID = str_c("BARIA_", ID),
#         DateTime = ymd_hms(Date),
#         Time = format(DateTime, "%H:%M:%S"),
#         Date = as_date(DateTime),
#         Year = format(Date, "%Y"),
#         Water = case_when(Water > 25000 ~ Water / 1000, .default = Water), # correct unit
#     )

# bariayear <- baria %>%
#     group_by(ID, Year, Date) %>%
#     summarise(across(TotalCal:PolyunsatFat, ~ sum(.x, na.rm = TRUE), .names = "daily_sum_{col}"), .groups = "drop") %>%
#     group_by(ID, Year) %>%
#     mutate(across(starts_with("daily_sum"), ~ slider::slide_dbl(.x, mean, .before = 2, .complete = FALSE), 
#                   .names = "rolling_avg_{col}"))
# head(bariayear)

bariabase <- bariayear %>% group_by(ID) %>% 
    filter(Visit == "V1") |> 
    slice(1) %>% 
    select(-starts_with("daily_sum_"), n_days) %>% 
    rename_with(~gsub("^mean_daily_sum_", "", .), starts_with("mean_")) |> 
    mutate(diary_tool = "Bonstat")
dim(bariabase)
head(bariabase)

bariallyears <- bariayear %>% group_by(ID, Visit) %>% 
    slice(1) %>% 
    select(-starts_with("daily_sum_"), n_days) %>% 
    rename_with(~gsub("^mean_daily_sum_", "", .), starts_with("mean_"))
dim(bariallyears)
head(bariallyears)

baria_eetmeter <- rio::import("data/BARIA_data_H1.xlsx") |> 
    select(ID = `Deelnemer ID + Dag`, Day = `Deelnemer ID + Dag`, Visit = `Deelnemer ID + Dag`, 
           Meal = Eetmoment,
           TotalCal_kcal = `Energie (kcal)`, 
           Protein_g = `Eiwit (g)`,
           Carbs_g = `Kolhydr (g)`, 
           Fat_g = `Vet (g)`,
           SatFat_g = `Verz. vet (g)`, 
           Alcohol_g = `Alcohl (g)`,
           Fibers_g = `Vezels (g)`, 
           Water_ml = `Water (g)`,
           Salt_g = `Zout (g)`,
           Sodium_mg = `Natrium (g)`,
           Potassium_mg = `Kalium (mg)`,
           Calcium_mg = `Calcium (mg)`,
           Magnesium_mg = `Magnesium (mg)`,
           Iron_mg = `IJzer (mg)`,
           Selenium_ug = `Selenium (µg)`,
           Zinc_mg = `Zink (mg)`,
           VitaminA_ug = `Vit. A (µg)`,
           VitaminD_ug = `Vit. D (µg)`,
           VitaminE_mg = `Vit. E (mg)`,
           VitaminB1_mg = `Vit. B1 (mg)`,
           VitaminB2_mg = `Vit. B2 (mg)`,
           VitaminB6_mg = `Vit. B6 (mg)`,
           FolicAcid_ug = `Foliumzuur (µg)`,
           VitaminB12_ug = `Vit. B12 (µg)`,
           Niacin_mg = `Niacine (mg)`,
           VitaminC_mg = `Vit. C (mg)`,
           Iodine_ug = `Jodium (µg)`,
           Phosphorus_mg = `Fosfor (mg)`,
           Sugars_g = `Suikers (g)`
        ) |> 
    mutate(
        ID = str_c("BARIA_", str_extract(ID, "(?<=\\s)\\d+(?=\\s)")),
        Visit = str_extract(Visit, "^\\S+"),
        Visit = as.factor(Visit),
        Day = str_extract(Day, "D\\d+"),
        Day = as.factor(Day)
    ) |> 
    mutate(across(TotalCal_kcal:Sugars_g, as.numeric))
# names(baria)
# dim(baria)
# baria$Visit[1:50]

bariayear_eetmeter <- baria_eetmeter %>%
    group_by(ID, Visit, Day) %>%
    summarise(across(TotalCal_kcal:Sugars_g, 
        ~ sum(.x, na.rm = TRUE), .names = "daily_sum_{col}"), .groups = "drop") %>%
    ungroup() %>% 
    group_by(ID, Visit) %>%
    mutate(across(starts_with("daily_sum"), 
        ~ mean(.x, na.rm = TRUE), .names = "mean_{col}")) |> 
    mutate(n_days = n())
head(bariayear_eetmeter)
dim(bariayear_eetmeter)

bariabase_eetmeter <- bariayear_eetmeter %>% group_by(ID) %>% 
    filter(Visit == "V1") %>% 
    slice(1) %>% 
    select(-starts_with("daily_sum_"), n_days) %>% 
    rename_with(~gsub("^mean_daily_sum_", "", .), starts_with("mean_")) |> 
    mutate(diary_tool = "Eetmeter")
dim(bariabase_eetmeter)
head(bariabase_eetmeter)



# bariallyears <- bariayear %>% group_by(ID, Visit) %>% 
#     slice(1) %>% 
#     select(-starts_with("daily_sum_")) %>% 
#     rename_with(~gsub("^mean_daily_sum_", "", .), starts_with("mean_"))
# dim(bariallyears)
# head(bariallyears)

diabar <- rio::import("data/diabar_dieet.xlsx") %>% 
    select(ID = naam, Date = datum, Meal = moment,
           TotalCal_kcal = `energie (kcal)`, Protein_g = `eiwit totaal (g)`,
           Carbs_g = `koolhydraten totaal (g)`, Fat_g = `vet totaal (g)`,
           SatFat_g = `vetzuren verzadigd (g)`, UnsatFat_g = `vetzuren onverzadigd (g)`,
           TransFat_g = `trans vetzuren (g)`, LinolicAcid_g = `linolzuur (g)`,
           Cholesterol_mg = `cholesterol (mg)`, Alcohol_g = `alcohol (g)`,
           Fibers_g = `voedingsvezels (g)`, Water_ml = `water (ml)`, MonoDiSacch_g = `tot mono disach (g)`,
           PolySacch_g = `polysacchariden (g)`, PolyunsatFat_g = `vetzuren meerv onverzadigd (g)`) %>% 
    mutate(
        ID = str_c("BARIA_", ID),
        DateTime = ymd_hms(Date),
        Time = format(DateTime, "%H:%M:%S"),
        Date = as_date(DateTime),
        Year = format(Date, "%Y")
    )

diabaryear <- diabar %>%
    group_by(ID, Year, Date) %>%
    summarise(across(TotalCal_kcal:PolyunsatFat_g, ~ sum(.x, na.rm = TRUE), .names = "daily_sum_{col}"), .groups = "drop") %>%
    group_by(ID, Year) %>%
    mutate(across(starts_with("daily_sum"), ~ slider::slide_dbl(.x, mean, .before = 2, .complete = FALSE), 
                  .names = "rolling_avg_{col}")) |> 
    mutate(n_days = n())
head(diabaryear)

diabarbase <- diabaryear %>% group_by(ID) %>% 
    filter(Year == min(Year)) %>%
    slice_tail(n = 1) %>% 
    ungroup(.) %>% 
    select(-starts_with("daily_sum_"), n_days) %>% 
    rename_with(~gsub("^rolling_avg_daily_sum_", "", .), starts_with("rolling_avg_")) |> 
    mutate(diary_tool = "Bonstat")

diettot <- rbind(bariabase, diabarbase) |> full_join(bariabase_eetmeter)
bariatot <- left_join(diettot, clin1) %>% 
    filter(TotalCal_kcal < 5000 | TotalCal_kcal > 500) %>% # one subject has >8000 kcal
    mutate(
        baselineYear = format(dmy(v1_date), "%Y"),
        diffDate = time_length(interval(start = Date, end = dmy(v1_date)), unit = "days"),
        # TotalCalBin = case_when(TotalCal < median(TotalCal_kcal, na.rm = TRUE) ~ 
        #                             str_c("<", round(median(TotalCal_kcal, na.rm = TRUE), 0), " kcal"),
        #                         TotalCal_kcal >= median(TotalCal_kcal, na.rm = TRUE) ~ 
        #                             str_c(">=", round(median(TotalCal, na.rm = TRUE), 0), " kcal")),
        # ProteinBin = case_when(Protein_g < median(Protein_g, na.rm = TRUE) ~ 
        #                            str_c("<", round(median(Protein, na.rm = TRUE), 0), " g"),
        #                        Protein_g >= median(Protein_g, na.rm = TRUE) ~ 
        #                            str_c(">=", round(median(Protein_g, na.rm = TRUE), 0), " g")),
        # FibersBin = case_when(Fibers < median(Fibers, na.rm = TRUE) ~ 
        #                           str_c("<", round(median(Fibers, na.rm = TRUE),0), " g"),
        #                       Fibers >= median(Fibers, na.rm = TRUE) ~ 
        #                           str_c(">=", round(median(Fibers, na.rm = TRUE), 0), " g")),
        # CarbsBin = case_when(Carbs < median(Carbs, na.rm = TRUE) ~ 
        #                          str_c("<", round(median(Carbs, na.rm = TRUE), 0), " g"),
        #                      Carbs >= median(Carbs, na.rm = TRUE) ~ 
        #                          str_c(">=", round(median(Carbs, na.rm = TRUE),0), " g")),
        # FatBin = case_when(Carbs < median(Fat, na.rm = TRUE) ~ 
        #                        str_c("<", round(median(Fat, na.rm = TRUE), 0), " g"),
        #                    Carbs >= median(Fat, na.rm = TRUE) ~ 
        #                        str_c(">=", round(median(Fat, na.rm = TRUE), 0), " g")),
        # SatFatBin = case_when(SatFat < median(SatFat, na.rm = TRUE) ~ 
        #                           str_c("<", round(median(SatFat, na.rm = TRUE),0), " g"),
        #                       SatFat >= median(SatFat, na.rm = TRUE) ~ 
        #                           str_c(">=", round(median(SatFat, na.rm = TRUE),0), " g")),
        # UnsatFatBin = case_when(UnsatFat < median(UnsatFat, na.rm = TRUE) ~ 
        #                             str_c("<", round(median(UnsatFat, na.rm = TRUE),0), " g"),
        #                         UnsatFat >= median(UnsatFat, na.rm = TRUE) ~ 
        #                             str_c(">=", round(median(UnsatFat, na.rm = TRUE),0), " g")),
        # TransFatBin = case_when(TransFat < median(TransFat, na.rm = TRUE) ~ 
        #                             str_c("<", round(median(TransFat, na.rm = TRUE),0), " g"),
        #                         TransFat >= median(TransFat, na.rm = TRUE) ~ 
        #                             str_c(">=", round(median(TransFat, na.rm = TRUE),0), " g")),
        # LinolicAcidBin = case_when(LinolicAcid < median(LinolicAcid, na.rm = TRUE) ~ 
        #                                str_c("<", round(median(LinolicAcid, na.rm = TRUE),0), " g"),
        #                            LinolicAcid >= median(LinolicAcid, na.rm = TRUE) ~ 
        #                                str_c(">=", round(median(LinolicAcid, na.rm = TRUE),0), " g")),
        # CholesterolBin = case_when(Cholesterol < median(Cholesterol, na.rm = TRUE) ~ 
        #                                str_c("<", round(median(Cholesterol, na.rm = TRUE),0), " mg"),
        #                            Cholesterol >= median(Cholesterol, na.rm = TRUE) ~ 
        #                                str_c(">=", round(median(Cholesterol, na.rm = TRUE),0), " mg")),
        # AlcoholBin = case_when(Alcohol == 0 ~ "0 g", Alcohol > 0 ~ ">0 g"),
        # MonoDiSacchBin = case_when(MonoDiSacch < median(MonoDiSacch, na.rm = TRUE) ~ 
        #                                str_c("<", round(median(MonoDiSacch, na.rm = TRUE),0), " g"),
        #                            MonoDiSacch >= median(MonoDiSacch, na.rm = TRUE) ~ 
        #                                str_c(">=", round(median(MonoDiSacch, na.rm = TRUE),0), " g")),
        # PolySacchBin = case_when(PolySacch < median(PolySacch, na.rm = TRUE) ~ 
        #                              str_c("<", round(median(PolySacch, na.rm = TRUE),0), " g"),
        #                          PolySacch >= median(PolySacch, na.rm = TRUE) ~ 
        #                              str_c(">=", round(median(PolySacch, na.rm = TRUE),0), " g")),
        # PolyunsatFatBin = case_when(PolyunsatFat < median(PolyunsatFat, na.rm = TRUE) ~ 
        #                                 str_c("<", round(median(PolyunsatFat, na.rm = TRUE),0), " g"),
        #                             PolyunsatFat >= median(PolyunsatFat, na.rm = TRUE) ~ 
        #                                 str_c(">=", round(median(PolyunsatFat, na.rm = TRUE),0), " g")),
        # across(contains("Bin"), as.factor)
    ) # %>% filter(diffDate >= 0 & diffDate < 365)
#bariadiet <- bariatot %>% filter(diffDate >= 0 & diffDate < 365)
any(duplicated(bariatot$ID)) # no duplicates
saveRDS(bariatot, "data/bariatot.RDS")

Amelia::missmap(diettot |> filter(diary_tool == "Eetmeter"))

bariatot_dii <- bariatot |> 
    filter(TotalCal_kcal < 5000) |> # one subject has >8000 kcal
    filter(TotalCal_kcal > 500) |> 
    mutate(
        TotalCal_DII = TotalCal_kcal * 0.180,
        Protein_DII = Protein_g * 0.021,
        Carbs_DII = Carbs_g * 0.097,
        Fat_DII = Fat_g * 0.298,
        SatFat_DII = SatFat_g * 0.373,
        TransFat_DII = TransFat_g * 0.229,    
        Cholesterol_DII = Cholesterol_mg * 0.110,
        Alcohol_DII = Alcohol_g * -0.278,
        Fiber_DII = Fibers_g * -0.663,
        PUFA_DII = PolyunsatFat_g * -0.337,
        Magnesium_DII = Magnesium_mg * -0.484,
        Iron_DII = Iron_mg * 0.032,
        Selenium_DII = Selenium_ug * -0.191,
        Zinc_DII = Zinc_mg * -0.313,
        VitaminA_DII = VitaminA_ug * -0.401,
        VitaminD_DII = VitaminD_ug * -0.446,
        VitaminE_DII = VitaminE_mg * -0.419,
        VitaminB1_DII = VitaminB1_mg * -0.098,
        VitaminB2_DII = VitaminB2_mg * -0.068,
        VitaminB6_DII = VitaminB6_mg * -0.365,
        FolicAcid_DII = FolicAcid_ug * -0.190,
        VitaminB12_DII = VitaminB12_ug * 0.106,
        Niacin_DII = Niacin_mg * -0.246,
        VitaminC_DII = VitaminC_mg * -0.424        
    ) |> 
    mutate(
        DII_Eetmeter = TotalCal_DII + Protein_DII + Carbs_DII + Fat_DII + SatFat_DII + Alcohol_DII + Fiber_DII + Magnesium_DII +
            Iron_DII + Selenium_DII + Zinc_DII + VitaminA_DII + VitaminD_DII + VitaminE_DII + VitaminB1_DII + VitaminB2_DII +
            VitaminB6_DII + FolicAcid_DII + VitaminB12_DII + Niacin_DII + VitaminC_DII,
        DII_Bonstat = TotalCal_DII + Protein_DII + Carbs_DII + Fat_DII + SatFat_DII + TransFat_DII + Cholesterol_DII + Alcohol_DII +
            Fiber_DII + PUFA_DII,
        DII_BonstatEetmeter = TotalCal_DII + Protein_DII + Carbs_DII + Fat_DII + SatFat_DII + Alcohol_DII + Fiber_DII 
    )
gghistogram(bariatot_dii$DII_Bonstat)
gghistogram(bariatot_dii$DII_Eetmeter)
gghistogram(bariatot_dii$DII_BonstatEetmeter)

## Gene list
gene_list <- read.csv("data/ICP_list.csv", sep = ';')

gene_list$ICP_symbol <- trimws(gene_list$ICP_symbol)
head(gene_list)
saveRDS(gene_list, "data/gene_list.RDS")

# Serum samples
serum <- readxl::read_xlsx("data/20250320_BARIA_300_Paris_Serum_Citraat_Nieuwdorp_MW.xlsx")
head(serum)
dim(serum)
serum$ID <- str_c("BARIA_", serum$ID)
serum %>% filter(Citraat < 1)
saveRDS(serum$ID, "data/serum_ids.RDS")
serum$ID
bariatot$ID
summary(serum$ID %in% bariatot$ID) 
missingids <- serum$ID[!serum$ID %in% bariatot$ID]
clin1 %>% filter(ID %in% missingids) %>%  pull(ID) %>% length() 
# 143 are all in clinical data
idsrnaseq <- readRDS("data/idsrnaseq.RDS")
summary(serum$ID %in% idsrnaseq) 
# 3 do not overlap with RNAseq data (overlap w missings above)

# concl: of 282 serum samples, 
# 53 do not have dietary data within 1 year, 
# 3 do not have RNAseq data

# filter diet data
diettot |> filter(TotalCal_kcal < 500) |> select(ID, TotalCal_kcal, diary_tool, n_days)
diettot2 <- diettot |> 
    filter(TotalCal_kcal > 500) |> 
    filter(ID %in% serum$ID) |> 
    mutate(diary_tool = as.factor(diary_tool))
dim(diettot2)
Amelia::missmap(diettot2)
diettot2 <- diettot2 |> select(ID, Visit, diary_tool, n_days, TotalCal_kcal, Carbs_g, Protein_g, Fat_g,
                                SatFat_g, Fibers_g, Alcohol_g, Water_ml)
gghistogram(diettot2$TotalCal_kcal)
gghistogram(diettot2$Carbs_g)
gghistogram(diettot2$Protein_g)
gghistogram(diettot2$Fat_g)
gghistogram(diettot2$SatFat_g)
gghistogram(diettot2$Fibers_g)
gghistogram(diettot2$Alcohol_g)
gghistogram(diettot2$Water_ml)
dim(diettot2)
write.csv(diettot2, "data/251002_BARIA_macronutrients.csv")
saveRDS(diettot2, "data/251002_BARIA_macronutrients.RDS")
