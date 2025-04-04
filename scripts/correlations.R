# Correlations with .. gene
# Barbara Verhaar

# Libraries
library(tidyverse)
library(ggpubr)

# Theme
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
meta <- readRDS("data/baria_metadata.RDS")
df <- readRDS("data/RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS")
names(meta)

# Get genes of interest
ensembls <- rownames(df)
which(str_detect(ensembls, "ENSG00000188389")) # PD1L #ENSG00000276977
which(str_detect(ensembls, "ENSG00000120217")) # CD274
	
pd1_liver <- as.data.frame(t(as.matrix(df[which(str_detect(ensembls, "ENSG00000188389")),str_detect(colnames(df), "Liver.V1")])))
pd1_vfat <- as.data.frame(t(as.matrix(df[which(str_detect(ensembls, "ENSG00000188389")),str_detect(colnames(df), "vFat")])))
pd1_sfat <- as.data.frame(t(as.matrix(df[which(str_detect(ensembls, "ENSG00000188389")),str_detect(colnames(df), "subFat")])))

cd274_liver <- as.data.frame(t(as.matrix(df[which(str_detect(ensembls, "ENSG00000120217")),str_detect(colnames(df), "Liver.V1")])))
cd274_vfat <- as.data.frame(t(as.matrix(df[which(str_detect(ensembls, "ENSG00000120217")),str_detect(colnames(df), "vFat.V1")])))
cd274_sfat <- as.data.frame(t(as.matrix(df[which(str_detect(ensembls, "ENSG00000120217")),str_detect(colnames(df), "subFat.V1")])))

pd1_liver <- pd1_liver %>% mutate(ID = rownames(.)) %>% rename("PD1_liver" = "ENSG00000188389.11")
pd1_vfat <- pd1_vfat %>% mutate(ID = rownames(.)) %>% rename("PD1_viscfat" = "ENSG00000188389.11")
pd1_sfat <- pd1_sfat %>% mutate(ID = rownames(.)) %>% rename("PD1_subcfat" = "ENSG00000188389.11")
pd1 <- full_join(pd1_liver, pd1_sfat, by = "ID") %>% 
            full_join(., pd1_vfat, by = "ID") %>% 
            mutate(ID = str_extract(ID, "([0-9]*)"),
                    ID = str_c("BARIA_", ID),
                    across(contains("PD1"), as.numeric)) %>% 
            left_join(., meta, by = "ID")

cd274_liver <- cd274_liver %>% mutate(ID = rownames(.)) %>% rename("CD274_liver" = "ENSG00000120217.14")
cd274_vfat <- cd274_vfat %>% mutate(ID = rownames(.)) %>% rename("CD274_viscfat" = "ENSG00000120217.14")
cd274_sfat <- cd274_sfat %>% mutate(ID = rownames(.)) %>% rename("CD274_subcfat" = "ENSG00000120217.14")
cd274 <- full_join(cd274_liver, cd274_sfat, by = "ID") %>% 
            full_join(., cd274_vfat, by = "ID") %>% 
            mutate(ID = str_extract(ID, "([0-9]*)"),
                    ID = str_c("BARIA_", ID),
                    across(contains("CD274"), as.numeric)) %>% 
                    left_join(., meta, by = "ID")

# Counts
gghistogram(pd1$PD1_liver, fill = "firebrick") + 
    theme_Publication() + 
    ggtitle("Distribution counts PD1 liver") +
    labs(x = "gene counts", y = "subjects")
ggsave("results/correlations/pd1/distribution_liver.pdf", width = 5, height = 5)
gghistogram(pd1$PD1_viscfat, fill = "gold") + 
    theme_Publication() +
    ggtitle("Distribution counts PD1 visceral fat") +
    labs(x = "gene counts", y = "subjects")
ggsave("results/correlations/pd1/distribution_viscfat.pdf", width = 5, height = 5)
gghistogram(pd1$PD1_subcfat, fill = "royalblue") + 
    theme_Publication() +
    ggtitle("Distribution counts PD1 subcutaneous fat") +
    labs(x = "gene counts", y = "subjects")
ggsave("results/correlations/pd1/distribution_subcfat.pdf", width = 5, height = 5)

gghistogram(cd274$CD274_liver, fill = "firebrick") + 
    theme_Publication() + 
    ggtitle("Distribution counts CD274/PDL1 liver") +
    labs(x = "gene counts", y = "subjects")
ggsave("results/correlations/cd274/distribution_liver.pdf", width = 5, height = 5)
gghistogram(cd274$CD274_viscfat, fill = "gold")  + 
    theme_Publication() + 
    ggtitle("Distribution counts CD274/PDL1 liver") +
    labs(x = "gene counts", y = "subjects")
ggsave("results/correlations/cd274/distribution_viscfat.pdf", width = 5, height = 5)
gghistogram(cd274$CD274_subcfat, fill = "royalblue")  + 
    theme_Publication() + 
    ggtitle("Distribution counts CD274/PDL1 liver") +
    labs(x = "gene counts", y = "subjects")
ggsave("results/correlations/cd274/distribution_subcfat.pdf", width = 5, height = 5)

# Correlations CRP
ggplot(data = pd1, aes(x = PD1_liver, y = log10(v0_crp+1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "log10(CRP (mg/l))") +
    theme_Publication()
ggsave("results/correlations/pd1/crp/crp_pd1_liver.pdf", width = 5, height = 5)

ggplot(data = pd1, aes(x = PD1_subcfat, y = log10(v0_crp+1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "log10(CRP (mg/l))") +
    theme_Publication()
ggsave("results/correlations/pd1/crp/crp_pd1_subcfat.pdf", width = 5, height = 5)

ggplot(data = pd1, aes(x = PD1_viscfat, y = log10(v0_crp+1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "log10(CRP (mg/l))") +
    theme_Publication()
ggsave("results/correlations/pd1/crp/crp_pd1_viscfat.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_liver, y = log10(v0_crp+1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "log10(CRP (mg/l))") +
    theme_Publication()
ggsave("results/correlations/cd274/crp/crp_cd274_liver.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_subcfat, y = log10(v0_crp+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "log10(CRP (mg/l))") +
    theme_Publication()
ggsave("results/correlations/cd274/crp/crp_cd274_subcfat.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_viscfat, y = log10(v0_crp+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "log10(CRP (mg/l))") +
    theme_Publication()
ggsave("results/correlations/cd274/crp/crp_cd274_viscfat.pdf", width = 5, height = 5)

# Correlations leukocytes
ggplot(data = pd1, aes(x = PD1_liver, y = log10(v0_leuko+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "Leukocytes x10E9/L") +
    theme_Publication()
ggsave("results/correlations/pd1/leuko/leuko_pd1_liver.pdf", width = 5, height = 5)

ggplot(data = pd1, aes(x = PD1_subcfat, y = log10(v0_leuko+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "Leukocytes x10E9/L") +
    theme_Publication()
ggsave("results/correlations/pd1/leuko/leuko_pd1_subcfat.pdf", width = 5, height = 5)

ggplot(data = pd1, aes(x = PD1_viscfat, y = log10(v0_leuko+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "Leukocytes x10E9/L") +
    theme_Publication()
ggsave("results/correlations/pd1/leuko/leuko_pd1_viscfat.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_liver, y = log10(v0_leuko+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "Leukocytes x10E9/L") +
    theme_Publication()
ggsave("results/correlations/cd274/leuko/leuko_cd274_liver.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_subcfat, y = log10(v0_leuko+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "Leukocytes x10E9/L") +
    theme_Publication()
ggsave("results/correlations/cd274/leuko/leuko_cd274_subcfat.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_viscfat, y = log10(v0_leuko+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "Leukocytes x10E9/L") +
    theme_Publication()
ggsave("results/correlations/cd274/leuko/leuko_cd274_viscfat.pdf", width = 5, height = 5)

# Correlations lipids
ggplot(data = pd1, aes(x = PD1_liver, y = v0_ldl)) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "LDL (mmol/L)") +
    theme_Publication()
ggsave("results/correlations/pd1/lipids/ldl_pd1_liver.pdf", width = 5, height = 5)

ggplot(data = pd1, aes(x = PD1_subcfat, y = v0_ldl)) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "LDL (mmol/L)") +
    theme_Publication()
ggsave("results/correlations/pd1/lipids/ldl_pd1_subfat.pdf", width = 5, height = 5)

ggplot(data = pd1, aes(x = PD1_viscfat, y = v0_ldl)) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "LDL (mmol/L)") +
    theme_Publication()
ggsave("results/correlations/pd1/lipids/ldl_pd1_viscfat.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_liver, y = v0_ldl)) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "LDL (mmol/L)") +
    theme_Publication()
ggsave("results/correlations/cd274/lipids/ldl_cd274_liver.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_subcfat, y = v0_ldl)) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "LDL (mmol/L)") +
    theme_Publication()
ggsave("results/correlations/cd274/lipids/ldl_cd274_subfat.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_viscfat, y = v0_ldl)) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "LDL (mmol/L)") +
    theme_Publication()
ggsave("results/correlations/cd274/lipids/ldl_cd274_viscfat.pdf", width = 5, height = 5)

# Correlations HOMA-IR
ggplot(data = pd1, aes(x = PD1_liver, y = log10(homaIR_v1+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "HOMA-IR") +
    theme_Publication()
ggsave("results/correlations/pd1/homair/homair_pd1_liver.pdf", width = 5, height = 5)

ggplot(data = pd1, aes(x = PD1_subcfat, y = log10(homaIR_v1+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "HOMA-IR") +
    theme_Publication()
ggsave("results/correlations/pd1/homair/homair_pd1_subfat.pdf", width = 5, height = 5)

ggplot(data = pd1, aes(x = PD1_viscfat, y = log10(homaIR_v1+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "PD1 (counts)", y = "HOMA-IR") +
    theme_Publication()
ggsave("results/correlations/pd1/homair/homair_pd1_viscfat.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_liver, y = log10(homaIR_v1+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "HOMA-IR") +
    theme_Publication()
ggsave("results/correlations/cd274/homair/homair_cd274_liver.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_subcfat, y = log10(homaIR_v1+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "HOMA-IR") +
    theme_Publication()
ggsave("results/correlations/cd274/homair/homair_cd274_subfat.pdf", width = 5, height = 5)

ggplot(data = cd274, aes(x = CD274_viscfat, y = log10(homaIR_v1+0.1))) +
    geom_point(alpha = 0.9, color = "royalblue") +
    geom_smooth(method = "lm", color = "firebrick4") +
    stat_cor(method = "spearman") +
    labs(x = "CD274 (counts)", y = "HOMA-IR") +
    theme_Publication()
ggsave("results/correlations/cd274/homair/homair_cd274_viscfat.pdf", width = 5, height = 5)
