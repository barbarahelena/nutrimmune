## Violin plots
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)

# Open data
baria <- readRDS("data/bariatot.RDS")
crp <- bariatot %>% filter(v0_crp < 100 & v4_crp < 100 & v5_crp < 60) # remove obvious outliers

#### CRP ####
plist <- c()
for(a in names(crp)[152:165]) {
    crp$var <- crp[[a]]
    pl <- ggplot(data = crp, aes(x = var, y = crp_v1v4)) +
        geom_violin(fill = ggsci::pal_bmj()(1)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        stat_compare_means() +
        theme_Publication() +
        labs(title = paste0(a), y = "crp change at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/violins/crpv1v4_violin_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(crp)[152:165]) {
    crp$var <- crp[[a]]
    pl <- ggplot(data = crp, aes(x = var, y = crp_v1v5)) +
        geom_violin(fill = ggsci::pal_bmj()(1)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        stat_compare_means() +
        theme_Publication() +
        labs(title = paste0(a), y = "crp change at 2 years",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/violins/crpv1v5_violin_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(crp)[152:165]) {
    crp$var <- crp[[a]]
    pl <- ggplot(data = crp, aes(x = var, y = v4_crp)) +
        geom_violin(fill = ggsci::pal_bmj()(1)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        stat_compare_means() +
        theme_Publication() +
        labs(title = paste0(a), y = "crp at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/violins/crpv4_diet.pdf", width = 10, height = 20)

#### BMI ####
plist <- c()
for(a in names(baria)[152:165]) {
    baria$var <- baria[[a]]
    pl <- ggplot(data = baria, aes(x = var, y = bmi_v1v4)) +
        geom_violin(fill = ggsci::pal_bmj()(1)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        stat_compare_means() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi change at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/violins/bmiv1v4_violin_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(baria)[152:165]) {
    baria$var <- baria[[a]]
    pl <- ggplot(data = baria, aes(x = var, y = bmi_v1v5)) +
        geom_violin(fill = ggsci::pal_bmj()(1)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        stat_compare_means() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi change at 2 years",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/violins/bmiv1v5_violin_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(baria)[152:165]) {
    baria$var <- baria[[a]]
    pl <- ggplot(data = baria, aes(x = var, y = bmiperc_v1v4)) +
        geom_violin(fill = ggsci::pal_bmj()(1)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        stat_compare_means() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi % change at 2 years",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/violins/bmipercv1v4_violin_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(baria)[152:165]) {
    baria$var <- baria[[a]]
    pl <- ggplot(data = baria, aes(x = var, y = bmiperc_v1v5)) +
        geom_violin(fill = ggsci::pal_bmj()(1)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        stat_compare_means() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi % change at 2 years",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/violins/bmipercv1v5_violin_diet.pdf", width = 10, height = 20)


#### HOMA IR ####
homair <- baria %>% filter(is.na(homaIR_v4) | homaIR_v4 < 25)
plist <- c()
for(a in names(homair)[152:165]) {
    homair$var <- homair[[a]]
    pl <- ggplot(data = homair, aes(x = var, y = homaIR_v1v4)) +
        geom_violin(fill = ggsci::pal_bmj()(1)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        stat_compare_means() +
        theme_Publication() +
        labs(title = paste0(a), y = "HOMA-IR change at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/violins/homairv1v4_violin_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(homair)[152:165]) {
    homair$var <- homair[[a]]
    pl <- ggplot(data = homair, aes(x = var, y = homaIR_v1v5)) +
        geom_violin(fill = ggsci::pal_bmj()(1)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        stat_compare_means() +
        theme_Publication() +
        labs(title = paste0(a), y = "HOMA-IR change at 2 years",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/violins/homairv1v5_violin_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(homair)[152:165]) {
    homair$var <- homair[[a]]
    pl <- ggplot(data = homair, aes(x = var, y = homaIR_v4)) +
        geom_violin(fill = ggsci::pal_bmj()(1)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        stat_compare_means() +
        theme_Publication() +
        labs(title = paste0(a), y = "HOMA-IR at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/violins/homairv4_violin_diet.pdf", width = 10, height = 20)
