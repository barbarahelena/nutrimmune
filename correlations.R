## Correlation plots
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)

# Open data
baria <- readRDS("data/bariatot.RDS")
crp <- bariatot %>% filter(v0_crp < 100 & v4_crp < 100)

#### CRP ####
plist <- c()
for(a in names(crp)[4:18]) {
    crp$var <- crp[[a]]
    print(gghistogram(crp$var, title = a, fill = ggsci::pal_bmj()(1)))
    pl <- ggplot(data = crp, aes(x = log10(var), y = crp_v1v4)) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "crp change at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/correlations/crpv1v4_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(crp)[4:18]) {
    crp$var <- crp[[a]]
    # print(gghistogram(crp$var, title = a, fill = ggsci::pal_bmj()(1)))
    pl <- ggplot(data = crp, aes(x = log10(var), y = crp_v1v5)) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "crp change at 2 years",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/correlations/crpv1v5_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(crp)[4:18]) {
    crp$var <- crp[[a]]
    pl <- ggplot(data = crp, aes(x = log10(var), y = v4_crp)) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "crp at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/correlations/crpv4_diet.pdf", width = 10, height = 20)

#### BMI ####
plist <- c()
for(a in names(baria)[4:18]) {
    baria$var <- baria[[a]]
    pl <- ggplot(data = baria, aes(x = log10(var), y = bmi_v1v4)) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi change at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/correlations/bmiv1v4_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(baria)[4:18]) {
    baria$var <- baria[[a]]
    pl <- ggplot(data = baria, aes(x = log10(var), y = bmi_v1v5)) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi change at 2 years",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/correlations/bmiv1v5_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(baria)[4:18]) {
    baria$var <- baria[[a]]
    pl <- ggplot(data = baria, aes(x = log10(var), y = bmiperc_v1v4)) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi % change at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/correlations/bmipercv1v4_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(baria)[4:18]) {
    baria$var <- baria[[a]]
    pl <- ggplot(data = baria, aes(x = log10(var), y = bmiperc_v1v5)) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "bmi % change at 2 years",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/correlations/bmipercv1v5_diet.pdf", width = 10, height = 20)


#### HOMA IR ####
plist <- c()
for(a in names(baria)[4:18]) {
    baria$var <- baria[[a]]
    pl <- ggplot(data = baria, aes(x = log10(var), y = log10(homaIR_v4+0.01))) +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "HOMA-IR at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/correlations/homairv4_diet.pdf", width = 10, height = 20)

# plist <- c()
# for(a in names(baria)[4:18]) {
#     baria$var <- baria[[a]]
#     pl <- ggplot(data = baria, aes(x = log10(var), y = log10(homaIR_v5+0.01))) +
#         geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
#         geom_smooth(method = "lm") +
#         stat_cor() +
#         theme_Publication() +
#         labs(title = paste0(a), y = "HOMA-IR at 2 years",
#              x = paste0(a))
#     plist[[a]] <- pl
# }
# ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
# ggsave("results/correlations/homairv5_diet.pdf", width = 10, height = 20)

plist <- c()
for(a in names(baria)[4:18]) {
    baria$var <- baria[[a]]
    pl <- ggplot(data = baria, aes(x = log10(var), y = homaIR_v1v4*-1)) +
        scale_y_log10() +
        geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
        geom_smooth(method = "lm") +
        stat_cor() +
        theme_Publication() +
        labs(title = paste0(a), y = "Change HOMA-IR at 1 year",
             x = paste0(a))
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
ggsave("results/correlations/homairv1v4_diet.pdf", width = 10, height = 20)


# plist <- c()
# for(a in names(baria)[4:18]) {
#     baria$var <- baria[[a]]
#     pl <- ggplot(data = baria, aes(x = log10(var), y = homaIR_v1v5*-1)) +
#         scale_y_log10() +
#         geom_point(alpha = 0.7, color = ggsci::pal_bmj()(1)) +
#         geom_smooth(method = "lm") +
#         stat_cor() +
#         theme_Publication() +
#         labs(title = paste0(a), y = "Change HOMA-IR at 1 year",
#              x = paste0(a))
#     plist[[a]] <- pl
# }
# ggarrange(plotlist = plist, ncol = 3, nrow = 5, labels = LETTERS[1:15])
# ggsave("results/correlations/homairv1v5_diet.pdf", width = 10, height = 20)

