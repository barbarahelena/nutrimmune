## Diet descriptives
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(corrplot)
library(Hmisc)

## Functions
draw_corrplot <- function(df){
    rescor <-   rcorr(as.matrix(df), type="spearman")
    corplot <- corrplot(rescor$r,  type = "upper", tl.col = "black", tl.cex = 0.6, cl.cex = 0.6, number.cex = 0.5,
                        order = 'hclust', hclust.method="ward.D", tl.srt = 45, insig="blank", sig.level = 0.05,
                        p.mat=rescor$P, method="color", mar=c(0,0,1,0), addCoef.col = "black", diag = F,
                        addgrid.col = "grey")
}

## Open data
baria <- readRDS("data/bariatot.RDS")

## Correlation plot
baria <- baria %>% select(TotalCal, Protein, Carbs, Fat, SatFat, UnsatFat, TransFat,
                          LinolicAcid, Cholesterol, Alcohol, Fibers,
                          Water, MonoDiSacch, PolySacch, PolyunsatFat)
pdf("results/descriptives/correlationplot_diet.pdf", width = 6, height = 6)
draw_corrplot(baria)
dev.off()

# Table
table2 <- baria %>% 
    CreateTableOne(data=., test = FALSE) %>% 
    print(nonnormal=c("Alcohol"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)
write.csv(table1, "results/descriptives/table2.csv")

## Diets
barialong <- baria %>% select(Carbs, Fibers, Protein, Fat) %>% 
    pivot_longer(., cols = 1:4, names_to = "var", values_to = "value")
ggplot(data = barialong, aes(x = fct_reorder(var, value, .desc = TRUE), 
                             y = value, fill = var)) +
    geom_violin(bw = 15) +
    geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
    scale_fill_bmj(guide = "none") +
    theme_Publication() +
    labs(x = "", y = "gram", title = "Macronutrients")
ggsave("results/descriptives/macronutrients.pdf", width = 4, height = 5)

barialong <- baria %>% select(SatFat, UnsatFat, TransFat, PolyunsatFat) %>% 
    pivot_longer(., cols = 1:4, names_to = "var", values_to = "value")
ggplot(data = barialong, aes(x = fct_reorder(var, value, .desc = TRUE), 
                             y = value, fill = var)) +
    geom_violin(bw = 5) +
    geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
    scale_fill_manual(guide = "none", values = pal_bmj()(9)[c(5:9)]) +
    theme_Publication() +
    labs(x = "", y = "gram", title = "Fat")
ggsave("results/descriptives/fat.pdf", width = 4, height = 5)

barialong <- baria %>% select(MonoDiSacch, PolySacch) %>% 
    pivot_longer(., cols = 1:2, names_to = "var", values_to = "value")
ggplot(data = barialong, aes(x = fct_reorder(var, value, .desc = TRUE), 
                             y = value, fill = var)) +
    geom_violin(bw = 5) +
    geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
    scale_fill_manual(guide = "none", values = pal_bmj()(9)[c(1,9)]) +
    theme_Publication() +
    labs(x = "", y = "gram", title = "Saccharides")
ggsave("results/descriptives/saccharides.pdf", width = 3, height = 5)
