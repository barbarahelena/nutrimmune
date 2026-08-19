# Script to generate scatter plots for significant correlations from heatmap analysis
# Barbara Verhaar

# Libraries
library(tidyverse)
library(ggpubr)
library(readxl)
library(circlize)
library(ggsci)
library(grid)
library(rlang)  # Added for .data pronoun

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
                axis.title = element_text(face = "bold", size = rel(0.8)),
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
                legend.key.size= unit(0.2, "cm"),
                legend.spacing = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
}

# Function to extract expression data for a gene across tissues
extract_gene_data <- function(ensembl_id, gene_symbol = NULL) {
  base_id <- str_extract(ensembl_id, "ENSG[0-9]+")
  matches <- which(str_detect(global_ensembls, base_id))
  if(length(matches) == 0) {
    warning(paste0("No matches found for gene ID: ", ensembl_id))
    return(NULL)
  }
  match_to_use <- matches[1]
  ensembl_with_version <- global_ensembls[match_to_use]

  liver_cols <- colnames(global_expr)[str_detect(colnames(global_expr), "\\.Liver\\.V1$")]
  vfat_cols <- colnames(global_expr)[str_detect(colnames(global_expr), "\\.vFat\\.V1$")]
  sfat_cols <- colnames(global_expr)[str_detect(colnames(global_expr), "\\.subFat\\.V1$")]

  process_tissue_data <- function(tissue_data, tissue_type) {
    if (nrow(tissue_data) == 0) return(NULL)
    subject_ids <- gsub("^([0-9]+)\\..+$", "\\1", rownames(tissue_data))
    formatted_ids <- paste0("BARIA_", subject_ids)
    values <- gsub(",", ".", as.character(tissue_data[,1]))
    numeric_values <- log10(as.numeric(values) + 1)
    symbol_to_use <- ifelse(is.null(gene_symbol) || gene_symbol == "", ensembl_id, gene_symbol)
    
    colname <- paste0(symbol_to_use, "_", tissue_type)
    df <- data.frame(
      ID = formatted_ids,
      value = numeric_values,
      stringsAsFactors = FALSE
    )
    colnames(df)[2] <- colname
    return(df)
  }

  liver_df <- process_tissue_data(as.data.frame(t(as.matrix(global_expr[match_to_use, liver_cols, drop=FALSE]))), "liver")
  vfat_df <- process_tissue_data(as.data.frame(t(as.matrix(global_expr[match_to_use, vfat_cols, drop=FALSE]))), "viscfat")
  sfat_df <- process_tissue_data(as.data.frame(t(as.matrix(global_expr[match_to_use, sfat_cols, drop=FALSE]))), "subcfat")

  if (is.null(liver_df) || is.null(vfat_df) || is.null(sfat_df)) return(NULL)

  gene_df <- full_join(liver_df, sfat_df, by = "ID") %>% 
    full_join(., vfat_df, by = "ID") %>% 
    left_join(., global_meta, by = "ID")
  return(gene_df)
}

# Function to create correlation plot
create_correlation_plot <- function(data, x_var, y_var, gene_symbol, tissue, clinical_var) {

  y_expr <- switch(clinical_var,
                   "TotalCalories" = "TotalCal_kcal",
                   "Carbohydrates" = "Carbs_g",
                   "Fat" = "Fat_g",
                   "SaturatedFat" = "SatFat_g",
                   "Protein" = "Protein_g",
                   "Fibers" = "Fibers_g",
                   "Age" = "v0_age",
                   "BMI" = "v0_bmi",
                   clinical_var)

  y_label <- switch(clinical_var,
                    "TotalCalories" = "Total calorie intake (kcal)",
                    "Carbohydrates" = "Carbohydrates (g)",
                    "Fat" = "Fat (g)",
                    "SaturatedFat" = "Saturated fat (g)",
                    "Protein" = "Protein (g)",
                    "Fibers" = "Fibers (g)",
                    "Age" = "Age (years)",
                    "BMI" = "BMI (kg/m²)",
                    clinical_var)

  dir_path <- file.path("results", "correlations_macronutrients", "significant_plots", gene_symbol, clinical_var)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  p <- ggplot(data = data) + 
    geom_point(aes(x = .data[[x_var]], y = .data[[y_expr]]), 
              alpha = 0.9, color = "royalblue") +
    geom_smooth(aes(x = .data[[x_var]], y = .data[[y_expr]]), 
               method = "lm", color = "firebrick4") +
    stat_cor(aes(x = .data[[x_var]], y = .data[[y_expr]]), 
            method = "spearman") +
    labs(
      x = paste0("log10(", gene_symbol, " counts + 1)"),
      y = y_label,
      title = paste0(gene_symbol, " in ", tissue, " vs ", clinical_var)
    ) +
    theme_Publication()

  ggsave(file.path(dir_path, paste0(clinical_var, "_", gene_symbol, "_", tissue, ".pdf")),
         plot = p, width = 5, height = 5)
}

# Main
main <- function() {
  cat("Loading data...\n")
  global_meta <<- readRDS("data/bariatot.RDS")
  global_expr <<- readRDS("data/RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS")
  global_ensembls <<- rownames(global_expr)
  gene_df <- readRDS("data/gene_list.RDS")
  gene_df$ensembl_id <- gene_df$Ensembl_ID
  gene_df$symbol <- gene_df$Gene

  tissues <- c("liver", "viscfat", "subcfat")
  clinical_vars <- c("TotalCalories", "Carbohydrates", "Fat", "SaturatedFat", "Protein", "Fibers", "Age", "BMI")

  all_correlations <- list()
  all_pvalues <- list()

  # Only load baseline (v0) correlations
  for (tissue in tissues) {
    corr_file <- paste0("results/correlations_macronutrients/heatmaps_v0/correlations_", tissue, ".csv")
    pval_file <- paste0("results/correlations_macronutrients/heatmaps_v0/pvalues_", tissue, ".csv")

    if (file.exists(corr_file) && file.exists(pval_file)) {
      all_correlations[[tissue]] <- read.csv(corr_file, row.names = 1)
      all_pvalues[[tissue]] <- read.csv(pval_file, row.names = 1)
      cat(sprintf("Loaded correlations for %s\n", tissue))
    }
  }

  significant_correlations <- data.frame()
  for (tissue in names(all_correlations)) {
    corr_matrix <- all_correlations[[tissue]]
    pval_matrix <- all_pvalues[[tissue]]
    for (i in 1:nrow(corr_matrix)) {
      for (j in 1:ncol(corr_matrix)) {
        if (!is.na(pval_matrix[i, j]) && pval_matrix[i, j] < 0.05) {
          gene_symbol <- rownames(corr_matrix)[i]
          clinical_var <- colnames(corr_matrix)[j]
          ensembl_id <- gene_df$ensembl_id[gene_df$symbol == gene_symbol][1]
          significant_correlations <- rbind(significant_correlations,
                                             data.frame(gene = gene_symbol,
                                                        ensembl_id = ensembl_id,
                                                        tissue = tissue,
                                                        clinical_var = clinical_var,
                                                        correlation = corr_matrix[i, j],
                                                        pvalue = pval_matrix[i, j]))
        }
      }
    }
  }

  write.csv(significant_correlations,
            "results/correlations_macronutrients/significant_correlations_baseline.csv",
            row.names = FALSE)
  cat(sprintf("Found %d significant baseline correlations\n", nrow(significant_correlations)))

  plots_created <- 0
  for (i in 1:nrow(significant_correlations)) {
    gene_symbol <- significant_correlations$gene[i]
    ensembl_id <- significant_correlations$ensembl_id[i]
    tissue <- significant_correlations$tissue[i]
    clinical_var <- significant_correlations$clinical_var[i]

    gene_data <- extract_gene_data(ensembl_id, gene_symbol)
    if (!is.null(gene_data)) {
      y_var <- switch(clinical_var,
                      "TotalCalories" = "TotalCal_kcal",
                      "Carbohydrates" = "Carbs_g",
                      "Fat" = "Fat_g",
                      "SaturatedFat" = "SatFat_g",
                      "Protein" = "Protein_g",
                      "Fibers" = "Fibers_g",
                      "Age" = "v0_age",
                      "BMI" = "v0_bmi",
                      clinical_var)
      x_var <- paste0(gene_symbol, "_", tissue)
      if (x_var %in% colnames(gene_data) && y_var %in% colnames(gene_data)) {
        create_correlation_plot(gene_data, x_var, y_var, gene_symbol, tissue, clinical_var)
        plots_created <- plots_created + 1
      }
    }
  }
  cat(sprintf("Created %d baseline correlation plots\n", plots_created))
}

main()