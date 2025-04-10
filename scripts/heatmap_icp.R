# Horizontal Heatmap of Median Gene Expression in Different Tissues (ICP Genes)
# Barbara Verhaar

# Libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggsci)

# Load data
cat("Loading data...\n")
global_expr <- readRDS("data/RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS")
gene_list <- readRDS("data/gene_list.RDS") %>% slice(1:45)

# Clean gene symbols
gene_list <- gene_list %>%
  mutate(
    symbol = gsub("\\s*\\([^)]*\\)", "", ICP_symbol),  # Remove parenthetical text
    symbol = gsub("[,;].*", "", symbol)               # Remove text after commas/semicolons
  )

# Extract genes of interest
genes_of_interest <- gene_list$symbol
ensembl_ids <- gene_list$ENSEMBL

# Remove version numbers from Ensembl IDs in global_expr
cat("Removing version numbers from Ensembl IDs in global_expr...\n")
rownames(global_expr) <- gsub("\\..*", "", rownames(global_expr))

# Handle European decimal format (replace commas with periods)
cat("Converting numbers with commas to numeric format...\n")
global_expr <- global_expr %>%
  mutate(across(everything(), ~ as.numeric(gsub(",", ".", as.character(.)))))

# Check matching rows
cat("Number of matching rows in expression data after cleaning Ensembl IDs:", sum(rownames(global_expr) %in% ensembl_ids), "\n")
if (sum(rownames(global_expr) %in% ensembl_ids) == 0) {
  stop("No matching rows found even after cleaning Ensembl IDs. Check data consistency.")
}

# Subset expression data for the genes of interest
cat("Subsetting expression data...\n")
expression_data <- global_expr[rownames(global_expr) %in% ensembl_ids, ] %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_id")

# Reshape data for easier manipulation
expression_long <- expression_data %>%
  pivot_longer(
    cols = -ensembl_id,
    names_to = "sample",
    values_to = "expression"
  ) %>%
  mutate(
    tissue = case_when(
      str_detect(sample, "\\.Liver\\.V1$") ~ "Liver",
      str_detect(sample, "\\.vFat\\.V1$") ~ "VisceralFat",
      str_detect(sample, "\\.subFat\\.V1$") ~ "SubcutaneousFat",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(tissue))  # Keep only relevant tissues

# Aggregate median expression per gene per tissue
cat("Aggregating expression data by tissue (median)...\n")
median_expression <- expression_long %>%
  group_by(ensembl_id, tissue) %>%
  summarize(median_expression = median(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = tissue,
    values_from = median_expression
  )

# Add gene symbols to the aggregated data
heatmap_data <- median_expression %>%
  left_join(gene_list, by = c("ensembl_id" = "ENSEMBL")) %>%
  select(symbol, Liver, VisceralFat, SubcutaneousFat) %>%
  column_to_rownames("symbol") %>%
  as.matrix()

# Remove rows with all NA or zero values
heatmap_data <- heatmap_data[rowSums(is.na(heatmap_data)) < ncol(heatmap_data), ]
heatmap_data <- heatmap_data[rowSums(heatmap_data != 0, na.rm = TRUE) > 0, ]

# Normalize data for better visualization
if (nrow(heatmap_data) > 0) {
  heatmap_data <- scale(t(heatmap_data))  # Z-score normalization
} else {
  stop("No valid data available for heatmap after filtering.")
}

# Define color palette
heatmap_colors <- colorRamp2(c(-2, 0, 2), c(pal_nejm()(6)[6], "white", pal_nejm()(3)[3]))

# Create horizontal heatmap
cat("Creating horizontal heatmap...\n")
heatmap <- Heatmap(
            heatmap_data,
            name = "Expression (Z-score)",
            col = heatmap_colors,
            row_names_side = "left",
            column_names_side = "top",
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 12, fontface = "bold"),
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            column_title = "Median ICP gene expression",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
                title = "Z-score",
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                labels_gp = gpar(fontsize = 9)
            )
            )

# Save heatmap as PDF
cat("Saving heatmap...\n")
pdf("results/heatmap_median_gene_expression_icp_horizontal.pdf", width = 12, height = 6)
draw(heatmap)
dev.off()