# Check RNA Availability per Tissue for Serum/Diet Overlap
# Barbara Verhaar

# Libraries
library(tidyverse)

# Load data
cat("Loading data...\n")

# Serum IDs
serum <- readxl::read_xlsx("data/20250320_BARIA_300_Paris_Serum_Citraat_Nieuwdorp_MW.xlsx")
serum$ID <- str_c("BARIA_", serum$ID)
cat("Number of serum samples:", length(unique(serum$ID)), "\n")

# Diet data (filtered for serum IDs with TotalCal > 500)
diettot <- readRDS("data/bariatot.RDS")
diettot2 <- diettot |>
    filter(TotalCal_kcal > 500) |>
    filter(ID %in% serum$ID)
cat("Number of diet samples (filtered, overlapping with serum):", nrow(diettot2), "\n")

# RNA seq data
global_expr <- readRDS("data/RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS")
rownames(global_expr) <- gsub("\\..*", "", rownames(global_expr))

# Get all column names (samples)
all_samples <- colnames(global_expr)

# Extract sample ID and tissue information
sample_info <- tibble(
  full_sample = all_samples
) %>%
  mutate(
    # Extract the ID part (everything before the tissue identifier)
    sample_id = str_extract(full_sample, "^[^.]+"),
    # Add BARIA_ prefix to match other datasets
    ID = str_c("BARIA_", sample_id),
    # Extract tissue type
    tissue = case_when(
      str_detect(full_sample, "\\.Liver\\.V1$") ~ "Liver",
      str_detect(full_sample, "\\.vFat\\.V1$") ~ "VisceralFat",
      str_detect(full_sample, "\\.subFat\\.V1$") ~ "SubcutaneousFat",
      TRUE ~ "Other"
    )
  )

# Filter to only relevant tissues
relevant_samples <- sample_info %>%
  filter(tissue %in% c("Liver", "VisceralFat", "SubcutaneousFat"))

# Get unique RNAseq IDs
idsrnaseq <- unique(relevant_samples$ID)
cat("Number of unique IDs in RNAseq data:", length(idsrnaseq), "\n\n")

# Check overlap between serum and RNAseq
cat("=== OVERLAP: SERUM & RNAseq ===\n")
overlap_serum_rna <- serum$ID[serum$ID %in% idsrnaseq]
cat("Serum IDs with RNAseq data:", length(overlap_serum_rna), "\n")
cat("Serum IDs without RNAseq data:", sum(!serum$ID %in% idsrnaseq), "\n\n")

# Check overlap between diet and RNAseq
cat("=== OVERLAP: DIET & RNAseq ===\n")
overlap_diet_rna <- diettot2$ID[diettot2$ID %in% idsrnaseq]
cat("Diet IDs with RNAseq data:", length(overlap_diet_rna), "\n")
cat("Diet IDs without RNAseq data:", sum(!diettot2$ID %in% idsrnaseq), "\n\n")

# Get the triple overlap: Serum + Diet + RNAseq
cat("=== TRIPLE OVERLAP: SERUM + DIET + RNAseq ===\n")
overlap_all <- intersect(intersect(serum$ID, diettot2$ID), idsrnaseq)
cat("IDs with all three (Serum + Diet + RNAseq):", length(overlap_all), "\n\n")

# Now check tissue availability for the triple overlap
tissue_availability <- relevant_samples %>%
  filter(ID %in% overlap_all) %>%
  group_by(ID, tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = tissue,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    Liver = ifelse(Liver > 0, 1, 0),
    VisceralFat = ifelse(VisceralFat > 0, 1, 0),
    SubcutaneousFat = ifelse(SubcutaneousFat > 0, 1, 0)
  )

cat("=== TISSUE AVAILABILITY FOR TRIPLE OVERLAP (Serum + Diet + RNAseq) ===\n")
cat("Total IDs in overlap:", nrow(tissue_availability), "\n\n")

# Count different patterns of data availability
pattern_summary <- tissue_availability %>%
  count(Liver, VisceralFat, SubcutaneousFat) %>%
  arrange(desc(n))

cat("Tissue availability patterns:\n")
print(pattern_summary, n = Inf)

# IDs with all three tissues
all_three_tissues <- tissue_availability %>%
  filter(Liver == 1 & VisceralFat == 1 & SubcutaneousFat == 1)

cat("\n=== COMPLETE DATA (Serum + Diet + All 3 Tissues) ===\n")
cat("Number of IDs with complete data:", nrow(all_three_tissues), "\n")
if(nrow(all_three_tissues) > 0) {
  cat("\nFirst 20 IDs with complete data:\n")
  print(head(all_three_tissues$ID, 20))
}

# IDs missing at least one tissue
incomplete <- tissue_availability %>%
  filter(!(Liver == 1 & VisceralFat == 1 & SubcutaneousFat == 1))

cat("\n=== INCOMPLETE TISSUE DATA ===\n")
cat("Number of IDs with incomplete tissue data:", nrow(incomplete), "\n")

if(nrow(incomplete) > 0) {
  cat("\nBreakdown by missing tissue:\n")
  missing_liver <- sum(tissue_availability$Liver == 0)
  missing_visceral <- sum(tissue_availability$VisceralFat == 0)
  missing_subcutaneous <- sum(tissue_availability$SubcutaneousFat == 0)

  cat("  Missing Liver:", missing_liver, "\n")
  cat("  Missing Visceral Fat:", missing_visceral, "\n")
  cat("  Missing Subcutaneous Fat:", missing_subcutaneous, "\n")

  cat("\nIDs with incomplete tissue data:\n")
  print(incomplete, n = Inf)
}

# Summary counts per tissue type
cat("\n=== TISSUE-SPECIFIC COUNTS (Serum + Diet + RNAseq overlap) ===\n")
has_liver <- sum(tissue_availability$Liver == 1)
has_visceral <- sum(tissue_availability$VisceralFat == 1)
has_subcutaneous <- sum(tissue_availability$SubcutaneousFat == 1)

cat("Samples with Liver RNAseq:", has_liver, "\n")
cat("Samples with Visceral Fat RNAseq:", has_visceral, "\n")
cat("Samples with Subcutaneous Fat RNAseq:", has_subcutaneous, "\n")
cat("Samples with ALL 3 tissues:", nrow(all_three_tissues), "\n")

# Save results
write_csv(tissue_availability, "results/tissue_availability_serum_diet_overlap.csv")
cat("\n=== SAVED ===\n")
cat("Results saved to: results/tissue_availability_serum_diet_overlap.csv\n")

# Save list of complete IDs for easy filtering
complete_ids <- all_three_tissues$ID
saveRDS(complete_ids, "data/complete_ids_serum_diet_all3tissues.RDS")
cat("Complete IDs saved to: data/complete_ids_serum_diet_all3tissues.RDS\n")

# Create summary visualization
library(ggplot2)

# Summary bar plot
summary_data <- tibble(
  Category = c("Serum only", "Serum + RNAseq", "Serum + Diet + RNAseq", "Complete (All 3 tissues)"),
  Count = c(
    length(unique(serum$ID)),
    length(overlap_serum_rna),
    length(overlap_all),
    nrow(all_three_tissues)
  )
) %>%
  mutate(Category = factor(Category, levels = Category))

p <- ggplot(summary_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Data Availability Cascade: Serum → RNAseq → Complete Tissues",
    x = "",
    y = "Number of Samples"
  ) +
  ylim(0, max(summary_data$Count) * 1.1)

ggsave("results/data_overlap_summary.pdf", p, width = 10, height = 6)
cat("Summary plot saved to: results/data_overlap_summary.pdf\n")

# Detailed tissue heatmap for overlap
availability_long <- tissue_availability %>%
  pivot_longer(
    cols = c(Liver, VisceralFat, SubcutaneousFat),
    names_to = "Tissue",
    values_to = "Available"
  ) %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Liver", "VisceralFat", "SubcutaneousFat")),
    Available = factor(Available, levels = c(0, 1), labels = c("Missing", "Available")),
    ID = str_remove(ID, "BARIA_")  # Shorten for visualization
  )

p2 <- ggplot(availability_long, aes(x = ID, y = Tissue, fill = Available)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("Missing" = "#E64B35FF", "Available" = "#4DBBD5FF")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 10),
    legend.position = "top",
    panel.grid = element_blank()
  ) +
  labs(
    title = "Tissue Availability for Samples with Serum + Diet + RNAseq Data",
    x = "Sample ID",
    y = "Tissue Type",
    fill = "Status"
  )

ggsave("results/tissue_availability_overlap_heatmap.pdf", p2, width = 14, height = 6)
cat("Tissue heatmap saved to: results/tissue_availability_overlap_heatmap.pdf\n")
