# Check Data Availability Across Samples and Tissues
# Barbara Verhaar

# Libraries
library(tidyverse)

# Load data
cat("Loading data...\n")
global_expr <- readRDS("data/RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS")

# Remove version numbers from Ensembl IDs in global_expr
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

cat("\n=== SUMMARY ===\n")
cat("Total samples in dataset:", length(all_samples), "\n")
cat("Relevant tissue samples:", nrow(relevant_samples), "\n")
cat("  - Liver:", sum(relevant_samples$tissue == "Liver"), "\n")
cat("  - Visceral Fat:", sum(relevant_samples$tissue == "VisceralFat"), "\n")
cat("  - Subcutaneous Fat:", sum(relevant_samples$tissue == "SubcutaneousFat"), "\n")

# Check which IDs have which tissues
tissue_availability <- relevant_samples %>%
  group_by(sample_id, tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = tissue,
    values_from = n,
    values_fill = 0
  )

cat("\n=== UNIQUE SAMPLE IDs ===\n")
cat("Total unique IDs:", nrow(tissue_availability), "\n\n")

# Convert to binary (0/1) for presence/absence
tissue_availability_binary <- tissue_availability %>%
  mutate(
    Liver = ifelse(Liver > 0, 1, 0),
    VisceralFat = ifelse(VisceralFat > 0, 1, 0),
    SubcutaneousFat = ifelse(SubcutaneousFat > 0, 1, 0)
  )

# Count different patterns of data availability
pattern_summary <- tissue_availability_binary %>%
  count(Liver, VisceralFat, SubcutaneousFat) %>%
  arrange(desc(n))

cat("=== DATA AVAILABILITY PATTERNS ===\n")
print(pattern_summary, n = Inf)

# IDs with all three tissues
all_three <- tissue_availability_binary %>%
  filter(Liver == 1 & VisceralFat == 1 & SubcutaneousFat == 1)

cat("\n=== COMPLETE DATA (All 3 Tissues) ===\n")
cat("Number of IDs with all three tissues:", nrow(all_three), "\n")
if(nrow(all_three) > 0) {
  cat("IDs with complete data:\n")
  print(all_three$sample_id)
}

# IDs missing at least one tissue
incomplete <- tissue_availability_binary %>%
  filter(!(Liver == 1 & VisceralFat == 1 & SubcutaneousFat == 1))

cat("\n=== INCOMPLETE DATA (Missing at least 1 tissue) ===\n")
cat("Number of IDs with incomplete data:", nrow(incomplete), "\n")
if(nrow(incomplete) > 0) {
  cat("\nDetails of incomplete samples:\n")
  print(incomplete, n = Inf)
}

# Check for missing data patterns
cat("\n=== MISSING DATA PATTERNS ===\n")
missing_liver <- sum(tissue_availability_binary$Liver == 0)
missing_visceral <- sum(tissue_availability_binary$VisceralFat == 0)
missing_subcutaneous <- sum(tissue_availability_binary$SubcutaneousFat == 0)

cat("IDs missing Liver:", missing_liver, "\n")
cat("IDs missing Visceral Fat:", missing_visceral, "\n")
cat("IDs missing Subcutaneous Fat:", missing_subcutaneous, "\n")

# Save detailed results
write_csv(tissue_availability_binary, "results/tissue_availability_by_id.csv")
cat("\n=== SAVED ===\n")
cat("Detailed availability table saved to: results/tissue_availability_by_id.csv\n")

# Optional: Create a visual summary
library(ggplot2)

availability_long <- tissue_availability_binary %>%
  pivot_longer(
    cols = c(Liver, VisceralFat, SubcutaneousFat),
    names_to = "Tissue",
    values_to = "Available"
  ) %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Liver", "VisceralFat", "SubcutaneousFat")),
    Available = factor(Available, levels = c(0, 1), labels = c("Missing", "Available"))
  )

p <- ggplot(availability_long, aes(x = sample_id, y = Tissue, fill = Available)) +
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
    title = "Tissue Data Availability by Sample ID",
    x = "Sample ID",
    y = "Tissue Type",
    fill = "Status"
  )

ggsave("results/tissue_availability_heatmap.pdf", p, width = 14, height = 6)
cat("Visualization saved to: results/tissue_availability_heatmap.pdf\n")