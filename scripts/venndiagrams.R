# Venn Diagram for Overlap Between Significant Correlations
# Barbara Verhaar

# Libraries
library(tidyverse)
library(ggVennDiagram)
library(RColorBrewer)
library(showtext)

# Function to create and save Venn diagrams
create_venn <- function(data_list, title, output_file) {
  cat(sprintf("Creating Venn diagram for %s...\n", title))
  
  venn_plot <- ggVennDiagram(data_list, label_alpha = 0, set_name_gp = gpar(fontface = "bold", fontsize = 14)) +
    scale_color_manual(values = rep("black", 3)) +
    scale_fill_gradient(low = "#D6EAF8", high = "#2E86C1", guide = "none") +  # Shades of blue
    theme(
      text = element_text(family = "roboto", size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "roboto"),
      plot.margin = margin(20, 20, 20, 20)  # Add padding around the graph
    ) +
    ggtitle(title)
  
  # Save the Venn diagram as a PDF
  ggsave(output_file, venn_plot, width = 10, height = 10)  # Increased dimensions for better fit
  cat(sprintf("Saved Venn diagram for %s to %s\n", title, output_file))
}

# Load custom fonts
font_add_google("Roboto", "roboto")  # Add Google font "Roboto"
showtext_auto()  # Enable showtext for custom fonts

# Load significant correlations data
cat("Loading significant correlations data...\n")
sig_liver <- read.csv("results/correlations/heatmaps_v0v4/pvalues_liver.csv")
sig_sfat <- read.csv("results/correlations/heatmaps_v0v4/pvalues_subcfat.csv")
sig_vfat <- read.csv("results/correlations/heatmaps_v0v4/pvalues_viscfat.csv")

siglivercrp <- sig_liver %>% filter(deltaCRP < 0.05) %>% pull(X)
sigsfatcrp <- sig_sfat %>% filter(deltaCRP < 0.05) %>% pull(X)
sigvfatcrp <- sig_vfat %>% filter(deltaCRP < 0.05) %>% pull(X)
sigliverbmi <- sig_liver %>% filter(deltaBMI < 0.05) %>% pull(X)
sigsfatbmi <- sig_sfat %>% filter(deltaBMI < 0.05) %>% pull(X)
sigvfatbmi <- sig_vfat %>% filter(deltaBMI < 0.05) %>% pull(X)
sigliverhomair <- sig_liver %>% filter(deltaHOMA.IR < 0.05) %>% pull(X)
sigsfathomair <- sig_sfat %>% filter(deltaHOMA.IR < 0.05) %>% pull(X)
sigvfathomair <- sig_vfat %>% filter(deltaHOMA.IR < 0.05) %>% pull(X)

# Create Venn diagrams for CRP, BMI, and HOMA-IR
create_venn(
  data_list = list(
    Liver = siglivercrp,
    SubcutaneousFat = sigsfatcrp,
    VisceralFat = sigvfatcrp
  ),
  title = "Overlap CRP correlations",
  output_file = "results/venn_diagram_crp.pdf"
)

create_venn(
  data_list = list(
    Liver = sigliverbmi,
    SubcutaneousFat = sigsfatbmi,
    VisceralFat = sigvfatbmi
  ),
  title = "Overlap BMI correlations",
  output_file = "results/venn_diagram_bmi.pdf"
)

create_venn(
  data_list = list(
    Liver = sigliverhomair,
    SubcutaneousFat = sigsfathomair,
    VisceralFat = sigvfathomair
  ),
  title = "Overlap HOMA-IR correlations",
  output_file = "results/venn_diagram_homair.pdf"
)


# Load significant correlations data
cat("Loading significant correlations data...\n")
sig_liver <- read.csv("results/correlations/heatmaps_v0/pvalues_liver.csv")
sig_sfat <- read.csv("results/correlations/heatmaps_v0/pvalues_subcfat.csv")
sig_vfat <- read.csv("results/correlations/heatmaps_v0/pvalues_viscfat.csv")

siglivercrp <- sig_liver %>% filter(CRP < 0.05) %>% pull(X)
sigsfatcrp <- sig_sfat %>% filter(CRP < 0.05) %>% pull(X)
sigvfatcrp <- sig_vfat %>% filter(CRP < 0.05) %>% pull(X)
sigliverbmi <- sig_liver %>% filter(BMI < 0.05) %>% pull(X)
sigsfatbmi <- sig_sfat %>% filter(BMI < 0.05) %>% pull(X)
sigvfatbmi <- sig_vfat %>% filter(BMI < 0.05) %>% pull(X)
sigliverhomair <- sig_liver %>% filter(HOMA.IR < 0.05) %>% pull(X)
sigsfathomair <- sig_sfat %>% filter(HOMA.IR < 0.05) %>% pull(X)
sigvfathomair <- sig_vfat %>% filter(HOMA.IR < 0.05) %>% pull(X)

# Create Venn diagrams for CRP, BMI, and HOMA-IR
create_venn(
  data_list = list(
    Liver = siglivercrp,
    SubcutaneousFat = sigsfatcrp,
    VisceralFat = sigvfatcrp
  ),
  title = "Overlap CRP correlations",
  output_file = "results/venn_diagram_crp_baseline.pdf"
)

create_venn(
  data_list = list(
    Liver = sigliverbmi,
    SubcutaneousFat = sigsfatbmi,
    VisceralFat = sigvfatbmi
  ),
  title = "Overlap BMI correlations",
  output_file = "results/venn_diagram_bmi_baseline.pdf"
)

create_venn(
  data_list = list(
    Liver = sigliverhomair,
    SubcutaneousFat = sigsfathomair,
    VisceralFat = sigvfathomair
  ),
  title = "Overlap HOMA-IR correlations",
  output_file = "results/venn_diagram_baseline.pdf"
)
