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
  # Find the base Ensembl ID (without version number)
  base_id <- str_extract(ensembl_id, "ENSG[0-9]+")
  
  # Find all rows that match the ensembl ID pattern
  matches <- which(str_detect(global_ensembls, base_id))
  
  if(length(matches) == 0) {
    warning(paste0("No matches found for gene ID: ", ensembl_id))
    return(NULL)
  }
  
  # Use the first match if multiple are found
  match_to_use <- matches[1]
  ensembl_with_version <- global_ensembls[match_to_use]
  
  # Extract expression from different tissues - ONLY using .V1 columns
  liver_cols <- colnames(global_expr)[str_detect(colnames(global_expr), "\\.Liver\\.V1$")]
  vfat_cols <- colnames(global_expr)[str_detect(colnames(global_expr), "\\.vFat\\.V1$")]
  sfat_cols <- colnames(global_expr)[str_detect(colnames(global_expr), "\\.subFat\\.V1$")]
  
  # Extract data for each tissue
  liver_data <- as.data.frame(t(as.matrix(global_expr[match_to_use, liver_cols, drop=FALSE])))
  vfat_data <- as.data.frame(t(as.matrix(global_expr[match_to_use, vfat_cols, drop=FALSE])))
  sfat_data <- as.data.frame(t(as.matrix(global_expr[match_to_use, sfat_cols, drop=FALSE])))
  
  # Function to process tissue data with correct ID extraction
  process_tissue_data <- function(tissue_data, tissue_type) {
    if (nrow(tissue_data) == 0) {
      warning(paste0("No data found for tissue: ", tissue_type))
      return(NULL)
    }
    
    # Extract the subject ID from column names (format: "subjectID.Tissue.Visit")
    subject_ids <- gsub("^([0-9]+)\\..+$", "\\1", rownames(tissue_data))
    
    # Format as "BARIA_X" to match clinical data
    formatted_ids <- paste0("BARIA_", subject_ids)
    
    # Handle European decimal format
    values <- as.character(tissue_data[,1])
    values <- gsub(",", ".", values)  # Replace comma with period
    numeric_values <- as.numeric(values)
    
    # Apply log10 transformation to gene counts
    numeric_values <- log10(numeric_values + 1)
    
    # Create data frame with correctly formatted IDs
    result_df <- data.frame(
      value = numeric_values,
      ID = formatted_ids,
      stringsAsFactors = FALSE
    )
    
    # Use gene symbol in column names if provided
    symbol_to_use <- ifelse(is.null(gene_symbol) || gene_symbol == "", ensembl_id, gene_symbol)
    colnames(result_df)[1] <- paste0(symbol_to_use, "_", tissue_type)
    
    return(result_df)
  }
  
  # Process each tissue dataset with correct tissue type names
  liver_df <- process_tissue_data(liver_data, "liver")
  vfat_df <- process_tissue_data(vfat_data, "viscfat") 
  sfat_df <- process_tissue_data(sfat_data, "subcfat")
  
  # Skip if any tissue has no data
  if (is.null(liver_df) || is.null(vfat_df) || is.null(sfat_df)) {
    warning("One or more tissues have no data")
    return(NULL)
  }
  
  # Join data with metadata - IDs are already in the correct format
  gene_df <- full_join(liver_df, sfat_df, by = "ID") %>% 
    full_join(., vfat_df, by = "ID") %>% 
    left_join(., global_meta, by = "ID")
  
  return(gene_df)
}

# Add this function to normalize delta variable names
normalize_delta_name <- function(var_name) {
  # Fix "delta_delta" pattern
  if (str_detect(tolower(var_name), "delta_delta")) {
    var_name <- str_replace(var_name, "(?i)delta_delta", "delta_")
  }
  
  # Fix other variations like "delta_delta_"
  if (str_detect(tolower(var_name), "delta_delta_")) {
    var_name <- str_replace(var_name, "(?i)delta_delta_", "delta_")
  }
  
  # Fix case where "delta" appears twice without underscore
  if (str_detect(tolower(var_name), "delta.*delta")) {
    var_name <- str_replace(var_name, "(?i)delta(.*?)delta", "delta\\1")
  }
  
  # Standard formatting
  var_name <- gsub("\\.", "-", var_name)
  
  return(var_name)
}

# Function to create a correlation plot
create_correlation_plot <- function(data, x_var, y_var, x_lab, y_lab, gene_symbol, tissue, clinical_var, is_delta = FALSE) {
  # Fix delta naming issues
  display_var <- normalize_delta_name(clinical_var)
  
  # Define clinical variable transformations as in the heatmap script
  if (str_detect(clinical_var, "CRP") && !str_detect(clinical_var, "delta")) {
    y_expr <- "log10(v0_crp+1)"
    y_label <- "log10(CRP (mg/l))"
  } else if (str_detect(clinical_var, "Leukocytes") && !str_detect(clinical_var, "delta")) {
    y_expr <- "log10(v0_leuko+0.1)"
    y_label <- "Leukocytes x10E9/L"
  } else if (str_detect(clinical_var, "LDL") && !str_detect(clinical_var, "delta")) {
    y_expr <- "v0_ldl"
    y_label <- "LDL (mmol/L)"
  } else if ((str_detect(clinical_var, "HOMA-IR") || str_detect(clinical_var, "HOMA.IR")) && !str_detect(clinical_var, "delta")) {
    y_expr <- "log10(homaIR_v1+0.1)"
    y_label <- "HOMA-IR"
  } else if (str_detect(clinical_var, "Age") && !str_detect(clinical_var, "delta")) {
    y_expr <- "v0_age"
    y_label <- "Age (years)"
  } else if (str_detect(clinical_var, "BMI") && !str_detect(clinical_var, "delta")) {
    y_expr <- "v0_bmi"
    y_label <- "BMI (kg/m²)"
  } else {
    y_expr <- y_var
    y_label <- display_var
  }
  
  # Create directory path based on whether it's a delta correlation or not
  if (is_delta) {
    dir_path <- file.path("results", "correlations", "significant_delta_plots", gene_symbol, display_var)
  } else {
    dir_path <- file.path("results", "correlations", "significant_plots", gene_symbol, display_var)
  }
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  # Create plot
  p <- ggplot(data = data) + 
    geom_point(aes(x = .data[[x_var]], y = eval(parse(text = y_expr))), 
              alpha = 0.9, color = "royalblue") +
    geom_smooth(aes(x = .data[[x_var]], y = eval(parse(text = y_expr))), 
               method = "lm", color = "firebrick4") +
    stat_cor(aes(x = .data[[x_var]], y = eval(parse(text = y_expr))), 
            method = "spearman") +
    labs(
      x = paste0("log10(", gene_symbol, " counts + 1)"),
      y = y_label,
      title = paste0(gene_symbol, " in ", tissue, " vs ", display_var) # Use the normalized name
    ) +
    theme_Publication()
  
  # Save plot with normalized name
  filename <- paste0(display_var, "_", gene_symbol, "_", tissue, ".pdf")
  file_path <- file.path(dir_path, filename)
  ggsave(file_path, plot = p, width = 5, height = 5)
  
  return(p)
}

# Add this function to standardize clinical variables
standardize_clinical_var <- function(var_name, is_delta = FALSE) {
  # Clean up variable name
  var_name <- gsub("\\.", "-", var_name)  # Replace dots with hyphens
  
  # Return properly formatted name for display
  if (is_delta && !str_detect(tolower(var_name), "delta")) {
    return(paste0("delta_", var_name))
  } else {
    return(var_name)
  }
}

# Add this function to map clinical variables to column names
map_clinical_to_column <- function(var_name, is_delta = FALSE) {
  # First standardize the variable name
  std_name <- tolower(var_name)
  std_name <- gsub("delta_", "", std_name)
  std_name <- gsub("delta", "", std_name)
  std_name <- gsub("-", "", std_name)
  std_name <- gsub("\\.", "", std_name)
  
  # Map to actual column names
  if (str_detect(std_name, "crp")) {
    return(ifelse(is_delta, "crp_v1v4", "v0_crp"))
  } else if (str_detect(std_name, "homair")) {
    return(ifelse(is_delta, "homaIR_v1v4", "homaIR_v1"))
  } else if (str_detect(std_name, "bmi")) {
    return(ifelse(is_delta, "bmi_v1v4", "v0_bmi"))
  } else if (str_detect(std_name, "leuko")) {
    return(ifelse(is_delta, "leuko_v1v4", "v0_leuko"))
  } else if (str_detect(std_name, "ldl")) {
    return(ifelse(is_delta, "ldl_v1v4", "v0_ldl"))
  } else if (str_detect(std_name, "age")) {
    return("v0_age")
  } else {
    # For other variables, construct a sensible column name
    base_name <- gsub("[^a-z0-9]", "", std_name)
    return(ifelse(is_delta, paste0(base_name, "_v1v4"), paste0("v0_", base_name)))
  }
}

# Main function
main <- function() {
  # Load data
  cat("Loading RDS data...\n")
  global_meta <<- readRDS("data/baria_metadata.RDS")
  global_expr <<- readRDS("data/RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS")
  global_ensembls <<- rownames(global_expr)
  
  # Load gene list
  cat("Reading gene list...\n")
  gene_df <- readRDS("data/gene_list.RDS") %>% slice(1:45)
  gene_df$ensembl_id <- gene_df$ENSEMBL
  gene_df$symbol <- gsub("\\,.*", "", gsub(" \\(.*\\)", "", gene_df$ICP_symbol))
  
  # Define tissues and clinical variables
  tissues <- c("liver", "viscfat", "subcfat")
  tissue_labels <- c("Liver", "Visceral Fat", "Subcutaneous Fat")
  clinical_vars <- c("CRP", "Leukocytes", "LDL", "HOMA-IR", "Age", "BMI")
  
  # Load correlation and p-value matrices for baseline (v0) and delta (v0v4)
  all_correlations <- list()
  all_pvalues <- list()
  
  # Load baseline correlations
  for (tissue in tissues) {
    corr_file <- paste0("results/correlations/heatmaps_v0/correlations_", tissue, ".csv")
    pval_file <- paste0("results/correlations/heatmaps_v0/pvalues_", tissue, ".csv")
    
    if (file.exists(corr_file) && file.exists(pval_file)) {
      all_correlations[[paste0(tissue, "_v0")]] <- read.csv(corr_file, row.names = 1)
      all_pvalues[[paste0(tissue, "_v0")]] <- read.csv(pval_file, row.names = 1)
      cat(sprintf("Loaded baseline correlations for tissue %s\n", tissue))
    } else {
      cat(sprintf("Files for baseline tissue %s not found\n", tissue))
    }
  }
  
  # Load delta correlations
  for (tissue in tissues) {
    corr_file <- paste0("results/correlations/heatmaps_v0v4/correlations_", tissue, ".csv")
    pval_file <- paste0("results/correlations/heatmaps_v0v4/pvalues_", tissue, ".csv")
    
    if (file.exists(corr_file) && file.exists(pval_file)) {
      all_correlations[[paste0(tissue, "_v0v4")]] <- read.csv(corr_file, row.names = 1)
      all_pvalues[[paste0(tissue, "_v0v4")]] <- read.csv(pval_file, row.names = 1)
      cat(sprintf("Loaded delta correlations for tissue %s\n", tissue))
    } else {
      cat(sprintf("Files for delta tissue %s not found\n", tissue))
    }
  }
  
  # Create a data frame to store all significant correlations
  significant_correlations <- data.frame(
    gene = character(),
    ensembl_id = character(),
    tissue = character(),
    clinical_var = character(),
    correlation = numeric(),
    pvalue = numeric(),
    is_delta = logical(),
    stringsAsFactors = FALSE
  )
  
  # Find significant correlations across all tissues (baseline and delta)
  for (key in names(all_correlations)) {
    parts <- strsplit(key, "_")[[1]]
    tissue <- parts[1]
    type <- parts[2]  # v0 or v0v4
    is_delta <- (type == "v0v4")
    
    corr_matrix <- all_correlations[[key]]
    pval_matrix <- all_pvalues[[key]]
    
    for (i in 1:nrow(corr_matrix)) {
      for (j in 1:ncol(corr_matrix)) {
        if (!is.na(pval_matrix[i, j]) && pval_matrix[i, j] < 0.05) {
          gene_symbol <- rownames(corr_matrix)[i]
          clinical_var <- colnames(corr_matrix)[j]
          
          # For delta correlations, handle the prefix more carefully
          if (is_delta) {
            # Only add delta_ if it doesn't already exist in the name
            if (!str_detect(tolower(clinical_var), "delta")) {
              clinical_var <- paste0("delta_", clinical_var)
            }
            # Make sure it's normalized even if it already had delta
            clinical_var <- normalize_delta_name(clinical_var)
          }
          
          # Find corresponding ensembl ID
          ensembl_id <- gene_df$ensembl_id[gene_df$symbol == gene_symbol]
          if (length(ensembl_id) == 0) ensembl_id <- NA
          else ensembl_id <- ensembl_id[1]
          
          # Add to significant correlations
          significant_correlations <- rbind(
            significant_correlations,
            data.frame(
              gene = gene_symbol,
              ensembl_id = ensembl_id,
              tissue = tissue,
              clinical_var = clinical_var,
              correlation = corr_matrix[i, j],
              pvalue = pval_matrix[i, j],
              is_delta = is_delta,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }
  
  # Sort by significance
  significant_correlations <- significant_correlations %>%
    arrange(pvalue)
  
  # Save the list of significant correlations
  write.csv(
    significant_correlations,
    "results/correlations/significant_correlations_all.csv",
    row.names = FALSE
  )
  
  cat(sprintf("Found %d significant correlations (baseline and delta)\n", nrow(significant_correlations)))
  
  # Create plots for each significant correlation
  plots_created <- 0
  
  for (i in 1:nrow(significant_correlations)) {
    gene_symbol <- significant_correlations$gene[i]
    ensembl_id <- significant_correlations$ensembl_id[i]
    tissue <- significant_correlations$tissue[i]
    clinical_var <- significant_correlations$clinical_var[i]
    is_delta <- significant_correlations$is_delta[i]
    
    if (is.na(ensembl_id)) {
      cat(sprintf("Skipping %s (no ensembl ID found)\n", gene_symbol))
      next
    }

    if (is_delta) {
      # Check if clinical_var already contains "delta" in some form
      if (!str_detect(tolower(clinical_var), "delta")) {
        clinical_var <- paste0("delta_", clinical_var)
      }
    }
    
    # Extract gene data
    gene_data <- extract_gene_data(ensembl_id, gene_symbol)
    
    if (!is.null(gene_data)) {
      # Create variable names
      x_var <- paste0(gene_symbol, "_", tissue)
      
      # Map clinical variables to their actual column names in the data
      if (clinical_var == "CRP") {
        y_var <- "v0_crp"
      } else if (clinical_var == "Leukocytes") {
        y_var <- "v0_leuko"
      } else if (clinical_var == "LDL") {
        y_var <- "v0_ldl"
      } else if (clinical_var == "HOMA-IR" || clinical_var == "HOMA.IR") { # Handle both versions
        y_var <- "homaIR_v1"
      } else if (clinical_var == "Age") {
        y_var <- "v0_age"
      } else if (clinical_var == "BMI") {
        y_var <- "v0_bmi"
      } else if (str_detect(clinical_var, "delta_") || str_detect(tolower(clinical_var), "delta")) {
        # Handle delta variables - normalize the name
        clinical_var_clean <- str_replace_all(clinical_var, "(?i)delta_|delta", "")
        
        # Map to actual column names
        if (str_detect(clinical_var_clean, "CRP")) {
          y_var <- "crp_v1v4" 
        } else if (str_detect(clinical_var_clean, "HOMA.IR|HOMA-IR")) {
          y_var <- "homaIR_v1v4"
        } else if (str_detect(clinical_var_clean, "BMI")) {
          y_var <- "bmi_v1v4"
        } else if (str_detect(clinical_var_clean, "Leuko")) {
          y_var <- "leuko_v1v4"
        } else if (str_detect(clinical_var_clean, "LDL")) {
          y_var <- "ldl_v1v4"
        } else {
          # For other variables, construct a column name
          y_var <- paste0(tolower(clinical_var_clean), "_v1v4")
          cat(sprintf("Delta variable: %s mapped to %s\n", clinical_var, y_var))
        }
      }
      else {
        # Add this default case to handle unexpected clinical variable names
        cat(sprintf("Unknown clinical variable: %s - attempting to standardize\n", clinical_var))
        
        # Try to standardize common variations in clinical variable names
        standardized_var <- clinical_var
        standardized_var <- gsub("\\.", "-", standardized_var)  # Replace dots with hyphens
        
        if (standardized_var == "HOMA-IR") {
          y_var <- "homaIR_v1"
        } else {
          # If still unknown, use lowercase version as a fallback
          y_var <- tolower(clinical_var)
          cat(sprintf("Using column name: %s\n", y_var))
        }
      }
      # Check if both variables exist in the data
      if (x_var %in% colnames(gene_data) && y_var %in% colnames(gene_data)) {
        # Create plot
        tissue_label <- tissue_labels[which(tissues == tissue)]
        create_correlation_plot(
          gene_data, 
          x_var, 
          y_var,
          paste0(gene_symbol, " (counts)"),
          clinical_var,
          gene_symbol,
          tissue,
          clinical_var,
          is_delta
        )
        plots_created <- plots_created + 1
      } else {
        cat(sprintf("Variables not found for %s in %s vs %s\n", gene_symbol, tissue, clinical_var))
      }
    } else {
      cat(sprintf("Could not extract data for %s\n", gene_symbol))
    }
  }
  
  cat(sprintf("Created %d correlation plots\n", plots_created))
}

# Run the main function
main()