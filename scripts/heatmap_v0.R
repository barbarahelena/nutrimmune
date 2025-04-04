# Heatmap of genes of interest (baseline)
# Barbara Verhaar

# Libraries
library(tidyverse)
library(ggpubr)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(grid)

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
  
  # Print information about found columns for debugging
  cat(sprintf("Found %d Liver.V1, %d vFat.V1, and %d subFat.V1 columns\n", 
              length(liver_cols), length(vfat_cols), length(sfat_cols)))
  
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
  
  # Use gene symbol in column names if available
  col_prefix <- ifelse(is.null(gene_symbol) || gene_symbol == "", ensembl_id, gene_symbol)
  
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
  
  return(list(
    data = gene_df,
    ensembl_id = ensembl_id,
    symbol = col_prefix,
    ensembl_with_version = ensembl_with_version
  ))
}

# Function to calculate correlations with clinical variables
calculate_correlations <- function(gene_data) {
  gene_df <- gene_data$data
  gene_symbol <- gene_data$symbol
  
  # Define clinical variables and their expressions
  clinical_vars <- list(
    "CRP" = "log10(v0_crp+1)",
    "Leukocytes" = "log10(v0_leuko+0.1)",
    "LDL" = "v0_ldl",
    "HOMA-IR" = "log10(homaIR_v1+0.1)",
    "Age" = "v0_age",
    "BMI" = "v0_bmi"
  )
  
  tissues <- c("liver", "viscfat", "subcfat")
  
  # Initialize results matrices for correlations and p-values
  results <- matrix(NA, nrow = length(tissues), ncol = length(clinical_vars))
  rownames(results) <- tissues
  colnames(results) <- names(clinical_vars)
  
  p_values <- matrix(NA, nrow = length(tissues), ncol = length(clinical_vars))
  rownames(p_values) <- tissues
  colnames(p_values) <- names(clinical_vars)
  
  # Calculate correlations for each tissue and clinical variable
  for (i in seq_along(tissues)) {
    tissue <- tissues[i]
    gene_col <- paste0(gene_symbol, "_", tissue)
    
    # Skip if column doesn't exist
    if (!gene_col %in% colnames(gene_df)) {
      cat(sprintf("Warning: Column %s not found for gene %s\n", gene_col, gene_symbol))
      next
    }
    
    for (j in seq_along(clinical_vars)) {
      var_name <- names(clinical_vars)[j]
      var_expr <- clinical_vars[[j]]
      
      # Calculate Spearman correlation
      tryCatch({
        # Parse and evaluate the expression
        y_values <- eval(parse(text = var_expr), gene_df)
        
        # Remove rows with NA values
        valid_idx <- !is.na(gene_df[[gene_col]]) & !is.na(y_values)
        
        if (sum(valid_idx) > 5) {  # Ensure enough data points
          cor_test <- cor.test(
            x = gene_df[[gene_col]][valid_idx],
            y = y_values[valid_idx],
            method = "spearman",
            exact = FALSE
          )
          
          # Store correlation coefficient and p-value
          results[i, j] <- cor_test$estimate
          p_values[i, j] <- cor_test$p.value
        } else {
          cat(sprintf("Not enough valid data points for %s and %s (only %d points)\n", 
                    gene_col, var_name, sum(valid_idx)))
        }
      }, error = function(e) {
        cat(sprintf("Error calculating correlation for %s and %s: %s\n", 
                  gene_col, var_name, e$message))
      })
    }
  }

  # Check for valid results
  valid_corr_rows <- sum(!is.na(results))
  if (valid_corr_rows == 0) {
    cat(sprintf("No valid correlations found for gene %s\n", gene_symbol))
    return(NULL)
  }

  # Check for missing tissues in results
  missing_tissues <- tissues[!tissues %in% rownames(results)]
  if (length(missing_tissues) > 0) {
    cat(sprintf("Warning: Missing tissues for gene %s: %s\n", 
                gene_symbol, paste(missing_tissues, collapse=", ")))
  }

  return(list(
    correlations = results,
    p_values = p_values,
    gene_id = gene_data$ensembl_id,
    gene_symbol = gene_symbol
  ))
}

# Function to create significance symbols
get_significance_symbols <- function(p_values) {
  symbols <- matrix("", nrow = nrow(p_values), ncol = ncol(p_values))
  symbols[p_values < 0.05] <- "*"
  symbols[p_values < 0.01] <- "**"
  symbols[p_values < 0.001] <- "***"
  symbols[p_values < 0.0001] <- "****"
  return(symbols)
}

# Function to create heatmaps using ComplexHeatmap
create_complex_heatmaps <- function(all_results) {
  # Create output directory
  dir.create("results/correlations/heatmaps_v0", showWarnings = FALSE, recursive = TRUE)
  
  # Extract gene information
  gene_ids <- sapply(all_results, function(x) x$gene_id)
  gene_symbols <- sapply(all_results, function(x) x$gene_symbol)
  
  # Create unique labels
  gene_labels <- gene_symbols
  
  tissues <- c("liver", "viscfat", "subcfat")
  tissue_labels <- c("Liver", "Visceral Fat", "Subcutaneous Fat")
  clinical_vars <- c("CRP", "Leukocytes", "LDL", "HOMA-IR", "Age", "BMI")
  
  # Define color schemes
  corr_col_fun <- colorRamp2(c(-0.5, 0, 0.5), c(pal_nejm()(6)[6], "white", pal_nejm()(3)[3]))
  pvalue_col_fun <- colorRamp2(c(0, 1, 3), c("lightgrey", "white", pal_nejm()(5)[5]))
  
  # Create separate heatmap for each tissue
  for (i in seq_along(tissues)) {
    tissue <- tissues[i]
    tissue_label <- tissue_labels[i]
    
    # Create correlation and p-value matrices
    corr_matrix <- matrix(NA, nrow = length(all_results), ncol = length(clinical_vars))
    rownames(corr_matrix) <- gene_labels
    colnames(corr_matrix) <- clinical_vars
    
    p_matrix <- matrix(NA, nrow = length(all_results), ncol = length(clinical_vars))
    rownames(p_matrix) <- gene_labels
    colnames(p_matrix) <- clinical_vars
    
    # Fill matrices
    for (j in seq_along(all_results)) {
      # Add error checking to handle missing tissues or dimensions mismatch
      tryCatch({
        # Check if tissue exists in correlations
        if (tissue %in% rownames(all_results[[j]]$correlations)) {
          # Check if dimensions match
          corr_vals <- all_results[[j]]$correlations[tissue, ]
          p_vals <- all_results[[j]]$p_values[tissue, ]
          
          if (length(corr_vals) == ncol(corr_matrix)) {
            corr_matrix[j, ] <- corr_vals
            p_matrix[j, ] <- p_vals
          } else {
            cat(sprintf("Dimension mismatch for gene %s, tissue %s: expected %d values, got %d\n", 
                      gene_labels[j], tissue, ncol(corr_matrix), length(corr_vals)))
          }
        } else {
          cat(sprintf("Missing tissue '%s' for gene %s\n", tissue, gene_labels[j]))
        }
      }, error = function(e) {
        cat(sprintf("Error processing gene %s for tissue %s: %s\n", 
                  gene_labels[j], tissue, e$message))
      })
    }
    
    # Check if we have enough data for the heatmap
    if (all(is.na(corr_matrix))) {
      cat(sprintf("No valid correlations for tissue %s - skipping heatmap\n", tissue_label))
      next
    }
    
    # Get significance symbols
    sig_symbols <- get_significance_symbols(p_matrix)
    
    # Create cell function to add significance symbols
    cell_fun <- function(j, i, x, y, width, height, fill) {
      if (!is.na(p_matrix[i, j]) && p_matrix[i, j] < 0.05) {
        # Center the star horizontally and vertically in each cell
        grid.text(sig_symbols[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold"), 
                  just = "center", vjust = 0.75, hjust = 0.5)
      }
    }
    
    # Determine if we can cluster (if there are enough non-NA values)
    can_cluster_rows <- sum(complete.cases(corr_matrix)) >= 2
    can_cluster_cols <- sum(complete.cases(t(corr_matrix))) >= 2
    
    # Create heatmap
    ht <- Heatmap(
      corr_matrix,
      name = "Correlation",
      col = corr_col_fun,
      rect_gp = gpar(col = "white", lwd = 1),
      cell_fun = cell_fun,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 10),
      column_names_gp = gpar(fontsize = 12, fontface = "bold"),
      cluster_rows = TRUE,
      show_row_dend = FALSE,
      cluster_columns = TRUE,
      show_column_dend = FALSE,
      column_title = paste0("Correlations in ", tissue_label),
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      # width = unit(6, "cm"),
      # height = unit(max(8, length(all_results) * 0.3), "cm"),
      heatmap_legend_param = list(
        title = "Spearman\nCorrelation",
        at = c(-0.5, -0.25, 0, 0.25, 0.5),
        labels = c("-0.5", "-0.25", "0", "0.25", "0.5"),
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9)
      )
    )
    
    # Create legends 
    lgd_sig <- Legend(
      pch = c("*", "**", "***", "****"), 
      type = "points",
      labels = c("p<0.05", "p<0.01", "p<0.001", "p<0.0001"),
      title = "Significance",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      grid_height = unit(6, "mm"),
      grid_width = unit(10, "mm")
    )
    
    # Save the plot
    pdf(paste0("results/correlations/heatmaps_v0/heatmap_", tissue, ".pdf"), 
        width = 5, height = 9)
    draw(ht, annotation_legend_list = list(lgd_sig))
    dev.off()
    
    # Save correlation and p-value matrices
    write.csv(
      corr_matrix,
      paste0("results/correlations/heatmaps_v0/correlations_", tissue, ".csv")
    )
    
    write.csv(
      p_matrix,
      paste0("results/correlations/heatmaps_v0/pvalues_", tissue, ".csv")
    )
  }
  
  # Create combined heatmap with all tissues
  all_corrs <- matrix(NA, nrow = length(all_results) * length(tissues), ncol = length(clinical_vars))
  all_pvals <- matrix(NA, nrow = length(all_results) * length(tissues), ncol = length(clinical_vars))
  row_labels <- character(length(all_results) * length(tissues))
  tissue_annot <- character(length(all_results) * length(tissues))
  
  # Fill the combined matrices - some corrs may be missing so we need to check
  row_idx <- 1
  for (i in seq_along(all_results)) {
    for (j in seq_along(tissues)) {
      tissue <- tissues[j]
      tissue_label <- tissue_labels[j]
      
      # Check if tissue exists in the correlations
      if (tissue %in% rownames(all_results[[i]]$correlations)) {
        tryCatch({
          all_corrs[row_idx, ] <- all_results[[i]]$correlations[tissue, ]
          all_pvals[row_idx, ] <- all_results[[i]]$p_values[tissue, ]
          row_labels[row_idx] <- gene_labels[i]
          tissue_annot[row_idx] <- tissue_label
          row_idx <- row_idx + 1
        }, error = function(e) {
          cat(sprintf("Error processing gene %s for tissue %s in combined heatmap: %s\n",
                    gene_labels[i], tissue_label, e$message))
        })
      } else {
        cat(sprintf("Skipping missing tissue %s for gene %s in combined heatmap\n", 
                  tissue, gene_labels[i]))
      }
    }
  }
  
  # Trim matrices if not all rows were filled
  if (row_idx <= nrow(all_corrs)) {
    all_corrs <- all_corrs[1:(row_idx-1), , drop=FALSE]
    all_pvals <- all_pvals[1:(row_idx-1), , drop=FALSE]
    row_labels <- row_labels[1:(row_idx-1)]
    tissue_annot <- tissue_annot[1:(row_idx-1)]
  }
  
  # Check if we have enough data for the combined heatmap
  if (nrow(all_corrs) == 0 || all(is.na(all_corrs))) {
    cat("Not enough data for combined heatmap - skipping\n")
    return()
  }
  
  # Add row and column names
  rownames(all_corrs) <- row_labels
  colnames(all_corrs) <- clinical_vars
  rownames(all_pvals) <- row_labels
  colnames(all_pvals) <- clinical_vars
  
  # Define tissue colors
  tissue_colors <- c(
    "Liver" = "#E41A1C",
    "Visceral Fat" = "#377EB8", 
    "Subcutaneous Fat" = "#4DAF4A"
  )
  
  # Create annotation
  ha <- HeatmapAnnotation(
    "Tissue" = anno_block(
      gp = gpar(fill = NA, col = "black"),
      labels = unique(tissue_annot),
      labels_gp = gpar(fontsize = 12, fontface = "bold")
    ),
    show_annotation_name = FALSE,
    height = unit(0.4, "cm")
  )
  
  row_ha <- rowAnnotation(
    "Tissue" = tissue_annot,
    col = list(Tissue = tissue_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_rot = 90,
    width = unit(0.6, "cm")
  )
  
  # Get significance symbols
  sig_symbols <- get_significance_symbols(all_pvals)
  
  # Create cell function to add significance symbols
  cell_fun <- function(j, i, x, y, width, height, fill) {
    if (!is.na(all_pvals[i, j]) && all_pvals[i, j] < 0.05) {
      # Center the star horizontally and vertically in each cell
      grid.text(sig_symbols[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold"), 
                just = "center", vjust = 0.75, hjust = 0.5)
    }
  }
  
  # Create combined heatmap
  ht_all <- Heatmap(
    all_corrs,
    name = "Correlation",
    col = corr_col_fun,
    rect_gp = gpar(col = "white", lwd = 1),
    cell_fun = cell_fun,
    row_split = tissue_annot,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    show_row_dend = FALSE,
    cluster_columns = TRUE,
    show_column_dend = FALSE,
    column_title = "Gene Correlations with Clinical Variables Across Tissues",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    left_annotation = row_ha,
    # width = unit(6, "cm"),
    # height = unit(max(15, nrow(all_corrs) * 0.25), "cm"),
    heatmap_legend_param = list(
      title = "Spearman\nCorrelation",
      at = c(-0.5, -0.25, 0, 0.25, 0.5),
      labels = c("-0.5", "-0.25", "0", "0.25", "0.5"),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )
  
  # Create legends for combined heatmap
  lgd_sig <- Legend(
    pch = c("*", "**", "***", "****"), 
    type = "points",
    labels = c("p<0.05", "p<0.01", "p<0.001", "p<0.0001"),
    title = "Significance",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    grid_height = unit(6, "mm"),
    grid_width = unit(10, "mm")
  )
  
  lgd_tissue <- Legend(
    title = "Tissue",
    legend_gp = gpar(fill = tissue_colors),
    labels = names(tissue_colors),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
  
  # Save the combined plot
  pdf("results/correlations/heatmaps_v0/heatmap_all_tissues.pdf", 
      width = 9, 
      height = max(15, nrow(all_corrs) * 0.25) + 2)
  draw(ht_all, annotation_legend_list = list(lgd_sig))
  dev.off()
  
  # Save as SVG for better editing possibilities
  svg("results/correlations/heatmaps_v0/heatmap_all_tissues.svg", 
      width = 6, 
      height = 25)
  draw(ht_all, annotation_legend_list = list(lgd_sig))
  dev.off()
  
  # Save combined data
  write.csv(
    all_corrs,
    "results/correlations/heatmaps_v0/correlations_all_tissues.csv"
  )
  
  write.csv(
    all_pvals,
    "results/correlations/heatmaps_v0/pvalues_all_tissues.csv"
  )
}

# Main function
main <- function() {
  # Load data
  cat("Loading RDS data...\n")
  global_meta <<- readRDS("data/baria_metadata.RDS")
  global_expr <<- readRDS("data/RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS")
  global_ensembls <<- rownames(global_expr)
  
  # Gene list
  cat("Reading gene list...\n")
  gene_df <- readRDS("data/gene_list.RDS") %>% slice(1:45)
  gene_df$ensembl_id <- gene_df$ENSEMBL
  gene_df$symbol <- gsub("\\s*\\([^)]*\\)", "", gene_df$ICP_symbol)  # First remove all parenthetical text 
  gene_df$symbol <- gsub("[,;].*", "", gene_df$symbol)  # Then remove everything after comma or semicolon
  cat("Found", nrow(gene_df), "genes in gene list\n") # Display gene list info
  
  # Process each gene
  cat("Processing genes...\n")
  all_results <- list()
  valid_count <- 0
  missing_genes <- c("ENSG00000273936", "ENSG00000276977", "ENSG00000236315") # to skip
  
  for (i in 1:nrow(gene_df)) {
    ensembl_id <- gene_df$ensembl_id[i]
    gene_symbol <- gene_df$symbol[i]
    
    cat(sprintf("Processing gene %d/%d: %s (%s)\n", 
                i, nrow(gene_df), 
                ifelse(is.null(gene_symbol) || gene_symbol == "", "Unknown", gene_symbol), 
                ensembl_id))
    
    # Skip genes that are known to be missing
    if (ensembl_id %in% missing_genes) {
      cat("Skipping known missing gene:", ensembl_id, "\n")
      next
    }
    
    # Extract gene data and calculate correlations
    gene_data <- extract_gene_data(ensembl_id, gene_symbol)
    
    if (!is.null(gene_data)) {
      results <- calculate_correlations(gene_data)
      
      # Only add valid results
      if (!is.null(results)) {
        valid_count <- valid_count + 1
        all_results[[valid_count]] <- results
      }
    }
  }
  
  # Check if we have any valid results before creating heatmaps
  if (length(all_results) > 0) {
    cat("Creating heatmaps with", length(all_results), "genes...\n")
    create_complex_heatmaps(all_results)
    cat("Heatmaps saved in results/correlations/heatmaps_v0/\n")
  } else {
    cat("No valid gene results to create heatmaps.\n")
  }
  
  cat("Analysis complete.\n")
}

# Run the main function
main()
