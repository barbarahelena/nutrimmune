# Nutrimmune Project

## Overview
This repository contains all scripts I use for the Nutraimmune consortium project. My analyses in this project were focused on understanding the relationships between gene expression (liver, subcutaneous fat, visceral fat), clinical variables, and dietary factors in the BARIA cohort. The project involves generating heatmaps, correlation analyses, and visualizations to explore gene expression patterns and their associations with clinical and dietary data.

## Project Structure
The project is organized into the following scripts and directories:

Scripts
1. `heatmap_v0.R`:
    - Generates heatmaps for baseline gene expression data.
    - Includes correlation analyses with clinical variables.
2. `heatmap_icp.R`:
    - Focuses on ICP genes and creates horizontal heatmaps for median expression across tissues.
3. `heatmap_v0v4.R`:
    - Compares baseline and one-year post-op gene expression data.
    - Includes combined heatmaps for all tissues.
4. `correlations.R`:
    - Generates scatter plots for significant correlations.
    - Handles clinical variable transformations and delta correlations.
5. `datacleaning.R`:
    - Prepares and cleans clinical and dietary data.
    - Handles missing values, date formatting, and unit conversions.
6. `violinplots.R`:
    - Creates violin plots for visualizing changes in clinical variables (e.g., CRP, BMI, HOMA-IR).
7. `diet_descriptives.R`:
    - Provides descriptive statistics and visualizations for dietary data.
8. `venndiagrams.R`
    - Check overlap in correlations derived from heatmap scripts and draw Venn diagrams.

## How to Run
1. Install Dependencies:
    - Ensure the following R packages are installed:
    `tidyverse, ComplexHeatmap, circlize, ggsci, ggpubr, readxl, grid, tableone, corrplot, Hmisc`
2. Prepare Data:
    - Place input files in the data/ directory.
    - Ensure RNA-Seq data and clinical metadata are formatted correctly.
3. Run Scripts: execute the datacleaning script first, before running other scripts.

## Contact
For questions or issues, contact Barbara Verhaar at b.j.verhaar@amsterdamumc.nl.