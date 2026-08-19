# Nutrimmune Project

## Overview
This repository contains the analysis scripts for the Nutrimmune consortium project. The analyses focus on the BARIA cohort (bariatric surgery patients) and examine gene expression of immune checkpoint (ICP) genes across liver, subcutaneous fat and visceral fat, and how this expression relates to patient characteristics, metabolic changes and diet.

## Objectives
1. The BARIA cohort: patient characteristics
2. Gene expression of ICP in liver, subcutaneous fat and visceral fat
3. Correlations with baseline characteristics
4. Correlations with changes in BMI, CRP and HOMA-IR
5. Correlations with macronutrients (baseline)

## Project structure
```
data/           raw and processed input data (gitignored)
scripts/        analysis scripts, numbered by pipeline stage
scripts/archive/ older/exploratory scripts kept for reference
results/        script outputs (gitignored)
results_2025/   outputs from an earlier round of analyses (gitignored)
```

## Pipeline
Scripts are numbered by the stage of the analysis they belong to; run them in this order.

**1. Setup & descriptives**
- `1_datacleaning.R` — cleans and derives variables from the raw BARIA clinical and dietary data (medication use, dates, unit conversions, etc.). Run this first, as later scripts depend on its output.
- `1_table1.R` — generates the cohort characteristics table (objective 1).
- `1_heatmap_icp.R` — heatmap of median ICP gene expression across liver, subcutaneous fat and visceral fat (objective 2).
- `1_diet_descriptives.R` — descriptive statistics and correlation plots for dietary intake.

**2. Sample selection**
- `2_matching.R` — selects a sex-matched subgroup of BARIA participants for downstream analyses. Eligible participants are those with dietary data, a serum sample, and a follow-up (visit 4) BMI measurement, excluding recent antibiotic use/diarrhoea at baseline. Females are optimally matched to all eligible males (`MatchIt`) on age, lipid-lowering drug use and metformin use, so that sex comparisons aren't confounded by these variables. Outputs a balance table and the list of matched IDs to `results/matching/`.

**3. Correlations with baseline characteristics & metabolic change**
- `3_heatmap_v0.R` — heatmap of ICP gene expression correlated with baseline clinical characteristics (objective 3).
- `3_heatmap_v0v4.R` — same, but for the baseline-to-one-year change in BMI, CRP and HOMA-IR (objective 4).
- `3_correlations.R` — scatter plots for the significant correlations identified in the heatmap scripts above.

**4. Correlations with macronutrients**
- `4_heatmap_v0.R` — heatmap of ICP gene expression correlated with baseline macronutrient intake (objective 5).
- `4_correlations.R` — scatter plots for the significant macronutrient correlations.

**Archived** (`scripts/archive/`): earlier or exploratory scripts (data availability checks, RNA/tissue overlap checks, Venn diagrams, violin plots, additional correlation variants) kept for reference but not part of the current pipeline.

## How to run
Dependencies are managed with [pixi](https://pixi.sh); the environment is defined in `pixi.toml`.

```bash
pixi install                 # set up the R environment
pixi run datacleaning        # run scripts/1_datacleaning.R (run this first)
pixi run table1
pixi run heatmap-icp
pixi run diet-descriptives
pixi run matching
pixi run heatmap-v0-3
pixi run heatmap-v0v4-3
pixi run correlations-3
pixi run heatmap-v0-4
pixi run correlations-4
```
See the `[tasks]` section of `pixi.toml` for the full list, including archived scripts.

Input files are expected in `data/` (not tracked in git); outputs are written to `results/`.

## Contact
Barbara Verhaar — bjh.verhaar@gmail.com
