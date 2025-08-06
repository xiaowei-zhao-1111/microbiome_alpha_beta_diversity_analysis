# Alpha and Beta Diversity Analysis in Microbial Data

This repository contains R code and example data for performing alpha and beta diversity analysis on microbial relative abundance profiles derived from tools such as MetaPhlAn.

## Repository Structure

- `scripts/`  
  Contains the core R script for the diversity analysis pipeline.

- `data/`  
  Example input data:
  - `species_abundance.tsv`: Species-level relative abundance table.
  - `metadata.csv`: Sample metadata including variables such as age, sex, BMI, condition, and household ID.

- `results/`  
  Output directory containing:
  - `rank_plot_species_cutoff.pdf`: Rank plot to determine abundance threshold.
  - `alpha_diversity_results/`: CSV and plots for diversity metrics and statistical results.

## Analysis Overview

The pipeline includes:

1. **Preprocessing**
   - Cleaning species names and sample identifiers
   - Converting relative abundance values to proportions

2. **Filtering**
   - Applying abundance threshold to remove noise
   - Transforming data for diversity calculations

3. **Alpha Diversity**
   - Shannon, Simpson, Inverse Simpson, Richness
   - Statistical comparison (Wilcoxon, mixed-effects linear model)

4. **Beta Diversity**
   - Bray–Curtis dissimilarity
   - PCoA visualization
   - PERMANOVA for group differences

## Dependencies

Make sure you have the following R packages installed:

```r
install.packages(c("ggplot2", "reshape2", "dplyr", "ggpubr", "tidyr", "RColorBrewer"))
install.packages("vegan")
install.packages("ade4")
install.packages("lmerTest")
