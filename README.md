# Alpha and Beta Diversity Analysis of Microbiome Data

This repository provides an R-based pipeline to analyze microbial community diversity using **species-level relative abundance data** and associated **sample metadata**. The workflow includes data preprocessing, alpha diversity analysis, beta diversity (PERMANOVA, PCoA), and statistical testing using both **mixed-effects models** and **non-parametric tests**.

---

## 📁 Repository Structure



---

## 🧬 Data Requirements

- **species_abundance.tsv**: Output from MetaPhlAn, containing species-level relative abundance for each sample (rows = taxa, columns = samples).
- **metadata.csv**: Sample metadata with rows as sample IDs and columns including:
  - `condition` (e.g., case/control)
  - `age`, `sex`, `BMI`
  - `household_id` (for stratified analysis)

---

## 🔧 Dependencies

Install required R packages:

```r
install.packages(c("lmerTest", "vegan", "ade4", "ggplot2", "reshape2", "dplyr", "ggpubr", "tidyr", "RColorBrewer"))
🧪 Analysis Steps
1. Data Preprocessing
Filters species-level clades (s__) from MetaPhlAn output
Converts relative abundances to proportions
Applies a rank-based cutoff for low-abundance filtering
Merges species table with metadata
2. Alpha Diversity
Calculates: Shannon, Simpson, Inverse Simpson, Richness
Performs:
Mixed-effects modeling (adjusting for household_id)
Wilcoxon rank-sum tests between conditions
3. Beta Diversity
Transforms data (arcsine square root)
Computes Bray–Curtis dissimilarity
Runs PERMANOVA:
Marginal model (adjusting for covariates)
Univariate models (condition, age, sex, BMI)
Performs Principal Coordinate Analysis (PCoA)
Visualizes results with group centroids and confidence ellipses
