# MR-fatigue-to-strokeoutcome
Scripts
1. MR-fatigue to strokeoutcome.R
Description: This script is designed to assess the causal effect of fatigue on stroke outcomes using Mendelian Randomization. The analysis includes data pre-processing, instrumental variable selection, and the application of MR methods to evaluate the potential causal relationship.
Key Functions: Data pre-processing and quality control. Instrumental variable selection based on genome-wide association study (GWAS) data. MR analysis using various statistical techniques.
2. MR-mediation-lipid.R
Description: This script investigates the mediating role of lipid levels in the relationship between genetic variants and stroke outcomes. The mediation analysis helps to understand whether lipid levels act as a mediator in the pathway from genetic variants to stroke.
Key Functions: Mediation analysis framework implementation. Estimation of direct and indirect effects. Statistical testing of the mediation effects.
3. MR-Data pre-processing.R
Description: This script focuses on the pre-processing steps required for MR analysis. It includes data cleaning, harmonization of exposure and outcome datasets, and the preparation of genetic instruments.
Key Functions: Data cleaning and harmonization. Preparation of genetic instruments. Quality control measures to ensure robust analysis. 
4. MR-Sensitivity analysis.R
Description and key Functions: This script focuses on Cochran's Q-test and retention analysis to test the reliability and robustness of the evaluation results. Using MR Egger intercept method for genetic pleiotropy detection.
Data sources 
The data called by these scripts comes from a public database, so there is no need to upload data in this GitHub repository. You can download according to the data source link provided in the Methods section of the manuscript.
“The GWAS study was based on a UK Biobank sample which was 449,019 individuals in total (ID: ukb-b-929, https://gwas.mrcieu.ac.uk/datasets/ukb-b-929/ ) in Table 1.”
“The data can be accessed via the following website: https://cd.hugeamp.org/dinspector.html?dataset=GWAS_GISCOME_eu.”
Usage
To run these scripts, ensure you have the necessary R packages installed. Each script can be executed in an R environment, and it is recommended to follow the sequence of scripts for a comprehensive analysis starting from data pre-processing to causal inference and mediation analysis.
1. Installation of Required Packages
install.packages(c("MendelianRandomization", "TwoSampleMR", "dplyr", "tidyverse"))
2. Running the Scripts
# Pre-processing
source("MR-Data pre-processing.R")
# Fatigue to Stroke Outcome Analysis
source("MR-fatigue to strokeoutcome.R")
# MR-Sensitivity analysis
source("MR-Sensitivity analysis.R")
# Mediation Analysis
source("MR-mediation-lipid.R")
This repository aims to provide a structured approach to conducting MR analyses, ensuring reproducibility and clarity in investigating the genetic underpinnings of health outcomes.
