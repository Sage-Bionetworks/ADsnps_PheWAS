# ADsnps_PheWAS
Code to perform PheWAS analysis in Heath et al., manuscript


The R code to run analyses in the manuscript entitled "Manifestations of genetic risk for Alzheimer Disease in the blood: a cross-sectional multi-mic analysis in healthy adults aged 18-90+" is posted here. You will need to install the package “PheWAS” and all dependencies to run the main analysis. See the manuscript for Arivale data acquisition details.

PheWAS_Rcode.r
Perform phenome-wide association study using R package "PheWAS" with Arivale data Running the R package “PheWAS”
Directions from Vanderbilt: https://www.vumc.org/cpm/cpm-blog/phewas-r-package
Github: https://github.com/PheWAS/PheWAS
Paper to be cited: https://www.nature.com/articles/nbt.2749
This package performs a phenome-wide association study, automatically performs linear or logistic regression based on input.
Input requires:
	1.	Genotype file, containing subject ID and genotype, coded numerically 0,1,2 for typical SNP analyses (but this value can be any predictor).
	2.	Covariate file, containing subject ID and any covariates (sex, age, BMI, comoborbidities, medication use, etc).
	3.	Phenotype file, containing subject ID and each outcome (ICD-9 code status or quantitative analyte).

PheWAS output for a single SNP includes (for each analyte): the FDR-adjusted p-value after rejoining the two dataframes (chemistries and proteins/metabolites, run separately to account for a difference in covariates), phenotype (analyte), snp, beta coefficient, standard error, p-value, N, Hardy-Weingberg Equilibrium p-value, allele frequency, and TRUE/FALSE flag for whether that analyte passed the FDR threshold in the initial analysis.

Gene by sex interaction analysis code: fit a linear regression model for each SNP with a sex-by-SNP interaction term, plus all covariates (chemistries run separately from proteins/metabolites). FDR-adjustment is performed on the interaction terms only. Results output as csv files for each SNP.

Code to generate box plot figures is also included.
