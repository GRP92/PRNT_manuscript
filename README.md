# PRNT_manuscript
R-Code for published paper: "Virological and Immunological Features of Sars-CoV-2 Infected Children Who Develop Neutralizing Antibodies"

Correlation analysis and Multi‐Omics Factor Analysis (MOFA) on FACS, virological and proteomic data of Sars-CoV-2 infected children.

Mofa.R: performs Multi‐Omics Factor Analysis (MOFA) on FACS and proteomic data in order to deconvolute the main sources of variation. MOFA is publicly accessible here: https://github.com/bioFAM/MOFA2

AUC_NP_swabs.R: calculate the area under the curve (AUC) for each patient. The curves were drawn plotting on the y-axis the viral copies and on the x-axis the number of days between the date of swab collection and the date of onset of symptoms for symptomatic patients and the hospital admission date for non-symptomatic ones.  

Correlations.R performs spearman correlation analysis between FACS, proteomic and virological data. The results are plotted as correlation heatmap and single correlation plot.

Project is created with RStudio (version 3.6.2).
