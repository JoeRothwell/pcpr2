# The Principal Component Partial R-squared method (PC-PR2)


#### Background

The PC-PR2 is a statistical method, developed by Fages *et al*. (2014) (1), used to investigate sources of variability in metabolomics or other omics data. It combines features of principal component analysis and multivariable linear regression analyses. The input is a complete X-matrix of omics data and a corresponding set of descriptive Y-data (subject metadata). The output is the proportion of variation in the omics data attributed to each Y-variable.

The original version of the PCPR-2 code is stored in this repository, as well as the latest script, which is in continuous development.

Test data, consisting of a sample of a transcriptomics dataset, is also included. This consists of five descriptive variables for the 124 subjects (two categorical, three numeric) and 3000 corresponding transcriptomics intensities.

Reference

(1) Fages et al (2014) Investigating sources of variability in metabolomic data in the EPIC study: the Principal Component Partial R-square (PC-PR2) method. *Metabolomics* 10(6): 1074-1083, DOI: 10.1007/s11306-014-0647-9

#### Instructions

PC-PR2 is performed using the function `runPCPR2` which outputs partial R2 values for each covariate as a named vector. The variability in the omics data desired to be explained can be set with the argument `pct_threshold`, which defaults to 0.8.

A subset of a transcriptomics dataset is provided as an example.

````r
output <- runPCPR2(transcripts, Y_metadata)
output
           sex         height         weight smoking.status     age.sample             R2 
    1.24647643     2.48569520     0.10218837     2.94946793     0.03072886     4.91513509 
````
The function `plotProp` is a wrapper for `barplot` that plots this output.

````r
par(mar=c(6,5,4,2)) # To fit on x-axis labels
plotProp(output)
title("Variability in transcriptomics data explained by covariates")
````
<p align="center">
<img src="example_plot.png">
</p>
