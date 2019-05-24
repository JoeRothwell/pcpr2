# The Principal Component Partial R-squared method 
# (PC-PR2)


#### Background

The PC-PR2 is a statistical method, developed by Fages *et al*. (2014) (1, 2), used to investigate sources of variability in metabolomics or other omics data. It combines features of principal component and multivariable linear regression analyses. The input is a complete X-matrix of omics data and a corresponding set of descriptive Y-data (subject metadata). The output is the proportion of variation in the omics data attributed to each Y-variable, expressed as Rpartial2.

Test data, consisting of a sample of a transcriptomics dataset, is also included. This consists of five descriptive variables for the 124 subjects (two categorical, three numeric) and 3000 corresponding transcriptomics intensities.

#### Installation and usage

It's a good idea to first update all installed R packages.

````r
update.packages()
````

Make sure `devtools` is installed and loaded on your system. Now install `pcpr2` as follows:

````r
install_github("JoeRothwell/pcpr2")
````

You may be prompted to update other packages.

PC-PR2 is performed using the function `runPCPR2` which outputs partial R2 values for each covariate as a named vector. The variability in the omics data desired to be explained can be set with the argument `pct_threshold`, which defaults to 0.8.

A sample of transcriptomics data is provided as an example.

````r
library(pcpr2)
output <- runPCPR2(transcripts, Y_metadata)
output
           sex         height         weight smoking.status     age.sample             R2 
    1.24647643     2.48569520     0.10218837     2.94946793     0.03072886     4.91513509 
````
The function `plotProp` is a wrapper for `barplot` that quickly plots this output.

````r
par(mar=c(6,5,4,2)) # Adjust plot area margins to fit on x-axis labels
plotProp(output, main = "Variability in transcriptomics data explained by covariates")
````
<p align="center">
<img src="example_plot.png">
</p>

#### Reference

(1) Fages et al (2014) Investigating sources of variability in metabolomic data in the EPIC study: the Principal Component Partial R-square (PC-PR2) method. *Metabolomics* 10(6): 1074-1083, DOI: 10.1007/s11306-014-0647-9

(2) Perrier et al (2018) Identifying and correction epigenetics measurements for systematic sources of variation.
*Clin Epigenetics* 21(10): 38, DOI: 10.1186/s13148-018-0471-6
