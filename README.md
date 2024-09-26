# utility
## GREAT.Annotation
Using `rGREAT` package to map a list of CpGs (a data frame contains at least three columns: "chr", "start", "end") to the clostest Transcription Start Site (TSS). For more details see [rGREAT](https://jokergoo.github.io/rGREAT/index.html). 
## Biomart.FindGenes.R
It's an R function to find all genes annotated around a list of CpGs (CpG location +- a distance) using Biomart package. 
## CategoricalBoxPlot.R
The following example describes the `categorical.box.plot` function in this script. 

```R
expr.data <- as.data.frame(matrix(data = rnorm(n = 800 , mean = 0 , sd = 0.6), nrow = 20))
names(expr.data) <- paste0("Sample" , c(1:40))
rownames(expr.data) <- paste0("Gene" , c(1:20))
head(expr.data)

phenotype <- as.data.frame(matrix(data = NA , nrow = 40 , ncol = 2))
rownames(phenotype) <- colnames(expr.data)
names(phenotype) <- c("category" , "facet")
phenotype$category[c(1:10 , 30:40)] <- "Case"
phenotype$category[11:29] <- "Control"
phenotype$facet[1:20] <- "Cohort1"
phenotype$facet[21:40] <- "Cohort2"
head(phenotype)

categorical.box.plot(expression = expr.data,phenotype = phenotype , category.col = "category" , facet.col = "facet" , method = "mean")
```

## CpG.Enrichment.R
It's an R function to run an enrichment analysis on two list of CpGs (or gene symbols) to find whether list A is overrepresented in list B using a Fisher's exact test. If you set the `use.annotation` parameter to `TRUE`, it will use `gsameth` function in `missMethyl` package in R to run the enrichment analysis in gene level. 
If one list contains CpG IDs and the other list gene symbol, it will use illumina annotation to map the CpG IDs in the second list and then run the Fisher's exact test. 
