# utility
## Biomart.FindGenes.R
It's an R function to find all genes annotated around a list of CpGs (CpG location +- a distance) using Biomart package. 
## CategoricalBoxPlot.R
The following example will describe the `categorical.box.plot` function. 

`
expr.data <- as.data.frame(matrix(data = rnorm(n = 800 , mean = 0 , sd = 0.6), nrow = 20))
names(expr.data) <- paste0("Sample" , c(1:40))
rownames(expr.data) <- paste0("Gene" , c(1:20))

phenotype <- as.data.frame(matrix(data = NA , nrow = 40 , ncol = 2))
rownames(phenotype) <- colnames(expr.data)
names(phenotype) <- c("category" , "facet")
phenotype$category[c(1:10 , 30:40)] <- "Case"
phenotype$category[11:29] <- "Control"
table(phenotype$category)
phenotype$facet[1:20] <- "Cohort1"
phenotype$facet[21:40] <- "Cohort2"

categorical.box.plot(expression = expr.data,phenotype = phenotype , category.col = "category" , facet.col = "facet" , method = "mean")
`

## CpG.Enrichment.R
It's an R function to run an enrichment analysis on two list of CpGs to find whether list A is overrepresented in list B. If you set the `use.annotation` parameter to `TRUE`, it will use `gsameth` function in `missMethyl` package in R to run the enrichment analysis in gene level.
