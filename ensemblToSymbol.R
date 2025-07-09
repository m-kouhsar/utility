ensembleToSymbol <- function(ensemblIDs){
  library(biomaRt)
  library(stringr)
  
  ensemblIDs <- str_split(ensemblIDs , pattern = "[.]" , simplify = T)[,1]
  
  mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  gene.symbols = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","entrezgene_id" ,
                                                                 "chromosome_name","start_position","end_position"),
                       values=ensemblIDs,mart= mart)
  
  return(gene.symbols)
}