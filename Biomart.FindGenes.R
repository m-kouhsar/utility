find.gene.list <- function(cpg.list , manifest, distance, genome){
  
  
  if(genome=="hg19"){
    ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",version = "GRCh37")
  }else
    if(genome=="hg38"){
      ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    }else 
      print("Please specify genome version: hg19 or hg38")
  
  filters <- c("chromosome_name","start","end")
  
  find.gene <- function(x){
    chr=as.character(x[2])
    start=as.character(as.numeric(x[3])-distance)
    end=as.character(as.numeric(x[4])+distance)
    values <- list(chromosome=chr,start=start,end=end)
    tryCatch({
      all.genes <- biomaRt::getBM(attributes=c("external_gene_name","ensembl_gene_id"), filters=filters, values=values, mart=ensembl)
      if(length(all.genes$external_gene_name)>0){
        genes <- paste(all.genes$external_gene_name,collapse =  ";")
        ens.id <- paste(all.genes$ensembl_gene_id,collapse =  ";")
        return(list(genes=genes,ens.id=ens.id))
      }
      else{
        return(list(genes="",ens.id=""))
      }
    },
    error = function(e) {
      message("The following error has been occured during fetching the results of",x[1],"(chr",chr,":",start,"-",end,") from biomart")
      message(conditionMessage(e))
      return(list(genes="Error",ens.id="Error"))})
  }
  
  cpg <- data.frame(ID=cpg.list)
  cpg$chr = NA
  cpg$start = NA
  cpg$end = NA
  cpg$genes = NA
  cpg$ensembl = NA
  
  index = match(cpg$ID , manifest$IlmnID)
  cpg$chr = manifest$CHR[index]
  cpg$start = manifest$MAPINFO[index]
  cpg$end = manifest$MAPINFO[index]
  
  for (i in 1:nrow(cpg)) {
    results = find.gene(cpg[i,])
    cpg$genes[i] <- results$genes
    cpg$ensembl[i] <- results$ens.id
  }
  
  return(cpg)
}

