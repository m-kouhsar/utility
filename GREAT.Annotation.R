great.annotation <- function(CpG,genome="hg19",distance=1000){
  suppressMessages(library(rGREAT))
  suppressMessages(library(GenomicRanges))
  
  if((!is.data.frame(CpG))|(!all(c("chr","start","end") %in% colnames(CpG))))
    stop("CpG must be a data frame which contains at least three columns (chr, start and end)")
  
  CpG.range = makeGRangesFromDataFrame(CpG,ignore.strand = T)
  job = submitGreatJob(gr = CpG.range , genome = genome, rule = "oneClosest")
  results = as.data.frame(getRegionGeneAssociations(job))
  results$ID = paste(results$seqnames,results$start,results$end,sep = ":")
  
  index <- match(results$ID , paste(CpG$chr,CpG$start,CpG$end,sep = ":"))
  
  CpG$GREAT.Annotation = NA
  CpG$dist.to.TSS = NA
  CpG$GREAT.Annotation[index] = as.character(results$annotated_genes)
  CpG$dist.to.TSS[index] = results$dist_to_TSS
  
  return(CpG)
}
