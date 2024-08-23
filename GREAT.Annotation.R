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

#################################################################
setwd("Brain EWAS Paper/WGCNA.Paper/Revision.June2024/")
cpg.list <- read.csv("darkgreen.CpGs.Annot.csv")

cpg.data <- data.frame(CpG.ID=cpg.list$CpG.ID,
                       chr=paste0("chr",cpg.list$Chr),
                       start=as.numeric(cpg.list$Genomic.position),
                       end=as.numeric(cpg.list$Genomic.position))
cpg.annot = great.annotation(CpG = cpg.data , genome = "hg19",distance = 1000)

identical(cpg.list$CpG.ID , cpg.annot$CpG.ID)
cpg.list$GREAT.annotation = cpg.annot$GREAT.Annotation
write.csv(cpg.list , file = "darkgreen.CpGs.Annot1.csv", row.names = F)
