CpG.Enrichment <- function(listA , listB ,listA.type = "illumina", listB.type="illumina", Background.Size = NA, 
                           Background.list = NA, use.missMethyl=F , arrayTypeA="EPIC", arrayTypeB = "EPIC"){
  
  #######################################################################################
  # listA and listB: A character vector of CpG IDs (illumina) or Gene Symbol
  # listA.type and listB.type: "illumina" if the corresponding list contains CpG IDs. "symbol" if corresponding list contains gene symbol
  # Background.list: A character vector of CpGs or NA. 
  # Background.Size: An intager value or NA. It will be ignored if use.annotation=T
  # use.missMethyl: If TRUE the gsameth function in missMethyl package will be used. In this case both listA and listB must contain CpG IDs
  
  #######################################################################################
  suppressMessages(library(missMethyl))
  listA <- unique(listA[!is.na(listA)])
  listB <- unique(listB[!is.na(listB)])
  
  if(use.missMethyl){
    
    arrayTypeA <- match.arg(arrayTypeA, choices = c("450K" , "EPIC"))
    arrayTypeB <- match.arg(arrayTypeB, choices = c("450K" , "EPIC"))
    arrayB <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeB)
    
    if(arrayTypeA != arrayTypeB){
      warning("Two different array types!")
      array1 <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeA)
      if(all(is.na(Background.list))){
        all_cpg <- intersect(array1$cpg , arrayB$cpg)
      }else{
        all_cpg <- Background.list
      }
    }else{
      if(all(is.na(Background.list))){
        all_cpg <- unique(arrayB$cpg)
      }else{
        all_cpg <- Background.list
      }
    }
    
    listB = unique(arrayB$entrezid[arrayB$cpg %in% listB])
    
    results <- suppressMessages(gsameth(sig.cpg = listA, collection = listB ,all.cpg = all_cpg, 
                       array.type = arrayTypeA, sig.genes=T))[,-4]
    results <- cbind.data.frame(length(listA) , results)
    names(results) <- c("N.listA","N.listB" , "N.Shared" , "P.value" , "Shared.Genes")
    
    }else{
      
      listA.type = match.arg(listA.type , choices = c("illumina" , "symbol"))
      listB.type = match.arg(listB.type , choices = c("illumina" , "symbol"))
      
      if(listA.type != listB.type){
        if(listA.type == "symbol"){
          message("Mapping listB to gene symbols...")
          arrayB = missMethyl:::.getFlatAnnotation(array.type = arrayTypeB)
          listB =  unique(arrayB$symbol[arrayB$cpg %in% listB])
        }
        if(listB.type == "symbol"){
          message("Mapping listA to gene symbols...")
          arrayA = missMethyl:::.getFlatAnnotation(array.type = arrayTypeA)
          listA =  unique(arrayA$symbol[arrayA$cpg %in% listA])
        }
      }
      
      if(all(is.na(Background.list))){
        if(is.na(Background.Size)){
          stop("In case of use.annotation=F , one of the Background.list or Background.Size must be specified.")
        }else{
          N <- Background.Size
        }
      }else{
        N <- length(unique(Background.list))
      }
      
      n <- length(listA[!is.na(listA)])
      m <- length(listB[!is.na(listB)])
      shared.CpGs <- intersect(listA , listB)
      shared.CpGs <- shared.CpGs[!is.na(shared.CpGs)]
      k <- length(shared.CpGs)
      
      con.table <- matrix(c(k, m-k, n-k, N-n-m+k),
                          nrow = 2,
                          dimnames = list(B=c("inB","notB"),
                                          A=c("inA","notA")))
      f.test <- fisher.test(con.table, alternative = "greater")
      results <- list(N.listA = n, N.listB=m , N.Shared=k , P.value = f.test$p.value , Shared.ID = paste(shared.CpGs,collapse = ";"))
    }
  
  return(results)
}
