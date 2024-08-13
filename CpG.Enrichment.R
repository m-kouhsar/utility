CpG.Enrichment <- function(listA , listB ,Background.Size = NA, Background.list = NA, use.annotation=F , arrayTypeA="EPIC", arrayTypeB = "EPIC"){
  
  # listA: A character vector of CpGs
  # listB: A character vector of CpGs
  # Background.list: A character vector of CpGs or NA. 
  # Background.Size: An intager value or NA. It will be ignored if use.annotation=T
  listA <- unique(listA[!is.na(listA)])
  listB <- unique(listB[!is.na(listB)])
  
  if(use.annotation){
    suppressMessages(library(missMethyl))
    arrayTypeA <- match.arg(arrayTypeA, choices = c("450K" , "EPIC"))
    arrayTypeB <- match.arg(arrayTypeB, choices = c("450K" , "EPIC"))
    array2 <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeB)
    if(arrayTypeA != arrayTypeB){
      array1 <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeA)
      if(all(is.na(Background.list))){
        all_cpg <- intersect(array1$cpg , array2$cpg)
      }else{
        all_cpg <- Background.list
      }
    }else{
      if(all(is.na(Background.list))){
        all_cpg <- array2$cpg
      }else{
        all_cpg <- Background.list
      }
    }
    
    listB = unique(array2$entrezid[array2$cpg %in% listB])
    
    results <- suppressMessages(gsameth(sig.cpg = listA, collection = listB ,all.cpg = all_cpg, 
                       array.type = arrayTypeA, sig.genes=T))[,-4]
    names(results) <- c("N.listB" , "N.Shared" , "P.value" , "Shared.Genes")
    }else{
      if(all(is.na(Background.list))){
        if(is.na(Background.Size)){
          stop("In case of use.annotation=F , one of the Background.list or Background.Size must be specified.")
        }else{
          N <- Background.Size
        }
      }else{
        message("Background.list is specefied. Background.Size will be ignored.")
        N <- length(Background.list)
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
      results <- list(N.listB=m , N.Shared=k , P.value = f.test$p.value , Shared.CpGs = paste(shared.CpGs,collapse = ";"))
    }
  
  return(results)
}
