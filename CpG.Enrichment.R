CpG.Enrichment <- function(listA , listB ,Background.Size = NA, Background.list = NA, use.annotation=F , arrayTypeA="EPIC", arrayTypeB = "EPIC"){
  
  # listA: A character vector of CpGs
  # listB: A character vector of CpGs
  # Background.list: A character vector of CpGs or NA. 
  # Background.Size: An intager value or NA. It will be ignored if use.annotation=T
  
  if(use.annotation){
    suppressMessages(library(missMethyl))
    arrayTypeA <- match.arg(arrayTypeA, choices = c("450K" , "EPIC"))
    arrayTypeB <- match.arg(arrayTypeB, choices = c("450K" , "EPIC"))
    
    array1 <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeA)
    array2 <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeB)
    
    
    listB = array2$entrezid[array2$cpg %in% listB]
    if(all(is.na(Background.list))){
      all_cpg <- intersect(array1$cpg , array2$cpg)
    }else{
      all_cpg <- Background.list
    }
    results <- gsameth(sig.cpg = listA, all.cpg = all_cpg, collection = listB , 
                       array.type = arrayTypeA, sig.genes=T)
    
    }else{
      if(all(is.na(Background.list))){
        if(is.na(Background.Size)){
          stop("In case of use.annotation=F , one of the Background.list or Background.Size must be specified.")
        }else{
          N <- Background.Size
        }
      }else{
        warning("Background.list is specefied. Background.Size will be ignored.")
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
      results <- list(N=m , DE=k , P.DE = f.test$p.value , Shared.CpGs = paste(intersect(CpGs.Region , Input.CpGs),collapse = ";"))
    }
  
  return(results)
}
