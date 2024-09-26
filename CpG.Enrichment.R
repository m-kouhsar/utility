CpG.Enrichment <- function(listA , listB , Background.Size = NA, 
                           Background.list = NA, use.missMethyl=F , arrayTypeA="EPIC", arrayTypeB = "EPIC"){
  
  #######################################################################################
  # listA and listB: A character vector of CpG IDs (illumina) or Gene Symbol
  
  # arrayTypeA and arrayTypeB: 
  #                            "450K" if the corresponding list contains CpG IDs from 450K illumina array. 
  #                            "EPIC" if corresponding list contains CpG IDs from EPIC array.
  #                            "symbol" if corresponding list contains gene symbols.
  
  # Background.list: A character vector of CpGs/gene symbols or NA. 
  
  # Background.Size: An intager value or NA. It will be ignored if use.annotation=T
  
  # use.missMethyl: If TRUE the gsameth function in missMethyl package will be used. 
  #                 In this case both listA and listB must contain CpG IDs
  
  #######################################################################################
  suppressMessages(library(missMethyl))
  listA <- unique(listA[!is.na(listA)])
  listB <- unique(listB[!is.na(listB)])
  arrayTypeA <- match.arg(arrayTypeA, choices = c("450K" , "EPIC", "symbol"))
  arrayTypeB <- match.arg(arrayTypeB, choices = c("450K" , "EPIC", "symbol"))
  
  
  if(use.missMethyl){
    if((arrayTypeA == "symbol")|(arrayTypeB == "symbol")){
      stop("Only illumina IDs allowed if use.missMethyl=TRUE")
    }
    if(!is.na(Background.Size)){
      warning("If use.missMethyl=TRUE Background.Size will be ignored.")
    }
  
    arrayB <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeB)
    
    if(!all(is.na(Background.list))){ #background list is specified 
      
      all_cpg <- unique(Background.list)
      
    }else{ # No background list: We need to define it
      
      if(arrayTypeA != arrayTypeB){
        
        warning("Two different array types!")
        arrayA <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeA)
        all_cpg <- intersect(arrayA$cpg , arrayB$cpg)
      }else{
        all_cpg <- unique(arrayB$cpg)
      }
    }
    
    listB = unique(arrayB$entrezid[arrayB$cpg %in% listB])
    
    results <- suppressMessages(gsameth(sig.cpg = listA, collection = listB ,all.cpg = all_cpg, 
                       array.type = arrayTypeA, sig.genes=T))[,-4]
    results <- cbind.data.frame(length(listA) , results)
    names(results) <- c("N.listA","N.listB" , "N.Shared" , "P.value" , "Shared.Genes")
    
    }else{ # don't use missMethyl
      
      if(arrayTypeA != arrayTypeB){ # Dirrefent Array type
        if(arrayTypeA == "symbol"){
          message("Mapping listB to gene symbols...")
          arrayB = missMethyl:::.getFlatAnnotation(array.type = arrayTypeB)
          listB =  unique(arrayB$symbol[arrayB$cpg %in% listB])
          if(all(is.na(Background.list))){
            if(is.na(Background.Size))
              N <- length(unique(arrayB$symbol))
            else
              N <- Background.Size
          }else{
            N <- length(unique(Background.list))
          }
        }
        if(arrayTypeB == "symbol"){
          message("Mapping listA to gene symbols...")
          arrayA = missMethyl:::.getFlatAnnotation(array.type = arrayTypeA)
          listA =  unique(arrayA$symbol[arrayA$cpg %in% listA])
          if(all(is.na(Background.list))){
            if(is.na(Background.Size))
              N <- length(unique(arrayA$symbol))
            else
              N <- Background.Size
          }else{
            N <- length(unique(Background.list))
          }
        }
        
      }else{ #I: simillar array type
        
        if(!all(is.na(Background.list))){ #background list is specified 
          
          N <- length(unique(Background.list))
          
        }else{ #II: 
          if(!is.na(Background.Size)){ #background size is specified
            
            N <- Background.Size
            
          }else{ #III: No background list, No background size: We need to define it
            
            arrayA <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeA)
            
            if((arrayTypeA == "450K")|(arrayTypeA == "EPIC")){ # both lists are illumina IDs
              
              if(arrayTypeA != arrayTypeB){
                
                warning("Two different array types!")
                arrayB <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeB)
                N <- length(unique(intersect(arrayA$cpg , arrayB$cpg)))
                
              }else{
                
                N <- length(unique(arrayA$cpg))
                
              }
            }else{ # both lists are gene symbols
              
              stop("Both list are gene symbols. Background.list or Background.Size must be specified.")
              
            }
          }# End of III
        } # End of II
      } # End of I
      
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
