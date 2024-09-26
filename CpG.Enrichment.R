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
  arrayTypeA <- match.arg(arrayTypeA, choices = c("450K" , "EPIC"))
  arrayTypeB <- match.arg(arrayTypeB, choices = c("450K" , "EPIC"))
  
  
  if(use.missMethyl){
    if((listA.type != "illumina")|(listB.type != "illumina")){
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
      
      listA.type = match.arg(listA.type , choices = c("illumina" , "symbol"))
      listB.type = match.arg(listB.type , choices = c("illumina" , "symbol"))
      
      if(listA.type != listB.type){ # Different ID lists (CpG IDs and Gene Symbols): CpGs in one list needs to be converted to gene symbol
        if(listA.type == "symbol"){
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
        if(listB.type == "symbol"){
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
        
      }else{ #I: both lists have the same IDs
        
        if(!all(is.na(Background.list))){ #background list is specified 
          
          N <- length(unique(Background.list))
          
        }else{ #II: 
          if(!is.na(Background.Size)){ #background size is specified
            
            N <- Background.Size
            
          }else{ #III: No background list, No background size: We need to define it
            
            arrayA <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeA)
            
            if(listA.type == "illumina"){ # both lists are illumina IDs
              
              if(arrayTypeA != arrayTypeB){
                
                warning("Two different array types!")
                arrayB <- missMethyl:::.getFlatAnnotation(array.type = arrayTypeB)
                N <- length(unique(intersect(arrayA$cpg , arrayB$cpg)))
                
              }else{
                
                N <- length(unique(arrayA$cpg))
                
              }
            }else{ # both lists are gene symbols
              
              N <- length(unique(arrayA$symbol))
              
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
