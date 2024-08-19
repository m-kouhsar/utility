
categorical.box.plot <- function(expression , phenotype, category.col, genes = NA, method, y.lab=NA, title="", facet.col = NA,category.labels){
  
  suppressMessages(library(ggplot2))
  if(!identical(rownames(phenotype) , colnames(expression))){
    stop("Row names of phenotype data should be matched with colnames of expression data")
  }
  if(length(genes) == 1){
    if(all(is.na(genes)){
      genes <- rownames(expression)
    }
  } 
  
  method <- match.arg(method , c("pca","mean"))
  expression <- expression[genes , ]
  index <- apply(expression, 1, function(x){return(!all(is.na(x)))})
  expression <- expression[index , ]
  phenotype$x_ <- factor(phenotype[,category.col])

  if(method == "pca"){
    pca <- prcomp(t(expression), center = T , scale. = T)
    phenotype$y_ <- pca$x[,1]
    
    if(is.na(y.lab))
      y.lab <- "PC1"
  }
  if(method == "mean"){
    phenotype$y_ <- colMeans(expression)
    if(is.na(y.lab))
      y.lab <- "Mean expression"
  }
  
  if(!is.na(facet.col)){
    phenotype$f_ <- phenotype[,facet.col]
  }
  p <- ggplot(data = phenotype,aes(x=x_,y=y_, fill=x_))+ 
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_boxplot() + 
    scale_fill_manual(values = c("#01BEC3","#F8756D"),labels=category.labels)+
    scale_x_discrete(labels=category.labels)+
    geom_jitter(color="black", size=0.8, alpha=0.9) +
    ylab(y.lab) + 
    xlab(category.col) + 
    ggtitle(title) + 
    guides(fill=guide_legend(title=category.col)) + 
    theme_bw()
   
   if(!is.na(facet.col)){
     p <- p + facet_grid(~f_)
   }
   
   return(p)
}

