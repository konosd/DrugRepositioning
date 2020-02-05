finding.degs <- function( cluster) {
  set.seed(1)
  ngenes = replicate( 1000,  {
    
    
    index = c( which( pheno.global$Status %in% c("Control", "Immortalized")),
               sample( which(pheno.global$Status == "Melanoma" &
                               df.pca.corrected$Cluster == cluster), 
                       6, replace = FALSE))
    
    tmp = limma.top.table(name = paste("Melanomas",
                                       paste(as.character(tail(index,6)), 
                                             collapse = "-") ),
                          dfactor = factor(gsub("Immortalized", "Control", 
                                                pheno.global$Status[index])),
                          expression = global.exp.corrected[,c(1,index+1)] ,
                          contrast.vector = "Melanoma-Control")
    
    tmp = tmp[[1]] %>% filter(adj.P.Val < 0.05, abs(logFC) > 1) %>% #nrow()
      dplyr::select( ID, logFC, t, adj.P.Val) %>% 
      mutate(MelCells = paste(as.character(tail(index,6)), collapse = "-"))
    return(tmp)
  })
  ngenes <- apply(X = ngenes, MARGIN = 1, FUN = unlist) %>% 
    as.data.frame(., stringsAsFactors = FALSE) %>% 
    mutate(logFC = as.double(logFC), 
           t = as.double(t),
           adj.P.Val = as.double(adj.P.Val),
           Set = paste("From cluster",cluster)) 
  return(ngenes)  
}