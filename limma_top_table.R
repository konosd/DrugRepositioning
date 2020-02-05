limma.top.table <- function( name, dfactor, expression, contrast.vector){
  design.matrix = model.matrix(~0+dfactor)
  colnames(design.matrix) <- levels(dfactor)
  
  #making contrasts
  contrast.matrix = makeContrasts(contrasts = contrast.vector, 
                                  levels = design.matrix)
  fit = lmFit( object = expression[,-1], 
               design = design.matrix)
  fit = contrasts.fit( fit, 
                       contrast.matrix)
  fit = eBayes( fit ) 
  
  limma.top.list =lapply(c(1:length(contrast.vector)), 
                         FUN =  topTable, 
                         fit = fit, 
                         genelist = expression$SYMBOL, 
                         number = Inf, adjust.method = "BH",sort.by = "p")
  names(limma.top.list)  = paste(name, contrast.vector)
  return(limma.top.list)
}