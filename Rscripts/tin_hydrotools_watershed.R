### Delineate the watershed
watershed <- function(cr.arc, dfr){
  
  from.node <- paste(round(flowpaths$from.node.x,4),round(flowpaths$from.node.y,4),sep="_")
  to.node <- paste(round(flowpaths$to.node.x,4),round(flowpaths$to.node.y,4),sep="_")
  
  dffp <- data.frame(from.node, 
                     from.node.z=flowpaths$from.node.z, 
                     to.node, 
                     to.node.z=flowpaths$to.node.z,
                     stringsAsFactors=FALSE); head(dffp)

  # Find the corresponding triangles and colour them 
 
  src.arcs <- NULL
  arcs.tr <- NULL
  x=0
  while(!is.null(cr.arc)){x=x+1; print(x); 
    
    src.arc <- NULL
    tmp <- NULL
   
    for(i in 1:length(cr.arc)){#i=1
    
      up.arc <- which(dffp$from.node[cr.arc[i]]==dffp$to.node); up.arc
      
      if(length(up.arc)!=0){
        tmp <- c(tmp,as.numeric(up.arc)); 
      } else{
        {src.arc <- c(src.arc,cr.arc[i])}; 
      }
      
    }#for  

    arcs.tr <- c(arcs.tr, tmp)
    cr.arc <- tmp
   
    if(!is.null(src.arc)){src.arcs <- c(src.arcs, src.arc)}; print(src.arcs)
 
  }#while 
  
  tri.rel <- src.arcs #optionally: unique(src.arcs)
  plot(tin.tri[tri.rel,], add=TRUE, col='palegreen')
 
  return(tri.rel)
}  
  


