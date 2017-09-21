### Calculation of the relative contributing area
rel_contr_area <- function(dfr=dftin.ws){
  
  from.node <- paste(round(flowpaths$from.node.x,4),round(flowpaths$from.node.y,4),sep="_")
  to.node <- paste(round(flowpaths$to.node.x,4),round(flowpaths$to.node.y,4),sep="_")
  
  dffp <- data.frame(from.node, 
                     from.node.z=flowpaths$from.node.z, 
                     to.node, 
                     to.node.z=flowpaths$to.node.z,
                     stringsAsFactors=FALSE); head(dffp)
  
  pts.gwn.xy <- paste(round(coordinates(pts.gwn.uni)[,1],4),round(coordinates(pts.gwn.uni)[,2],4),sep="_")
  
  cen.xy <- paste(round(dfr$cen_x,4),round(dfr$cen_y,4),sep="_")
  poi.xy <- paste(round(dfr$poi_x,4),round(dfr$poi_y,4),sep="_")
  
  for(i in 1:nrow(dfr)){#i=1
   
    print(i)
    cr.tri <- i
    cr.arc <- which(cen.xy[i]==dffp$from.node)
    cd.values <- NULL
    while(!length(cr.arc)==0 && all(pts.gwn.xy!=dffp$to.node[cr.arc])){
      
      # Length of the section
      v.fp <- c(flowpaths$to.node.x[cr.arc]-flowpaths$from.node.x[cr.arc],
                flowpaths$to.node.y[cr.arc]-flowpaths$from.node.y[cr.arc],
                flowpaths$to.node.z[cr.arc]-flowpaths$from.node.z[cr.arc])
      dist <- sqrt(v.fp[1]^2+v.fp[2]^2) 
      
      # Contribution coefficient of the section
      if(flowpaths$section[cr.arc]==3 || flowpaths$section[cr.arc]==5 || flowpaths$section[cr.arc]==7){
        
        tan.alpha <- atan(abs(v.fp[3])/dist)
        slope <- tan.alpha*180/pi
        k <- sqrt(90/slope)
      
      } else{ 
        k <- sqrt(90/dfr$slope[cr.tri])
      }   
        
      # Contribution distance value
      cd <- k*dist
      
      if(is.null(cd)){cd.values <- cd} else{cd.values <- c(cd.values,cd)}; cd.values
      
      cr.arc <- which(dffp$to.node[cr.arc]==dffp$from.node); cr.arc
      
      if(!length(cr.arc)==0){
        if(flowpaths$section[cr.arc]==2 || flowpaths$section[cr.arc]==4 || flowpaths$section[cr.arc]==6){
          cr.tri <- which(poi.xy==dffp$to.node[cr.arc])
        }
      }
       
    }#while
  
    if(!length(cr.arc)==0){
     
      # Length of the section
      v.fp <- c(flowpaths$to.node.x[cr.arc]-flowpaths$from.node.x[cr.arc],
                flowpaths$to.node.y[cr.arc]-flowpaths$from.node.y[cr.arc],
                flowpaths$to.node.z[cr.arc]-flowpaths$from.node.z[cr.arc])
      dist <- sqrt(v.fp[1]^2+v.fp[2]^2) 
      
      # Contribution coefficient of the section
      if(flowpaths$section[cr.arc]==3 || flowpaths$section[cr.arc]==5 || flowpaths$section[cr.arc]==7){
  
        tan.alpha <- atan(abs(v.fp[3])/dist)
        slope <- tan.alpha*180/pi
        k <- sqrt(90/slope)
        
      } else{ 
        k <- sqrt(90/dfr$slope[cr.tri])
      }   
      
      # Contribution distance value
      cd <- k*dist
      
      if(is.null(cd)){cd.values <- cd} else{cd.values <- c(cd.values,cd)}; cd.values
      
      g <- 1/sum(cd.values)*1000; g
    } else{g <- NA}
    
    print(g)
    dfr$rel_contrar[i] <- g
      
  }#for
  
  return(dfr)
}