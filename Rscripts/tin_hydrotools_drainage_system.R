# Definition of all flow sections (drainage system)

### Flow section from centroid to point of intersection (overland)
flowsec_cen_poi <- function(dfr=dftin){ 
  data.frame(from.node.x=dfr$cen_x, 
             from.node.y=dfr$cen_y, 
             from.node.z=dfr$cen_z,
             to.node.x=dfr$poi_x, 
             to.node.y=dfr$poi_y,
             to.node.z=dfr$poi_z,
             section=1)
}  


### Flow section from point of intersection to point of intersection (overland)
flowsec_poi_poi <- function(dfr=dftin){
  ssi <- which(dfr$flowtype_poi=='overland');ssi
  sdfr <- dfr[ssi,];sdfr
    
  ov <- data.frame(from.node.x=sdfr$poi_x,
                   from.node.y=sdfr$poi_y,
                   from.node.z=sdfr$poi_z,
                   to.node.x=NA,
                   to.node.y=NA,
                   to.node.z=NA)    

  for(i in 1:nrow(ov)){#i=1
    ov$to.node.x[i] <- dfr$poi_x[which(dfr$TID==sdfr$nb_poi[i])]
    ov$to.node.y[i] <- dfr$poi_y[which(dfr$TID==sdfr$nb_poi[i])]
    ov$to.node.z[i] <- dfr$poi_z[which(dfr$TID==sdfr$nb_poi[i])]
    print(i)
  }
  
  ov$section <- 2
  
  return(ov)
}


### Flow section from point of intersection to vertex (channel)
flowsec_poi_ver <- function(dfr=dftin){
  
  ssi <- which(dfr$flowtype_poi=='channel');ssi
  sdfr <- dfr[ssi,];sdfr
  
  cha <- data.frame(from.node.x=sdfr$poi_x,
                    from.node.y=sdfr$poi_y,
                    from.node.z=sdfr$poi_z,
                    to.node.x=NA,
                    to.node.y=NA,
                    to.node.z=NA)   
  
  for(i in 1:nrow(cha)){#i=1
    cha$to.node.x[i] <- sdfr[i,c(sdfr$vertex[i]*3)-1]
    cha$to.node.y[i] <- sdfr[i,c(sdfr$vertex[i]*3)]
    cha$to.node.z[i] <- sdfr[i,c(sdfr$vertex[i]*3)+1]
  }
  
  cha$section <- 3
  
  section <- paste(round(cha$from.node.x,4),round(cha$from.node.y,4),round(cha$to.node.x,4),round(cha$to.node.y,4),sep="_")
  dupl <- which(duplicated(section))
  if(!length(dupl)==0){
    cha <- cha[-dupl,]
  }
  
  return(cha)
}


### Flow section from vertex to point of intersection (overland)
flowsec_ver_poi <- function(dfr=dftin){
  
  ssi <- which(dfr$flowtype_vertex=='overland');ssi
  sdfr <- dfr[ssi,];sdfr
  
  ver.poi <- data.frame(from.node.x=rep(NA,nrow(sdfr)),
                        from.node.y=rep(NA,nrow(sdfr)),
                        from.node.z=rep(NA,nrow(sdfr)),
                        to.node.x=rep(NA,nrow(sdfr)),
                        to.node.y=rep(NA,nrow(sdfr)),
                        to.node.z=rep(NA,nrow(sdfr)))   
  
  for(i in 1:nrow(sdfr)){#i=1
    ver.poi$from.node.x[i] <- sdfr[i,c(sdfr$vertex[i]*3)-1]
    ver.poi$from.node.y[i] <- sdfr[i,c(sdfr$vertex[i]*3)]
    ver.poi$from.node.z[i] <- sdfr[i,c(sdfr$vertex[i]*3)+1]
    ver.poi$to.node.x[i] <- dfr$poi_x[sdfr$nb_vertex[i]]
    ver.poi$to.node.y[i] <- dfr$poi_y[sdfr$nb_vertex[i]]
    ver.poi$to.node.z[i] <- dfr$poi_z[sdfr$nb_vertex[i]]
  }
  
  ver.poi$section <- 4
  
  ver.poi <- unique(ver.poi)
  return(ver.poi)
}


### Flow section from vertex to vertex (channel)
flowsec_ver_ver <- function(dfr=dftin){
  
  ssi <- c(which(dfr$flowtype_vertex=='channel'));ssi
  sdfr <- dfr[ssi,];sdfr
  
  ver.ver <- data.frame(from.node.x=rep(NA,nrow(sdfr)),
                        from.node.y=rep(NA,nrow(sdfr)),
                        from.node.z=rep(NA,nrow(sdfr)),
                        to.node.x=rep(NA,nrow(sdfr)),
                        to.node.y=rep(NA,nrow(sdfr)),
                        to.node.z=rep(NA,nrow(sdfr)))   
  
  for(i in 1:nrow(sdfr)){#i=1
    ver.ver$from.node.x[i] <- sdfr[i,c(sdfr$vertex[i]*3)-1]
    ver.ver$from.node.y[i] <- sdfr[i,c(sdfr$vertex[i]*3)]
    ver.ver$from.node.z[i] <- sdfr[i,c(sdfr$vertex[i]*3)+1]
    ver.ver$to.node.x[i] <- dfr[sdfr$nb_vertex[i],c(sdfr$vertex_2[i]*3)-1]
    ver.ver$to.node.y[i] <- dfr[sdfr$nb_vertex[i],c(sdfr$vertex_2[i]*3)]
    ver.ver$to.node.z[i] <- dfr[sdfr$nb_vertex[i],c(sdfr$vertex_2[i]*3)+1]
  }
  
  ver.ver$section <- 5
  
  ver.ver <- unique(ver.ver)
  return(ver.ver)
}


### Determine missing flow sections
missing_flowsections <- function(dfr=dftin, flowpaths, channel){
  
  pt.a<-paste(dfr$ax,dfr$ay,round(dfr$az,0),sep='_');pt.a
  pt.b<-paste(dfr$bx,dfr$by,round(dfr$bz,0),sep='_');pt.b
  pt.c<-paste(dfr$cx,dfr$cy,round(dfr$cz,0),sep='_');pt.c
  
  ssi <- c(which(dfr$flowtype_poi=='channel'),which(dfr$flowtype_poi=='no.flowtype'));ssi
  sdfr <- dfr[ssi,];sdfr 
  
  vertex.x <- NULL
  vertex.y <- NULL
  for(i in 1:nrow(sdfr)){#i=1
    
    vertex.x[i] <- sdfr[i,c(sdfr$vertex[i]*3-1)] 
    vertex.y[i] <- sdfr[i,c(sdfr$vertex[i]*3)] 
    
    vertices <- paste(vertex.x,vertex.y,sep='_')
    vertices <- unique(vertices);vertices
  }
  
  ssi2 <- which(dfr$flowtype_vertex=='channel');ssi2
  sdfr2 <- dfr[ssi2,];sdfr2
  
  tmp <- paste(sdfr2$nb_vertex,sdfr2$vertex_2,sep='_')
  dupl <- which(duplicated(tmp))
  ssi2 <- ssi2[-dupl]
  sdfr2 <- sdfr2[-dupl,];sdfr2
  
  for(i in 1:nrow(sdfr2)){#i=1
    
    ver.cr <- paste(dfr[sdfr2$nb_vertex[i],c(sdfr2$vertex_2[i]*3-1)],
                    dfr[sdfr2$nb_vertex[i],c(sdfr2$vertex_2[i]*3)],sep='_'); ver.cr
    
    first.run <- TRUE
    overland <- FALSE
    while(all(vertices!=ver.cr) && overland==FALSE){
      
      # Define the vertex
      if(first.run==TRUE){
        pt.rel <- sdfr2$vertex_2[i];pt.rel
        ver.x <- dfr[sdfr2$nb_vertex[i],c(pt.rel*3)-1]; ver.x
        ver.y <- dfr[sdfr2$nb_vertex[i],c(pt.rel*3)]; ver.y
        ver.z <- dfr[sdfr2$nb_vertex[i],c(pt.rel*3)+1]; ver.z
        first.run <- FALSE
      }
      vertex <- paste(ver.x,ver.y,round(ver.z,0),sep='_'); vertex
      
      
      # Find adjacent triangles
      nbpts <- c(which(pt.a==vertex),which(pt.b==vertex),which(pt.c==vertex));nbpts
      TID.nb <- dfr$TID[nbpts[which(nbpts!=ssi2[i])]];TID.nb
      
      ai <- NULL
      bi <- NULL
      ci <- NULL
      if(any(dfr$az[nbpts]<ver.z)){
        ai <- nbpts[which(dfr$az[nbpts]<ver.z)];ai
      }
      
      if(any(dfr$bz[nbpts]<ver.z)){
        bi <- nbpts[which(dfr$bz[nbpts]<ver.z)];bi
      }
      
      if(any(dfr$cz[nbpts]<ver.z)){
        ci <- nbpts[which(dfr$cz[nbpts]<ver.z)];ci
      }
      
      nb.rel <- unique(c(ai,bi,ci));nb.rel
      
      if(is.null(nb.rel) || length(nb.rel)==1){
        overland <- TRUE
        print(i)
        print("sink/end of TIN")
      } else{
        
        # Order relevant neighbour triangles in clockwise direction
        nb.rel.cw <- nb.rel[circ_sort(rpt2d = data.frame(x=ver.x,y=ver.y), 
                                      pts2d = data.frame(x=dfr$cen_x[nb.rel],y=dfr$cen_y[nb.rel]),
                                      clockwise=TRUE)];nb.rel.cw
        
        dfv <- data.frame(flowtype=rep(NA,length(nb.rel.cw)),
                          nb.ver=rep(NA,length(nb.rel.cw)),
                          ver.rel=rep(NA,length(nb.rel.cw)));dfv
        
        for(j in 1:length(nb.rel.cw)){#j=1
          
          # Determine relevant points and their heights  
          pts.rel <- which(c(pt.a[nb.rel.cw[j]],pt.b[nb.rel.cw[j]],pt.c[nb.rel.cw[j]])!=vertex);pts.rel
          if(pts.rel[1]==1 && pts.rel[2]==3){pts.rel <- c(pts.rel[2],pts.rel[1])}
          z.pt <- as.numeric(dfr[nb.rel.cw[j],c(pts.rel*3)+1]);z.pt
          
          v.verpt.x <- NA
          v.verpt.y <- NA
          v.ptver.x <- NA
          v.ptver.y <- NA
          v.verpt2.x <- NA
          v.verpt2.y <- NA
          
          # Determine vectors left edge
          if(z.pt[1]<ver.z){
            v.verpt.x <- dfr[nb.rel.cw[j],c(pts.rel[1]*3)-1]-ver.x;
            v.verpt.y <- dfr[nb.rel.cw[j],c(pts.rel[1]*3)]-ver.y;
            v.ptver.x <- ver.x-dfr[nb.rel.cw[j],c(pts.rel[1]*3)-1];
            v.ptver.y <- ver.y-dfr[nb.rel.cw[j],c(pts.rel[1]*3)];
          }  
          
          # Determine gradient vectors
          if(j==1){
            vg.le.x <- dfr$grdv_x[nb.rel.cw[length(nb.rel.cw)]]
            vg.le.y <- dfr$grdv_y[nb.rel.cw[length(nb.rel.cw)]]
          } else{
            vg.le.x <- dfr$grdv_x[nb.rel.cw[j-1]]
            vg.le.y <- dfr$grdv_y[nb.rel.cw[j-1]]
          }
          vg.cr.x <- dfr$grdv_x[nb.rel.cw[j]]
          vg.cr.y <- dfr$grdv_y[nb.rel.cw[j]]
          
          # Calculate vector products to the test the left edge
          matrixvp1 <- matrix(c(vg.le.x,vg.le.y,v.verpt.x,v.verpt.y),nrow=2,byrow=TRUE); matrixvp1
          det.le <- det(matrixvp1);det.le
          matrixvp2 <- matrix(c(vg.cr.x,vg.cr.y,v.ptver.x,v.ptver.y),nrow=2,byrow=TRUE); matrixvp2
          det.cr <- det(matrixvp2);det.cr
          
          # Determine vectors edges
          v.verpt2.x <- dfr[nb.rel.cw[j],c(pts.rel[2]*3)-1]-ver.x
          v.verpt2.y <- dfr[nb.rel.cw[j],c(pts.rel[2]*3)]-ver.y
          
          if(z.pt[1]>=ver.z){
            v.verpt.x <- dfr[nb.rel.cw[j],c(pts.rel[1]*3)-1]-ver.x;
            v.verpt.y <- dfr[nb.rel.cw[j],c(pts.rel[1]*3)]-ver.y;
          }
          
          # Calculate vector products to test the triangle
          matrixvp3 <- matrix(c(vg.cr.x,vg.cr.y,v.verpt2.x,v.verpt2.y),nrow=2,byrow=TRUE); matrixvp3
          det.eri <- det(matrixvp3);det.eri
          matrixvp4 <- matrix(c(vg.cr.x,vg.cr.y,v.verpt.x,v.verpt.y),nrow=2,byrow=TRUE); matrixvp4
          det.ecr <- det(matrixvp4);det.ecr
          
          # Determine flow type
          if(!is.nan(det.le) && det.le>0 && det.cr>0){
            flowtype  <- "channel"} else if(det.eri<=0 && det.ecr>=0){
              flowtype <- "overland"} else{flowtype <- NA}                
          
          dfv$flowtype[j] <- flowtype
          if(!is.na(flowtype) && flowtype=="overland"){
            dfv$nb.ver[j] <- nb.rel.cw[j]}
          if(!is.na(flowtype) && flowtype=="channel"){
            dfv$nb.ver[j] <- nb.rel.cw[j]
            dfv$ver.rel[j] <- pts.rel[1]
          }
          
        }#j 
        
        ft.rel <- which(!is.na(dfv$flowtype)==TRUE); ft.rel 
        print(i)     
        print(dfv$flowtype[ft.rel])     
        
        if(!length(ft.rel)){
          overland <- TRUE
          print(i)
          print("sink/end of TIN")
        }
        if(length(ft.rel)==1){
          if(dfv$flowtype[ft.rel]=="overland"){
            from.node.x <- ver.x
            from.node.y <- ver.y
            from.node.z <- ver.z
            to.node.x <- dfr$poi_x[dfv$nb.ver[ft.rel]]
            to.node.y <- dfr$poi_y[dfv$nb.ver[ft.rel]]
            to.node.z <- dfr$poi_z[dfv$nb.ver[ft.rel]]
            flowpaths <- rbind(flowpaths, c(from.node.x, from.node.y, from.node.z, to.node.x, to.node.y, to.node.z,6))
            
            overland <- TRUE
          }
          if(dfv$flowtype[ft.rel]=="channel"){
            from.node.x <- ver.x
            from.node.y <- ver.y
            from.node.z <- ver.z
            ver.x <- dfr[dfv$nb.ver[ft.rel],c(dfv$ver.rel[ft.rel]*3)-1]
            ver.y <- dfr[dfv$nb.ver[ft.rel],c(dfv$ver.rel[ft.rel]*3)]
            ver.z <- dfr[dfv$nb.ver[ft.rel],c(dfv$ver.rel[ft.rel]*3)+1]
            channel <- rbind(channel, c(from.node.x, from.node.y, from.node.z, ver.x, ver.y, ver.z,7))
            
            ver.cr <- paste(ver.x,ver.y,sep='_'); ver.cr
          }
        }
        if(length(ft.rel)>1){
          pts.z <- NULL
          for (k in 1:length(ft.rel)){#k=1
            if(dfv$flowtype[ft.rel[k]]=="overland"){pt.z <- dfr$poi_z[dfv$nb.ver[ft.rel[k]]]} else{
              pt.z <- dfr[dfv$nb.ver[ft.rel[k]], dfv$ver.rel[ft.rel[k]]*3+1]}
            pts.z[k] <- pt.z
          }#k
          ft.rel <- ft.rel[which.min(pts.z)]
          
          if(dfv$flowtype[ft.rel]=="overland"){
            from.node.x <- ver.x
            from.node.y <- ver.y
            from.node.z <- ver.z
            to.node.x <- dfr$poi_x[dfv$nb.ver[ft.rel]]
            to.node.y <- dfr$poi_y[dfv$nb.ver[ft.rel]]
            to.node.z <- dfr$poi_z[dfv$nb.ver[ft.rel]]
            flowpaths <- rbind(flowpaths, c(from.node.x, from.node.y, from.node.z, to.node.x, to.node.y, to.node.z,6))
            
            overland <- TRUE
          }
          if(dfv$flowtype[ft.rel]=="channel"){
            from.node.x <- ver.x
            from.node.y <- ver.y
            from.node.z <- ver.z
            ver.x <- dfr[dfv$nb.ver[ft.rel],c(dfv$ver.rel[ft.rel]*3)-1]
            ver.y <- dfr[dfv$nb.ver[ft.rel],c(dfv$ver.rel[ft.rel]*3)]
            ver.z <- dfr[dfv$nb.ver[ft.rel],c(dfv$ver.rel[ft.rel]*3)+1]
            channel <- rbind(channel, c(from.node.x, from.node.y, from.node.z, ver.x, ver.y, ver.z,7))
            
            ver.cr <- paste(ver.x,ver.y,sep='_'); ver.cr
          }
        }
        
      }#if
    }#while
  }#i 
  
  return(list(flowpaths, channel))
}


### Subdivison of mulitple channel sections 
arc_split <- function(channel){
  
  # Calculation of the slope and the length of the sections
  v.fp.x <- channel$from.node.x-channel$to.node.x
  v.fp.y <- channel$from.node.y-channel$to.node.y
  
  slope <- round(v.fp.y/v.fp.x,4)
  
  dist <- round(sqrt(v.fp.x^2+v.fp.y^2),4)
  
  to.node <- paste(round(channel$to.node.x,4),round(channel$to.node.y,4),sep="_")
  
  dftn <- data.frame(to.node, slope, dist); dftn
  
  dupl <- which(duplicated(dftn$to.node))
  
  dftn.uni <- dftn[-dupl,]; dftn.uni
  
  # Find identical to nodes with same slopes
  for(i in 1:nrow(dftn.uni)){#i=1
    
    tns.sim <- which(dftn$to.node==dftn.uni$to.node[i]); tns.sim
    if(length(tns.sim)>1){
      
      diff.slopes <- unique(dftn$slope[tns.sim]); diff.slopes
      
      for(j in 1:length(diff.slopes)){#j=1
      
        slct <- which(dftn$slope[tns.sim]==diff.slopes[j])
        
        if(length(slct)>1){
          tns.rel <- tns.sim[slct]
          dupl <- which(duplicated(dftn$dist[tns.rel]))
          if(!length(dupl)==0){
            tns.rel <- tns.rel[-dupl]
          }
          tns.ord <- tns.rel[order(dftn$dist[tns.rel])]
          print(i)
          print(tns.ord)
          
          for(k in 1:length(tns.ord)){#k=1
            
            if(k>1){
              channel$to.node.x[tns.ord[k]] <- channel$from.node.x[tns.ord[k-1]] 
              channel$to.node.y[tns.ord[k]] <- channel$from.node.y[tns.ord[k-1]]
              channel$to.node.z[tns.ord[k]] <- channel$from.node.z[tns.ord[k-1]]
            }
          }#k
        }#if
      }#j 
      
    }#if
  }#i  

  flowpaths.cha <- channel
  return(flowpaths.cha)
}#fun


### Plot drainage system 
flowpaths_lines <- function(fp=flowpaths){
  
  #points(x=dftin$cen_x, y=dftin$cen_y, col="skyblue", add=TRUE)
  
  fls.coords <- NULL
  fls.lines <- vector("list", nrow(fp))
  for(i in 1:nrow(fp)){#i=1
    
    fls.coords <- Line(cbind(c(fp$from.node.x[i],fp$to.node.x[i]),c(fp$from.node.y[i],fp$to.node.y[i])))
   
    fls.lines[i] <- Lines(list(fls.coords), ID=i)
  }
    
  fls.SL <- SpatialLines(fls.lines)
  fls.SL$section <- fp$section
  
  #shapefile(x=fls.SL, filename='flowsections_R.shp')
  return(fls.SL)
}



