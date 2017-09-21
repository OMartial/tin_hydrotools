### Calculation of the path of steepest ascent on a triangle 
steepest_ascent_poi <- function(dfr,crpt,i){
  
  # Variables of the line equation of the vector gradient
  pv.crpt <- c(crpt[1],crpt[2]); pv.crpt
  vg <- c(dfr$grdv_x[i],dfr$grdv_y[i])
  
  # Parameters of the equation of the plane
  A <- dfr$ay[i]*(dfr$bz[i]-dfr$cz[i])+dfr$by[i]*(dfr$cz[i]-dfr$az[i])+dfr$cy[i]*(dfr$az[i]-dfr$bz[i]);A
  B <- dfr$az[i]*(dfr$bx[i]-dfr$cx[i])+dfr$bz[i]*(dfr$cx[i]-dfr$ax[i])+dfr$cz[i]*(dfr$ax[i]-dfr$bx[i]);B
  C <- dfr$ax[i]*(dfr$by[i]-dfr$cy[i])+dfr$bx[i]*(dfr$cy[i]-dfr$ay[i])+dfr$cx[i]*(dfr$ay[i]-dfr$by[i]);C
  D <- -A*dfr$ax[i]-B*dfr$ay[i]-C*dfr$az[i];D
  
  # Determine the relevant edges
  # Highest vertex of the triangle
  p.max <- which.max(c(dfr$az[i],dfr$bz[i],dfr$cz[i])); p.max
  
  # Position vector of the vertex a, b or c
  pv.edge <- c(dfr[i,c(p.max[1]*3)-1],dfr[i,c(p.max[1]*3)])
  
  # Selection of the other two vertices
  pts.tri <- c(1,2,3)
  p.rem <- which(pts.tri!=p.max); p.rem
  
  # First direction vector
  dv1 <- c(dfr[i,c(p.rem[1]*3)-1],dfr[i,c(p.rem[1]*3)])-pv.edge; dv1
  # First system of equations
  par.pq.1 <- par_eqli(dv1,vg,pv.edge,pv.crpt); par.pq.1
  # Coordinates of the POI
  poi.1 <- eqli(pv.edge,par.pq.1[1],dv1); poi.1
  poi.z.1 <- z(x=poi.1[1],y=poi.1[2],A,B,C,D); poi.z.1
  
  # Erster direction vector
  dv2 <- c(dfr[i,c(p.rem[2]*3)-1],dfr[i,c(p.rem[2]*3)])-pv.edge; dv2
  # Second system of equations
  par.pq.2 <- par_eqli(dv2,vg,pv.edge,pv.crpt); par.pq.2
  # Coordinates of the POI
  poi.2 <- eqli(pv.edge,par.pq.2[1],dv2); poi.2
  poi.z.2 <- z(x=poi.2[1],y=poi.2[2],A,B,C,D); poi.z.2 
  
  # Selection of the correct POI
  if(!is.na(par.pq.1[1]) && round(par.pq.1[1]>=0,4) && round(par.pq.1[1]<=1,4) && 
     round(poi.z.1,4)>round(crpt[3],4)){
    poi <- c(poi.1,poi.z.1) 
    edge.poi <- sort(c(p.rem[1],p.max))
    p.sep <- p.rem[2]
  } else if(!is.na(par.pq.2[1]) && round(par.pq.2[1]>=0,4) && round(par.pq.2[1]<=1,4) && 
            round(poi.z.2,4)>round(crpt[3],4)){
    poi <- c(poi.2,poi.z.2) 
    edge.poi <- sort(c(p.rem[2],p.max))
    p.sep <- p.rem[1]
  } else {print(i); print('There is a problem!')}
  
  return(list(poi=poi,edge.poi=edge.poi,p.sep=p.sep))
}


### Subdivison of the corresponding triangles
tri_subdiv <- function(dfr=dftin, dfr.up=dftin.up, crpt, i){
  
  pt.a <- paste(dfr$ax,dfr$ay,round(dfr$az,0),sep='_');pt.a
  pt.b <- paste(dfr$bx,dfr$by,round(dfr$bz,0),sep='_');pt.b
  pt.c <- paste(dfr$cx,dfr$cy,round(dfr$cz,0),sep='_');pt.c
  pt.lst <- list(pt.a,pt.b,pt.c);pt.lst
  
  stas <- steepest_ascent_poi(dfr,crpt,i)
  points(x=stas$poi[1],y=stas$poi[2],col='blue', pch=16, cex=2.5, add=TRUE)
  
  end <- FALSE
  quadr <- FALSE
  ver <- FALSE
  while(end==FALSE){    
    
    ### Subdivide the triangles
    TID.max <- max(dfr.up$TID)
    tri.sub1 <- data.frame(TID=TID.max+1,
                           ax=crpt[1], ay=crpt[2], az=crpt[3],
                           bx=stas$poi[1], by=stas$poi[2], bz=stas$poi[3],
                           cx=dfr[i,c(stas$edge.poi[1]*3)-1], cy=dfr[i,c(stas$edge.poi[1]*3)], cz=dfr[i,c(stas$edge.poi[1]*3)+1]); tri.sub1
    dfr.up <- rbind(dfr.up,tri.sub1)
    
    tri.sub2 <- data.frame(TID=TID.max+2,
                           ax=crpt[1], ay=crpt[2], az=crpt[3],
                           bx=stas$poi[1], by=stas$poi[2], bz=stas$poi[3],
                           cx=dfr[i,c(stas$edge.poi[2]*3)-1], cy=dfr[i,c(stas$edge.poi[2]*3)], cz=dfr[i,c(stas$edge.poi[2]*3)+1]); tri.sub2
    dfr.up <- rbind(dfr.up,tri.sub2)
    
    # Plot new edge of the triangles
    seg1 <- cbind(c(tri.sub1$ax,tri.sub1$bx),c(tri.sub1$ay,tri.sub1$by));seg1
    seg1L <- Line(seg1);seg1L
    seg1Ls <- Lines(list(seg1L), ID=1);seg1Ls
    seg1SL <- SpatialLines(list(seg1Ls));seg1SL
    plot(seg1SL, col="magenta", lwd=4, add=TRUE)
    
    # Subdivide the quadrilateral
    if(quadr==TRUE){
      
      crds.edge.1 <- c(paste(dfr[i.pre,c(edge.poi.pre[1]*3)-1], dfr[i.pre,c(edge.poi.pre[1]*3)],sep="_"),
                       paste(dfr[i.pre,c(edge.poi.pre[2]*3)-1], dfr[i.pre,c(edge.poi.pre[2]*3)],sep="_")) 
      
      crds.edge.2 <- c(paste(dfr[i,c(stas$edge.poi[1]*3)-1], dfr[i,c(stas$edge.poi[1]*3)], sep="_"),
                       paste(dfr[i,c(stas$edge.poi[2]*3)-1], dfr[i,c(stas$edge.poi[2]*3)], sep="_"))
      
      pt.3 <- stas$edge.poi[-which(crds.edge.2%in%crds.edge.1)]; pt.3
      
      tri.sub3 <- data.frame(TID=TID.max+3,
                             ax=crpt[1], ay=crpt[2], az=crpt[3],
                             bx=dfr[i,c(stas$p.sep*3)-1], by=dfr[i,c(stas$p.sep*3)], bz=dfr[i,c(stas$p.sep*3)+1],
                             cx=dfr[i,c(pt.3*3)-1], cy=dfr[i,c(pt.3*3)], cz=dfr[i,c(pt.3*3)+1]); tri.sub3
      dfr.up <- rbind(dfr.up,tri.sub3)
      
      # Plot new edge
      seg1 <- cbind(c(tri.sub3$ax,tri.sub3$cx),c(tri.sub3$ay,tri.sub3$cy));seg1
      seg1L <- Line(seg1);seg1L
      seg1Ls <- Lines(list(seg1L), ID=1);seg1Ls
      seg1SL <- SpatialLines(list(seg1Ls));seg1SL
      plot(seg1SL, col="orange", lwd=4, add=TRUE)
      
      # Testing if the resulting triangle is too thin
      v1 <- stas$poi-c(crpt[1],crpt[2],crpt[3])
      v2 <- c(dfr[i,c(pt.3*3)-1],dfr[i,c(pt.3*3)],dfr[i,c(pt.3*3)+1])-c(crpt[1],crpt[2],crpt[3])
     
      cos.phi <- abs(v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3])/
        (sqrt(v1[1]^2+v1[2]^2+v1[3]^2)*sqrt(v2[1]^2+v2[2]^2+v2[3]^2))
      phi <- acos(cos.phi)*180/pi; phi
      
      if(phi<3){#Adjust the triangle
        
        # tri.sub1
        dfr.up$bx[nrow(dfr.up)-2] <- dfr[i,c(stas$edge.poi[2]*3)-1]
        dfr.up$by[nrow(dfr.up)-2] <- dfr[i,c(stas$edge.poi[2]*3)]
        dfr.up$bz[nrow(dfr.up)-2] <- dfr[i,c(stas$edge.poi[2]*3)+1]
        
        # tri.sub2
        dfr.up <- dfr.up[-c(nrow(dfr.up)-1),]
        
        plot(seg1SL, col="magenta", lwd=4, add=TRUE)
        
        ver <- TRUE
      }  
    }    
            
    dfr.up <- dfr.up[-which(dfr.up$TID==dfr$TID[i]),]; tail(dfr.up) 
    quadr <- TRUE
    
    ### Find the relevant neighbour  
    na <- c(which(pt.a==pt.lst[[stas$edge.poi[1]]][i]),which(pt.b==pt.lst[[stas$edge.poi[1]]][i]),which(pt.c==pt.lst[[stas$edge.poi[1]]][i]));na
    nb <- c(which(pt.a==pt.lst[[stas$edge.poi[2]]][i]),which(pt.b==pt.lst[[stas$edge.poi[2]]][i]),which(pt.c==pt.lst[[stas$edge.poi[2]]][i]));nb
    nab <- na[na%in%nb];nab
        
    if(length(nab)>1){
      TID.nb <- dfr$TID[nab[which(nab!=i)]]
    } else{end <- TRUE}; which(dfr$TID==TID.nb)     
    
    if(end==FALSE){
      
      # Determine where the path continues at an edge
      # Vector of the common edge of the triangles
      x <- dfr[i,c(stas$edge.poi*3)-1];x  
      y <- dfr[i,c(stas$edge.poi*3)];y  
          
      if(stas$edge.poi[1]==1 && stas$edge.poi[2]==3){
        ve.x <- as.numeric(x[1]-x[2]); 
        ve.y <- as.numeric(y[1]-y[2]); 
      } else{
        ve.x <- as.numeric(x[2]-x[1]);
        ve.y <- as.numeric(y[2]-y[1]); 
      }
          
      # Gradient vector of the neighbour triangle
      vg.x <- dfr$grdv_x[which(dfr$TID==TID.nb)]
      vg.y <- dfr$grdv_y[which(dfr$TID==TID.nb)]
      
      # Calculate vector product
      matrixvp <- matrix(c(vg.x,vg.y,ve.x,ve.y),nrow=2,byrow=TRUE); matrixvp
      det <- det(matrixvp); det
      
      if(ver==TRUE){
        det <- -1
        ver <- FALSE
      }
        
        ### Determine where the path continues at a vertex
        if(det<0){
                             
          first.run <- TRUE
          while(quadr==TRUE){
          
            # Define the vertex
            if(first.run==TRUE){
            pt.rel <- stas$edge.poi[which.max(c(dfr[i,c(stas$edge.poi[1]*3)+1],dfr[i,c(stas$edge.poi[2]*3)+1]))];pt.rel
            first.run <- FALSE
            }
            ver.x <- dfr[i,c(pt.rel*3)-1]; ver.x
            ver.y <- dfr[i,c(pt.rel*3)]; ver.y
            ver.z <- dfr[i,c(pt.rel*3)+1]; ver.z
            vertex <- paste(ver.x,ver.y,round(ver.z,0),sep='_');vertex
            
            points(x=ver.x,y=ver.y,col='green', pch=16, cex=2.5, add=TRUE)
  
            # Find adjacent triangles
            nbpts <- c(which(pt.a==vertex),which(pt.b==vertex),which(pt.c==vertex));nbpts
            nbs.TID <- dfr$TID[nbpts[which(nbpts!=i)]]; nbs.TID
            
            ai <- NULL
            bi <- NULL
            ci <- NULL
            if(any(dfr$az[nbpts]>ver.z)){
              ai <- nbpts[which(dfr$az[nbpts]>ver.z)];ai
            }
            
            if(any(dfr$bz[nbpts]>ver.z)){
              bi <- nbpts[which(dfr$bz[nbpts]>ver.z)];bi
            }
            
            if(any(dfr$cz[nbpts]>ver.z)){
              ci <- nbpts[which(dfr$cz[nbpts]>ver.z)];ci
            }
            
            nb.rel <- unique(c(ai,bi,ci)); nb.rel
                      
            if(is.null(nb.rel) || length(nb.rel)==1){end <- TRUE; quadr <- FALSE; print(i);print("end")} else{
              
              # Order relevant neighbour triangles in clockwise direction
              nb.rel.cw <- nb.rel[circ_sort(rpt2d = data.frame(x=ver.x,y=ver.y), 
                                            pts2d = data.frame(x=dfr$cen_x[nb.rel],y=dfr$cen_y[nb.rel]),
                                            clockwise=TRUE)];nb.rel.cw
              
              dfv <- data.frame(path=rep(NA,length(nb.rel.cw)),
                                nb.ver=rep(NA,length(nb.rel.cw)),
                                ver.rel=rep(NA,length(nb.rel.cw))); dfv
                               
              for(j in 1:length(nb.rel.cw)){#j=1
                
                #  Determine relevant points and their heights    
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
                if(z.pt[1]>ver.z){
                  v.verpt.x <- dfr[nb.rel.cw[j],c(pts.rel[1]*3)-1]-ver.x
                  v.verpt.y <- dfr[nb.rel.cw[j],c(pts.rel[1]*3)]-ver.y
                  v.ptver.x <- ver.x-dfr[nb.rel.cw[j],c(pts.rel[1]*3)-1]
                  v.ptver.y <- ver.y-dfr[nb.rel.cw[j],c(pts.rel[1]*3)]
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
                
                if(z.pt[1]<=ver.z){
                  v.verpt.x <- dfr[nb.rel.cw[j],c(pts.rel[1]*3)-1]-ver.x;
                  v.verpt.y <- dfr[nb.rel.cw[j],c(pts.rel[1]*3)]-ver.y;
                }
                
                # Calculate vector products to test the triangle
                matrixvp3 <- matrix(c(vg.cr.x,vg.cr.y,v.verpt2.x,v.verpt2.y),nrow=2,byrow=TRUE); matrixvp3
                det.eri <- det(matrixvp3);det.eri
                matrixvp4 <- matrix(c(vg.cr.x,vg.cr.y,v.verpt.x,v.verpt.y),nrow=2,byrow=TRUE); matrixvp4
                det.ecr <- det(matrixvp4);det.ecr
                
                # Determine the course
                if(!is.nan(det.le) && det.le<0 && det.cr<0){
                  path  <- "ridge"} else if(det.eri>=0 && det.ecr<=0){
                  path <- "surface"} else{path <- NA}
                
                dfv$path[j] <- path
                if(!is.na(path) && path=="surface"){
                  dfv$nb.ver[j] <- nb.rel.cw[j]}
                if(!is.na(path) && path=="ridge"){
                  dfv$nb.ver[j] <- nb.rel.cw[j]
                  dfv$ver.rel[j] <- pts.rel[1]
                }                    
                
              }#j
              
              pa.rel <- which(!is.na(dfv$path)==TRUE); pa.rel 
              print(i)     
              print(dfv$path[pa.rel])
              
              if(!length(pa.rel)){
                end <- TRUE
                quadr <- FALSE
              }
              if(length(pa.rel)==1){
                if(dfv$path[pa.rel]=="surface"){
                  i <- dfv$nb.ver[pa.rel]
                  crpt[1] <- ver.x
                  crpt[2] <- ver.y
                  crpt[3] <- ver.z
                  quadr <- FALSE
                  stas <- steepest_ascent_poi(dfr,crpt,i)
                 
                  points(x=stas$poi[1],y=stas$poi[2],col='red', pch=16, cex=2.5,add=TRUE)
                  
                  # Testing if the resulting triangle is too thin
                  v1 <- stas$poi-c(crpt[1],crpt[2],crpt[3])
                  v2 <- c(dfr[i,c(stas$edge.poi[1]*3)-1],dfr[i,c(stas$edge.poi[1]*3)],dfr[i,c(stas$edge.poi[1]*3)+1])-
                    c(crpt[1],crpt[2],crpt[3])
                  v3 <- c(dfr[i,c(stas$edge.poi[2]*3)-1],dfr[i,c(stas$edge.poi[2]*3)],dfr[i,c(stas$edge.poi[2]*3)+1])-
                    c(crpt[1],crpt[2],crpt[3])
                  
                  cos.phi <- abs(v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3])/
                    (sqrt(v1[1]^2+v1[2]^2+v1[3]^2)*sqrt(v2[1]^2+v2[2]^2+v2[3]^2))
                  phi1 <- acos(cos.phi)*180/pi; phi1
                  
                  cos.phi <- abs(v1[1]*v3[1]+v1[2]*v3[2]+v1[3]*v3[3])/
                    (sqrt(v1[1]^2+v1[2]^2+v1[3]^2)*sqrt(v3[1]^2+v3[2]^2+v3[3]^2))
                  phi2 <- acos(cos.phi)*180/pi; phi2
                                                    
                  if(phi1<3){
                    pt.rel <- stas$edge.poi[1]
                    quadr <- TRUE
                  }
                  if(phi2<3){
                    pt.rel <- stas$edge.poi[2]
                    quadr <- TRUE
                  }      
                }
                if(dfv$path[pa.rel]=="ridge"){
                  i <- dfv$nb.ver[pa.rel]
                  pt.rel <- dfv$ver.rel[pa.rel]
                }
              }
              if(length(pa.rel)>1){
                pts.z <- NULL
                for (k in 1:length(pa.rel)){#k=1
                  if(dfv$path[pa.rel[k]]=="surface"){
                    stas <- steepest_ascent_poi(dfr,crpt=c(ver.x,ver.y,ver.z),i=dfv$nb.ver[pa.rel[k]])
                    pt.z <- stas$poi[3]
                    } else{pt.z <- dfr[dfv$nb.ver[pa.rel[k]], dfv$ver.rel[pa.rel[k]]*3+1]}
                  pts.z[k] <- pt.z
                }#for
                pa.rel <- pa.rel[which.max(pts.z)]
                
                if(dfv$path[pa.rel]=="surface"){
                  i <- dfv$nb.ver[pa.rel]
                  crpt[1] <- ver.x
                  crpt[2] <- ver.y
                  crpt[3] <- ver.z
                  quadr <- FALSE
                  stas <- steepest_ascent_poi(dfr,crpt,i)
                  
                  points(x=stas$poi[1],y=stas$poi[2],col='red',pch=16, cex=2.5,add=TRUE)
                  
                  # Testing if the resulting triangle is too thin
                  v1 <- stas$poi-c(crpt[1],crpt[2],crpt[3])
                  v2 <- c(dfr[i,c(stas$edge.poi[1]*3)-1],dfr[i,c(stas$edge.poi[1]*3)],dfr[i,c(stas$edge.poi[1]*3)+1])-
                    c(crpt[1],crpt[2],crpt[3])
                  v3 <- c(dfr[i,c(stas$edge.poi[2]*3)-1],dfr[i,c(stas$edge.poi[2]*3)],dfr[i,c(stas$edge.poi[2]*3)+1])-
                    c(crpt[1],crpt[2],crpt[3])
                  
                  cos.phi <- abs(v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3])/
                    (sqrt(v1[1]^2+v1[2]^2+v1[3]^2)*sqrt(v2[1]^2+v2[2]^2+v2[3]^2))
                  phi1 <- acos(cos.phi)*180/pi; phi1
                  
                  cos.phi <- abs(v1[1]*v3[1]+v1[2]*v3[2]+v1[3]*v3[3])/
                    (sqrt(v1[1]^2+v1[2]^2+v1[3]^2)*sqrt(v3[1]^2+v3[2]^2+v3[3]^2))
                  phi2 <- acos(cos.phi)*180/pi; phi2
                  
                  if(phi1<3){
                    pt.rel <- stas$edge.poi[1]
                    quadr <- TRUE
                  }
                  if(phi2<3){
                    pt.rel <- stas$edge.poi[2]
                    quadr <- TRUE
                  }
                }
                if(dfv$path[pa.rel]=="ridge"){
                  i <- dfv$nb.ver[pa.rel]
                  pt.rel <- dfv$ver.rel[pa.rel]
                }
              }
            
            }#if 
          }#while  
        }#if
    
       if(quadr==TRUE){
         edge.poi.pre <- stas$edge.poi 
         i.pre <- i
         
         i <- which(dfr$TID==TID.nb); i #i triangle
         crpt <- stas$poi
         stas <- steepest_ascent_poi(dfr,crpt,i)
         
         points(x=stas$poi[1],y=stas$poi[2],col='red', pch=16, cex=2.5,add=TRUE)        
      }
    }#if
  }#while 
  
  return(dfr.up)
}#fun 
