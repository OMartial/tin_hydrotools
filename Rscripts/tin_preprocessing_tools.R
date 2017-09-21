# Preprocessing of data

### Removal of flat triangles
mod_tri_flat <- function(dfr=dftin, tpts){
  
  pt.a<-paste(dfr$ax,dfr$ay,round(dfr$az,0),sep='_');pt.a
  pt.b<-paste(dfr$bx,dfr$by,round(dfr$bz,0),sep='_');pt.b
  pt.c<-paste(dfr$cx,dfr$cy,round(dfr$cz,0),sep='_');pt.c
  
  pts.xyz <- paste(tp$P[,1],tp$P[,2],round(tp$PA,0),sep='_'); tail(pts.xyz)
  
  # Find flat triangles
  tri.faulty1 <- NULL
  for(i in 1:nrow(dfr)){
    if(round(dfr$az[i],4)==round(dfr$bz[i],4) && round(dfr$az[i],4)==round(dfr$cz[i],4) && 
         round(dfr$bz[i],4)==round(dfr$cz[i],4)){
      if(is.null(tri.faulty1)){tri.faulty1 <- i} else{tri.faulty1 <- c(tri.faulty1,i)}
    }
  }
  
  if(!is.null(tri.faulty1)){
    for(i in 1:length(tri.faulty1)){#i=1
      
      if(round(dfr$az[tri.faulty1[i]],4)==round(dfr$bz[tri.faulty1[i]],4) && round(dfr$az[tri.faulty1[i]],4)==round(dfr$cz[tri.faulty1[i]],4) && 
           round(dfr$bz[tri.faulty1[i]],4)==round(dfr$cz[tri.faulty1[i]],4)){
      
        pt.abc.xyz <- c(paste(dfr$ax[tri.faulty1[i]],dfr$ay[tri.faulty1[i]],round(dfr$az[tri.faulty1[i]],0),sep='_'),
                        paste(dfr$bx[tri.faulty1[i]],dfr$by[tri.faulty1[i]],round(dfr$bz[tri.faulty1[i]],0),sep='_'),
                        paste(dfr$cx[tri.faulty1[i]],dfr$cy[tri.faulty1[i]],round(dfr$cz[tri.faulty1[i]],0),sep='_'))
        
        ptz.means <- NULL 
        for(j in 1:3){#j=1
          
          nbpts <- c(which(pt.a==pt.abc.xyz[j]),which(pt.b==pt.abc.xyz[j]),which(pt.c==pt.abc.xyz[j]));nbpts
        
          pt.a.coord <- as.matrix(dfr[nbpts,c(2:4)])
          pt.b.coord <- as.matrix(dfr[nbpts,c(5:7)])
          pt.c.coord <- as.matrix(dfr[nbpts,c(8:10)])
          pt.coord <- unique(rbind(pt.a.coord,pt.b.coord,pt.c.coord))
          
          ptz.mean <- mean(pt.coord[,3]); ptz.mean
               
          ptz.means[j] <- ptz.mean
        }
        
        ptz.means.ord <- order(ptz.means)
        
        # Modify the heigths in tpts
        print(i); print(dfr$az[tri.faulty1[i]])
        pti <- which(pts.xyz==pt.abc.xyz[ptz.means.ord[2]]); tpts[pti,]
        tpts[pti,3] <- tpts[pti,3]+0.01; print(tpts[pti,3])
        
        pti <- which(pts.xyz==pt.abc.xyz[ptz.means.ord[3]]); tpts[pti,]
        tpts[pti,3] <- tpts[pti,3]+0.02; print(tpts[pti,3])
        
        dfr <- data.frame(TID=1:nrow(tp$T),
                          ax=tpts[tp$T[,1],1],
                          ay=tpts[tp$T[,1],2],
                          az=tpts[tp$T[,1],3],
                          bx=tpts[tp$T[,2],1],
                          by=tpts[tp$T[,2],2],
                          bz=tpts[tp$T[,2],3],
                          cx=tpts[tp$T[,3],1],
                          cy=tpts[tp$T[,3],2],
                          cz=tpts[tp$T[,3],3])
        
        pt.a<-paste(dfr$ax,dfr$ay,round(dfr$az,0),sep='_')
        pt.b<-paste(dfr$bx,dfr$by,round(dfr$bz,0),sep='_')
        pt.c<-paste(dfr$cx,dfr$cy,round(dfr$cz,0),sep='_')
        
        pts.xyz <- paste(tpts[,1],tpts[,2],round(tpts[,3],0),sep='_')
      
      }#if(flat triangle)
    }#for
  } else{print('There are no flat triangles!')}
  
  return(tpts)
}

### Removal of pits
rm_pits <- function(dfr=dftin, tpts, pts.brd){

  pt.a<-paste(dfr$ax,dfr$ay,round(dfr$az,0),sep='_');pt.a
  pt.b<-paste(dfr$bx,dfr$by,round(dfr$bz,0),sep='_');pt.b
  pt.c<-paste(dfr$cx,dfr$cy,round(dfr$cz,0),sep='_');pt.c
  
  brd.coords <- NULL
  brd.lines <- vector("list", nrow(pts.brd))
  for(i in 1:nrow(pts.brd)){#i=1
    
    if(i<nrow(pts.brd)){
      brd.coords <- Line(cbind(c(pts.brd[i,1],pts.brd[i+1,1]),c(pts.brd[i,2],pts.brd[i+1,2])))
    } else {
      brd.coords <- Line(cbind(c(pts.brd[i,1],pts.brd[1,1]),c(pts.brd[i,2],pts.brd[1,2])))
    }
    brd.lines[i] <- Lines(list(brd.coords), ID=i)
  }
  
  brd.SL <- SpatialLines(brd.lines)
  #plot(brd.SL, col="yellow", add=TRUE)
  
  pts.xyz <- paste(tpts[,1],tpts[,2],round(tpts[,3],0),sep='_')
  
  # Find pits
  pits <- NULL
  for(i in 1:length(pts.xyz)){#i=1
    
    pt.z <- tpts[i,3]; pt.z
  
    # Find adjacent vertices
    nbpts <- c(which(pt.a==pts.xyz[i]),which(pt.b==pts.xyz[i]),which(pt.c==pts.xyz[i]));nbpts
    #TID.nb <- dfr$TID[nbpts[which(nbpts!=ssi[i])]];TID.nb
    
    if(!any(dfr$az[nbpts]<pt.z) && !any(dfr$bz[nbpts]<pt.z) && !any(dfr$cz[nbpts]<pt.z)){
      
      pit <- matrix(tpts[i,], nrow=1)
      pit.SP <- SpatialPoints(pit)
      pit.SP.buf <- buffer(pit.SP, 100)
      
      if(rgeos::gDisjoint(brd.SL, pit.SP.buf)){
        plot(pit.SP, col='green', pch=16 , add=TRUE)
      
        if(is.null(pits)){pits <- i} else{pits <- c(pits,i)}
      }
    }
  }#for 
  pits.ord <- pits[order(tpts[pits,3],decreasing=TRUE)]
  print(pits.ord)
  
  # Find breaching path
  for(i in 1:length(pits.ord)){#i=1
    
    pit.ori <- c(tpts[pits.ord[i],1],tpts[pits.ord[i],2],(tpts[pits.ord[i],3]))
    pit.xyz <- paste(pit.ori[1],pit.ori[2],round(pit.ori[3],0),sep='_')
    nbpts <- c(which(pt.a==pit.xyz),which(pt.b==pit.xyz),which(pt.c==pit.xyz));nbpts
    
    if(!any(dfr$az[nbpts]<pit.ori[3]) && !any(dfr$bz[nbpts]<pit.ori[3]) && !any(dfr$cz[nbpts]<pit.ori[3])){
    
      points(pit.ori[1],pit.ori[2], pch=16, col='red', add=TRUE)
  
      pt.a.coord <- as.matrix(dfr[nbpts,c(2:4)])
      pt.b.coord <- as.matrix(dfr[nbpts,c(5:7)])
      pt.c.coord <- as.matrix(dfr[nbpts,c(8:10)])
      pt.coord.ori <- unique(rbind(pt.a.coord,pt.b.coord,pt.c.coord))
      pt.coord.ord <- pt.coord.ori[order(pt.coord.ori[,3]),]
      xslct <- which(pt.coord.ord[,1]==pit.ori[1]);xslct
      yslct <- which(pt.coord.ord[,2]==pit.ori[2]);yslct
      pt.coord.ord <- pt.coord.ord[-xslct[which(xslct%in%yslct)],]
      
      x=0
      repeat{x=x+1
        pt.coord.pre <- pt.coord.ori
        
        pt.min <- pt.coord.ord[x,]
        pt.xyz <- paste(pt.min[1],pt.min[2],round(pt.min[3],0),sep='_')
        #points(pt.min[1],pt.min[2], col='orange', add=TRUE, pch=16)
        
        brch.path <- c(pits.ord[i], which(pts.xyz==pt.xyz))
        while(pit.ori[3]<=pt.min[3] && length(brch.path)<9){
          
          nbpts <- c(which(pt.a==pt.xyz),which(pt.b==pt.xyz),which(pt.c==pt.xyz));nbpts
          
          pt.a.coord <- as.matrix(dfr[nbpts,c(2:4)])
          pt.b.coord <- as.matrix(dfr[nbpts,c(5:7)])
          pt.c.coord <- as.matrix(dfr[nbpts,c(8:10)])
          pt.coord <- unique(rbind(pt.a.coord,pt.b.coord,pt.c.coord))
      
          pt.coord.xy <- paste(pt.coord[,1],pt.coord[,2],round(pt.coord[,3],0),sep='_') 
          pt.coord.pre.xy <- paste(pt.coord.pre[,1],pt.coord.pre[,2],round(pt.coord.pre[,3],0),sep='_') 
          pts.rel <-  matrix(pt.coord[which(!(pt.coord.xy%in%pt.coord.pre.xy)),],ncol=3)
         
          if(length(pts.rel)!=0){
            pt.min <- pts.rel[which.min(pts.rel[,3]),]
            pt.xyz <- paste(pt.min[1],pt.min[2],round(pt.min[3],0),sep='_')
            #points(pt.min[1],pt.min[2], col='yellow', add=TRUE, pch=16)
            
            pt.coord.pre <- rbind(pt.coord.pre, pts.rel)
            brch.path <- c(brch.path,which(pts.xyz==pt.xyz)); brch.path 
          } else{
            brch.path <- c(brch.path,1:7)
          }
        }#while
        
        if(pit.ori[3]>pt.min[3]){print(i);print(brch.path);break}
        if(x==nrow(pt.coord.ord)){print(i);print('There is no lower point in the vicinity!');break}
      }#repeat
      
      # Modify the heigths in tpts
      if(pit.ori[3]>pt.min[3]){
        pts.z <- tpts[brch.path,3]
        z.diff <- pts.z[1]-pts.z[length(pts.z)]
        z.int <- z.diff/(length(brch.path)-1)
        
        y <- length(brch.path)-2
        for(j in 1:y){#j=1
          tpts[brch.path[j+1],3] <- tpts[brch.path[j],3]-z.int
          print(pts.z[j+1])
          print(tpts[brch.path[j],3]-z.int)
        }
        
        dfr <- data.frame(TID=1:nrow(tp$T),
                          ax=tpts[tp$T[,1],1],
                          ay=tpts[tp$T[,1],2],
                          az=tpts[tp$T[,1],3],
                          bx=tpts[tp$T[,2],1],
                          by=tpts[tp$T[,2],2],
                          bz=tpts[tp$T[,2],3],
                          cx=tpts[tp$T[,3],1],
                          cy=tpts[tp$T[,3],2],
                          cz=tpts[tp$T[,3],3])
        
        pt.a<-paste(dfr$ax,dfr$ay,round(dfr$az,0),sep='_');pt.a
        pt.b<-paste(dfr$bx,dfr$by,round(dfr$bz,0),sep='_');pt.b
        pt.c<-paste(dfr$cx,dfr$cy,round(dfr$cz,0),sep='_');pt.c
        
        pts.xyz <- paste(tpts[,1],tpts[,2],round(tpts[,3],0),sep='_')
       
      }
    }#if(pit)
  }#for

  return(tpts)
}  
  
### Find faulty triangles
find_tri_faulty <- function(dfr=dftin){#i=1
  tri.faulty <- NULL
  for(i in 1:nrow(dfr)){
    if(round(dfr$ax[i],4)==round(dfr$bx[i],4) && round(dfr$ax[i],4)==round(dfr$cx[i],4) && 
         round(dfr$bx[i],4)==round(dfr$cx[i],4)){
      if(is.null(tri.faulty)){tri.faulty <- i} else{tri.faulty <- c(tri.faulty,i)}
    }
  }
  
  for(i in 1:nrow(dfr)){
    if(round(dfr$ay[i],4)==round(dfr$by[i],4) && round(dfr$ay[i],4)==round(dfr$cy[i],4) && 
         round(dfr$by[i],4)==round(dfr$cy[i],4)){
      if(is.null(tri.faulty)){tri.faulty <- i} else{tri.faulty <- c(tri.faulty,i)}
    }
  }
  
  for(i in 1:nrow(dfr)){
    if(round(dfr$ax[i],2)==round(dfr$bx[i],2) && round(dfr$ay[i],2)==round(dfr$by[i],2) ||
         round(dfr$bx[i],2)==round(dfr$cx[i],2) && round(dfr$by[i],2)==round(dfr$cy[i],2) ||
         round(dfr$ax[i],2)==round(dfr$cx[i],2) && round(dfr$ay[i],2)==round(dfr$cy[i],2)){
      if(is.null(tri.faulty)){tri.faulty <- i} else{tri.faulty <- c(tri.faulty,i)}
    }
  }
  
  if(!is.null(tri.faulty)){
    tri.faulty <- unique(tri.faulty)
  }
  
  return(tri.faulty)
} 

### Create SpatialPolygons from triangulated points
tri_polygons <- function(dfr=dftin){
  tri.coords <- NULL
  tri.poly <- vector("list", nrow(dfr))
  for(i in 1:nrow(dfr)){#i=1
    tri.coords <- Polygon(cbind(c(dfr$ax[i],dfr$bx[i],dfr$cx[i]),c(dfr$ay[i],dfr$by[i],dfr$cy[i])))
    
    tri.poly[i] <- Polygons(list(tri.coords), ID=dfr$TID[i])
  }  
  
  tri.SP <- SpatialPolygons(tri.poly)
  
  #shapefile(x=tri.SP, filename='triangles_R.shp')
  return(tri.SP)
}    

### Order the vertices of the triangle in clockwise direction
sort_ccw <- function(dfr=dftin){
  
  for (i in 1:nrow(dfr)){#i=1
    matrixpts <- matrix(c(dfr$ax[i],dfr$ay[i],1,
                          dfr$bx[i],dfr$by[i],1,
                          dfr$cx[i],dfr$cy[i],1),
                        nrow=3, byrow=TRUE);matrixpts
    det <- det(matrixpts);det
    
    if(det>0){print('Coordinates of the points B and C are changed.')}
    if(det<0){print('Coordinates of the points B und C are kept.')}
    
    if(det>0){ 
      tmp <- dfr$bx[i]
      dfr$bx[i] <- dfr$cx[i] 
      dfr$cx[i] <- tmp
      tmp <- dfr$by[i]
      dfr$by[i] <- dfr$cy[i]
      dfr$cy[i] <- tmp
      tmp <- dfr$bz[i]
      dfr$bz[i] <- dfr$cz[i] 
      dfr$cz[i] <- tmp
    }
  }
  
  return(dfr)
}
 