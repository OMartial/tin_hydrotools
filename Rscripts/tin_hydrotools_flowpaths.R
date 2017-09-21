# 

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

### Determine the centroid of the triangle
centroid_triangle <- function(dfr=dftin){
  
  for (i in 1:nrow(dfr)){#i=1
    dfp <- data.frame(x=c(dfr$ax[i],dfr$bx[i],dfr$cx[i]),
                      y=c(dfr$ay[i],dfr$by[i],dfr$cy[i]),
                      z=c(dfr$az[i],dfr$bz[i],dfr$cz[i]));dfp  
  
    cen <- c(sum(dfp[,1])/3,sum(dfp[,2])/3,sum(dfp[,3])/3)
    
    dfr$cen_x[i] <- cen[1]
    dfr$cen_y[i] <- cen[2]
    dfr$cen_z[i] <- cen[3]
  }
    
  return(dfr)
}

### Linear equation in parametric form
eqli <- function(pv,p,dv){
  pv + p*dv
}

### Equation of the plane
z <- function(x,y,A,B,C,D){
  -(A/C*x+B/C*y+D/C)
}

### Solve the system of equations with two unknown
par_eqli <- function(dv1, dv2, pv1, pv2){
  
  matrix.seq <- matrix(c(dv1,-dv2), nrow=2)
  vector.seq <- c(pv2-pv1)
  
  if(matrix.seq[1,1]==0 && matrix.seq[1,2]==0 || matrix.seq[2,1]==0 && matrix.seq[2,2]==0 || 
     matrix.seq[1,1]+matrix.seq[2,1]==0 && matrix.seq[1,2]+matrix.seq[2,2]==0 ||
     round(matrix.seq[1,2]*matrix.seq[2,1],4)==round(matrix.seq[2,2]*matrix.seq[1,1],4)){sol <- NA} else{
    sol <- solve(a=matrix.seq, b=vector.seq)
  }
  
  return(sol)
}


### Calculation of the path of steepest descent on a triangle 
steepest_descent_poi <- function(dfr=dftin){
  
  dfr$grdv_x <- NA
  dfr$grdv_y <- NA
  dfr$poi_x <- NA
  dfr$poi_y <- NA
  dfr$poi_z <- NA
  dfr$edge_poi <- NA
  dfr$vertex <- NA
  
  for(i in 1:nrow(dfr)){#i=1
    poi <- NULL
    
    # Position vector of the centroid
    pv.cen  <- c(dfr$cen_x[i],dfr$cen_y[i])
    
    # Parameters of the equation of the plane 
    A <- dfr$ay[i]*(dfr$bz[i]-dfr$cz[i])+dfr$by[i]*(dfr$cz[i]-dfr$az[i])+dfr$cy[i]*(dfr$az[i]-dfr$bz[i]);A
    B <- dfr$az[i]*(dfr$bx[i]-dfr$cx[i])+dfr$bz[i]*(dfr$cx[i]-dfr$ax[i])+dfr$cz[i]*(dfr$ax[i]-dfr$bx[i]);B
    C <- dfr$ax[i]*(dfr$by[i]-dfr$cy[i])+dfr$bx[i]*(dfr$cy[i]-dfr$ay[i])+dfr$cx[i]*(dfr$ay[i]-dfr$by[i]);C
    D <- -A*dfr$ax[i]-B*dfr$ay[i]-C*dfr$az[i];D
    
    # Gradient vector
    delta.x <- A/C
    delta.y <- B/C
    vg <-  c(delta.x,delta.y); vg
    
    # Position vector of the vertex a
    pva <- c(dfr$ax[i],dfr$ay[i])
    
    # Position vector of the vertex b
    pvb <- c(dfr$bx[i],dfr$by[i])
    
    # Direction vectors
    # Line ab
    dv.ab <- pvb-pva
    # Line ac
    dv.ac <- c(dfr$cx[i],dfr$cy[i])-pva
    # Line bc
    dv.bc <- c(dfr$cx[i],dfr$cy[i])-pvb
    
    print(i)
    # First system of equations
    par.pq.1 <- par_eqli(dv.ab,vg,pva,pv.cen); par.pq.1
    # Coordinates of the POI
    poi.1 <- eqli(pva,par.pq.1[1],dv.ab); poi.1
    poi.z.1 <- z(x=poi.1[1],y=poi.1[2],A,B,C,D)
    
    # Second system of equations
    par.pq.2 <- par_eqli(dv.ac,vg,pva,pv.cen); par.pq.2
    # Coordinates of the POI
    poi.2 <- eqli(pva,par.pq.2[1],dv.ac); poi.2
    poi.z.2 <- z(x=poi.2[1],y=poi.2[2],A,B,C,D) 
    
    # Third system of equations
    par.pq.3 <- par_eqli(dv.bc,vg,pvb,pv.cen); par.pq.3
    # Coordinates of the POI
    poi.3 <- eqli(pvb,par.pq.3[1],dv.bc); poi.3
    poi.z.3 <- z(x=poi.3[1],y=poi.3[2],A,B,C,D) 
    
    # Selection of the correct POI
    if(!is.na(par.pq.1[1]) && round(par.pq.1[1],4)>=0 && round(par.pq.1[1],4)<=1 && 
       poi.z.1<dfr$cen_z[i]){
      poi <- c(poi.1,poi.z.1) 
      edge.poi <- paste(1,2,sep='_')
    } else if(!is.na(par.pq.2[1]) && round(par.pq.2[1],4)>=0 && round(par.pq.2[1],4)<=1 && 
              poi.z.2<dfr$cen_z[i]){
      poi <- c(poi.2,poi.z.2) 
      edge.poi <- paste(1,3,sep='_')
    } else if(!is.na(par.pq.3[1]) && round(par.pq.3[1],4)>=0 && round(par.pq.3[1],4)<=1 && 
              poi.z.3<dfr$cen_z[i]){
      poi <- c(poi.3,poi.z.3) 
      edge.poi <- paste(2,3,sep='_')
    } else {print('There is a problem!')} 
    
    dfr$grdv_x[i] <- vg[1]
    dfr$grdv_y[i] <- vg[2]
    dfr$poi_x[i] <- poi[1]
    dfr$poi_y[i] <- poi[2]
    dfr$poi_z[i] <- poi[3]
    dfr$edge_poi[i] <- edge.poi
    
    poi.xyz <- paste(poi[1],poi[2],round(poi[3],0),sep='_')
    pt.a.xyz <- paste(dfr$ax[i],dfr$ay[i],round(dfr$az[i],0),sep='_')
    pt.b.xyz <- paste(dfr$bx[i],dfr$by[i],round(dfr$bz[i],0),sep='_')
    pt.c.xyz <- paste(dfr$cx[i],dfr$cy[i],round(dfr$cz[i],0),sep='_')
    if(poi.xyz==pt.a.xyz){dfr$edge_poi[i] <- 'no.edge'; dfr$vertex[i] <- 1}
    if(poi.xyz==pt.b.xyz){dfr$edge_poi[i] <- 'no.edge'; dfr$vertex[i] <- 2}
    if(poi.xyz==pt.c.xyz){dfr$edge_poi[i] <- 'no.edge'; dfr$vertex[i] <- 3}
  }
  
  return(dfr)  
}


### Find the relevant neighbour triangle
find_neigb<-function(dfr=dftin){

  pt.a <- paste(dfr$ax,dfr$ay,round(dfr$az,0),sep='_');pt.a
  pt.b <- paste(dfr$bx,dfr$by,round(dfr$bz,0),sep='_');pt.b
  pt.c <- paste(dfr$cx,dfr$cy,round(dfr$cz,0),sep='_');pt.c
  pt.lst <- list(pt.a,pt.b,pt.c);pt.lst
  
  for(i in 1:nrow(dfr)){
    if(dfr$edge[i]!='no.edge'){
      
      pts.rel <- as.numeric(unlist(strsplit(dfr$edge[i],split="_")))
   
      na <- c(which(pt.a==pt.lst[[pts.rel[1]]][i]),which(pt.b==pt.lst[[pts.rel[1]]][i]),which(pt.c==pt.lst[[pts.rel[1]]][i]));na
      nb <- c(which(pt.a==pt.lst[[pts.rel[2]]][i]),which(pt.b==pt.lst[[pts.rel[2]]][i]),which(pt.c==pt.lst[[pts.rel[2]]][i]));nb
      nab <- na[na%in%nb];nab
      
      if(length(nab)>1){
      TID.nb <- dfr$TID[nab[which(nab!=i)]]
      } else{TID.nb <- NA}
      dfr$nb_poi[i] <- TID.nb 
      
    } else{dfr$nb_poi[i] <- 'no.nb'}
  }
  
  return(dfr)  
}  


### Determine where the flow path continues at an edge
determine_flowtype<-function(dfr=dftin){
       
  for(i in 1:nrow(dfr)){#i=1
    if(!is.na(dfr$nb_poi[i]) && dfr$nb_poi[i]!='no.nb'){  
      
      pts.rel <- as.numeric(unlist(strsplit(dfr$edge_poi[i],split="_")));pts.rel
      
      # Vector of the common edge of the triangles
      x <- dfr[i,c(pts.rel*3)-1];x  
      y <- dfr[i,c(pts.rel*3)];y  
      
      if(pts.rel[1]==1 && pts.rel[2]==3){
        ve.x <- as.numeric(x[1]-x[2]); 
        ve.y <- as.numeric(y[1]-y[2]); 
      } else{
        ve.x <- as.numeric(x[2]-x[1]);
        ve.y <- as.numeric(y[2]-y[1]); 
      }
      
      # Gradient vector of the neighbour triangle
      if(is.na(dfr$nb_poi[i])){vg.x <- NA; vg.y <- NA} else{
        vg.x <- dfr$grdv_x[which(dfr$TID==dfr$nb_poi[i])]
        vg.y <- dfr$grdv_y[which(dfr$TID==dfr$nb_poi[i])]
      }
      
      # Calculate vector product
      a <- c(vg.x,vg.y)
      b <- c(ve.x,ve.y)
      det <- a[1]*b[2]-a[2]*b[1]; det
      
      # Determine flow type
      if(det>0){flowtype  <- "channel"} else{flowtype  <- "overland"}
      
      dfr$flowtype_poi[i] <- flowtype
      
    } else if(is.na(dfr$nb_poi[i])){dfr$flowtype_poi[i] <- NA} else{
      dfr$flowtype_poi[i] <- 'no.flowtype'
    }    
  }
   
  return(dfr)
}  


### Order points around a reference point
circ_sort <- function(rpt2d,pts2d,clockwise=TRUE){
  a <- atan2(rpt2d$y-pts2d$y,rpt2d$x-pts2d$x);a
  if(clockwise==FALSE){return(order(a))}
  if(clockwise==TRUE){return(order(a,decreasing = TRUE))}
}


### Determine where the flow path continues at a vertex
flowpath_vertex <- function(dfr=dftin){
  
  dfr$flowtype_vertex <- NA
  dfr$nb_vertex <- NA
  dfr$vertex_2 <- NA
    
  pt.a<-paste(dfr$ax,dfr$ay,round(dfr$az,0),sep='_');pt.a
  pt.b<-paste(dfr$bx,dfr$by,round(dfr$bz,0),sep='_');pt.b
  pt.c<-paste(dfr$cx,dfr$cy,round(dfr$cz,0),sep='_');pt.c
  
  ssi <- c(which(dfr$flowtype_poi=='channel'),which(dfr$flowtype_poi=='no.flowtype'));ssi
  sdfr <- dfr[ssi,];sdfr
  for(i in 1:nrow(sdfr)){#i=1
    if(sdfr$nb_poi[i]!='no.nb'){  
    
      # Define the vertex
      pts.edge <- as.numeric(unlist(strsplit(sdfr$edge[i],split="_"))); pts.edge
      pts.edge.z <- c(sdfr[i,c(pts.edge[1]*3)+1],sdfr[i,c(pts.edge[2]*3)+1])
      pt.rel <- pts.edge[which.min(pts.edge.z)];pt.rel
      if(pts.edge.z[1]==pts.edge.z[2]){
        ptz.means <- NULL
        for(j in 1:length(pts.edge)){#j=1
          vertex <- paste(sdfr[i,c(pts.edge[j]*3)-1],sdfr[i,c(pts.edge[j]*3)],round(pts.edge.z[j],0),sep='_')
          
          nbpts <- c(which(pt.a==vertex),which(pt.b==vertex),which(pt.c==vertex));nbpts
          
          pt.az <- dfr$az[nbpts]
          pt.bz <- dfr$bz[nbpts]  
          pt.cz <- dfr$cz[nbpts]
          ptz <- unique(c(pt.az,pt.bz,pt.cz))
          ptz.rel <- ptz[-which(ptz==pts.edge.z[j])]
          
          ptz.mean <- mean(ptz); ptz.mean
          ptz.means[j] <- ptz.mean
        }#for
        pt.rel <- pts.edge[which.min(ptz.means)]
      }#if
    } else{pt.rel <- sdfr$vertex[i]} 
    
    ver.x <- sdfr[i,c(pt.rel*3)-1]; ver.x
    ver.y <- sdfr[i,c(pt.rel*3)]; ver.y
    ver.z <- sdfr[i,c(pt.rel*3)+1]; ver.z
    vertex <- paste(ver.x,ver.y,round(ver.z,0),sep='_');vertex
    
    # Find adjacent triangles
    nbpts <- c(which(pt.a==vertex),which(pt.b==vertex),which(pt.c==vertex));nbpts
    TID.nb <- dfr$TID[nbpts[which(nbpts!=ssi[i])]];TID.nb
    
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
      dfr$flowtype_vertex[ssi][i] <- "sink/end of TIN"
      dfr$vertex[ssi][i] <- pt.rel
      print(i)
      print("sink/end of TIN")
    } else{
      
      # Order relevant neighbour triangles in clockwise direction
      nb.rel.cw <- nb.rel[circ_sort(rpt2d = data.frame(x=ver.x,y=ver.y), 
                                    pts2d = data.frame(x=dfr$cen_x[nb.rel],y=dfr$cen_y[nb.rel]),
                                    clockwise=TRUE)];nb.rel.cw
      
      dfv <- data.frame(flowtype=rep(NA,length(nb.rel.cw)),
                        nb.ver=rep(NA,length(nb.rel.cw)),
                        ver.rel=rep(NA,length(nb.rel.cw)))
                             
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
        
        if(z.pt[1]>=ver.z){
          v.verpt.x <- dfr[nb.rel.cw[j],c(pts.rel[1]*3)-1]-ver.x
          v.verpt.y <- dfr[nb.rel.cw[j],c(pts.rel[1]*3)]-ver.y
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
      
      dfr$vertex[ssi][i] <- pt.rel
      if(!length(ft.rel)){
        dfr$flowtype_vertex[ssi][i] <- "end?"}
      if(length(ft.rel)==1){
        dfr$flowtype_vertex[ssi][i] <- dfv$flowtype[ft.rel] 
        dfr$nb_vertex[ssi][i] <- dfv$nb.ver[ft.rel]
        dfr$vertex_2[ssi][i] <- dfv$ver.rel[ft.rel]
      }
      if(length(ft.rel)>1){
        pts.z <- NULL
        for (k in 1:length(ft.rel)){#k=1
          if(dfv$flowtype[ft.rel[k]]=="overland"){pt.z <- dfr$poi_z[dfv$nb.ver[ft.rel[k]]]} else{
             pt.z <- dfr[dfv$nb.ver[ft.rel[k]], dfv$ver.rel[ft.rel[k]]*3+1]}
          pts.z[k] <- pt.z
        }#k
        ft.rel <- ft.rel[which.min(pts.z)]
        
        dfr$flowtype_vertex[ssi][i] <- dfv$flowtype[ft.rel] 
        dfr$nb_vertex[ssi][i] <- dfv$nb.ver[ft.rel]
        dfr$vertex_2[ssi][i] <- dfv$ver.rel[ft.rel]
      }
    
    }#if 
  }#i
  
  return(dfr)    
}#fun    




  
  