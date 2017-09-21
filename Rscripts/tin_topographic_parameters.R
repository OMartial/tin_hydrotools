# Topographic parameters

# Calculate the area of a triangle
tri_area <- function(dfr=dftin.ws){
  
  for(i in 1:nrow(dfr)){#i=1
    
    # 2D
    matrixpts <- matrix(c(1,1,1,
                          dfr$ax[i],dfr$bx[i],dfr$cx[i],
                          dfr$ay[i],dfr$by[i],dfr$cy[i]),
                        nrow=3, byrow=TRUE);matrixpts
    
    area2d <- abs(det(matrixpts))/2; area2d
      
    # 3D
    v1 <- c(dfr$bx[i]-dfr$ax[i],dfr$by[i]-dfr$ay[i],dfr$bz[i]-dfr$az[i])
    
    v2 <- c(dfr$cx[i]-dfr$ax[i],dfr$cy[i]-dfr$ay[i],dfr$cz[i]-dfr$az[i])
    
    area3d <- sqrt((v1[2]*v2[3]-v1[3]*v2[2])^2
                  +(v1[3]*v2[1]-v1[1]*v2[3])^2
                  +(v1[1]*v2[2]-v1[2]*v2[1])^2)/2; area3d
    
    print(i); print(paste('area2d:',area2d)); print(paste('area3d:',area3d))
    dfr$area_2d[i] <- area2d 
    dfr$area_3d[i] <- area3d
  }  
    
  return(dfr)
}

# Calculate the slope and the aspect of a triangle
tri_slope_asp <- function(dfr=dftin.ws,n2=c(0,0,1)){
  
  for(i in 1:nrow(dfr)){#i=1
  
    # Direction vectors of the plane
    dv1.x <- dfr$bx[i]-dfr$ax[i]
    dv1.y <- dfr$by[i]-dfr$ay[i]
    dv1.z <- dfr$bz[i]-dfr$az[i]
    
    dv2.x <- dfr$cx[i]-dfr$ax[i]
    dv2.y <- dfr$cy[i]-dfr$ay[i]
    dv2.z <- dfr$cz[i]-dfr$az[i]
    
    # Cross product (normal vector)
    nx <- dv1.y * dv2.z - dv1.z * dv2.y
    ny <- dv1.z * dv2.x - dv1.x * dv2.z
    nz <- dv1.x * dv2.y - dv1.y * dv2.x
    n1 <- c(nx,ny,nz)
    
    # calculate angle between two normal vectors n1 and n2
    # result: angle in degree
    cos.alpha <- abs(n1[1]*n2[1]+n1[2]*n2[2]+n1[3]*n2[3])/
      (sqrt(n1[1]^2+n1[2]^2+n1[3]^2)*sqrt(n2[1]^2+n2[2]^2+n2[3]^2))
    slope <- acos(cos.alpha)*180/pi

    # calculate the aspect
    a <- n1[3]/n1[1];a
    b <- n1[3]/n1[2];b
    
    if(b>0 & a>0){asp <- 180/pi*atan2(b,a);asp}
    if(b>0 & a<0){asp <- 180/pi*atan2(b,a)+180;asp}
    if(b<0 & a>0){asp <- (180/pi*atan2(b,a)) %% 180;asp}
    if(b<0 & a<0){asp <- (180/pi*atan2(b,a)) %% 360;asp}
    if(b==0 & a<0){asp <- 90}
    if(b==0 & a>0){asp <- 270}
    if(a==0 & b==0){asp <- NA}
    
    print(i); print(paste('slope:',slope)); print(paste('aspect:',asp))
    dfr$slope[i] <- slope 
    dfr$aspect[i] <- asp
  }
    
  return(dfr)
}






