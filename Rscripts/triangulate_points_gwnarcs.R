library(raster)
library(RTriangle)
setwd("C:/TEMP/moggier/Flow_Path_Project/TIN_tools")


# Read input
# Data points of the terrain surface
fn <- 'Flow_Path_Project/TIN_tools/input_8_final/points_turt_z20.shp'
pts <- shapefile(fn)
plot(pts)

# Points of the river network
fn <- 'GIS/TLM_Fliessgewässer_nodes_turt.shp'
pts.gwn <- shapefile(fn)
plot(pts.gwn)
dupl <- which(duplicated(pts.gwn@coords))
pts.gwn <- pts.gwn[-dupl,]

# Arcs of the river network
fn <- 'GIS/TLM_Fliessgewässer_arcs_turt_clean.shp'
gwn <- shapefile(fn)
plot(gwn, col='blue')

pts.gwn.xy <- paste(round(pts.gwn@coords[,1],1),round(pts.gwn@coords[,2],1),sep='_')

# Identify the relevant points of the river network
fntn <- matrix(nrow=nrow(gwn), ncol=2)
for(i in 1:nrow(gwn)){#i=1
  fn.str <- gwn$fnode_wkt[i]
  fn <- wktpz2df(x=fn.str)
  fn.xy <- paste(round(fn$x,1),round(fn$y,1),sep='_')
  fni <- which(fn.xy==pts.gwn.xy); fni
  if(!length(fni)){fni <- NA}
  fntn[i,1] <- fni
  
  tn.str <- gwn$tnode_wkt[i]
  tn <- wktpz2df(x=tn.str)
  tn.xy <- paste(round(tn$x,1),round(tn$y,1),sep='_')
  tni <- which(tn.xy==pts.gwn.xy); tni
  if(!length(tni)){tni <- NA}
  fntn[i,2] <- tni
} 

nodes <-  na.omit(unique(c(fntn[,1],fntn[,2])))
pts.gwn.rel <- pts.gwn[nodes,]

dem <- raster('C:/Users/Martial Oggier/Dropbox/Eigene Dateien/Masterarbeit/swissalti_2011_turt_10m.tif')
pts.gwn.rel <- extract(dem, pts.gwn.rel, sp=TRUE)

pts.coords <- rbind(coordinates(pts)[,1:2],coordinates(pts.gwn.rel)[,1:2])
pts.z  <- c(coordinates(pts)[,3],pts.gwn.rel$swissalti_2011_turt_10m)

pts.xy <- paste(round(pts.coords[,1],1),round(pts.coords[,2],1),sep='_')

# Identify the segments for the PSLG 
arcs <- matrix(nrow=nrow(gwn), ncol=2)
for(i in 1:nrow(gwn)){#i=1
  fn.str <- gwn$fnode_wkt[i]
  fn <- wktpz2df(x=fn.str)
  fn.xy <- paste(round(fn$x,1),round(fn$y,1),sep='_')
  fni <- which(fn.xy==pts.xy); fni
  if(!length(fni)){fni <- NA}
  arcs[i,1] <- fni
  
  tn.str <- gwn$tnode_wkt[i]
  tn <- wktpz2df(x=tn.str)
  tn.xy <- paste(round(tn$x,1),round(tn$y,1),sep='_')
  tni <- which(tn.xy==pts.xy); tni
  if(!length(tni)){tni <- NA}
  arcs[i,2] <- tni
} 

noarc <- c(which(is.na(arcs[,1])),which(is.na(arcs[,2])))
arcs <- arcs[-noarc,]

pts.brd <- rbind(pts.all[8405,],pts.all[9955,],pts.all[30660,],pts.all[25694,])
points(x=pts.brd[,1], y=pts.brd[,2], col='yellow', pch=16, add=TRUE)

# Triangulation
p <- RTriangle::pslg(P = pts.coords, PA=pts.z,
                     S=rbind(c(8405, 9955), c(9955, 30660), c(30660, 25694), c(25694, 8405), arcs))

tp <- RTriangle::triangulate(p, D=TRUE)
plot(tp)
str(tp)

# Assign heights to the added points
pts.xy.ori <- paste(pts.coords[,1],pts.coords[,2],sep='_')
pts.xy.all <- paste(tp$P[,1],tp$P[,2],sep='_')

pts.new <- which(!(pts.xy.all%in%pts.xy.ori)); pts.new
pts.new.z <- extract(dem.turt, tp$P[pts.new,])
tp$PA[pts.new] <- pts.new.z

tpts <- cbind(tp$P,tp$PA);head(tpts)

dftin <- data.frame(TID=1:nrow(tp$T),
                    ax=tpts[tp$T[,1],1],
                    ay=tpts[tp$T[,1],2],
                    az=tpts[tp$T[,1],3],
                    bx=tpts[tp$T[,2],1],
                    by=tpts[tp$T[,2],2],
                    bz=tpts[tp$T[,2],3],
                    cx=tpts[tp$T[,3],1],
                    cy=tpts[tp$T[,3],2],
                    cz=tpts[tp$T[,3],3])
tail(dftin)
