library(spatstat)
library(raster) # read tif
library(rgdal) # read shp
library(maptools) # read shp
require(ggplot2) # historgram
library(RColorBrewer) # color bar
library(classInt)

# Read observation data and study area
mydata <- readShapeSpatial("DATA/Obspoint_Proj.shp")
mydata_window <- readShapeSpatial("DATA/Obs_bound_Proj.shp")

obs_bound <- as(mydata_window, "owin")
obspoint <- as(mydata,"ppp")

mypattern <- ppp(obspoint$x, obspoint$y, obs_bound)


# density 
dens <- nrow(mydata)/area(mydata_window) 
print(paste("The density is:", dens)) # 1.01177249175275e-09

#Quadrat counting tests for CSR
M <- quadrat.test(mypattern, nx = 10, ny = 15)
plot(mypattern)
plot(M, add=T,cex = 1.25)
M$p.value
#K-S test
ks <- cdf.test(mypattern, "x")
plot(ks)
pval <- ks$p.value
pval

# Load covariate 
imgpop <- raster('DATA/Pop.tif')
imgroad <- raster('DATA/Roaddensity.tif')
imglc_11 <-raster('DATA/cover_11.tif') # open water  11
imglc_21 <-raster('DATA/cover_21.tif') # Developed, open space  21
imglc_22 <-raster('DATA/cover_22.tif') # Developed, Low intensity  22
imglc_23 <-raster('DATA/cover_23.tif') # Developed, Medium Intensity  23
imglc_24 <-raster('DATA/cover_24.tif') # Developed, High Intensity 24
imglc_31 <-raster('DATA/cover_31.tif') # Barren 31 
imglc_41 <-raster('DATA/cover_41.tif') # Deciduous 41
imglc_42 <-raster('DATA/cover_42.tif') # Evergrenn 42
imglc_43 <-raster('DATA/cover_43.tif') # Mixed 43
imglc_52 <-raster('DATA/cover_52.tif') # Shrub 52
imglc_71 <-raster('DATA/cover_71.tif') # Grassland 71
imglc_81 <-raster('DATA/cover_81.tif') # Pasture 81
imglc_82 <-raster('DATA/cover_82.tif') # Cultivated crops 82
imglc_90 <-raster('DATA/cover_90.tif') # Woody Wetlands 90
imglc_95 <-raster('DATA/cover_95.tif') # Emergent Herbaceous Wetlands 95

c1<-as.im(imgpop)
c2<-as.im(imgroad)
c3 <- as.im(imglc_11)
c4 <- as.im(imglc_21)
c5<- as.im(imglc_22)
c6<- as.im(imglc_23)
c7<- as.im(imglc_24)
c8<- as.im(imglc_31)
c9<- as.im(imglc_41) 
c10<- as.im(imglc_42)
c11<- as.im(imglc_43)
c12<- as.im(imglc_52)
c13 <- as.im(imglc_71)
c14 <- as.im(imglc_81)
c15 <- as.im(imglc_82)
c16 <- as.im(imglc_90)
c17 <- as.im(imglc_95)  
clist<-list(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17)

# Spatial distribution
par(mfrow=c(3,2), mar=c(2,1,2,1))
cl <- rev(terrain.colors(100))
for (i in clist){
  plot(i, las=1, col=cl, main=NULL)
}

# density on a covaraite
par(mfrow=c(3,2), mar=c(3,2,3,2))
cl2 <- brewer.pal(4, "YlOrRd")

for (i in clist){
  brk  <- c( -Inf, 0.25, 0.5, 0.75 , Inf)  # Define the breaks
  Zcut <- cut(i, breaks=brk, labels=1:4)  # Classify the raster
  E    <- tess(image=Zcut)  # Create a tesselated surface
  Q   <- quadratcount(mypattern, tess = E)  # Tally counts
  Q.d <- intensity(Q) 
  plot(intensity(Q, image=TRUE), las=1, col=cl2, main=NULL)
}

# MODEL

R <- ppm(mypattern,~pop+rd+lc11+lc21+lc22+lc23+lc24+lc31+lc41+lc42+lc43+lc43+lc52+lc71+lc81+lc82+lc90+lc95,
          covariates = list(pop=c1,rd=c2,
                            lc11=c3,lc21=c4,lc22=c5,lc23=c6,lc24=c7,
                            lc31=c8,lc41=c9,lc42=c10,lc43=c11,lc52=c12,
                            lc71=c13,lc81=c14,lc82=c15,lc90=c16,lc95=c17))


fitnull<-ppm(mypattern,~1)
AIC(fitnull)
AIC(R)
step(R) # Automatic model selection

# Plot fitted intensity, trend - two different angles
plot.ppm(R,cif=TRUE,ngrid=100,how="persp")
plot.ppm(R,cif=TRUE,ngrid=100,how="persp",theta=-30,phi=40,d=4,ticktype="detailed",zlab="z")


# extract quadrature points in corresponding order
quadpoints <- union.quad(quad.ppm(R))

# plot conditional intensity values
# as circles centred on the quadrature points 
quadmarked <- setmarks(quadpoints, fitted(R))
plot(quadmarked)
# plot effect size of each covariate
coef<- coef(R)
fitcif<-coef[2:18]
#hist(fitcif,freq = FALSE,labels=TRUE,right=FALSE,breaks=17)
df.cif<- data.frame(fitcif)
rowname<-c("Pop","Road","Open water","Open space","Low intensity","Medium intensity","High intensity","Barren","Deciduous","Evergreen","Mixed","Shrub","Grassland","Pasture","Cultivated crops","Woody wetlands","Herbaceous Wetlands")

dfcif<-data.frame(cif=fitcif,name=rowname)
ggplot(dfcif,aes(x=cif,y=name)) + 
  theme(text = element_text(size=35)) + geom_bar(stat = "identity")

quad.ppm(R)
# Quadrature scheme (Berman-Turner)
# 2916 data points, 5652 dummy points # 8568 in total
# 110 x 110 grid of dummy points, plus 4 corner points
# dummy spacing: 23056.27 x 23345.00 units
# Total weight 2873006897997.52


# smooth residual
diagnose.ppm(R, which = "smooth")
# Model diagnostics (raw residuals)
# Diagnostics available:
#   smoothed residual field
# range of smoothed field =  [-5.273e-10, 1.318e-09]

# extract residuals
rr <- residuals(R, type="raw")
rr.value <-rr$val # extract value # hist(rr.value)
rr.x<- rr$loc$x
rr.y<- rr$loc$y
res2<- data.frame(x = rr.x, y = rr.y)
res <- data.frame(x = rr.x, y = rr.y,Residuals = rr.value)
# rasterize resigual
#imgbd <- raster('E:/Thesis/DATA/Road/Obs_boundary_raster1.tif')
ras.rr <- rasterize(res2,imgpop,res$Residuals) # imgtree 
hist(ras.rr, maxpixels = ncell(ras.rr))
tmp <- tempdir()
writeRaster(ras.rr, filename=file.path(tmp, "residual.tif"), format="GTiff", overwrite=TRUE)

h = hist(rr.value) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE)

