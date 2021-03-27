# library
library(tidyverse)
library(geoR)
library(gstat)
library(sp)
library(maps)
ussoils <- read_csv("soils.csv")

### 2.1 ######################################################
par(mfrow=c(1,3))
boxplot(ussoils$Se_ppm, main = "Boxplot of Se_ppm")
hist(ussoils$Se_ppm, xlab = "Se_ppm",main = "Histogram of Se_ppm")
plot(ecdf(ussoils$Se_ppm), xlab = "Se_ppm", ylab = "Cumulative Probability", 
     main = "Empirical Cumulative Distribution of Se_ppm")
curve(pnorm(x,mean(ussoils$Se_ppm),sd(ussoils$Se_ppm)),-1, 4, add=T, col =2)
grid(NULL,NULL)
legend("bottomright",legend=c("ecdf","normal cdf"),col=1:2,lty=1,cex=1,box.lty = 0)

### 2.2 ######################################################
par(mfrow=c(1,3))
log_Se <- log(ussoils$Se_ppm)
boxplot(log_Se, main = "Boxplot of log_Se")
hist(log_Se,freq = F, xlab = "log_Se",main = "Histogram of log_Se")
lines(density(log_Se))
plot(ecdf(log_Se),xlab = "log_Se", ylab = "Cumulative Probability", 
     main = "ECDF of log_Se")
curve(pnorm(x,mean(log_Se),sd(log_Se)),-3, 2, add=T, col=2)
grid(NULL,NULL)
legend("bottomright",legend=c("ecdf","normal cdf"),col=1:2,lty=1,cex=1,box.lty = 0)

### 2.3 ######################################################
par(mfrow=c(1,2))
plot(ussoils$Longitude,log_Se,xlab="Longitude",ylab = "log_Se", 
     main="Plot of log_Se against Longitude")
abline(lm(log_Se~Longitude,data=ussoils),col="red")
plot(ussoils$Latitude,log_Se,xlab="Latitude",ylab = "log_Se", 
     main="Plot of log_Se against Latitude")
abline(lm(log_Se~Latitude,data=ussoils),col="red")

### 3.1 ######################################################
hscat<-cbind(ussoils[,c(2,3)],log_Se)
coordinates(hscat) <- ~Longitude+Latitude 
hscat(log_Se~1,hscat, seq(0,3,0.5)) 

### 3.2 ######################################################
plot(ussoils$Longitude,ussoils$Latitude, xlim=c(-125,-65),ylim=c(25,50), xlab="Longitude",
     ylab="Latitude", main="Log_Selemium Concentration in the United States", "n")
map("usa",add=TRUE)
points(ussoils$Longitude,ussoils$Latitude, 
       cex=log_Se/mean(log_Se),pch=19,col="orange")


geodata<-as.geodata(cbind(ussoils[,c(2,3)],log_Se))
v_classic<-variog(geodata,max.dist = 20)
v_robust<-variog(geodata, max.dist=20,estimator.type="modulus")
plot(v_classic,main = "Classical Estimator")
plot(v_classic,main = "Robust Estimator")

cloud_classic <- variog(geodata,max.dist=20,option = "cloud")
cloud_robust <- variog(geodata,estimator.type="modulus", max.dist=20,option = "cloud")
par(mfrow=c(2,2))
plot(cloud_classic, main="classical estimator")
plot(cloud_robust, main = "robust estimator")
plot(variog(geodata, bin.cloud=T, max.dist=20), bin.cloud=T)
plot(variog(geodata, bin.cloud=T, estimator.type="modulus", max.dist=20), bin.cloud=T)

par(mfrow=c(1,2))
v_30 <- variog(geodata, dir=30*pi/180, max.dist=20)
plot(v_30, main = "dir = 30 degrees")
v_60 <- variog(geodata, dir=60*pi/180, max.dist=20)
plot(v_60, main = "dir = 60 degrees")

v_4dir <- variog4(geodata, max.dist=20)
plot(v_4dir, main = "dir = 0, 45, 90, 135 degrees")

fit_exp<-variofit(v_robust,cov.model = "exp",ini.cov.pars = c(0.1,10),
                  fix.nugget = F, nugget = 0.35)
fit_exp_cressie<-variofit(v_robust,cov.model = "exp",ini.cov.pars = c(0.1,10),
                          fix.nugget = F, nugget = 0.35, weights = "cressie")
fit_sph<-variofit(v_classic,cov.model = "sph",ini.cov.pars = c(0.1,10),
                  fix.nugget = F, nugget = 0.35)
fit_sph_cressie<-variofit(v_classic,cov.model = "exp",ini.cov.pars = c(0.1,10),
                          fix.nugget = F, nugget = 0.35, weights = "cressie")

plot(v_classic, main = "Omnidirectional Variogram")
lines(fit_exp,col=2)
lines(fit_exp_cressie,col=3)
lines(fit_sph,col=4)
lines(fit_sph_cressie,col=5)
legend("bottomright",
       legend=c("exp and default weights","exp and Cressie's weights",
                "sph","sph and Cressie's weights"),col=2:5,lty=1,cex=0.6,box.lty = 0)


cv_exp <- xvalid(geodata, model=fit_exp)
cv_exp_cressie <- xvalid(geodata, model=fit_exp_cressie)
cv_sph <- xvalid(geodata, model=fit_sph)
cv_sph_cressie <- xvalid(geodata, model=fit_sph_cressie)

dif_exp <- log(ussoils$Se_ppm) - cv_exp$predicted
dif_exp_cressie <- log(ussoils$Se_ppm) - cv_exp_cressie$predicted
dif_sph <- log(ussoils$Se_ppm) - cv_sph$predicted
dif_sph_cressie <- log(ussoils$Se_ppm) - cv_exp_cressie$predicted
paste("PRESS of exponential model is:",sum(dif_exp^2)) 
paste("PRESS of exponential model using cressie's weights is:",sum(dif_exp_cressie^2))
paste("PRESS of spherical model is:",sum(dif_sph^2)) 
sum(dif_sph_cressie^2)


### 4.1 ######################################################
g_ok<-gstat(id="log_Se",formula =log(Se_ppm)~1,
            locations=~Longitude+Latitude,data=ussoils)
fit_ok<-fit.variogram(variogram(g_ok), vgm(0.1,"Sph",10,0.35))
g_uk<-gstat(id="log_Se",formula =log(Se_ppm)~Longitude+Latitude,
            locations=~Longitude+Latitude,data=ussoils)
fit_uk<-fit.variogram(variogram(g_uk), vgm(0.1,"Sph",10,0.35))
g_co<-gstat(id="log_Se",formula =log(Se_ppm)~1,
            locations=~Longitude+Latitude,data=ussoils)
g_co<-gstat(g_co,id="log_Cu",formula =log(Cu_ppm)~1,
            locations=~Longitude+Latitude,data=ussoils)
g_co<-gstat(g_co,id="log_Ba",formula =log(Ba_ppm)~1,
            locations=~Longitude+Latitude,data=ussoils)
fit_co<-fit.lmc(variogram(g_co),g_co,vgm(0.1,"Sph",10,0.35))

plot(variogram(g_ok),fit_ok, main = "ordinary kriging")
plot(variogram(g_uk),fit_uk,main = "universal kriging")
plot(variogram(g_co),fit_co,main="cokriging")

cv_ok<- krige.cv(log(Se_ppm)~1,data=ussoils, locations=~Longitude+Latitude,
                 model=fit_ok,nfold=nrow(ussoils))
cv_uk<- krige.cv(log(Se_ppm)~Longitude+Latitude,data=ussoils,
                 locations=~Longitude+Latitude, model=fit_uk,nfold=nrow(ussoils))
cv_co <- gstat.cv(fit_co,nfold = nrow(ussoils))

paste("PRESS of ordinary kriging is:",sum(cv_ok$residual^2)/nrow(ussoils)) 
paste("PRESS of universal kriging is:",sum(cv_uk$residual^2)/nrow(ussoils)) 
paste("PRESS of cokriging is:",sum(cv_co$residual^2)/nrow(ussoils)) 

### 4.2 ######################################################
x.range <- as.integer(range(ussoils[,2]))
y.range <- as.integer(range(ussoils[,3]))
grd <- expand.grid(Longitude=seq(from=x.range[1], to=x.range[2], by=1),
                   Latitude=seq(from=y.range[1], to=y.range[2], by=1))

pr_ck <- predict(fit_co,grd)
pred_table <- matrix(pr_ck$log_Se.pred,
                     length(seq(from=x.range[1], to=x.range[2], by=1)),
                     length(seq(from=y.range[1], to=y.range[2], by=1)))

image(seq(from=x.range[1], to=x.range[2], by=1),
      seq(from=y.range[1], to=y.range[2], by=1), pred_table,
      xlab="Longitude", ylab="Latitude", 
      main="Raster map of the Predicted Values")
contour(seq(from=x.range[1], to=x.range[2], by=1), 
        seq(from=y.range[1],to=y.range[2], by=1), pred_table, add=TRUE, labcex=1)
points(cbind(ussoils$Longitude,ussoils$Latitude))

var_table <- matrix(pr_ck$log_Se.var,
                    length(seq(from=x.range[1], to=x.range[2], by=1)),
                    length(seq(from=y.range[1], to=y.range[2], by=1)))
image(seq(from=x.range[1], to=x.range[2], by=1),
      seq(from=y.range[1], to=y.range[2], by=1), var_table,
      xlab="Longitude", ylab="Latitude", 
      main="Raster map of Variances")
contour(seq(from=x.range[1], to=x.range[2], by=1), 
        seq(from=y.range[1],to=y.range[2], by=1), var_table, add=TRUE,
        col="black",labcex=1)