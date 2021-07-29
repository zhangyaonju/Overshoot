#### get average contribution, sensitivity, trend and per-biome statistics
args = commandArgs(trailingOnly=F)
print(args)

myargs1 <-sub('-','',args[length(args)-3])
myargs2 <-sub('-','',args[length(args)-2])
myargs3 <-sub('-','',args[length(args)-1])
myargs4 <-sub('-','',args[length(args)])

spinup<-as.numeric(myargs1)
spinrep<-as.numeric(myargs2)
delta<-c(0.95,0.98,0.99,0.995,0.999)[as.numeric(myargs3)+1]
ind_rand<-as.numeric(myargs4)
indicators<-c("prec",'prec_temp','prec_temp_solr','randomndvi')

if (ind_rand>3){
  ind <-4
  rand<-(ind_rand-4)%%5+1
  window_ind<-ceiling((ind_rand-3)/5)-1
  indicator<-indicators[ind]
  out_label<-paste0(indicator,window_ind,rand)
}else{
  indicator<-indicators[ind_rand]
  out_label<-indicator
}



setwd('/global/scratch/yaozhang/Project/overshooting/')

library(raster)
library(ncdf4)
library(abind)

## deseasonlized and detrended variables are included
scenario<-paste0(out_label,".spinup",spinup,".rep",spinrep,".delta",delta)
ncin<-nc_open(paste0("./analysis/gimms/gimms.DLM.results.",scenario,".nc"))
predictmean<-ncvar_get(ncin,varid="Predictmean")
Xpre<-ncvar_get(ncin,varid="Regrevariable")
nc_close(ncin)
xdim = dim(Xpre)[4]

trend<-predictmean[,,1:(414+spinup*12*spinrep),1]
seas1<-predictmean[,,1:(414+spinup*12*spinrep),3+xdim]
seas2<-predictmean[,,1:(414+spinup*12*spinrep),5+xdim]
seas3<-predictmean[,,1:(414+spinup*12*spinrep),7+xdim]
#seas4<-predictmean[,,1:(414+spinup*12*spinrep),9+xdim]
y_pre1<-predictmean[,,1:(414+spinup*12*spinrep),3]*Xpre[,,1:(414+spinup*12*spinrep),1]
y_pre23<-predictmean[,,1:(414+spinup*12*spinrep),4]*Xpre[,,1:(414+spinup*12*spinrep),2]
y_pre46<-predictmean[,,1:(414+spinup*12*spinrep),5]*Xpre[,,1:(414+spinup*12*spinrep),3]
y_pre712<-predictmean[,,1:(414+spinup*12*spinrep),6]*Xpre[,,1:(414+spinup*12*spinrep),4]
y_pre1324<-predictmean[,,1:(414+spinup*12*spinrep),7]*Xpre[,,1:(414+spinup*12*spinrep),5]
y_pre13<-predictmean[,,1:(414+spinup*12*spinrep),8]*Xpre[,,1:(414+spinup*12*spinrep),6]
y_temp<-predictmean[,,1:(414+spinup*12*spinrep),9]*Xpre[,,1:(414+spinup*12*spinrep),7]
y_solr<-predictmean[,,1:(414+spinup*12*spinrep),10]*Xpre[,,1:(414+spinup*12*spinrep),8]
rm(predictmean)

lat<-ncdim_def(name = 'latitude',units="degree",vals=(-179:180-0.5)/2)
lon<-ncdim_def(name = "longitude",units="degree",vals=(-359:360-0.5)/2)

tim<-ncdim_def(name = 'time', units='months since 1981/07/01',vals=(0-spinup*12*spinrep):(413))
tre<-ncvar_def(name = "trend", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)
# contribution from the seasonal component
s1<-ncvar_def(name = "seasonal1", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)
s2<-ncvar_def(name = "seasonal2", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)
s3<-ncvar_def(name = "seasonal3", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)

# contribution from previous NDVI at different time scale.
y_p1<-ncvar_def(name = "pre_contri1", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)
y_p23<-ncvar_def(name = "pre_contri23", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)
y_p46<-ncvar_def(name = "pre_contri46", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)
y_p712<-ncvar_def(name = "pre_contri712", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)
y_p1324<-ncvar_def(name = "pre_contri1324", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)

# contribution from precipitation, temperature and radiation.
y_13<-ncvar_def(name = "pre_contri13", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)
y_t<-ncvar_def(name = "tmp_contri", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)
y_r<-ncvar_def(name = "rad_contri", units = 'NA (NDVI)', dim = list(lon,lat,tim),missval = -9999,compression = 9)

outfile<-paste0("./analysis/gimms/gimms.component_pred.",scenario,".nc")
if(file.exists(outfile)){
  file.remove(outfile)
}

ncout<-nc_create(outfile,vars = list(tre,s1,s2,s3,y_p1,y_p23,y_p46,y_p712,y_p1324,y_13,y_t,y_r))

ncvar_put(ncout,varid=tre,vals=trend)
ncvar_put(ncout,varid=s1,vals=seas1)
ncvar_put(ncout,varid=s2,vals=seas2)
ncvar_put(ncout,varid=s3,vals=seas3)

ncvar_put(ncout,varid=y_p1,vals=y_pre1)
ncvar_put(ncout,varid=y_p23,vals=y_pre23)
ncvar_put(ncout,varid=y_p46,vals=y_pre46)
ncvar_put(ncout,varid=y_p712,vals=y_pre712)
ncvar_put(ncout,varid=y_p1324,vals=y_pre1324)

ncvar_put(ncout,varid=y_13,vals=y_pre13)
ncvar_put(ncout,varid=y_t,vals=y_temp)
ncvar_put(ncout,varid=y_r,vals=y_solr)

nc_close(ncout)


