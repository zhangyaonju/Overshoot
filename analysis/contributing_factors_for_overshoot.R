#### compare the fraction of numbers of events and fraction of variations.

###  factors to be considered
#   timing of events.
#   temperature percentile


#   soil property
#   aridity
#   vegetation types
#   surface mositure residence time  (SMAP)
#   total water residence time    (GRACE)
#   mean VI
#   magnitude of VI
#   total precipitation
#   percipitation seasonality
#   MAT
#   maginitude of Ta
#   growing season length
#   

### prepare data

library(ncdf4)
#library(raster)
setwd("/global/scratch/yaozhang/Project/overshooting/")

#setwd("~/Documents/Project/overshooting/")
#load("~/Dropbox/YAOZHANG/code/R-code/tools/RCOLOR.RData")
#load("./analysis/gimms_overshoot/overshoot_number_conitrbution.RData")

#num_percent<-num_overshoot/total_drought
#comp_percent<-total_overshoot_compare[2,,]/total_overshoot_compare[1,,]

mean_magni<-function(vec){
  n_years<-length(vec)/12
  dim(vec)<-c(12,n_years)
  annual<-apply(vec,2,mean,na.rm=T)
  iav<-sd(annual,na.rm=T)
  msc<-apply(vec,1,mean,na.rm=T)
  if (sum(is.nan(msc))>10){
    return(rep(NA,4))
  }
  maxmon<-which.max(msc)
  mean_val<-mean(msc,na.rm=T)
  mag<-sd(msc,na.rm=T)
  return(c(mean_val,mag,iav,maxmon))
}


var_file<-"./analysis/other_data/climate_variables.RData"

if(!file.exists(var_file)){
  
  ncin<-nc_open("./Data/CRU_TS_403_pre.nc")
  prec<-ncvar_get(ncin,varid="pre")
  print(dim(prec))
  nc_close(ncin)
  pre<-apply(prec[,,19:438],c(1,2),mean_magni)
  meaPrec<-pre[1,,]
  magPrec<-pre[2,,]
  iavPrec<-pre[3,,]
  maxPrec<-pre[4,,]
  
  ncin<-nc_open("./Data/CRU_TS_403_tmp.nc")
  temp<-ncvar_get(ncin,varid="tmp")
  print(dim(temp))
  nc_close(ncin)
  tmp<-apply(temp[,,19:438],c(1,2),mean_magni)
  meaTemp<-tmp[1,,]
  magTemp<-tmp[2,,]
  iavTemp<-tmp[3,,]
  maxTemp<-tmp[4,,]
  
  ncin<-nc_open("./Data/GIMMS3g_v1.mon.hd.growingseason.nc")
  ndvi<-ncvar_get(ncin,varid="ndvi")
  print(dim(ndvi))
  nc_close(ncin)
  vi<-apply(ndvi[,,7:414],c(1,2),mean_magni)
  meaVI<-vi[1,,]
  magVI<-vi[2,,]
  iavVI<-vi[3,,]
  
  
  ncin<-nc_open("./Data/MOD13C2_NDVI_2000_2018_growingseason.nc")
  modisndvi<-ncvar_get(ncin,varid="ndvi")
  print(dim(modisndvi))
  nc_close(ncin)
  modisvi<-apply(modisndvi[,,12:227],c(1,2),mean_magni)
  modismeaVI<-modisvi[1,,]
  modismagVI<-modisvi[2,,]
  modisiavVI<-modisvi[3,,]
  
  ## vegetation types 
  lc_ras<-raster("./Data/IGBP_HD_2001.tif")
  lc<-t(as.matrix(lc_ras))[,360:1]
  #image(lc)
  
  ## drought recovery
  dat<-read.csv("./Data/41586_2017_BFnature23021_MOESM5_ESM.csv",stringsAsFactors = F)
  droughtRec<-t(as.matrix(dat))[,360:1]
  #image(drought_recover)
  
  ## glace residence time
  ncin<-nc_open("~/Documents/Data/GRACE_REC/GRACE_REC_v03_JPL_ERA5_residence_time.nc")
  graceTau<-ncvar_get(ncin,varid="median_residence_time")
  model_performance<-ncvar_get(ncin,varid="model_skill")
  graceTau[model_performance<0.5]<-NA
  
  ## smap residence time
  library(R.matlab)
  dat<-readMat("./Data/tau.mat")
  #image(dat$count.all)
  dat_com<-cbind(as.vector(dat$lon),as.vector(dat$lat),as.vector(dat$tau.median))
  e <- extent(c(-180,180,-90,90))
  r <- raster(e, ncol=720, nrow=360)
  x <- rasterize(dat_com[, 1:2], r, dat_com[,3], fun=mean)
  
  smapTau<-t(as.matrix(x))[,360:1]
  
  ## aridity index
  load("./Data/ai_HD.RData")
  ai<-ai_dat/1000
  ai_ras<-raster(t(ai[,360:1]))
  extent(ai_ras)<-c(-180,180,-90,90)
  
  ## elevation
  ncin<-nc_open("./Data/elev.0.5-deg.nc")
  ele<-ncvar_get(ncin,varid="data")[c(361:720,1:360),360:1]
  
  ### biodiversity
  biodiver<-shapefile("~/Documents/Data/plant_diversity_shapefile/ellis_2012_l8_dataset_2012_01_17.shp")
  biodras<-rasterize(biodiver,ai_ras,field="N")
  biod<-t(as.matrix(biodras))[,360:1]
  
  ## get growing season length
  lgs<-raster("~/Documents/Data/phenology/lgs_hd.tif")
  lgs[lgs==-1]<-NA
  lgsmat<-t(as.matrix(lgs)[360:1,])
  
  ### asynchronology
  #get the maximum difference between climate maximum and temperature maximum
  # get multi year average temperature and precipitation/
  asycho<-maxTemp-maxPrec
  asycho<-ifelse(asycho<= -6,asycho+12,asycho)
  asycho<-ifelse(asycho> 6,asycho-12,asycho)
  
  ### get the soil clay sand silt
  ncin<-nc_open("./Data/soil_data_HD.nc")
  dat<-ncvar_get(ncin,varid = "dat")
  pct_clay<-(dat[,,1]+dat[,,4])/2
  pct_sand<-(dat[,,2]+dat[,,5])/2
  pct_silt<-(dat[,,3]+dat[,,6])/2
  
  save(meaPrec,magPrec,iavPrec,meaTemp,magTemp,iavTemp,ai,maxTemp,maxPrec,asycho,ele,
       meaVI,magVI,iavVI,modismeaVI,modismagVI,modisiavVI,
       lc,droughtRec,graceTau,smapTau,biod,lgsmat,pct_clay,pct_sand,pct_silt,
       file=var_file)
  
}



# prepare data -----------------------------------------------------------
library(randomForest)
library(ranger)
library(dplyr)
#setwd("~/Documents/Project/overshooting/")
load("./analysis/other_data/climate_variables.RData")
msk<-meaVI/10000>0.15
msk[msk==0]<-NA
load("./analysis/gimms_overshoot_prec_temp/prec_temp_overshoot_num_impact_compare_1_414_spei05.RData")

set_range<-function(mat,mn,mx){
  mat[mat<mn]<-NA
  mat[mat>mx]<-NA
  return(mat)
}
num_frac<-set_range(num_overshoot/num_drought,0,1)   #num_overshoot/num_drought#
imp_frac<-set_range(overshoot_drought_ndvi/total_drought_ndvi_decrease,0,1)       #overshoot_drought_ndvi/total_drought_ndvi_decrease#
com_frac<-set_range(overshoot_component_ndvi_sum/overshoot_drought_ndvi,0,2)  #overshoot_component_ndvi_sum/predicted_drought_ndvi_sum#

data<-data.frame(as.vector(num_frac*msk),as.vector(imp_frac*msk),as.vector(com_frac*msk),
                 as.vector(num_overshoot*msk),as.vector(overshoot_drought_ndvi*msk),as.vector(overshoot_component_ndvi_sum*msk),
                 as.vector(num_drought*msk),as.vector(total_drought_ndvi_decrease*msk),as.vector(predicted_drought_ndvi_sum*msk),
                 as.vector(meaPrec*msk),as.vector(meaTemp*msk),as.vector(meaVI/10000*msk),as.vector(iavPrec*msk),as.vector(iavTemp*msk),
                 as.vector(iavVI/10000*msk),as.factor(lc*msk),as.vector(droughtRec*msk),as.vector(ele*msk),as.vector(log(graceTau*msk)),
                 as.vector(biod*msk),as.vector(lgsmat*msk),as.vector(abs(asycho*msk)),as.vector(ai*msk),
                 as.vector(pct_clay*msk),as.vector(pct_sand*msk),as.vector(pct_silt*msk),as.vector(magPrec*msk),as.vector(magTemp*msk),as.vector(magVI*msk))

var_names<-c("num_frac","imp_frac","com_frac",
               "overshoot_num","overshoot_ndvi","overshoot_comp",
               "drought_num","drought_ndvi","modelpred",
               "meanPrec","meanTemp","meanVI","iavPrec","iavTemp","iavVI","lc","droughtRec","ele","graceTau",
               "biodiversity","lgs","asychro","ai","pct_clay","pct_sand","pct_silt","magPrec","magTemp","magVI")

names(data)<-var_names

save(data, file="./analysis/other_data/overshoot_with_climate_variables.RData")

get_ranges<-function(vec){
  quantile(vec,c(0.025,0.05,0.95,0.975),na.rm=T)
}
outlier<-function(dat){
  for (i in 1:dim(dat)[2]){
    if (is.factor(dat[1,i])){
      next
    }
    sdvar<-sd(dat[,i],na.rm=T)
    medianval<-median(dat[,i],na.rm=T)
    dat[(dat[,i]<medianval-sdvar*3)|(dat[,i]>medianval+sdvar*3),i]<-NA
  }
  return(dat)
}


# for number of overshoot (use all cases) ------------------------------------------------------
load("./analysis/other_data/overshoot_with_climate_variables.RData")
var_names<-c("num_frac","imp_frac","com_frac",
             "overshoot_num","overshoot_ndvi","overshoot_comp",
             "drought_num","drought_ndvi","modelpred",
             "meanPrec","meanTemp","meanVI","iavPrec","iavTemp","iavVI","lc","droughtRec","ele","graceTau",
             "biodiversity","lgs","asychro","ai","pct_clay","pct_sand","pct_silt","magPrec","magTemp","magVI")
num_data<-subset(data,select=var_names[c(1,4,7:29)])
# num_data[,2]<-log(num_data[,2],10)
# num_data[is.infinite(num_data[,2]),2]<-NA
# num_data[,3]<-log(num_data[,3],10)
# num_data[is.infinite(num_data[,3]),3]<-NA
dat<-num_data[complete.cases(num_data),]
#dim(dat)


## get the ranges
range_stat<-as.data.frame(array(NA,dim=c(length(var_names),5)))
names(range_stat)<-c("var_name","p1",'p99',"p2_5","p97_5")
range_stat$var_name<-var_names
for (i in 1:length(var_names)){
  if(length(unique(data[,i]))<20){
    if (is.factor(data[,i])){
      range_stat[i,2:5]<-as.numeric(rep(range(levels(data[,i])),2))
    }else{
      range_stat[i,2:5]<-as.numeric(rep(range(data[,i],na.rm=T),2))
    }
  }else{
    range_stat[i,2:5]<-quantile(data[,i],c(0.01,0.99,0.025,0.975),na.rm=T)
  }
}
save(range_stat,file="./analysis/other_data/Random_forest_variables_ranges.RData")



# overshoot number  --------------------------------------------------------
set.seed(123)
#tune.res <- tuneRF(dat[,c(4,7:13,16:21)], dat[,1], stepFactor=1.5)
load("./analysis/other_data/overshoot_with_climate_variables.RData")
num_data<-subset(data,select=var_names[c(1,3,4,7:29)])
# num_data[,2]<-log(num_data[,2],10)
# num_data[is.infinite(num_data[,2]),2]<-NA
# num_data[,3]<-log(num_data[,3],10)
# num_data[is.infinite(num_data[,3]),3]<-NA
dat<-num_data[complete.cases(num_data),]
#dim(dat)

num_rf<-randomForest(overshoot_num~drought_num+meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+graceTau+biodiversity+lgs+asychro+ai+pct_clay,#+ele,
                     data=dat,
                     mtry=4,
                     min.node.size=5,
                     sample.fraction=0.8,
                     ntree=500)

save(num_rf,file="./analysis/randomforest/trained_num_rf.RData")


library(pdp)
load("./analysis/randomforest/trained_num_rf.RData")
varresponse<-list()
varImpPlot(num_rf)
var_importance<-data.frame(attributes(num_rf$importance)$dimnames[[1]],as.numeric(num_rf$importance))
names(var_importance)<-c("var","importance")
# decending
var_importance<-var_importance[order(var_importance[,2],decreasing = T),]

for (i in 1:dim(var_importance)[1]){
  # partplt<-partialPlot(num_rf,pred.data = dat, x.var = as.character(var_importance[i,1]), 
  #                      plot = F, rug = F)
  partplt<-partial(num_rf, pred.var = as.character(var_importance[i,1]), 
                   plot = F, rug = F)
  varresponse[[i]]<-partplt
}

save(varresponse,var_importance,file="./analysis/randomforest/num_rf_var_response.RData")



# overshoot num fraction --------------------------------------------------
set.seed(456)
#tune.res <- tuneRF(dat[,c(4,7:13,16:21)], dat[,1], stepFactor=1.5)
var_names<-c("num_frac","imp_frac","com_frac",
             "overshoot_num","overshoot_ndvi","overshoot_comp",
             "drought_num","drought_ndvi","modelpred",
             "meanPrec","meanTemp","meanVI","iavPrec","iavTemp","iavVI","lc","droughtRec","ele","graceTau",
             "biodiversity","lgs","asychro","ai","pct_clay","pct_sand","pct_silt","magPrec","magTemp","magVI")
num_data<-subset(data,select=var_names[c(1,3,4,7:29)])
# num_data[,2]<-log(num_data[,2],10)
# num_data[is.infinite(num_data[,2]),2]<-NA
# num_data[,3]<-log(num_data[,3],10)
# num_data[is.infinite(num_data[,3]),3]<-NA
num_frac_data<-num_data[complete.cases(num_data),]
num_rf<-randomForest(num_frac~meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+graceTau+biodiversity+lgs+asychro+ai+pct_clay,#+ele,
                     data=num_frac_data,
                     mtry=4,
                     min.node.size=5,
                     sample.fraction=0.8,
                     ntree=500)

save(num_rf,file="./analysis/randomforest/trained_num_frac_rf.RData")


library(pdp)
load("./analysis/randomforest/trained_num_frac_rf.RData")
varresponse<-list()
varImpPlot(num_rf)
var_importance<-data.frame(attributes(num_rf$importance)$dimnames[[1]],as.numeric(num_rf$importance))
names(var_importance)<-c("var","importance")
# decending
var_importance<-var_importance[order(var_importance[,2],decreasing = T),]

for (i in 1:dim(var_importance)[1]){
  # partplt<-partialPlot(num_rf,pred.data = dat, x.var = as.character(var_importance[i,1]), 
  #                      plot = F, rug = F)
  partplt<-partial(num_rf, pred.var = as.character(var_importance[i,1]), 
                   plot = F, rug = F)
  varresponse[[i]]<-partplt
}

save(varresponse,var_importance,file="./analysis/randomforest/num_frac_rf_var_response.RData")






# overshoot comp to overshoot total ---------------------------------------
var_names<-c("num_frac","imp_frac","com_frac",
             "overshoot_num","overshoot_ndvi","overshoot_comp",
             "drought_num","drought_ndvi","modelpred",
             "meanPrec","meanTemp","meanVI","iavPrec","iavTemp","iavVI","lc","droughtRec","ele","graceTau",
             "biodiversity","lgs","asychro","ai","pct_clay","pct_sand","pct_silt","magPrec","magTemp","magVI")
load("./analysis/other_data/overshoot_with_climate_variables.RData")
imp_data<-subset(data,select=var_names[c(3,5,6,8:29)])
imp_frac_data<-imp_data[complete.cases(imp_data),]

set.seed(746)
impact_rf<-randomForest(overshoot_comp~overshoot_ndvi+meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+graceTau+biodiversity+lgs+asychro+ai+pct_clay,
                        data=imp_frac_data,
                        mtry=4,
                        min.node.size=5,
                        sample.fraction=0.8,
                        ntree=500)

save(impact_rf,file="./analysis/randomforest/trained_impact_overshoot_rf.RData")


library(pdp)
load("./analysis/randomforest/trained_impact_overshoot_rf.RData")
varresponse<-list()
varImpPlot(impact_rf)
var_importance<-data.frame(attributes(impact_rf$importance)$dimnames[[1]],as.numeric(impact_rf$importance))
names(var_importance)<-c("var","importance")
# decending
var_importance<-var_importance[order(var_importance[,2],decreasing = T),]

for (i in 1:dim(var_importance)[1]){
  # partplt<-partialPlot(num_rf,pred.data = dat, x.var = as.character(var_importance[i,1]), 
  #                      plot = F, rug = F)
  partplt<-partial(impact_rf, pred.var = as.character(var_importance[i,1]), 
                   plot = F, rug = F)
  varresponse[[i]]<-partplt
}

save(varresponse,var_importance,file="./analysis/randomforest/impact_overshoot_rf_var_response.RData")






# overshoot comp to overshoot fraction ---------------------------------------
library(randomForest)
library(ranger)
var_names<-c("num_frac","imp_frac","com_frac",
             "overshoot_num","overshoot_ndvi","overshoot_comp",
             "drought_num","drought_ndvi","modelpred",
             "meanPrec","meanTemp","meanVI","iavPrec","iavTemp","iavVI","lc","droughtRec","ele","graceTau",
             "biodiversity","lgs","asychro","ai","pct_clay","pct_sand","pct_silt","magPrec","magTemp","magVI")
load("./analysis/other_data/overshoot_with_climate_variables.RData")
imp_data<-subset(data,select=var_names[c(3,5,6,8:29)])
imp_frac_data<-imp_data[complete.cases(imp_data),]

set.seed(3231)
impact_rf<-randomForest(com_frac~meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+graceTau+biodiversity+lgs+asychro+ai+pct_clay,
                        data=imp_frac_data,
                        mtry=4,
                        min.node.size=5,
                        sample.fraction=0.8,
                        ntree=500)

save(impact_rf,file="./analysis/randomforest/trained_impact_frac_overshoot_rf.RData")


library(pdp)
load("./analysis/randomforest/trained_impact_frac_overshoot_rf.RData")
varresponse<-list()
varImpPlot(impact_rf)
var_importance<-data.frame(attributes(impact_rf$importance)$dimnames[[1]],as.numeric(impact_rf$importance))
names(var_importance)<-c("var","importance")
# decending
var_importance<-var_importance[order(var_importance[,2],decreasing = T),]

for (i in 1:dim(var_importance)[1]){
  # partplt<-partialPlot(num_rf,pred.data = dat, x.var = as.character(var_importance[i,1]), 
  #                      plot = F, rug = F)
  partplt<-partial(impact_rf, pred.var = as.character(var_importance[i,1]), 
                   plot = F, rug = F)
  varresponse[[i]]<-partplt
}

save(varresponse,var_importance,file="./analysis/randomforest/impact_overshoot_frac_rf_var_response.RData")





#stopCluster(cl)


# pltdat<-var1$panel.args[[1]]
# plot(pltdat$x,pltdat$y)
# # set.seed(133)
# # num_rf<-randomForest(overshoot_num~drought_num+droughtRec+meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+droughtRec+graceTau+biodiversity+lgs+asychro+ai,
# #                      data=dat,ntree=2000)
# library(LSD)
# plot(num_rf)
# varImpPlot(num_rf)
# pred.f<-predict(num_rf,dat,type="response")
# heatscatter(pred.f,dat$overshoot_num)
# 
# lm_reg<-lm(dat$overshoot_num~pred.f)
# summary(lm_reg)
# abline(0,1)
# 
# frac_rf<-randomForest(num_frac~droughtRec+meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+droughtRec+graceTau+biodiversity+lgs+asychro+ai,
#                      data=dat)
# 
# plot(frac_rf)
# varImpPlot(frac_rf)
# pred.f<-predict(frac_rf,dat,type="response")
# heatscatter(pred.f,dat$num_frac)
# 
# lm_reg<-lm(dat$num_frac~pred.f)
# summary(lm_reg)
# abline(0,1)
# 
# partialPlot(frac_rf,pred.data=dat,x.var=c(1))
# 
# print(frac_rf)
# heatscatter(frac_rf$predicted,dat$num_frac)
# 
# lm_reg<-lm(dat$num_frac~frac_rf$predicted)
# abline(0,1)
# abline(lm_reg)
# 
# varUsed(frac_rf)
# 
# pdf("~/Dropbox/YAOZHANG/paper/2020_overshooting/figures/Figure1_num_frac_importance.pdf",width=5,height=15)
# contribution<-summary(num_gbm,n.trees=best.iter)
# dev.off()
# best.iter <- gbm.perf(num_gbm,method="cv")
# contribution<-summary(num_gbm,n.trees=best.iter)
# print(contribution)
# f.predict <- predict(num_gbm,dat,best.iter)
# lm_reg<-lm(dat[,1]~f.predict)
# summary(lm_reg)
# plot(f.predict,dat[,1])
# dim(dat)
# 
# library(LSD)
# pdf("~/Dropbox/YAOZHANG/paper/2020_overshooting/figures/Figure1_num_frac_model_perform.pdf",width=5,height=5)
# heatscatter(f.predict,dat[,1],xlim=c(0,1),ylim=c(0,1))
# abline(0,1)
# dev.off()
# #cor.test(dat[,1],f.predict)
# pdf("~/Dropbox/YAOZHANG/paper/2020_overshooting/figures/Figure_variation.pdf",width=5,height=5)
# par(fig=c(0,1,0,0.5))
# for(i in 1:14){
#   plot.gbm(num_gbm, i.var = 2)
# }
# dev.off()
# 
# pretty.gbm.tree(num_gbm)
# 
# 
# 
# imp_gbm<-gbm(imp_frac~meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+droughtRec+graceTau+smapTau+biodiversity+lgs+asychro+ai,
#              data=dat,
#              var.monotone=rep(0,14),#,0,0),
#              distribution="gaussian",
#              n.trees=2000,
#              shrinkage=0.05,
#              interaction.depth=2,
#              train.fraction = 0.8,
#              bag.fraction=0.5,
#              cv.folds = 10,
#              n.minobsinnode = 5,
#              n.cores=6)
# best.iter <- gbm.perf(imp_gbm,method="test")
# contribution<-summary(imp_gbm,n.trees=best.iter)
# print(contribution)
# f.predict <- predict(imp_gbm,dat[,-(1:3)],best.iter)
# lm_reg<-lm(dat[,2]~f.predict)
# summary(lm_reg)
# plot(f.predict,dat[,2])
# 
# 
# com_gbm<-gbm(com_frac~meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+droughtRec+graceTau+smapTau+biodiversity+lgs+asychro+ai,
#              data=dat,
#              var.monotone=rep(0,14),#,0,0),
#              distribution="gaussian",
#              n.trees=2000,
#              shrinkage=0.05,
#              interaction.depth=2,
#              train.fraction = 0.8,
#              bag.fraction=0.5,
#              cv.folds = 10,
#              n.minobsinnode = 5,
#              n.cores=6)
# best.iter <- gbm.perf(com_gbm,method="cv")
# contribution<-summary(com_gbm,n.trees=best.iter)
# print(contribution)
# f.predict <- predict(com_gbm,dat,best.iter)
# lm_reg<-lm(dat[,3]~f.predict)
# summary(lm_reg)
# plot(f.predict,dat[,3])
# 
# 
# BRT analysis ------------------------------------------------------------
# library("gbm")
# library("dismo")
# setwd("~/Documents/Project/overshooting/")
# load("./analysis/other_data/overshoot_with_climate_variables.RData")
# load("~/Dropbox/YAOZHANG/code/R-code/tools/RCOLOR.RData")
# 
# 
# outlier<-function(dat){
#   for (i in 1:dim(dat)[2]){
#     if (is.factor(dat[,i])){
#       next
#     }
#     sdvar<-sd(dat[,i],na.rm=T)
#     medianval<-median(dat[,i],na.rm=T)
#     dat[(dat[,i]<medianval-sdvar*3)|(dat[,i]>medianval+sdvar*3),i]<-NA
#   }
#   return(dat)
# }
# 
# dat<-subset(data,select=var_names[c(1,4,7:29)])
# dat[,2]<-log(dat[,2],10)
# dat[,3]<-log(dat[,3],10)
# dat<-dat[complete.cases(dat),]
# # dat<-outlier(dat)
# # dat<-dat[complete.cases(dat),]
# 
# num_gbm<-gbm(overshoot_num~drought_num+meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+droughtRec+graceTau+biodiversity+lgs+asychro+ai,
#             data=dat,
#             verbose = T,
#             var.monotone=rep(0,14),#,0,0),
#             distribution="poisson",
#             n.trees=1000,
#             shrinkage=0.02,
#             interaction.depth=2,
#             train.fraction = 0.9,
#             bag.fraction = 0.5,
#             cv.folds = 10,
#             n.minobsinnode = 2,
#             n.cores=4)
# 
# best.iter <- gbm.perf(num_gbm,method="cv")
# contribution<-summary(num_gbm,n.trees=best.iter)
# print(contribution)
# f.predict <- predict(num_gbm,dat,best.iter,type="response")
# lm_reg<-lm(dat$overshoot_num~f.predict)
# summary(lm_reg)
# dev.off()
# heatscatter(f.predict,dat$overshoot_num)
# abline(0,1)
# 


# 
# 
# # predict overshoot impact  -----------------------------------------------
# set.seed(123)
# # hyper_grid <- expand.grid(
# #   mtry       = seq(4, 10, by = 2),
# #   node_size  = seq(3, 9, by = 2),
# #   sampe_size = c(.55, .632, .70, .80),
# #   OOB_RMSE   = 0
# # )
# # 
# # for(i in 1:nrow(hyper_grid)) {
# #   
# #   # train model
# #   model <- ranger(
# #     formula         = overshoot_comp ~ drought_ndvi+ meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+droughtRec+graceTau+biodiversity+lgs+asychro+ai, 
# #     data            = dat, 
# #     num.trees       = 500,
# #     mtry            = hyper_grid$mtry[i],
# #     min.node.size   = hyper_grid$node_size[i],
# #     sample.fraction = hyper_grid$sampe_size[i],
# #     seed            = 123
# #   )
# #   
# #   # add OOB error to grid
# #   hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
# # }
# # 
# # hyper_grid %>% 
# #   dplyr::arrange(OOB_RMSE) %>%
# #   head(10)
# # 
# #graceTau+
# impact_rf<-randomForest(overshoot_comp~drought_ndvi+meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+graceTau+biodiversity+lgs+asychro+ai+pct_clay,
#                         data=dat,
#                         mtry=4,
#                         min.node.size=5,
#                         sample.fraction=0.8,
#                         ntree=500)
# 
# save(impact_rf,file="./analysis/trained_impact_rf.RData")
# 
# 
# library(pdp)
# load("./analysis/trained_impact_rf.RData")
# varresponse<-list()
# varImpPlot(impact_rf)
# var_importance<-data.frame(attributes(impact_rf$importance)$dimnames[[1]],as.numeric(impact_rf$importance))
# names(var_importance)<-c("var","importance")
# # decending
# var_importance<-var_importance[order(var_importance[,2],decreasing = T),]
# 
# for (i in 1:dim(var_importance)[1]){
#   # partplt<-partialPlot(num_rf,pred.data = dat, x.var = as.character(var_importance[i,1]), 
#   #                      plot = F, rug = F)
#   partplt<-partial(impact_rf, pred.var = as.character(var_importance[i,1]), 
#                    plot = F, rug = F)
#   varresponse[[i]]<-partplt
# }
# 
# save(varresponse,var_importance,file="./analysis/impact_rf_var_response.RData")
# 
# pred.f<-predict(impact_rf,dat,type="response")
# heatscatter(pred.f,dat$overshoot_comp)
# lmr<-lm(pred.f~dat$overshoot_comp)
# summary(lmr)
# 
# hist(dat$overshoot_ndvi)
# image(overshoot_drought_ndvi>0)
# # 

#dat<-outlier(dat)
#dat<-dat[complete.cases(dat),]



# variable selection ------------------------------------------------------
# image(asycho)
# 
# library(caret)
# library(doParallel)
# cl <- makePSOCKcluster(10)
# #stopCluster(cl)
# registerDoParallel(cl)
# # load the data
# # define the control using a random forest selection function
# control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# # run the RFE algorithm
# subsets <- c(1:17)
# results <- rfe(dat[,c(7,10:25)], dat[,4], sizes=subsets, rfeControl=control)
# # summarize the results
# print(results)
# # list the chosen features
# predictors(results)

#stopCluster(cl)
# 
# hyper_grid <- expand.grid(
#   mtry       = seq(4, 10, by = 2),
#   node_size  = seq(3, 9, by = 2),
#   sampe_size = c(.55, .632, .70, .80),
#   OOB_RMSE   = 0
# )
# 
# for(i in 1:nrow(hyper_grid)) {
#   
#   # train model
#   model <- ranger(
#     formula         = overshoot_num ~ drought_num+ meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+droughtRec+graceTau+biodiversity+lgs+asychro+ai, 
#     data            = dat, 
#     num.trees       = 500,
#     mtry            = hyper_grid$mtry[i],
#     min.node.size   = hyper_grid$node_size[i],
#     sample.fraction = hyper_grid$sampe_size[i],
#     seed            = 123
#   )
#   
#   # add OOB error to grid
#   hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
# }
# 
# hyper_grid %>% 
#   dplyr::arrange(OOB_RMSE) %>%
#   head(10)
#graceTau+
# 
# hyper_grid <- expand.grid(
#   mtry       = seq(4, 10, by = 2),
#   node_size  = seq(3, 9, by = 2),
#   sampe_size = c(.55, .70, .80, 0.9),
#   OOB_RMSE   = 0
# )
# 
# for(i in 1:nrow(hyper_grid)) {
#   
#   # train model
#   model <- ranger(
#     formula         = overshoot_num ~ drought_num+ meanPrec+meanTemp+meanVI+iavPrec+iavTemp+iavVI+lc+droughtRec+graceTau+biodiversity+lgs+asychro+ai,
#     data            = dat,
#     num.trees       = 500,
#     mtry            = hyper_grid$mtry[i],
#     min.node.size   = hyper_grid$node_size[i],
#     sample.fraction = hyper_grid$sampe_size[i],
#     seed            = 123
#   )
#   
#   # add OOB error to grid
#   hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
# }
# 
# hyper_grid %>%
#   dplyr::arrange(OOB_RMSE) %>%
#   head(10)