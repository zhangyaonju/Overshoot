args = commandArgs(trailingOnly=F)
print(args)
# args=c("-5","-2","-1","-2")
myargs1 <-sub('-','',args[length(args)-3])
myargs2 <-sub('-','',args[length(args)-2])
myargs3 <-sub('-','',args[length(args)-1])
myargs4 <-sub('-','',args[length(args)])

spinup<-as.numeric(myargs1)
spinrep<-as.numeric(myargs2)
delta<-c(0.95,0.98,0.99,0.995,0.999)[as.numeric(myargs3)+1]
indcators<-c('prec',"prec_temp",'randomndvi')
ind_rand<-as.numeric(myargs4)

if (ind_rand>=3){
  ind <-3
  rand<-(ind_rand-3)%%5+1
  window_ind<-ceiling((ind_rand-2)/5)-1
  indcator<-indcators[ind]
  out_label<-paste0(indcator,window_ind,rand)
}else{
  indcator<-indcators[ind_rand]
  out_label<-indcator
}


library(raster)
library(ncdf4)
library(abind)
library(zoo)
library(parallel)

setwd("/global/scratch/yaozhang/Project/overshooting/")
#setwd("~/Documents/Project/overshooting/")

mean_ind<-c('trend','seasonal1','seasonal2','seasonal3',
            'pre_contri1','pre_contri23','pre_contri46','pre_contri712','pre_contri1324','pre_contri13')

cl <- makeCluster(8)
clusterEvalQ(cl, {
  library(zoo)
})


ind_data_for_figure1<-T
ind_for_annual_maps<-F
ind_test_for_difference<-T

ind_speed_for_drought_dev<-T
ind_test_drought_impact_anova<-T
ind_test_drought_impact_anova_low<-T
ind_fit_drought_timing<-T
ind_os_nos_timing<-T
ind_os_temeprature<-T

if (as.numeric(myargs4)>=3){
  ind_for_annual_maps<-F
  ind_test_for_difference<-F
  
  ind_speed_for_drought_dev<-F
  ind_test_drought_impact_anova<-F
  ind_test_drought_impact_anova_low<-F
  ind_fit_drought_timing<-F
  ind_os_nos_timing<-F
  ind_os_temeprature<-F
}

drought_f<-paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/drought_with_spei_",out_label,"_spei05.RData")


# prepare data for first figure ---------------------------------------
##### dataset for figure 1. 6 maps, statistics along latitude for 4 legacy periods.
stat_for_overshoot_maps<-function(t_range){
  ## 1. total drought numbers.
  drought_f<-paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/drought_with_spei_",out_label,"_spei05.RData")
  load(drought_f)
  get_drought_number<-function(drought_event){
    num_drought_id<-unique(drought_event,na.rm=T)[-1]
    num_drought<-length(num_drought_id)
    return(num_drought)
  }
  num_drought<-parApply(cl,drought_all[t_range,,],c(2,3),get_drought_number)
  
  ## 2. total overshoot numbers
  file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_",
              t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
  load(file)
  num_overshoot<-parApply(cl,overshoot_boolean,c(2,3),sum,na.rm=T)
  #prepare data along lat
  # get the overshoot_comp_boolean and 
  dim(overshoot_comp_boolean)<-c(50,4,720,360)
  overshoot_comp_num<-parApply(cl,overshoot_comp_boolean,c(2,3,4),sum,na.rm=T)
  lat_stat_overshoot_drought_num<-array(NA,dim=c(90,5))
  d2num_drought<-num_drought
  dim(d2num_drought)<-c(2880,90)
  lat_num_drought<-parApply(cl,d2num_drought,2,sum,na.rm=T)
  
  d2num_overshoot_comp_num<-overshoot_comp_num
  dim(d2num_overshoot_comp_num)<-c(4,2880,90)
  lat_stat_overshoot_drought_num[,1:4]<-t(apply(d2num_overshoot_comp_num,c(1,3),sum,na.rm=T))/lat_num_drought
  d2_os_num<-num_overshoot
  dim(d2_os_num)<-c(2880,90)
  lat_os_num<-apply(d2_os_num,2,sum,na.rm=T)
  lat_stat_overshoot_drought_num[,5]<-lat_os_num/lat_num_drought
  
  ## 3. total drought caused NDVI decrease
  file2= paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_drought_impact_method2_",
                t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
  load(file2)
  total_drought_ndvi_decrease<-parApply(cl,sum_drought_ndvi,c(2,3),sum,na.rm=T)
  
  ## 4. total overshoot event caused NDVI decrease
  file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_",
              t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
  load(file)
  file2= paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_drought_impact_method2_",
                t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
  load(file2)
  dat<-abind(sum_drought_ndvi,overshoot_boolean,along=1)
  get_overshoot_ndvi_decline<-function(vec){
    n_obs<-length(vec)/2
    sum_drought_ndvi<-vec[1:n_obs]
    overshoot_boolean<-vec[(1+n_obs):(n_obs*2)]
    overshoot_drought_ndvi<-sum(sum_drought_ndvi[overshoot_boolean==1],na.rm=T)
    return(overshoot_drought_ndvi)
  }
  overshoot_drought_ndvi<-parApply(cl,dat,c(2,3),get_overshoot_ndvi_decline)
  
  d2_drought_ndvi<-total_drought_ndvi_decrease
  dim(d2_drought_ndvi)<-c(720*4,90)
  lat_drought_ndvi_impact<-apply(d2_drought_ndvi,2,sum,na.rm=T)
  
  d2_overshoot_ndvi<-overshoot_drought_ndvi
  dim(d2_overshoot_ndvi)<-c(720*4,90)
  lat_overshoot_ndvi_impact<-apply(d2_overshoot_ndvi,2,sum,na.rm=T)
  
  lat_stat_overshoot_drought_impact<-lat_overshoot_ndvi_impact/lat_drought_ndvi_impact
  
  ## 5. total model predicted NDVI decrease for the overshoot events
  file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_",
              t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
  load(file)
  comp_drought_ndvi_list<-list()
  var_id<-c("1","23","46","712","1324","13")
  for (i in 2:5){
    file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_pre",var_id[i],"_overshoot_method2_",
                t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
    load(file)
    comp_drought_ndvi_list[[i-1]]<-aperm(comp_drought_ndvi,c(2,3,1))
  }
  comp_drought_ndvi<-do.call(abind, comp_drought_ndvi_list)
  
  overshoot_comp_boolean_aperm<-aperm(overshoot_comp_boolean,c(2,3,1))
  dat<-abind(comp_drought_ndvi,overshoot_comp_boolean_aperm,along=3)
  
  sum_drought_pred_ndvi_aperm<-aperm(sum_drought_pred_ndvi,c(2,3,1))
  dat_all<-abind(dat,sum_drought_pred_ndvi_aperm,along=3)
  
  compare_model_predict_with_overshoot<-function(vec){   ### update the overshoot prediction
    comp_drought_ndvi<-vec[1:200]
    overshoot_comp_boolean<-vec[201:400]
    drought_pred_ndvi<-vec[401:450]
    dim(comp_drought_ndvi)<-c(50,4)
    dim(overshoot_comp_boolean)<-c(50,4)
    overshoot_boolean<-apply(overshoot_comp_boolean,1,sum,na.rm=T)
    
    drought_pred_ndvi[overshoot_boolean==0]<-NA
    comp_drought_ndvi_sum<-sum(drought_pred_ndvi,na.rm=T)
    comp_drought_ndvi[overshoot_comp_boolean==F]<-NA
    overshoot_component_ndvi_sum<-sum(comp_drought_ndvi,na.rm=T)
    return(c(comp_drought_ndvi_sum,overshoot_component_ndvi_sum))
  }
  res<-parApply(cl,dat_all,c(1,2), compare_model_predict_with_overshoot)
  predicted_drought_ndvi_sum<-res[1,,]
  ## 6. total overshoot component decrease during the overshoot events.
  overshoot_component_ndvi_sum<-res[2,,]
  
  ### get the fraction of overshoot component to total drought decline.
  overshoot_comp_drought<-function(vec){   ### update the overshoot prediction
    comp_drought_ndvi<-vec[1:200]
    overshoot_comp_boolean<-vec[201:400]
    dim(comp_drought_ndvi)<-c(50,4)
    dim(overshoot_comp_boolean)<-c(50,4)
    
    comp_drought_ndvi[overshoot_comp_boolean==F]<-NA
    comp_os_drought<-apply(comp_drought_ndvi,2,sum,na.rm=T)
    return(comp_os_drought)
  }
  os_comp_ndvi<-parApply(cl,dat,c(1,2), overshoot_comp_drought)
  
  total_ndvi_decrease<-total_drought_ndvi_decrease
  dim(total_ndvi_decrease)<-c(2880,90)
  lat_total_ndvi_decrease<-apply(total_ndvi_decrease,2,sum,na.rm=T)
  lat_stat_overshoot_comp_total_projected<-array(NA,dim=c(90,5))
  
  d2_os_comp_ndvi<-os_comp_ndvi
  dim(d2_os_comp_ndvi)<-c(4,2880,90)
  lat_stat_overshoot_comp_total_projected[,1:4]<-t(apply(d2_os_comp_ndvi,c(1,3),sum,na.rm=T))/lat_total_ndvi_decrease
  lat_stat_overshoot_comp_total_projected[,5]<-apply(lat_stat_overshoot_comp_total_projected[,1:4],1,sum,na.rm=T)
  
  #overshoot_drought_ndvi
  overshoot_ndvi_decrease<-overshoot_drought_ndvi
  dim(overshoot_ndvi_decrease)<-c(2880,90)
  lat_overshoot_ndvi_decrease<-apply(overshoot_ndvi_decrease,2,sum,na.rm=T)
  lat_stat_overshoot_comp_overshoot<-array(NA,dim=c(90,5))
  
  lat_stat_overshoot_comp_overshoot[,1:4]<-t(apply(d2_os_comp_ndvi,c(1,3),sum,na.rm=T))/lat_overshoot_ndvi_decrease
  lat_stat_overshoot_comp_overshoot[,5]<-apply(lat_stat_overshoot_comp_overshoot[,1:4],1,sum,na.rm=T)
  
  ##### get the pct decrease for 1. overshoot component to non_drought reference
  #  2. total drought induced ndvi decline to non-drought reference
  
  file1=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_drought_impact_method2_",
               t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
  load(file1)
  file2=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_",
               t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
  load(file2)
  dat_ref<-abind(sum_nondrought_ref_ndvi,overshoot_boolean,along=1)
  overshoot_ref<-apply(dat_ref,c(2,3), get_overshoot_ndvi_decline)
  drought_ref<-apply(sum_nondrought_ref_ndvi,c(2,3),sum,na.rm=T)
  all_drought_pct_decrease<-total_drought_ndvi_decrease/drought_ref
  overshoot_drought_pct_decrease<-overshoot_drought_ndvi/overshoot_ref
  overshoot_component_pct_decrease<-overshoot_component_ndvi_sum/overshoot_ref
  
  save(num_drought, num_overshoot, total_drought_ndvi_decrease, overshoot_drought_ndvi,
       overshoot_comp_num,os_comp_ndvi,
       predicted_drought_ndvi_sum, overshoot_component_ndvi_sum,
       lat_stat_overshoot_drought_num,lat_stat_overshoot_drought_impact,
       lat_stat_overshoot_comp_overshoot,lat_stat_overshoot_comp_total_projected,
       overshoot_component_pct_decrease,overshoot_drought_pct_decrease, all_drought_pct_decrease,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_num_impact_compare_",
                   t_range[1],"_",t_range[length(t_range)],"_spei05.RData"))
}

if(ind_data_for_figure1){
  stat_for_overshoot_maps(t_range = 1:414)
  
  stat_for_overshoot_maps(t_range = 7:210)
  stat_for_overshoot_maps(t_range = 211:414)
}


# Test overshoot impact --------------------------------------------

testsig<-function(vec){
  n_obs<-length(vec)/2
  withover<-vec[1:n_obs][vec[(n_obs+1):(n_obs*2)]==1]
  withoutover<-vec[1:n_obs][vec[(n_obs+1):(n_obs*2)]==0]
  if (sum(!is.na(withoutover))>=2&sum(!is.na(withover))>=2){
    if (sum(!is.na(unique(withover)))==1&sum(!is.na(unique(withoutover)))==1){
      return(c(mean(withover,na.rm=T),mean(withoutover,na.rm=T),0))
    }
    tst<-t.test(withover,withoutover,alternative = "less")
    return(c(tst$estimate,tst$p.value))
  }else{
    return(rep(NA,3))
  }
}

change_res<-function(mat){  # matsize = 50,720,360
  dim(mat)<-c(250,144,360)
  mat_aperm<-aperm(mat,c(1,3,2))
  dim(mat_aperm)<-c(1250,72,144)
  mat_lowres<-aperm(mat_aperm,c(1,3,2))
  return(mat_lowres)
}

if (ind_test_for_difference==T){
  #var_id<-c("1","23","46","712","1324","13")
  file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_drought_impact_method2_1_414_spei05.RData")
  load(file)
  
  load(paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_1_414_spei05.RData"))
  minimum_drought_spei[minimum_drought_spei==0]<-NA
  minimum_drought_spei[is.infinite(minimum_drought_spei)]<-NA
  sum_drought_spei[sum_drought_spei==0]<-NA
  
  ndvidat<-abind(sum_drought_ndvi,overshoot_boolean,along=1)
  test_ndvi<-parApply(cl,ndvidat,c(2,3),testsig)
  
  ndvimindat<-abind(minimum_drought_ndvi,overshoot_boolean,along=1)
  test_minndvi<-parApply(cl,ndvimindat,c(2,3),testsig)
  
  speidat<-abind(sum_drought_spei,overshoot_boolean,along=1)
  test_spei3<-parApply(cl,speidat,c(2,3),testsig)
  
  speimindat<-abind(minimum_drought_spei,overshoot_boolean,along=1)
  test_minspei3<-parApply(cl,speimindat,c(2,3),testsig)
  
  save(test_spei3,test_ndvi,test_minndvi,test_minspei3,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_test_overshoot_impact_method2_1_414_spei05.RData"))
  
  lowspeidat<-abind(change_res(sum_drought_spei),change_res(overshoot_boolean),along=1)
  test_low_spei3<-parApply(cl,lowspeidat,c(2,3),testsig)
  
  lowspeimindat<-abind(change_res(minimum_drought_spei),change_res(overshoot_boolean),along=1)
  test_low_minspei3<-parApply(cl,lowspeimindat,c(2,3),testsig)
  
  save(test_low_spei3,test_low_minspei3,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_test_low_overshoot_impact_method2_1_414_spei05.RData"))
  
}


# Development speed ---------------------------------------

### does events with overshoot has faster ndvi decrease speed?
## from the deseasonlized and detrended anomlies, together with the drought event and overshoot events.

#load the total drought.

analyze_speed_of_change<-function(vec){
  n_obs<-length(vec)/2
  drought_ts<-vec[1:n_obs]
  deseason_vi<-vec[(n_obs+1):(n_obs*2)]
  num_droughts<-min(max(drought_ts),50)
  if (sum(is.na(deseason_vi))>n_obs*0.75|sum(is.na(drought_ts))>0|num_droughts==0){
    return(rep(NA,50*7))
  }
  #### get the monotonic decrease during the drought period
  deseason_vi[is.na(deseason_vi)]<-0
  diff_ndvi<-c(diff(deseason_vi),0)
  decrease_vi<- (diff_ndvi< 0)
  rle_decrease<-rle(decrease_vi)
  ndvi_drought<-rle_decrease$values==T
  ndvi_drought[ndvi_drought==1]<-1:sum(ndvi_drought,na.rm=T)
  decrease_indi<-c(0,rep(ndvi_drought,rle_decrease$lengths))
  
  ## two things to be done here
  # 1. get the timing of the start of drought, and the timing of severe drought
  timing_of_start<-rep(NA,50)
  timing_of_drought<-rep(NA,50)
  
  ### for each drought event, search for the start, and the minimum value.
  decrease_slope1<-rep(NA,50)     ### this is used to get the median of the speed of change
  decrease_slope2<-rep(NA,50)     ### this is used to get the 75 quantile of the speed of change
  decrease_slope3<-rep(NA,50)     ### get the change of speed when cross zero.
  drought_length<-rep(NA,50)
  develop_length<-rep(NA,50)
  if (num_droughts==0){
    return(rep(NA,50*7))
  }
  for (j in 1:num_droughts){
    # get the date when minimum is reached
    drought_vi<-deseason_vi[drought_ts==j]
    drought_length[j]<-length(drought_vi)
    timing_of_start[j]<-which(drought_ts==j)[1]
    timing_of_drought[j]<-timing_of_start[j]-1+which.min(drought_vi)
    ## get the drought starting date
    decrease_event_id<-decrease_indi[timing_of_start[j]]
    decrease_start_date<-which(decrease_indi==decrease_event_id)[1]-1
    # get median of the slope
    res<-quantile(diff_ndvi[decrease_start_date:timing_of_drought[j]],c(0.25,0.5),na.rm=T)
    #zero_crossing<-which(deseason_vi[timing_of_start[j]:timing_of_drought[j]]<0)[1]
    decrease_slope1[j]<-res[1]
    decrease_slope2[j]<-res[2]
    decrease_slope3[j]<-diff_ndvi[max(timing_of_start[j]-1,1)]
    # develop_length
    develop_length[j]<-timing_of_drought[j]-timing_of_start[j]+1
  }
  return(c(decrease_slope1,decrease_slope2,decrease_slope3,develop_length,
           drought_length,timing_of_start,timing_of_drought))
}


get_ts_speed_of_change<-function(vec){
  n_obs<-(length(vec)-50)/2
  drought_ts<-vec[1:n_obs]
  num_droughts<-min(max(drought_ts),50)
  deseason_vi<-vec[(n_obs+1):(n_obs*2)]
  overshootid<-vec[(1+n_obs*2):(n_obs*2+50)]
  overshootid<-overshootid[!is.na(overshootid)]
  if (sum(is.na(deseason_vi))>n_obs*0.75|sum(is.na(drought_ts))>0|num_droughts==0){
    return(rep(NA,750))
  }
  #### get the monotonic decrease during the drought period
  deseason_vi[is.na(deseason_vi)]<-0
  diff_ndvi<-c(diff(deseason_vi),0)
  decrease_vi<- (diff_ndvi< 0)
  rle_decrease<-rle(decrease_vi)
  ndvi_drought<-rle_decrease$values==T
  ndvi_drought[ndvi_drought==1]<-1:sum(ndvi_drought,na.rm=T)
  decrease_indi<-c(0,rep(ndvi_drought,rle_decrease$lengths))
  # 2. get the decline time series of drought development here. array of 15/events, 50 events
  ts_drought_develop<-array(NA,dim=c(15,50))
  
  for (j in 1:num_droughts){
    # get the date when minimum is reached
    drought_vi<-deseason_vi[drought_ts==j]
    start_of_drought_date<-which(drought_ts==j)[1]
    minimum_date<-start_of_drought_date-1+which.min(drought_vi)
    ## get the drought starting date
    decrease_event_id<-decrease_indi[start_of_drought_date]
    decrease_start_date<-which(decrease_indi==decrease_event_id)[1]-1
    # get deseason_vi   1-7 8-15
    positive_per<-deseason_vi[decrease_start_date:start_of_drought_date]
    negative_per<-deseason_vi[start_of_drought_date:minimum_date]
    len_pos<-min(length(positive_per),8)
    len_neg<-min(length(negative_per),8)
    ts_drought_develop[(9-len_pos):8,j]<-positive_per[(length(positive_per)-len_pos+1):(length(positive_per))]
    ts_drought_develop[8:(7+len_neg),j]<-negative_per[1:len_neg]
  }
  ## overshoot related 
  ts_drought_develop_overshoot<-array(NA,dim=c(15,50))
  ## first 15 for the overshoot events, 16-50 for nonovershoot
  overshoot_num<-min(sum(overshootid,na.rm=T),15)
  overshoot_cum_sum<-cumsum(overshootid)
  nonovershootid<-overshootid==0
  overshootid[overshoot_cum_sum>15]<-0
  nonovershoot_cum_sum<-cumsum(nonovershootid)
  nonovershoot_num<-min(sum(nonovershootid,na.rm=T),35)
  nonovershootid[nonovershoot_cum_sum>35]<-0
  
  if(overshoot_num>=1){
    ts_drought_develop_overshoot[,1:overshoot_num]<-ts_drought_develop[,which(overshootid==1)]
  }
  if(nonovershoot_num>=1){
    ts_drought_develop_overshoot[,16:(15+nonovershoot_num)]<-ts_drought_develop[,which(nonovershootid==1)]
  }
  return(as.vector(ts_drought_develop_overshoot))
}

testsig<-function(vec){
  n_obs<-length(vec)/2
  withover<-vec[1:n_obs][vec[(n_obs+1):(n_obs*2)]==1]
  withoutover<-vec[1:n_obs][vec[(n_obs+1):(n_obs*2)]==0]
  if (sum(!is.na(withoutover))>=2&sum(!is.na(withover))>=2){
    if (sum(!is.na(unique(withover)))==1&sum(!is.na(unique(withoutover)))==1){
      return(c(mean(withover,na.rm=T),mean(withoutover,na.rm=T),0))
    }
    tst<-t.test(withover,withoutover,alternative = "less")
    return(c(tst$estimate,tst$p.value))
  }else{
    return(rep(NA,3))
  }
}
get_normalized<-function(vec){
  return(vec[1:750]/vec[751])
}


if (ind_speed_for_drought_dev){
  ### get the drought events and the deseasonalized detrended NDVI anomaly
  ncin<-nc_open(paste0("./analysis/gimms/gimms.component_pred.",out_label,".spinup5.rep2.delta0.98.nc"))
  pred_comp<-list()
  for (i in 1:length(mean_ind)){
    pred_comp[[i]]<-ncvar_get(ncin,varid=mean_ind[i])[,,121:534]
  }
  nc_close(ncin)
  ## read the coefficients
  ncin<-nc_open(paste0("./analysis/gimms/gimms.DLM.results.",out_label,".spinup5.rep2.delta0.98.nc"))
  Yano<-ncvar_get(ncin,varid="Yanomaly")[,,121:534]
  nc_close(ncin)
  load(drought_f)
  drought_ts<-aperm(drought_all,c(2,3,1))
  desea_detr<-Yano-pred_comp[[1]]-pred_comp[[2]]-pred_comp[[3]]-pred_comp[[4]]
  drought<-abind(drought_ts[,,1:414],desea_detr[,,1:414],along=3)
  
  results<-parApply(cl,drought,c(1,2),analyze_speed_of_change)
  
  decrease_speed1<-results[1:50,,]
  decrease_speed2<-results[51:100,,]
  decrease_speed3<-results[101:150,,]
  develop_length<-results[151:200,,]
  drought_length<-results[201:250,,]
  start_timing<-results[251:300,,]
  minimum_timing<-results[301:350,,]
  
  save(decrease_speed1,decrease_speed2,decrease_speed3,drought_length,develop_length,start_timing,minimum_timing,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_decreasespeed_droughtlength_developlength_method2_spei05.RData"))
  
  # get time series of drought development
  
  load(paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_1_414_spei05.RData"))
  #ndvidat<-abind(sum_drought_ndvi,overshoot_boolean,along=1)
  overshoot<-aperm(overshoot_boolean,c(2,3,1))
  drought_with_overshoot<-abind(drought,overshoot,along=3)
  
  time_series_change<-parApply(cl,drought_with_overshoot,c(1,2),get_ts_speed_of_change)
  
  #### normalized speed changes.
  ## get the std of variation
  std_desea_detr<-parApply(cl,desea_detr,c(1,2),sd,na.rm=T)
  comb_ts_change<-abind(time_series_change,std_desea_detr,along=1)
  normalized_time_series_changes<-parApply(cl, comb_ts_change,c(2,3),get_normalized)
  
  save(std_desea_detr,time_series_change,normalized_time_series_changes,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_decrease_time_series_method2_spei05.RData"))
  
  ### get area of interest
  
  cords<-t(array(c(63,70,-140,-150,
                   33,43,-98,-88,
                   -30,-20,-58,-48,
                   -25,-15,15,25,
                   20,30,75,85,
                   60,70,100,115),
                 dim=c(4,6)))
  
  cordsxy<-cbind(180+cords[,1:2]*2,360+cords[,3:4]*2)
  
  regions_ts<-list()
  regions_normalized_ts<-list()
  
  for (i in 1:6){
    sel_region1<-time_series_change[,(cordsxy[i,3]+1):cordsxy[i,4],(cordsxy[i,1]+1):cordsxy[i,2]]
    sel_region2<-normalized_time_series_changes[,(cordsxy[i,3]+1):cordsxy[i,4],(cordsxy[i,1]+1):cordsxy[i,2]]
    regions_ts[[i]]<-sel_region1
    regions_normalized_ts[[i]]<-sel_region2
  }
  
  save(regions_ts,regions_normalized_ts,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_6region_ts_changes_method2_spei05.RData"))
  
  # test the speed of drought development
  
  file1=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_decreasespeed_droughtlength_developlength_method2_spei05.RData")
  load(file1)
  file2=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_1_414_spei05.RData")
  load(file2)
  dat1<-abind(decrease_speed1,overshoot_boolean,along=1)
  test_decrease_speed1<-parApply(cl,dat1,c(2,3),testsig)
  dat2<-abind(decrease_speed2,overshoot_boolean,along=1)
  test_decrease_speed2<-parApply(cl,dat2,c(2,3),testsig)
  dat3<-abind(decrease_speed3,overshoot_boolean,along=1)
  test_decrease_speed3<-parApply(cl,dat3,c(2,3),testsig)
  
  dat<-abind(drought_length,overshoot_boolean,along=1)
  test_drought_length<-parApply(cl,dat,c(2,3),testsig)
  dat<-abind(develop_length,overshoot_boolean,along=1)
  test_develop_length<-parApply(cl,dat,c(2,3),testsig)
  
  save(test_decrease_speed1,test_decrease_speed2,test_decrease_speed3,test_drought_length,test_develop_length,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_test_decreasespeed_droughtlength_developlength_mehtod2_spei05.RData"))
  
}

# ANOVA overshoot impact 3 model-----------------------------------------------


change_res<-function(mat){  # matsize = 50,720,360
  dim(mat)<-c(250,144,360)
  mat_aperm<-aperm(mat,c(1,3,2))
  dim(mat_aperm)<-c(1250,72,144)
  mat_lowres<-aperm(mat_aperm,c(1,3,2))
  return(mat_lowres)
}
partial_f_test_low<-function(vec){
  ## build two models and compare the model with anova
  n_obs<-length(vec)/3
  ndvi<-vec[1:n_obs]
  spei3<-vec[(1+n_obs):(2*n_obs)]
  overshoot_com<-vec[(1+n_obs*2):(3*n_obs)]==1
  overshoot_com[is.na(ndvi)]<-NA
  if (sum(overshoot_com==T,na.rm=T)>=10&sum(overshoot_com==F,na.rm=T)>=10){
    overshoot<-as.factor(overshoot_com)
    model1<-lm(ndvi~spei3)
    model2<-lm(ndvi~spei3+overshoot)
    model3<-lm(ndvi~spei3*overshoot)
    av<-anova(model1,model2,model3)
    # compare the models whether
    sig<-av$`Pr(>F)`[2:3]
    ## for model3 it can provide additional information
    model1_coef<-as.vector(summary(model1)$coefficients[,c(1,4)])
    model2_coef<-as.vector(summary(model2)$coefficients[,c(1,4)])
    model3_coef<-as.vector(summary(model3)$coefficients[,c(1,4)])
    return(c(sig,as.vector(model1_coef),as.vector(model2_coef),as.vector(model3_coef),
             summary(model1)$r.squared,summary(model2)$r.squared,summary(model3)$r.squared))
  }else{
    return(rep(NA,23))
  }
}

get_normalized<-function(vec){
  return(vec[1:50]/vec[51])
}

if (ind_test_drought_impact_anova_low){
  ## compare two models, with different slopes and different intercepts,test conducted at 5 by 5 deg
  file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_drought_impact_method2_1_414_spei05.RData")
  load(file)
  
  load(paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_1_414_spei05.RData"))
  load(paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_decrease_time_series_method2_spei05.RData"))
  minimum_drought_ndvi[is.infinite(minimum_drought_ndvi)]<-NA
  
  all_drought_ndvi<-abind(sum_drought_ndvi,std_desea_detr,along=1)
  normalized_sum_drought_ndvi<-parApply(cl,all_drought_ndvi,c(2,3),get_normalized)
  
  ### compare NDVI and SPEI relationship
  #ndvispei<-abind(change_res(normalized_sum_drought_ndvi),change_res(sum_drought_spei),along=1)
  ndvispeidat<-abind(change_res(normalized_sum_drought_ndvi),change_res(sum_drought_spei),change_res(overshoot_boolean),along=1)
  test_anova_low<-parApply(cl,ndvispeidat,c(2,3),partial_f_test_low)
  
  anovam2_pval<-test_anova_low[1,,]
  anovam3_pval<-test_anova_low[2,,]
  model1_coef<-test_anova_low[3:4,,]
  model1_pval<-test_anova_low[5:6,,]
  model2_coef<-test_anova_low[7:9,,]
  model2_pval<-test_anova_low[10:12,,]
  model3_coef<-test_anova_low[13:16,,]
  model3_pval<-test_anova_low[17:20,,]
  model1_rsq<-test_anova_low[21,,]
  model2_rsq<-test_anova_low[22,,]
  model3_rsq<-test_anova_low[23,,]
  
  save(anovam2_pval,anovam3_pval,model1_coef,model1_pval,model2_coef,model2_pval, model3_coef,model3_pval,
       model1_rsq,model2_rsq,model3_rsq,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_ndvi_anova_test_lowres_1_414_spei05.RData"))
  
  
  all_drought_ndvi<-abind(minimum_drought_ndvi,std_desea_detr,along=1)
  normalized_minimum_drought_ndvi<-parApply(cl,all_drought_ndvi,c(2,3),get_normalized)
  
  ### compare minimum NDVI and SPEI relationship
  #ndvispei<-abind(change_res(normalized_minimum_drought_ndvi),change_res(minimum_drought_spei),along=1)
  ndvispeidat<-abind(change_res(normalized_minimum_drought_ndvi),change_res(minimum_drought_spei),change_res(overshoot_boolean),along=1)
  test_anova_low<-parApply(cl,ndvispeidat,c(2,3),partial_f_test_low)
  
  anovam2_pval<-test_anova_low[1,,]
  anovam3_pval<-test_anova_low[2,,]
  model1_coef<-test_anova_low[3:4,,]
  model1_pval<-test_anova_low[5:6,,]
  model2_coef<-test_anova_low[7:9,,]
  model2_pval<-test_anova_low[10:12,,]
  model3_coef<-test_anova_low[13:16,,]
  model3_pval<-test_anova_low[17:20,,]
  model1_rsq<-test_anova_low[21,,]
  model2_rsq<-test_anova_low[22,,]
  model3_rsq<-test_anova_low[23,,]
  
  save(anovam2_pval,anovam3_pval,model1_coef,model1_pval,model2_coef,model2_pval, model3_coef,model3_pval,
       model1_rsq,model2_rsq,model3_rsq,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_ndvi_minimum_anova_test_lowres_1_414_spei05.RData"))
  
  
  ### compare sm and SPEI relationship
  #smspei<-abind(change_res(sum_drought_sm),change_res(sum_drought_spei),along=1)
  smspeidat<-abind(change_res(sum_drought_sm),change_res(sum_drought_spei),change_res(overshoot_boolean),along=1)
  test_anova_low<-parApply(cl,smspeidat,c(2,3),partial_f_test_low)
  
  anovam2_pval<-test_anova_low[1,,]
  anovam3_pval<-test_anova_low[2,,]
  model1_coef<-test_anova_low[3:4,,]
  model1_pval<-test_anova_low[5:6,,]
  model2_coef<-test_anova_low[7:9,,]
  model2_pval<-test_anova_low[10:12,,]
  model3_coef<-test_anova_low[13:16,,]
  model3_pval<-test_anova_low[17:20,,]
  model1_rsq<-test_anova_low[21,,]
  model2_rsq<-test_anova_low[22,,]
  model3_rsq<-test_anova_low[23,,]
  
  save(anovam2_pval,anovam3_pval,model1_coef,model1_pval,model2_coef,model2_pval, model3_coef,model3_pval,
       model1_rsq,model2_rsq,model3_rsq,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_sm_anova_test_lowres_1_414_spei05.RData"))
  
}


# analyze the timing of overshoot -----------------------------------------

get_overshoot_timing_comparison<-function(vec){
  
  hist_month<-function(vec){
    # convert id to months.
    vec<-vec[!is.na(vec)]
    if(length(vec)==1){
      return((vec+5)%%12+1)
    }
    mon1<-(vec+5)%%12+1    #start from Jan
    mon2<-(vec+2)%%12+1    #start from  Apr
    mon3<-(vec+11)%%12+1    #start from Jul
    mon4<-(vec+8)%%12+1    #start from  Oct
    den1<-density(mon1)
    den2<-density(mon2)
    den3<-density(mon3)
    den4<-density(mon4)
    den<-list(den1,den2,den3,den4)
    den_maxy<-c(max(den1$y),max(den2$y),max(den3$y),max(den4$y))
    if (length(unique(den_maxy))==1){
      return(0)
    }
    ind<-which.max(den_maxy)[1]
    monadj<-round(den[[ind]]$x[which.max(den[[ind]]$y)]+3*(ind-1))
    return((monadj-1)%%12+1)
  }
  
  n_obs<-length(vec)/2
  timing<-vec[1:n_obs]
  overshoot<-vec[(1+n_obs):(2*n_obs)]
  if (sum(!is.na(timing))==0){
    return(rep(NA,3))
  }
  timing_overshoot<-timing[overshoot==1]
  if(sum(!is.na(timing_overshoot))>0){
    overshoot_mon<-hist_month(timing_overshoot)
  }else{
    overshoot_mon<-NA
  }
  timing_nonovershoot<-timing[overshoot==0]
  
  if(sum(!is.na(timing_nonovershoot))>0){
    nonovershoot_mon<-hist_month(timing_nonovershoot)
  }else{
    nonovershoot_mon<-NA
  }
  
  all_mon<-hist_month(timing)
  return(c(overshoot_mon,nonovershoot_mon,all_mon))
}


if (ind_fit_drought_timing){
  file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_decreasespeed_droughtlength_developlength_method2_spei05.RData")
  load(file)
  overshoot_file<-paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_1_414_spei05.RData")
  load(overshoot_file)
  ## get two different period
  start_timing_overshoot<-abind(start_timing,overshoot_boolean,along=1)
  usual_start_drought<-parApply(cl,start_timing_overshoot,c(2,3),get_overshoot_timing_comparison)
  
  minimum_timing_overshoot<-abind(minimum_timing,overshoot_boolean,along=1)
  usual_minimum_drought<-parApply(cl,minimum_timing_overshoot,c(2,3),get_overshoot_timing_comparison)                              
  
  save(usual_start_drought,usual_minimum_drought, 
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_timing_method2_1_414_spei05.RData"))
}



# get overshoot timing with phenology -------------------------------------

get_overshoot_timing_stat<-function(vec){
  adjustmon<-function(vec){
    adj<-ifelse(vec<= -6,vec+12,vec)
    adj2<-ifelse(adj>= 6,adj-12,adj)
    return(adj2)
  }
  
  n_obs<-(length(vec)-1)/2
  timing<-(vec[1:n_obs]-7)%%12+1
  overshoot<-vec[1:n_obs+n_obs]
  overshoot[is.na(timing)]<-NA
  pos<-vec[2*n_obs+1]
  overshoot_timing<-timing[overshoot==1]-0.5-pos/30.5
  nonovershoot_timing<-timing[overshoot==0]-0.5-pos/30.5
  os_timing<-rep(NA,50)
  nos_timing<-rep(NA,50)
  
  os_timing[1:length(overshoot_timing)]<-adjustmon(overshoot_timing)
  nos_timing[1:length(nonovershoot_timing)]<-adjustmon(nonovershoot_timing)
  
  if (sum(!is.na(overshoot_timing))>=1&sum(!is.na(nonovershoot_timing))>=1){
    if (sum(!is.na(unique(overshoot_timing)))==1|sum(!is.na(unique(nonovershoot_timing)))==1){
      return(c(os_timing,nos_timing,mean(overshoot_timing,na.rm=T),mean(nonovershoot_timing,na.rm=T),0))
    }
    tst<-t.test(overshoot_timing,nonovershoot_timing,alternative = "two.sided")
    return(c(os_timing,nos_timing,tst$estimate,tst$p.value))
  }else{
    return(c(os_timing,nos_timing,mean(overshoot_timing,na.rm=T),mean(nonovershoot_timing,na.rm=T),0))
  }
}

### get os timing and phenology

if (ind_os_nos_timing){
  load("./Data/pheno_hd.RData")
  file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_decreasespeed_droughtlength_developlength_method2_spei05.RData")
  load(file)
  overshoot_file<-paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_1_414_spei05.RData")
  load(overshoot_file)
  ## get two different periods
  overshoot_pos<-abind(overshoot_boolean,pos,along=1)
  start_timing_overshoot_pos<-abind(start_timing,overshoot_pos,along=1)
  temp_start<-parApply(cl,start_timing_overshoot_pos,c(2,3),get_overshoot_timing_stat)
  
  minimum_timing_overshoot_pos<-abind(minimum_timing,overshoot_pos,along=1)
  temp_min<-parApply(cl,minimum_timing_overshoot_pos,c(2,3),get_overshoot_timing_stat)
  
  os_nos_start<-temp_start[1:100,,]
  os_nos_minimum<-temp_min[1:100,,]
  start_timing_diff<-temp_start[101:103,,]
  minimum_timing_diff<-temp_min[101:103,,]
  
  save(os_nos_start,os_nos_minimum,start_timing_diff,minimum_timing_diff,
       file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_timing_all_events_method2_1_414_spei05.RData"))
}


# get overshoot time and corresponding temperature ------------------------
get_overshoot_tmp_stat<-function(vec){
  t_test<-function(a,b=rep(0,length(a))){
    if (sum(!is.na(a))>=1&sum(!is.na(b))>=1){
      if (sum(!is.na(unique(a)))==1|sum(!is.na(unique(b)))==1){
        return(c(mean(a,na.rm=T),mean(b,na.rm=T),0))
      }
      tst<-t.test(a,b,alternative = "two.sided")
      return(c(tst$estimate,tst$p.value))
    }else{
      return(rep(NA,3))
    }
  }
  
  #n_obs<-(length(vec)-100)
  tmp<-vec[1:414]
  os_boolean<-vec[1:50+414]
  timing<-vec[51:100+414]
  
  os_month<- ifelse(vec[515]==0,NA,vec[515])
  nos_month<-ifelse(vec[516]==0,NA,vec[516])
  
  if(sum(os_boolean,na.rm=T)<1){
    return(rep(NA,14))
  }
  timing_months<-(timing-7)%%12+1
  msc_temp<-tmp[7:414]
  dim(msc_temp)<-c(12,34)
  msc_tmp<-apply(msc_temp,1,mean,na.rm=T)
  ano_tmp<-tmp-c(msc_tmp[7:12],rep(msc_tmp,34))
  ### get the obs months
  msc_event<-msc_tmp[timing_months]
  os_bool<-os_boolean[!is.na(timing)]==1
  os_tmp<-msc_event[os_bool]
  nos_tmp<-msc_event[!os_bool]
  ttest<-t_test(os_tmp,nos_tmp)
  ### actual temp for drought compared to normal
  ano_tmp_event<-ano_tmp[timing]
  os_tmp_ano<-ano_tmp_event[os_bool]
  nos_tmp_ano<-ano_tmp_event[!os_bool]
  ttest1<-t_test(os_tmp_ano)
  ttest2<-t_test(nos_tmp_ano)
  ttest3<-t_test(os_tmp_ano,nos_tmp_ano)
  return(c(ttest,ttest3,ttest1,ttest2,ifelse(is.na(os_month),NA,msc_tmp[os_month]),ifelse(is.na(nos_month),NA,msc_tmp[nos_month])))
}


if (ind_os_temeprature){
  ## load overshoot
  overshoot_file<-paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_method2_1_414_spei05.RData")
  load(overshoot_file)
  file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_decreasespeed_droughtlength_developlength_method2_spei05.RData")
  load(file)
  file=paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_overshoot_timing_method2_1_414_spei05.RData")
  load(file)
  ## load temperature
  tmp_file<-'/global/scratch/yaozhang/Data/CRU/cru_ts4.04.1901.2019.tmp.dat.nc'
  ncin<-nc_open(tmp_file)
  tmp<-ncvar_get(ncin,varid="tmp")[,,967:1380]    #temp start from 1979/07/01
  
  tmp_aperm<-aperm(tmp,c(3,1,2))
  # dat<-abind(tmp_aperm,overshoot_boolean,along=1)
  # start_timing_tmp<-abind(dat,start_timing,along=1)
  start_timing_tmp_usual<-abind(tmp_aperm,overshoot_boolean,start_timing,usual_start_drought,along=1)
  start_temperature<-parApply(cl,start_timing_tmp_usual,c(2,3),get_overshoot_tmp_stat)
  
  #minimum_timing_tmp<-abind(dat,minimum_timing,along=1)
  minimum_timing_tmp_usual<-abind(tmp_aperm,overshoot_boolean,minimum_timing,usual_minimum_drought,along=1)
  minimum_temperature<-parApply(cl,minimum_timing_tmp_usual,c(2,3),get_overshoot_tmp_stat)
  
  fileout<-paste0("./analysis/gimms_overshoot_",substr(out_label,1,11),"/",out_label,"_tmp_diff_ano_method2_spei05.RData")
  save(start_temperature,minimum_temperature,file=fileout)
}



