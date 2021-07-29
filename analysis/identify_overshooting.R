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
#indcators<-c("prec",'prec_temp','prec_temp_solr','randomndvi')
indcator<-'prec_temp_solr'
out_label<-indcator


library(raster)
library(ncdf4)
library(abind)
library(zoo)
library(parallel)


ind_calculate_drought<-T
ind_evaluate_drought<-T

ind_evaluate_overshoot<-T


setwd("/global/scratch/yaozhang/Project/overshooting/")

mean_ind<-c('trend','seasonal1','seasonal2','seasonal3',
            'pre_contri1','pre_contri23','pre_contri46','pre_contri712','pre_contri1324',
            'pre_contri13','tmp_contri','rad_contri')

cl <- makeCluster(8)
clusterEvalQ(cl, {
  library(zoo)
})



# 1.1 identify drought for each period ----------------------------------------
#get model performance
get_r_square<-function(vec){
  obs_n<-length(vec)/2
  obs<-vec[1:obs_n]
  pred<-vec[1:obs_n+obs_n]
  if(sum(!is.na(obs))<100|sum(!is.na(pred))<100){
    return(c(NA,NA))
  }
  rval<-cor.test(obs,pred)
  ### get rrmse
  reg<-lm(obs~pred)
  rrmse<-sqrt(mean(reg$residual^2,na.rm=T))/diff(range(obs,na.rm=T))
  return(c(rval$estimate,rrmse))
}

#identify drought events
# a. negative anomalies were identified first
# b. drought starts from first negative until recover to 70% of the minimum.
# c. during drought development, temperature contribution should be smaller than 
#    precipitation or temperature sensitivity should be negative.
# d. SPEI value should be smaller than -0.5


# for test
# vec<-c(obs_de,spei_ts,predmean_ts[121:534,9],component_ts[121:534,10],component_ts[121:534,11])
drought_identify<-function(vec){
  n_obs <- (length(vec)-1)/6
  ### use both the threshold and the spei to determine the drought events
  vi<-vec[1:n_obs]
  spei<-vec[1:n_obs+n_obs]
  prec_sensi<-vec[1:n_obs+2*n_obs]
  temp_sensi<-vec[1:n_obs+3*n_obs]
  prec_contri<-vec[1:n_obs+4*n_obs]
  temp_contri<-vec[1:n_obs+5*n_obs]
  meanvi<-vec[n_obs*6+1]
  if (sum(is.na(vi))>n_obs*0.75|sum(is.na(spei))>n_obs*0.75|is.nan(meanvi)){
    return(rep(NA,n_obs))
  }
  vi[is.na(vi)]<-0
  drought<- (vi< 0)   ##### drought should be negative
  rle_dro<-rle(drought)
  ndvi_drought<-rle_dro$values==T&rle_dro$lengths>=2
  ndvi_drought[ndvi_drought==1]<-1:sum(ndvi_drought)
  drought_indi<-rep(ndvi_drought,rle_dro$lengths)
  drought_indi_new<-drought_indi
  ##### drought should start from below zero to get to the minimum
  if (max(drought_indi)>=1){
    for (i in 1:max(drought_indi)){
      vi_dro<-vi[drought_indi==i]
      min_drought<-min(vi_dro)
      drou_thr<-min_drought*0.7
      ind_min<-which.min(vi_dro)[1]
      if (ind_min<length(vi_dro)){
        ind_end<-which(vi_dro[ind_min:length(vi_dro)]> drou_thr)[1]+ind_min-2
        if (is.na(ind_end)){
          ind_end<-ind_min
        }else{
          drought_indi_new[drought_indi==i][(ind_end+1):length(vi_dro)]<-0
        }
      }else{
        ind_end<-ind_min
      }
      spei_dro<-spei[drought_indi_new==i]   ### spei smaller than p< -0.5
      temp_sensi_dro<-temp_sensi[drought_indi_new==i]   ### mean temperature sensi negative
      prec_sensi_dro<-prec_sensi[drought_indi_new==i]
      temp_contr_dro<-temp_contri[drought_indi_new==i]
      prec_contr_dro<-prec_contri[drought_indi_new==i]
      if (sum(!is.na(spei_dro))==0){
        drought_indi_new[drought_indi==i]<-0
      }else if((mean(spei_dro,na.rm=T)> -0.5)|(min(vi_dro,na.rm=T)> -0.1*meanvi)|
               #(mean(prec_contr_dro)> 0)|(mean(prec_sensi_dro)<=0)|
               ((mean(temp_sensi_dro)>0)&(mean(temp_contr_dro)<0)&(abs(mean(temp_contr_dro))>abs(mean(prec_contr_dro))))){
        drought_indi_new[drought_indi==i]<-0
      }
    }
  }
  drought_indi_new[drought_indi_new!=0]<-1
  rle_dro<-rle(drought_indi_new)
  ndvi_drought<-rle_dro$values==T&rle_dro$lengths>=2
  ndvi_drought[ndvi_drought==1]<-1:sum(ndvi_drought)
  drought_indi<-rep(ndvi_drought,rle_dro$lengths)
  return(drought_indi)
}

if (ind_calculate_drought){
  ### read the components of DLM 
  var_file<-"./analysis/other_data/climate_variables.RData"
  load(var_file)
  
  drought_f<-paste0("./analysis/gimms_overshoot_",out_label,"/drought_with_spei_",out_label,"_spei05.RData")
  ncin<-nc_open(paste0("./analysis/gimms/gimms.component_pred.",out_label,".spinup5.rep2.delta0.98.nc"))
  pred_comp<-list()
  for (i in 1:length(mean_ind)){
    pred_comp[[i]]<-ncvar_get(ncin,varid=mean_ind[i])[,,121:534]
  }
  nc_close(ncin)
  ## read the coefficients
  ncin<-nc_open(paste0("./analysis/gimms/gimms.DLM.results.",out_label,".spinup5.rep2.delta0.98.nc"))
  Yano<-ncvar_get(ncin,varid="Yanomaly")[,,121:534]
  sensi<-ncvar_get(ncin,varid="Predictmean")[,,121:534,8:9]
  nc_close(ncin)
  
  ## read spei dataset

  ncin<-nc_open("./Data/spei03.nc")
  spei3<-ncvar_get(ncin,varid='spei')[,,967:1380]
  nc_close(ncin)

  
  seasonal_trend<-pred_comp[[1]]+pred_comp[[2]]+pred_comp[[3]]+pred_comp[[4]]
  desea_detr<-Yano-seasonal_trend
  #comp<-abind(desea_detr,spei3,along=3)
  ### add temp sensi into consideration, temp contribution, prec contribution  
  comp3<-abind(desea_detr,spei3,sensi[,,,1],sensi[,,,2],pred_comp[[10]],pred_comp[[11]],meaVI/10000,along=3)
  #contri<-abind(pred_comp[[10]],pred_comp[[11]],along=3)
  #comp3<-abind(comp2,contri,along=3)
  comp3[is.nan(comp3)]<-NA
  drought_all<-parApply(cl,comp3,c(1,2),drought_identify)
  
  ### get the total drought numbers.
  total_drought<-parApply(cl,drought_all,c(2,3),max,na.rm=T)
  save(desea_detr,drought_all,total_drought,
       file=drought_f)
  
  pred_anomaly<-pred_comp[[5]]+pred_comp[[6]]+pred_comp[[7]]+pred_comp[[8]]+pred_comp[[9]]+pred_comp[[10]]+pred_comp[[11]]+pred_comp[[12]]
  dat_comb<-abind(desea_detr,pred_anomaly,along=3)
  anomaly_performance<-parApply(cl,dat_comb,c(1,2),get_r_square)
  
  pred_y<-pred_comp[[1]]+pred_comp[[2]]+pred_comp[[3]]+pred_comp[[4]]+pred_anomaly
  dat_comb<-abind(Yano,pred_y,along=3)
  obs_performance<-parApply(cl,dat_comb,c(1,2),get_r_square)
  save(anomaly_performance,obs_performance,
       file=paste0("./analysis/gimms_overshoot_",out_label,"/DLM_performance_",out_label,"_spei05.RData"))
}


# 1.2 evaluate each drought event ---------------------------------------------
### input: a. drought identifier 
###        b. drought NDVI deseason detrend 
###        c. drought SPEI3 
###        d. DLM avg coef 
###        e. DLM avg component
###
evaluate_overshoot<-function(vec){
  ## negative contribution from previous month greater than negative from current months.
  # combination of drought events, coef and contribution
  n_obs<-(length(vec)-1)/6
  if (sum(is.na(vec))>n_obs*0.9*2){
    return(rep(NA,n_obs+1+50*2))
  }
  drought_event<-vec[1:n_obs]
  drought_ndvi<-vec[(n_obs+1):(n_obs*2)]
  spei3<- vec[(n_obs*2+1):(n_obs*3)]
  coef<-vec[(n_obs*3+1):(n_obs*4)]
  coef_upper<-vec[(n_obs*4+1):(n_obs*5)]
  comp<-vec[(n_obs*5+1):(n_obs*6)]
  meanvi<-vec[n_obs*6+1]
  num_drought_id<-unique(drought_event,na.rm=T)
  num_drought_id<-setdiff(num_drought_id,0)
  num_drought<-length(num_drought_id)
  num_drought<-min(num_drought,50)
  if (num_drought>=1){
    over_count<-1
    overshoot_ts_count<-rep(0,n_obs)
    comp_drought_ndvi<-rep(0,50)
    overshoot_drought_boolean<-rep(0,50)
    for (i in 1:num_drought){
      drought_id<-which(drought_event==num_drought_id[i])
      #sum_ndvi_ano<-sum(drought_ndvi[drought_id],na.rm=T)
      mean_coef<-mean(coef[drought_id],na.rm=T)
      mean_coef_upper<-min(coef_upper[drought_id],na.rm=T)
      comp_drought_ndvi[i] <-sum(comp[drought_id],na.rm=T)
      overshoot_drought_boolean[i]<-(mean_coef<0)&(comp_drought_ndvi[i]<0)&(mean_coef_upper<0)&
        (min(comp[drought_id],na.rm=T)< -0.02*meanvi)    # this may still need some discussion
      if (overshoot_drought_boolean[i]==1){
        overshoot_ts_count[drought_id]<-over_count
        over_count <- over_count+1
      }
    }
    comp_avg_overshoot_contri<-sum(comp[overshoot_ts_count!=0],na.rm=T)/sum(drought_ndvi[overshoot_ts_count!=0],na.rm=T)
    
    return(c(overshoot_ts_count,comp_avg_overshoot_contri,comp_drought_ndvi,overshoot_drought_boolean))
  }else{
    return(rep(0,n_obs+1+50*2))
  }
}

get_drought_stat<-function(vec){
  n_obs<-length(vec)/6
  if (sum(is.na(vec))>n_obs*0.9*6){
    return(rep(NA,400))
  }
  drought_event<-vec[1:n_obs]
  drought_ndvi<-vec[(n_obs+1):(n_obs*2)]
  drought_pred_ndvi<-vec[(n_obs*2+1):(n_obs*3)]
  nondrought_ref_ndvi<-vec[(n_obs*3+1):(n_obs*4)]
  spei3<- vec[(n_obs*4+1):(n_obs*5)]
  sm<- vec[(n_obs*5+1):(n_obs*6)]
  
  num_drought_id<-unique(drought_event,na.rm=T)
  num_drought_id<-setdiff(num_drought_id,0)
  num_drought<-length(num_drought_id)
  num_drought<-min(num_drought,50)
  
  sum_drought_ndvi<-rep(NA,50)
  minimum_drought_ndvi<-rep(NA,50)
  sum_nondrought_ref_ndvi<-rep(NA,50)
  sum_drought_pred_ndvi<-rep(NA,50)
  sum_drought_spei<-rep(NA,50)
  minimum_drought_spei<-rep(NA,50)
  sum_drought_sm<-rep(NA,50)
  minimum_drought_sm<-rep(NA,50)
  if (num_drought>0){
    for (i in 1:num_drought){
      drought_id<-which(drought_event==num_drought_id[i])
      sum_drought_ndvi[i]<-sum(drought_ndvi[drought_id],na.rm=T)
      sum_drought_pred_ndvi[i]<-sum(drought_pred_ndvi[drought_id],na.rm=T)
      sum_nondrought_ref_ndvi[i]<-sum(nondrought_ref_ndvi[drought_id],na.rm=T)
      sum_drought_spei[i]<-sum(spei3[drought_id],na.rm=T)
      sum_drought_sm[i]<-mean(sm[drought_id],na.rm=T)*length(drought_id)
      minimum_drought_ndvi[i]<-min(drought_ndvi[drought_id],na.rm=T)
      minimum_drought_spei[i]<-min(spei3[drought_id],na.rm=T)
      minimum_drought_sm[i]<-min(sm[drought_id],na.rm=T)
    }
    return(c(sum_drought_ndvi,sum_drought_pred_ndvi,sum_nondrought_ref_ndvi,sum_drought_spei,sum_drought_sm,
             minimum_drought_ndvi,minimum_drought_spei,minimum_drought_sm))
  }else{
    return(rep(0,400))
  }
}

# 1.3 identify all overshoot events separately ----------------------------------

get_overshoot_for_each_component<-function(t_range,var_comp,outfile){   # t_range 7:210,  211:414
  ### get stat period and for each component
  var_file<-"./analysis/other_data/climate_variables.RData"
  load(var_file)
  
  comb_temp<-abind(drought_ts[,,t_range],desea_detr[,,t_range],spei3[,,t_range],
                   predmean[,,t_range,var_comp+2],upper_uncert[,,t_range,var_comp],
                   pred_comp[[var_comp+4]][,,t_range],meaVI/10000,along=3)
  overshoot_res<-parApply(cl,comb_temp,c(1,2),evaluate_overshoot)
  n_obs<-(dim(overshoot_res)[1]-1-100)
  overshoot_ts_count<-overshoot_res[1:n_obs,,]
  comp_avg_overshoot_contri<-overshoot_res[n_obs+1,,]
  comp_drought_ndvi<-overshoot_res[n_obs+1+1:50,,]
  overshoot_drought_boolean<-overshoot_res[n_obs+1+51:100,,]
  
  save(overshoot_ts_count,comp_avg_overshoot_contri,comp_drought_ndvi,overshoot_drought_boolean,
       file=outfile)
}

evaluate_drought_impact<-function(t_range,outfile){   # t_range 7:210,  211:414
  ### get stat period and for each component
  # #### todo add model predicted drought impact (summation of all component exclude trend and seasonal ones)
  drought_ndvi_spei_sm<-abind(drought_ts[,,t_range],desea_detr[,,t_range],
                              pred_desea_detr[,,t_range],seasonal_ref[,,t_range],
                              spei3[,,t_range],sm[,,t_range],along=3)
  drought_res<-parApply(cl,drought_ndvi_spei_sm,c(1,2),get_drought_stat)
  #n_obs<-(dim(drought_res)[1]-1-100)
  sum_drought_ndvi<-drought_res[1:50,,]
  sum_drought_pred_ndvi<-drought_res[51:100,,]
  sum_nondrought_ref_ndvi<-drought_res[101:150,,]
  sum_drought_spei<-drought_res[151:200,,]
  sum_drought_sm<-drought_res[201:250,,]
  minimum_drought_ndvi<-drought_res[251:300,,]
  minimum_drought_spei<-drought_res[301:350,,]
  minimum_drought_sm<-drought_res[351:400,,]
  
  save(sum_drought_ndvi,sum_drought_pred_ndvi,sum_nondrought_ref_ndvi,minimum_drought_spei,
       sum_drought_spei,minimum_drought_ndvi,sum_drought_sm,minimum_drought_sm,file=outfile)
}

deseason<-function(vec){
  n_yr<-length(vec)/12
  if (sum(!is.na(vec))<3*n_yr){
    return(rep(NA,length(vec)))
  }
  if (n_yr!=floor(n_yr)){
    n_yr<-floor(n_yr)
    vecf<-vec[1:(n_yr*12)]
    dim(vecf)<-c(12,ceiling(n_yr))
    msc<-apply(vecf,1,mean,na.rm=T)
    desea<-vec-rep(msc,n_yr+1)[1:length(vec)]
  }else{
    dim(vec)<-c(12,n_yr)
    msc<-apply(vec,1,mean,na.rm=T)
    desea<-vec-rep(msc,n_yr)
  }
  return(desea)
}


if(ind_evaluate_drought){
  ## get the raw ndvi data
  # ncin<-nc_open(".//Data/GIMMS3g_v1.mon.hd.growingseason.nc")
  # ndvi<-ncvar_get(ncin,varid="ndvi")
  # ndvi[ndvi<0]<-NA
  # nc_close(ncin)
  # ndvi <- ndvi/10000
  # mean_ndvi<- parApply(cl, ndvi,c(1,2), mean,na.rm=T)
  var_file<-"./analysis/other_data/climate_variables.RData"
  load(var_file)
  
  mean_ndvi_mat<-rep(meaVI/10000,414)
  dim(mean_ndvi_mat)<-c(720,360,414)
  
  ## get DLM results
  ncin<-nc_open(paste0("./analysis/gimms/gimms.DLM.results.",out_label,".spinup5.rep2.delta0.98.nc"))
  predmean<-ncvar_get(ncin,varid="Predictmean")[,,121:534,]
  Xpre<-ncvar_get(ncin,varid="Regrevariable")[,,121:534,]
  Yano<-ncvar_get(ncin,varid="Yanomaly")[,,121:534]
  snu<-ncvar_get(ncin,varid="Degreefreedom")[,,121:534]
  predvar<-ncvar_get(ncin,varid="Predictvariance")[,,121:534,]
  nc_close(ncin)
  ### get the uncertainty for each regression coefficient
  unce_std<-qt(0.9,df=snu)
  upper_uncert<-array(NA,dim=c(720,360,414,5))
  for (i in 1:5){
    upper_uncert[,,,i]<-predmean[,,,2+i]+sqrt(predvar[,,,2+i])*unce_std
  }
  
  ## load component 
  ncin<-nc_open(paste0("./analysis/gimms/gimms.component_pred.",out_label,".spinup5.rep2.delta0.98.nc"))
  pred_comp<-list()
  for (i in 1:length(mean_ind)){
    pred_comp[[i]]<-ncvar_get(ncin,varid=mean_ind[i])[,,121:534]
  }
  nc_close(ncin)
  
  ## load SM
  sm_file<-'/global/scratch/yaozhang/Data/ESACCI/monthly_SM_HD.nc'
  ncin<-nc_open(sm_file)
  sm_raw<-ncvar_get(ncin,varid="soil moisture")[,,43:456]    #sm start from 1978/01/01
  sm_raw[is.nan(sm_raw)]<-NA
  ### get deseasonalized sm
  desea_sm<-parApply(cl,sm_raw,c(1,2),deseason)
  sm<-aperm(desea_sm,c(2,3,1))
  
  drought_f<-paste0("./analysis/gimms_overshoot_",out_label,"/drought_with_spei_",out_label,"_spei05.RData")
  load(drought_f)
  drought_ts<-aperm(drought_all,c(2,3,1))
  pred_desea_detr<-pred_comp[[5]]+pred_comp[[6]]+pred_comp[[7]]+pred_comp[[8]]+pred_comp[[9]]+pred_comp[[10]]+pred_comp[[11]]+pred_comp[[12]]
  seasonal_ref<-Yano-desea_detr+mean_ndvi_mat
  

  ncin<-nc_open("./Data/spei03.nc")
  spei3<-ncvar_get(ncin,varid='spei')[,,967:1380]
  nc_close(ncin)

  
  var_id<-c("1","23","46","712","1324","13",'0')
  
  ## all combined

  t_range=1:414
  for (i in 1:5){
    file=paste0("./analysis/gimms_overshoot_",out_label,"/",out_label,"_pre",var_id[i],"_overshoot_method2_",
                t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
    #if (!file.exists(file)){
    get_overshoot_for_each_component(t_range=t_range,var_comp=i,outfile=file)
    #}
  }
  ### for drought
  file2= paste0("./analysis/gimms_overshoot_",out_label,"/",out_label,"_drought_impact_method2_",
                t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
  evaluate_drought_impact(t_range=t_range,outfile=file2)

}


# 1.4 overshoot identification ----------------------------------------
### this is the second method used to identify overshoot, the results from the first step is further filtered.
##  two algorithms are adopted here. 
#   1. the overshoot related component should have overall larger effect than any other lagged effects
#   2. the overshoot component should be smaller than the standard deviation 
##   results expected: 1. whether an event is overshoot event or not. 
##                     2. the preseason contribution (which time period contribute most) 

### drought total ndvi decrease : drought_ndvi_sum
### drought overshoot component : comp_drought_ndvi
### drought overshoot identification:overshoot_drought_boolean

identify_overshoot_first_method<-function(vec){
  #comp_drought_ndvi,overshoot_drought_boolean
  overshoot_drought_boolean<-vec[1:200]
  comp_drought_ndvi<-vec[201:400]
  dim(overshoot_drought_boolean)<-c(50,4)
  dim(comp_drought_ndvi)<-c(50,4)
  overshoot_comp<-comp_drought_ndvi
  overshoot_comp[overshoot_drought_boolean==0]<-NA
  total_overshoot_comp<-apply(overshoot_comp,1,sum,na.rm=T)
  
  nonovershoot_comp<-comp_drought_ndvi
  nonovershoot_comp[overshoot_drought_boolean==1]<-NA
  total_nonovershoot_comp<-apply(nonovershoot_comp,1,sum,na.rm=T)
  ### if total overshoot comp is greater than nonovershoot comp, it is considered as an overshoot event
  overshoot_boolean<-abs(total_overshoot_comp)>abs(total_nonovershoot_comp)
  overshoot_boolean[is.na(overshoot_boolean)]<-FALSE
  overshoot_drought_boolean[overshoot_boolean==FALSE,]<-FALSE
  return(c(overshoot_boolean,as.vector(overshoot_drought_boolean)))
}

get_overshoot<-function(t_range,infiles,outfile){
  overshoot_drought_boolean_list<-list()
  comp_drought_ndvi_list<-list()
  for (i in 2:5){
    load(infiles[i])
    overshoot_drought_boolean_list[[i-1]]<-aperm(overshoot_drought_boolean,c(2,3,1))
    comp_drought_ndvi_list[[i-1]]<-aperm(comp_drought_ndvi,c(2,3,1))
  }
  overshoot_boolean_all<-do.call(abind, overshoot_drought_boolean_list)
  comp_drought_ndvi<-do.call(abind, comp_drought_ndvi_list)
  
  overshoot_component<-abind(overshoot_boolean_all,comp_drought_ndvi,along=3)
  output<-parApply(cl,overshoot_component,c(1,2),identify_overshoot_first_method)
  overshoot_boolean<-output[1:50,,]
  overshoot_comp_boolean<-output[51:250,,]
  save(overshoot_boolean,overshoot_comp_boolean,
       file=outfile)
}

if(ind_evaluate_overshoot){
  var_id<-c("1","23","46","712","1324","13")

  t_range=1:414
  infiles=paste0("./analysis/gimms_overshoot_",out_label,"/",out_label,"_pre",var_id[1:5],"_overshoot_method2_",
                 paste0(t_range[1],"_",t_range[length(t_range)],"_spei05.RData"))
  outfile=paste0("./analysis/gimms_overshoot_",out_label,"/",out_label,"_overshoot_method2_",
                 t_range[1],"_",t_range[length(t_range)],"_spei05.RData")
  get_overshoot(t_range,infiles,outfile)
}


