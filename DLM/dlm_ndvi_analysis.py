### dlm_ndvi_analysis
import numpy as np
from netCDF4 import Dataset
from numpy import matlib
from datetime import datetime
from dlm_functions import forwardFilteringM, Model, perpixel_vi,perpixel_randomvi_vary,perpixel_randomspei,de_seasonalize,deseasonalize_dlm
import pandas as pd
import sys
import random
from multiprocessing import Pool
#import dask.dataframe as dd
#from dask.multiprocessing import get
#%%
spin_up_years = int(sys.argv[1])
spin_up_rep = int(sys.argv[2])
delta_v = int(sys.argv[3])
ind_random = int(sys.argv[4])   # for ind_v >=3 random_seed = ind_v -2

deltas = [0.95,0.98,0.99,0.995,0.999]
delta = deltas[delta_v]
indis = ['prec','prec_temp','prec_temp_solr','randomndvi']
window_sizes = [2,6,24]

if ind_random>2:
    ind = indis[3]
    window_ind = int(np.ceil((ind_random-2)/5.0)-1)
    window_size = window_sizes[window_ind]
    rand_seed = np.mod(ind_random-3,5)+1              
    out_label = ind+str(window_ind)+str(rand_seed)
else:
    ind = indis[ind_random]
    out_label = ind
print(out_label)


print("spin_up_years: ",spin_up_years)
print("spin_up_repeat: ",spin_up_rep)
print("delta: ",delta)
print("dataset: ",ind)

if ind_random>2:
    print("random_seed: ",rand_seed)

rseas = [1,2,3]

#### for first three experiments, only 'ind_process_dlm' is true, 
# for randomize experiments, the 'ind_random_spei' and 'ind_random_id' should also be true,
# these two additional processes are to provide the order of each scene after randomized, 
# so that they can be traced back to their orginal timestamp.

ind_process_dlm=True  #True

ind_random_spei=False  #True
ind_random_id=False

#%%
def parallelize_dataframe(df, fun, n_cores=4):
    df_split = np.array_split(df, n_cores,axis=1)
    pool = Pool(n_cores)
    if fun==perpixel_randomvi_vary:
        df = np.concatenate(pool.map(unpacking_apply_along_axis_randomvi, df_split),axis=1)
    elif fun==perpixel_vi:
        df = np.concatenate(pool.map(unpacking_apply_along_axis_vi, df_split),axis=1)
    elif fun==perpixel_randomspei:
        df = np.concatenate(pool.map(unpacking_apply_along_axis_randomspei, df_split),axis=1)
    pool.close()
    pool.join()
    return df


def unpacking_apply_along_axis_vi(arr):
    #print(lag_mon,spin_up_years)
    return np.apply_along_axis(perpixel_vi,0, arr,spin_up_years=spin_up_years,spin_up_rep=spin_up_rep,delta=delta,rseas=rseas)
def unpacking_apply_along_axis_randomvi(arr):
    #print(lag_mon,spin_up_years)
    return np.apply_along_axis(perpixel_randomvi_vary,0, arr,spin_up_years=spin_up_years,spin_up_rep=spin_up_rep,delta=delta,rseas=rseas,rand_seed = rand_seed, window_size=window_size)
def unpacking_apply_along_axis_randomspei(arr):
    #print(lag_mon,spin_up_years)
    return np.apply_along_axis(perpixel_randomspei,0, arr,rand_seed = rand_seed, window_size=window_size)


def exportnc(results, outfile):
    n_regvar = 8
    n_var=2+n_regvar+len(rseas)*2
    out_dim = int(results.shape[0]/(n_var*2+2+n_regvar+1))

    slik = results[0:out_dim,:]
    sm = results[out_dim:((n_var+1)*out_dim),:]
    sC = results[((n_var+1)*out_dim):((n_var*2+1)*out_dim),:]
    snu = results[((n_var*2+1)*out_dim):((n_var*2+2)*out_dim),:]
    Xpre = results[((n_var*2+2)*out_dim):((n_var*2+2+n_regvar)*out_dim),:]    
    Yout = results[((n_var*2+2+n_regvar)*out_dim):((n_var*2+3+n_regvar)*out_dim),:]
 
    slik = slik.reshape((out_dim,360,720))
    sm = sm.reshape((n_var,out_dim,360,720))
    sC = sC.reshape((n_var,out_dim,360,720))
    snu = snu.reshape((out_dim,360,720))
    Xpre = Xpre.reshape((n_regvar,out_dim,360,720))
    Yout = Yout.reshape((out_dim,360,720))

    f = Dataset(outfile,'w',format="NETCDF4")
    lon = np.arange(-180+0.25,180-0.25+0.01,0.5)
    lat = np.arange(-90+0.25,90-0.25+0.01,0.5)
    ti = np.arange(0,out_dim)
    va = np.arange(0,n_var)
    rvar = np.arange(0,n_regvar)

    f.createDimension('lon',len(lon))
    f.createDimension('lat',len(lat))
    f.createDimension('tim',len(ti))
    f.createDimension('var',n_var)
    f.createDimension('regvar',n_regvar)

    longitude = f.createVariable('Longitude','f4','lon')
    latitude  = f.createVariable('Latitude','f4','lat')
    variable = f.createVariable('Variable','i4','var')
    time = f.createVariable('Time','i4','tim')
    regvar = f.createVariable('Regvar','i4','regvar')

    likelihood = f.createVariable('Likelihood','f4',('tim','lat','lon'),zlib=True)
    predictmean = f.createVariable('Predictmean','f4',('var','tim','lat','lon'),zlib=True)
    predictvariance = f.createVariable('Predictvariance','f4',('var','tim','lat','lon'),zlib=True)
    degreefreedom = f.createVariable('Degreefreedom','f4',('tim','lat','lon'),zlib=True)
    regrevariable = f.createVariable('Regrevariable','f4',('regvar','tim','lat','lon'),zlib=True)
    yanomaly = f.createVariable('Yanomaly','f4',('tim','lat','lon'),zlib=True)
    longitude[:] = lon
    latitude[:] = lat
    variable[:] = va
    time[:] = ti
    regvar[:] = rvar
    likelihood[:,:,:] = slik
    predictmean[:,:,:,:] = sm
    predictvariance[:,:,:,:] = sC
    degreefreedom[:,:,:] = snu 
    regrevariable[:,:,:,:] = Xpre
    yanomaly[:,:,:] = Yout
    f.close()



def exportspeinc(results, outfile):
    out_dim = int(results.shape[0])

    f = Dataset(outfile,'w',format="NETCDF4")
    lon = np.arange(-180+0.25,180-0.25+0.01,0.5)
    lat = np.arange(-90+0.25,90-0.25+0.01,0.5)
    ti = np.arange(0,out_dim)

    f.createDimension('lon',len(lon))
    f.createDimension('lat',len(lat))
    f.createDimension('tim',len(ti))

    longitude = f.createVariable('Longitude','f4','lon')
    latitude  = f.createVariable('Latitude','f4','lat')
    time = f.createVariable('Time','i4','tim')

    spei = f.createVariable('spei','f4',('tim','lat','lon'),zlib=True)
    longitude[:] = lon
    latitude[:] = lat
    time[:] = ti
    spei[:,:,:] = results
    f.close()


#%%
### using previous month prec for current ndvi  ## 0 using current prec for current ndvi
# NDVI dataset 
#root = '/Users/YaoZhang/'
#root = '/rigel/glab/users/zy2309/'
root = "/global/scratch/yaozhang/Project/overshooting/"

if ind_process_dlm:
    ndvi_file = root+'/Data/GIMMS3g_v1.mon.hd.growingseason.nc'
    #prec_file = root+'/Data/CRU_TS_403_pre.nc'
    ### use GPCC precipitation instead of CRU prec. start from 1891/01/01
    prec_file = '/global/scratch/yaozhang/Data/GPCC_mon/full_data_monthly_v2018_05.nc'

    temp_file = '/global/scratch/yaozhang/Data/CRU/cru_ts4.04.1901.2019.tmp.dat.nc'

    solr_files = ['/global/scratch/yaozhang/Data/CRUNCEP_v7_Mon/CRU_NCEP_V7_1981_1990_solr.nc',
        '/global/scratch/yaozhang/Data/CRUNCEP_v7_Mon/CRU_NCEP_V7_1991_2000_solr.nc',
        '/global/scratch/yaozhang/Data/CRUNCEP_v7_Mon/CRU_NCEP_V7_2001_2010_solr.nc',
        '/global/scratch/yaozhang/Data/CRUNCEP_v7_Mon/CRU_NCEP_V7_2011_2016_solr.nc']

    with Dataset(ndvi_file,'r') as fin:   # starts from 1981/7
        ndvi = fin.variables["ndvi"][:]
    with Dataset(prec_file,'r') as fin:   # starts from 1979/7
        prec = np.flip(fin.variables["precip"][:][1062:1500,:,:],axis=1)
    with Dataset(temp_file,'r') as fin:
        temp = fin.variables['tmp'][:][966:1380,:,:]
    for i in range(len(solr_files)):
        with Dataset(solr_files[i],'r') as fin:
            solr_temp = fin.variables['solr'][:]
            print(solr_temp.shape)
            if i>0:
                solr = np.concatenate((solr,solr_temp),axis=2)
            else:
                solr = solr_temp
    
    solr = np.transpose(solr[:,:,6:420], axes=[2, 0, 1])
    print(ndvi.shape)
    print(prec.shape)
    print(solr.shape)

    fill_value = -999
    ndvi = np.where(np.logical_or(ndvi== -3000, ndvi== -9999), np.nan, ndvi)
    
    ndvi = ndvi/10000.0
    
    obs_ndvi = ndvi.shape[0]
    print(obs_ndvi)

    if ind!='prec_temp_solr':
        solr[:] = 0
    if ind=='prec':
        temp[:] = 0


    combined = np.concatenate((ndvi,prec[0:(obs_ndvi+24),:,:],temp,solr),axis=0)
    df = combined.reshape([obs_ndvi*4+24,360*720])
    print(df.shape)
    
    for i in range(4):
        if ind=='randomndvi':
            resultstemp =  parallelize_dataframe(df[:,(i*90*720):((i+1)*90*720)], perpixel_randomvi_vary, n_cores=20)
            if i==0:
                results = resultstemp
            else:
                results = np.concatenate((results,resultstemp),axis=1)
        else:
            resultstemp =  parallelize_dataframe(df[:,(i*90*720):((i+1)*90*720)], perpixel_vi, n_cores=20)            
            if i==0:
                results = resultstemp
            else:
                results = np.concatenate((results,resultstemp),axis=1)

    ### export to nc
    outfile = root+'/analysis/gimms/gimms.DLM.results.'+out_label+'.spinup'+str(spin_up_years)+'.rep'+str(spin_up_rep)+'.delta'+str(delta)+'.nc'
    results = np.where(results==-9999,np.nan,results)
    exportnc(results, outfile)


if ind_random_spei:
    ### also create the spei_random
    spei_file = root+'/Data/spei03.nc'
    with Dataset(spei_file,'r') as fin:   # starts from 1981/7
        spei = fin.variables["spei"][:][966:1380,:,:]
    spei = np.where(np.absolute(spei)> 1000, np.nan, spei)
    df = spei.reshape([414,360*720])
    
    spei_random = parallelize_dataframe(df, perpixel_randomspei, n_cores=20)
    spei_random = spei_random.reshape([414,360,720])
    spei_out_file = root+'/Data/spei03_'+out_label+'.nc'
    exportspeinc(spei_random,spei_out_file)

if ind_random_id:
    ### also create the spei_random
    ind = np.arange(1,415)
    rand_ind_out = perpixel_randomspei(ind, rand_seed=rand_seed, window_size=window_size)
    ind_out_file = root+'/Data/id_'+out_label+'.nc'

    out_dim = len(rand_ind_out)

    f = Dataset(ind_out_file,'w',format="NETCDF4")
    ti = np.arange(0,out_dim)

    f.createDimension('tim',len(ti))
    time = f.createVariable('Time','i4','tim')
    indi = f.createVariable('indi','f4',('tim'),zlib=True)
    time[:] = ti
    indi[:] = rand_ind_out
    f.close()

