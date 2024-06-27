import numpy as np
import matplotlib.pyplot as plt
import glob
import datetime as dt
from datetime import timedelta
import sys
import config
import utils
import os
import pickle

# THIS WILL GENEREATE ABOUT 3GB OF DATA!

p = config.params()

reference_date = utils.constants.reference_date
i_date         = p.i_date #dt.datetime(2015,1,1)                                                                                                                              
f_date         = p.f_date # dt.datetime(2020,12,31)                                                                                                                           
period         = [i_date, f_date]
variables_dic=p.variables_dic
varname_dic   = utils.variables.names
varunit_dic   = utils.variables.units
simulations=p.simulations
variables=variables_dic[0]
for c_var, variable in enumerate(variables) :
    if c_var == 0 :
        str_vars = variable
    else:
        str_vars = str_vars + '-' + variable

#path_in = p.path_in_r
path_data = p.path_data
path_dates = p.path_dates
path_stats = p.path_stats
path_z0 = p.path_z0
path_out=p.path_amf_corrected
print(path_out)
os.system('mkdir -p ' + path_out)
vars_ind={}
for ivar,var in enumerate(variables):
    vars_ind[var]=ivar

stations=np.load(path_stats+'stats-selected.npy')
lat_stats=np.load(path_stats+'stats-lats.npy')
lon_stats=np.load(path_stats+'stats-lons.npy')
height_stats=np.load(path_stats+'stats-heights.npy')
regimes=['all','neutral','stable','unstable']
regimes=['all','neutral','stable','unstable']
sources=['AMF','AMFc-BEL_withoutz0','AMFc-BEL_withz0','GEM']
z0s=['all','low','mid','hig','p01-p1']
seasons=['DJF','MAM','JJA','SON']
sea_months={'DJF':[12,1,2],'MAM':[3,4,5],'JJA':[6,7,8],'SON':[9,10,11]}

dic_z0={}
for sim in simulations:
    with open(path_z0+sim+'_mean_z0_pickle','rb') as s:
        temp=pickle.load(s)
        print(temp.keys())
        print(len(temp))
        z0_t=[]
        for st in stations:
            print(st)
            z0_t.append(temp[st])
    dic_z0[sim]=np.asarray(z0_t)

AMF_stats_z0=np.load(path_stats+'/stats-z0.npy')

"""
cond=((AMF_stats_z0>0.01) & (AMF_stats_z0<0.1))
for sim in simulations:
    print(np.sum(cond))
    cond=cond & ((dic_z0[sim]>0.01) & (dic_z0[sim]<0.1))
"""
for z0 in z0s:
    indices=np.argsort(AMF_stats_z0)
    addout=''
    if z0=='all':
        indtemp1=indices
    if z0=='low':
        indtemp1=indices[:9]
        addout=z0
    if z0=='mid':
        indtemp1=indices[9:19]
        addout=z0
    if z0=='hig':
        indtemp1=indices[19:]
        addout=z0
    if z0=='p01-p1':
        indtemp1=((AMF_stats_z0>0.01) & (AMF_stats_z0<0.1))
        addout=z0

    for source in sources:
        if source[:3]=='AMF':
            simulations_tmp=['']
        else:
            simulations_tmp=simulations
        for sim in simulations_tmp:
            paths=[]
            if source=='GEM':
                if z0=='p01-p1':
                    indtemp2=((dic_z0[sim]>0.01) & (dic_z0[sim]<0.1))
                else:
                    indtemp2=indtemp1
                for stat in stations[indtemp2]:
                    paths.append(glob.glob(path_data+source+'_data_*'+stat+'_2015-2020_'+str_vars+'_'+sim+'_1.npy')[0])
            else:
                for stat in stations[indtemp1]:
                    paths.append(glob.glob(path_data+source+'_data_*'+stat+'_2015-2020_'+str_vars+'_1.npy')[0])

            print('number of stations: ',len(paths))
            print('Processing simulation: ',source+' - '+sim)
            stat_data=[]
            for j,path in enumerate(paths):
                print(path)
                data_stat = np.load(path,allow_pickle=True)
                if source[:3]=='AMF':
                    if 'AMFc' in path:
                        # I don't really like how this is done.
                        station=path.split("_")[3]+"_"+path.split("_")[4]
                        dates_stat = np.load(path_dates+f"AMF_dates_{station}_2015-2020_WS-USTAR-ZL-TA-H-LE-PA-RH_1.npy")
                    else:
                         dates_stat = np.load(path.replace('data','dates'))
                else:
                    station=path.split("_")[2]+"_"+path.split("_")[3]
                    dates_stat = np.load(path_dates+f"AMF_dates_{station}_2015-2020_WS-USTAR-ZL-TA-H-LE-PA-RH_1.npy") #np.load(path.replace('_data_','_dates_')) Use amf dates

                stat_data.append(dates_stat.shape[0])
                if j==0:
                    print(data_stat.shape)
                    data_sim=data_stat[:,:]
                    dates_sim=dates_stat[:]
                else:
                    data_sim=np.concatenate((data_sim,data_stat[:,:]),axis=0)
                    dates_sim=np.concatenate((dates_sim,dates_stat[:]),axis=0)

            # Separate regimes
            for reg in regimes:
                select=utils.get_stability(data_sim[:,vars_ind['ZL']],reg)
                data_temp=data_sim[select,:]
                date_temp=dates_sim[select]
                np.save(path_out+source+'_data_'+sim+'_'+reg+'_z0'+addout,data_temp)
                np.save(path_out+source+'_dates_'+sim+'_'+reg+'_z0'+addout,date_temp)
                mins=[]
                hour=[]
                month=[]
                for secs in date_temp:
                    date=reference_date + timedelta(seconds=int(secs))
                    month.append(date.month)
                    hour.append(date.hour)
                    mins.append(date.minute)
                np.save(path_out+source+'_month_'+sim+'_'+reg+'_z0'+addout,month)
                np.save(path_out+source+'_hours_'+sim+'_'+reg+'_z0'+addout,hour)
                np.save(path_out+source+'_mins_'+sim+'_'+reg+'_z0'+addout,mins)

                #Seasonal cycle
                m_mean=[]
                for mm in np.arange(1,13):
                    indices=mm==month
                    m_mean.append(np.mean(data_temp[indices,:],axis=0))
                np.save(path_out+source+'_monthly-mean_'+sim+'_'+reg+'_z0'+addout,m_mean)
                #print(m_mean)

                # Diurnal cycle                                                                                                                                           
                m_mean=[]
                for mm in np.arange(0,24):
                    indices=mm==hour
                    m_mean.append(np.mean(data_temp[indices,:],axis=0))
                np.save(path_out+source+'_hourly-mean_'+sim+'_'+reg+'_z0'+addout,m_mean)
                print(path_out+source+'_hourly-mean_'+sim+'_'+reg+'_z0'+addout)
                print('Number of data: ',len(indices))

                # Diurnal cycle for different seasons                                                          
                mms=np.asarray(month)
                for sea in seasons:
                    m_mean=[]
                    s_mean=[]
                    m_indices=np.in1d(mms, sea_months[sea])
                    temp=data_temp[m_indices,:]
                    h_indices=np.asarray(hour)[m_indices]
                    for mm in np.arange(0,24):
                        indices=mm==h_indices
                        m_mean.append(np.mean(temp[indices,:],axis=0))
                        s_mean.append(np.std(temp[indices,:],axis=0))
                    np.save(path_out+source+'_hourly-mean_'+sea+'_'+sim+'_'+reg+'_z0'+addout,m_mean)
                    np.save(path_out+source+'_hourly-std_'+sea+'_'+sim+'_'+reg+'_z0'+addout,s_mean)
                    print(path_out+source+'_hourly-mean_'+sea+'_'+sim+'_'+reg+'_z0'+addout)
                    print('Number of data: ',m_indices.sum())

np.save(path_out+'station_data-number'+'_z0'+addout,stat_data)


