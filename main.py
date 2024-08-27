import config
import utils
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as colors
import glob
from tqdm import tqdm
from datetime import datetime, timedelta
import xarray as xr
import warnings
import copy as cp
import sys
import pdb

warnings.filterwarnings("ignore") # Ignore annoying matplotlib warnings


def main():
    # Setup con
    p = config.params()
    simulation_names = p.simulations  
    variables_dic  = p.variables_dic
    varname_dic   = utils.variables.names
    varunit_dic   = utils.variables.units
    simulations=p.simulations
    terms=[ 'tot', 'freq', 'int', 'residual' ]
    variables=variables_dic[0]
    for c_var, variable in enumerate(variables) :
        if c_var == 0 :
            str_vars = variable
        else:
            str_vars = str_vars + '-' + variable
    
    # Setup paths
    path_in_r = p.path_data # Raw data
    path_in = p.path_amf_corrected # This is AMF_corrected
    path_stats = p.path_stats
    path_z0 = p.path_z0
    path_dates = p.path_dates
    path_out=p.path_out
    os.system('mkdir -p ' + path_out)

    vars_ind={}
    for ivar,var in enumerate(variables):
        vars_ind[var]=ivar
    stations=np.load(path_stats+'stats-selected.npy')
    
    term_names=utils.constants.term_names
    term_error_names=utils.constants.term_error_names

    # Plotting params
    sim_name=p.sim_name
    sources=['AMFc-BEL_withz0','GEM']
    dpi=300
    lw=5
    plt.rcParams.update({'font.size':19})
    plt.rcParams.update({'legend.fontsize':15})
    plt.rcParams['figure.constrained_layout.use'] = True
    
    term_names=utils.constants.term_names
    term_error_names=utils.constants.term_error_names
    regimes=p.regimes
    sea_months=utils.constants.sea_months

    # Get roughness length
    # Roughness lengths                                                                                                                                                              
    dic_z0={}
    for sim in simulations:
        with open(path_z0+sim+'_mean_z0_pickle','rb') as s:
            temp=pickle.load(s)
            z0_t=[]
            for st in stations:
                z0_t.append(temp[st])
        dic_z0[sim]=np.asarray(z0_t)

    
    # Compute regime error
    months={}
    hours={}
    minutes={}
    data_gem={}
    m_mean={}
    h_mean={}
    z0='all'
    biased=['biased','unbiased']

    for bb in biased:
        print(bb," case")
        add=bb
        addout=''
        if z0!='all':
            addout=z0
        for source in sources:
            if source[:3]=='AMF':
                simulations_tmp=['']
            else:
                simulations_tmp=simulations
            for sim in simulations_tmp:
                for reg in regimes:
                    #print('Processing simulation: ',source+' - '+sim)
                    data_gem[source+'-'+reg+'-'+sim+'-'+z0]=np.load(path_in+source+'_data_'+sim+'_'+reg+'_z0'+addout+'.npy',allow_pickle=True)
                    m_mean[source+'-'+reg+'-'+sim+'-'+z0]=np.load(path_in+source+'_monthly-mean_'+sim+'_'+reg+'_z0'+addout+'.npy',allow_pickle=True)
                    h_mean[source+'-'+reg+'-'+sim+'-'+z0]=np.load(path_in+source+'_hourly-mean_'+sim+'_'+reg+'_z0'+addout+'.npy',allow_pickle=True)
                    months[source+'-'+reg+'-'+sim+'-'+z0]=np.load(path_in+source+'_month_'+sim+'_'+reg+'_z0'+addout+'.npy',allow_pickle=True)
                    hours[source+'-'+reg+'-'+sim+'-'+z0]=np.load(path_in+source+'_hours_'+sim+'_'+reg+'_z0'+addout+'.npy',allow_pickle=True)
                    minutes[source+'-'+reg+'-'+sim+'-'+z0]=np.load(path_in+source+'_mins_'+sim+'_'+reg+'_z0'+addout+'.npy',allow_pickle=True)

        bins_all={}
        percentiles=np.arange(0,120,20)
        #print ('CALCULATING RANGE FOR VARIABLE: ')
        for c_var,var in enumerate(variables):
            #print(var)
            all_dat_t=np.concatenate((data_gem['AMFc-BEL_withz0-all-'+'-'+z0][:,c_var],data_gem['GEM-all-NAM-11m_GEM511_ERA5_CLASS_USGS_Delage_SN8'+'-'+z0][:,c_var])).ravel()
            minv=np.floor(np.min(all_dat_t))
            bins=[]
            bins.append(minv)
            for c_per,per in enumerate(percentiles[1:]):
                bins.append(np.percentile(all_dat_t,per))
                cond= (all_dat_t>=bins[c_per]) & (all_dat_t<bins[c_per+1])
                #print('# events between '+str(bins[c_per])+' and '+str(bins[c_per+1])+': '+str(np.sum(cond)))
            bins_all[var]=np.asarray(bins)
            
        bins_all['ZL'][1]=-.4
        bins_all['ZL'][2]=-.02
        bins_all['ZL'][3]=.02
        bins_all['ZL'][4]=.4
        
        variables_plot = ['WS','ZL','USTAR']
        reg='all'
        add=z0+'_'+bb
        var_indices_plot={}
        add_met=[]
        for source in tqdm(sources[1:]):
            if source[:3]=='AMF':
                simulations_tmp=['']
            else:
                simulations_tmp=simulations
            for sim in simulations_tmp:
                str_plot=''
                var_amf=[]
                var_gem=[]
                for c_var, var in enumerate(variables_plot):
                    str_plot = str_plot + var
                    var_amf.append(data_gem['AMFc-BEL_withz0-'+reg+'-'+'-'+z0][:, variables.index(var)])
                    var_gem.append(data_gem[source+'-'+reg+'-'+sim+'-'+z0][:, variables.index(var)])
                    var_indices_plot[var] = c_var

                var_amf = np.asarray(var_amf)
                var_gem = np.asarray(var_gem)
                if bb=='unbiased':
                    var_amf[var_indices_plot['WS'],:]=var_amf[var_indices_plot['WS'],:]-np.mean(var_amf[var_indices_plot['WS'],:])
                    var_gem[var_indices_plot['WS'],:]=var_gem[var_indices_plot['WS'],:]-np.mean(var_gem[var_indices_plot['WS'],:])
            
                freq_bin,int_bin,max_bin,max5_bin,std_bin,per_bin=utils.calculate_binning_mean(var_gem,bins_all,variables_plot,percentiles,prob=True)
                freq_bin_a,int_bin_a,max_bin,max5_bin,std_bin,per_bin=utils.calculate_binning_mean(var_amf,bins_all,variables_plot,percentiles,prob=True)
                #print(int_bin_a)
                var_errors      = np.zeros((len(terms),freq_bin.shape[0],freq_bin.shape[1]))
                var_errors      = utils.error_decomposition_taylor(var_errors,terms,freq_bin,int_bin,freq_bin_a,int_bin_a)
                var_errors      = utils.error_illdefined(var_errors,terms,freq_bin,int_bin,freq_bin_a,int_bin_a)

                var_perc_errors      = np.zeros((len(terms),freq_bin.shape[0],freq_bin.shape[1]))
                var_perc_errors = utils.percent_error_decomposition_taylor(var_perc_errors,terms,freq_bin,int_bin,freq_bin_a,int_bin_a)
                var_perc_errors      = utils.error_illdefined(var_perc_errors,terms,freq_bin,int_bin,freq_bin_a,int_bin_a)

                factor          = 1
                var_errors      = var_errors/factor
                var_perc_errors      = var_perc_errors/factor
                
                for term in ['int','freq','tot'] :
                    if term=='int':
                        cmap = plt.get_cmap('viridis')
                    else:
                        cmap = plt.get_cmap('cubehelix_r')
                    if term=='tot':
                        cmap = plt.get_cmap('cool')
                    cmap.set_bad( "gray", alpha=0 )
                    if term=='int':
                        Z=int_bin_a.T
                        tt=8.
                        levels = MaxNLocator(nbins=11).tick_values(0,tt)
                        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                        plt.title('AMF (' + str(np.round(np.ma.mean(Z),1)) +' '+ varunit_dic[variables_plot[0]] + ')')
                    if term=='freq':
                        Z=freq_bin_a.T
                        norm=colors.LogNorm(vmin=0.0005, vmax=.5, clip=True)
                        plt.title('AMF (' + str(int(np.sum(Z))) +')')
                    if term=='tot':
                        Z=int_bin_a.T*freq_bin_a.T
                        tt=Z.max()
                        levels = MaxNLocator(nbins=11).tick_values(0,tt)
                        norm=colors.LogNorm(vmin=0.00001, vmax=.5, clip=True)
                        plt.title('AMF (' + str(np.round(np.ma.mean(Z),1)) +' '+ varunit_dic[variables_plot[0]] + ')')
                    #print(Z)
                    
                    plt.ylabel(varname_dic[variables_plot[2]] + ' [' + varunit_dic[variables_plot[2]] + ']')
                    plt.xlabel('Stability [-]')
                    plt.title('')
                    labelsx = []
                    aa = bins_all[variables_plot[1]]
                    for ii in np.arange(aa.shape[0]-1) :
                        aa1=str("{:.0e}".format(aa[ii]))
                        aa2=str("{:.0e}".format(aa[ii+1]))
                        labelsx.append('[' + aa1 + ',' + aa2+ ')')
                    labelsx=['very unstable','unstable','neutral','stable','very stable']
                    labelsx=['VU','U','N','S','VS']
                    labelsy = []
                    aa = bins_all[variables_plot[2]]
                    for ii in np.arange(aa.shape[0]-1):
                        labelsy.append(str(round(aa[ii+1],2)))

                    plt.xticks(np.arange(len(labelsx))+.5, labelsx, rotation=0)
                    plt.yticks(np.arange(len(labelsy))+.5, labelsy, rotation=0)
                    plt.pcolormesh(Z, cmap=cmap, norm=norm)
                    for m in range(5):
                        for h in range(5):
                            plt.text(m+0.5,h+0.5,f'{Z[h, m]:.2f}', ha='center', va='center', color='black',fontsize=11)
                    cbar = plt.colorbar(extend='max')
                    if term in ['int','tot'] :
                        cbar.set_label(term_names[term] + ' [' + varunit_dic[variables_plot[0]] + ']')
                    else :
                        cbar.set_label(term_names[term] + ' ' )
                    fig_name = path_out+'binning_mean_AMF_'+term+'_'+add+'.png'
                    plt.savefig(fig_name, bbox_inches='tight')
                    plt.close()
                
                for c_term,term in enumerate(terms):
                    Z = var_perc_errors[c_term,:].T*100
                    Z_v = var_errors[c_term,:].T
                    if term == 'tot' :
                        Zmax = np.ceil(np.abs(Z).max())
                    Zmax=100
                    levels = MaxNLocator(nbins=11, symmetric=True).tick_values(-1*Zmax, Zmax)
                    levels = levels[levels != 0]
                    cmap   = plt.get_cmap('PiYG')
                    norm   = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                    cmap.set_bad( "gray", alpha=0 )
                    plt.title(term_error_names[term] + r' ($\epsilon_{percent}=$' + str(np.round(np.ma.mean(abs(Z)),2)) + ' '+varunit_dic[variables_plot[0]] + ')')
                    plt.title(sim_name[sim])
                    plt.ylabel(varname_dic[variables_plot[2]] + ' [' + varunit_dic[variables_plot[2]] + ']')
                    plt.xticks(np.arange(len(labelsx))+.5, labelsx)
                    plt.yticks(np.arange(len(labelsy))+.5, labelsy)
                    for m in range(5):
                        for h in range(5):
                            plt.text(m+0.5,h+0.5,f'{Z_v[h, m]:.3f}', ha='center', va='center', color='black',fontsize=11)
                    plt.pcolormesh(Z, cmap=cmap, norm=norm)
                    #plt.legend()
                    cbar = plt.colorbar(extend='both')
                    cbar.set_label(term_error_names[term] + ' [' + varunit_dic[variables_plot[0]] + '] (%)')
                    fig_name = path_out+ 'binning_mean_perc_error_' + source + '_' + term + '_' + sim_name[sim] +'_'+add+ '.png'
                    plt.savefig(fig_name, bbox_inches='tight')
                    
                    plt.close()
                    if term == 'tot' :
                        Z = np.sum(abs(var_errors[1:,:,:]),axis=0).T
                        #Z_v = np.sum(abs(var_errors[1:,:,:]),axis=0).T
                        Zmax=0.5
                        levels = MaxNLocator(nbins=11).tick_values(0, Zmax)
                        cmap   = plt.get_cmap('cool')
                        norm=colors.LogNorm(vmin=0.0005, vmax=Zmax, clip=True)
                        cmap.set_bad( "gray", alpha=0 )
                        plt.title(term_error_names[term] + r' ($\epsilon_{reg}=$' + str(np.round(np.ma.mean(abs(Z)),2)) + ' '+varunit_dic[variables_plot[0]] + ')')
                        plt.title(sim_name[sim])
                        plt.ylabel(varname_dic[variables_plot[2]] + ' [' + varunit_dic[variables_plot[2]] + ']')
                        plt.xticks(np.arange(len(labelsx))+.5, labelsx)
                        plt.yticks(np.arange(len(labelsy))+.5, labelsy)
                        plt.pcolormesh(Z, cmap=cmap, norm=norm)
                        for m in range(5):
                            for h in range(5):
                                plt.text(m+0.5,h+0.5,f'{Z[h, m]:.2f}', ha='center', va='center', color='black',fontsize=11)
                        #plt.legend()
                        cbar = plt.colorbar(extend='max')
                        cbar.set_label(r'$\epsilon_{reg}$'+ ' [' + varunit_dic[variables_plot[0]] + ']')
                        fig_name = path_out+ 'binning_mean_perc_reg-error_' + source + '_' + term + '_' + sim_name[sim] +'_'+add+ '.png'
                        plt.savefig(fig_name, bbox_inches='tight')
                        plt.close()

                    
                for c_term,term in enumerate(terms):
                    Z = var_errors[c_term,:].T
                    if term == 'tot' :
                        Zmax = np.ceil(np.abs(Z).max())
                    Zmax=.3
                    levels = MaxNLocator(nbins=11, symmetric=True).tick_values(-1*Zmax, Zmax)
                    levels = levels[levels != 0]
                    cmap   = plt.get_cmap('PiYG')
                    norm   = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                    cmap.set_bad( "gray", alpha=0 )
                    plt.title(term_error_names[term] + r' ($\epsilon_{abs}=$' + str(np.round(np.ma.mean(abs(Z)),2)) + ' '+varunit_dic[variables_plot[0]] + ')')
                    plt.title(sim_name[sim])
                    plt.ylabel(varname_dic[variables_plot[2]] + ' [' + varunit_dic[variables_plot[2]] + ']')
                    plt.xticks(np.arange(len(labelsx))+.5, labelsx)
                    plt.yticks(np.arange(len(labelsy))+.5, labelsy)
                    for m in range(5):
                        for h in range(5):
                            plt.text(m+0.5,h+0.5,f'{Z[h, m]:.2f}', ha='center', va='center', color='black',fontsize=11)
                    plt.pcolormesh(Z, cmap=cmap, norm=norm)
                    #plt.legend()
                    cbar = plt.colorbar(extend='both')
                    cbar.set_label(term_error_names[term] + ' [' + varunit_dic[variables_plot[0]] + ']')
                    fig_name = path_out+ 'binning_mean_error_' + source + '_' + term + '_' + sim_name[sim] +'_'+add+ '.png'
                    plt.savefig(fig_name, bbox_inches='tight')
                    
                    plt.close()
                    if term == 'tot' :
                        Z = np.sum(abs(var_errors[1:,:,:]),axis=0).T
                        Zmax=.5
                        levels = MaxNLocator(nbins=11).tick_values(0, Zmax)
                        cmap   = plt.get_cmap('cool')
                        norm=colors.LogNorm(vmin=0.0005, vmax=Zmax, clip=True)
                        cmap.set_bad( "gray", alpha=0 )
                        plt.title(term_error_names[term] + r' ($\epsilon_{reg}=$' + str(np.round(np.ma.mean(abs(Z)),2)) + ' '+varunit_dic[variables_plot[0]] + ')')
                        plt.title(sim_name[sim])
                        plt.ylabel(varname_dic[variables_plot[2]] + ' [' + varunit_dic[variables_plot[2]] + ']')
                        plt.xticks(np.arange(len(labelsx))+.5, labelsx)
                        plt.yticks(np.arange(len(labelsy))+.5, labelsy)
                        plt.pcolormesh(Z, cmap=cmap, norm=norm)
                        for m in range(5):
                            for h in range(5):
                                plt.text(m+0.5,h+0.5,f'{Z[h, m]:.2f}', ha='center', va='center', color='black',fontsize=11)
                        #plt.legend()
                        cbar = plt.colorbar(extend='max')
                        cbar.set_label(r'$\epsilon_{reg}$'+ ' [' + varunit_dic[variables_plot[0]] + ']')
                        fig_name = path_out+ 'binning_mean_reg-error_' + source + '_' + term + '_' + sim_name[sim] +'_'+add+ '.png'
                        plt.savefig(fig_name, bbox_inches='tight')
                        plt.close()
                        add_met.append(np.mean(Z))

        np.save('add_reg_'+bb+'.npy',add_met) 

    # Computing the mean, diurnal, and annual cycle errors
    # Get amf paths
    amf_data_path = glob.glob(path_in_r+"/AMFc-BEL_withz0_data_*.npy")
    amf_date_path = glob.glob(path_dates+"/AMF_dates_*.npy")


    # Sort paths
    amf_data_path.sort()
    amf_date_path.sort()
    N = len(amf_data_path)

    gem_ws = {} # Store u to compute mean and show dist.
    gem_int_per_station = {} # to store u_{h,m} intensity
    gem_freq_per_station = {} # to store u_{h,m} freq
    gem_diu_per_station = {} # to store u_{h} 
    gem_ann_per_station = {} # to store u_{m} 

    for sim in tqdm(simulation_names, desc='Computing errors', unit='simulation'):
        amf_ws = []
        amf_int_per_station = []
        amf_freq_per_station = []
        amf_diu_per_station = []
        amf_ann_per_station = []
        counter=0
        for date_p,data_p in zip(amf_date_path, amf_data_path):
            # Load data
            data = np.load(data_p)
            dates = np.load(date_p)

            # Switch paths to gem
            data_p_gem = data_p.replace("AMFc-BEL_withz0_data","GEM_data").replace("1.npy",sim+"_1.npy")
            #date_p_gem = date_p.replace("AMF_dates","GEM_dates").replace("1.npy",sim+"_1.npy")
            data_gem = np.load(data_p_gem)
            dates_gem =  dates #np.load(date_p_gem) # Use AMF dates

            # Correct dates
            reference_date=datetime(1971,1,1)
            corrected_dates = []
            corrected_dates_gem = []

            for i in range(len(dates)):
                first_date_timestamp = dates[i]
                first_date_timestamp_gem = dates_gem[i]
                timestamp_delta = timedelta(seconds=first_date_timestamp)
                timestamp_delta_gem = timedelta(seconds=first_date_timestamp_gem)
                corrected_dates.append(reference_date + timestamp_delta)
                corrected_dates_gem.append(reference_date + timestamp_delta_gem)

            WS = data[:,0]
            WS_gem = data_gem[:,0]

            # Save windspeeds
            amf_ws.append(WS)
            if sim in gem_ws:
                gem_ws[sim].append(WS_gem)
            else:
                gem_ws[sim] = [WS_gem]
                
            # Make dataframes to make quick diurnal and yearly cycles
            ds = xr.Dataset({
                'WS': (('time',), WS),
                }, coords={'time': corrected_dates})


            ds_gem = xr.Dataset({
                    'WS': (('time',), WS_gem),
                    }, coords={'time': corrected_dates_gem})


            # Convert corrected_dates to numpy array
            corrected_dates = np.array(corrected_dates)
            corrected_dates_gem = np.array(corrected_dates_gem)

            # Extract hour and month components from corrected_dates
            hours = np.array([date.hour for date in corrected_dates], dtype=int)
            months = np.array([date.month for date in corrected_dates], dtype=int)
            hours_gem = np.array([date.hour for date in corrected_dates_gem], dtype=int)
            months_gem = np.array([date.month for date in corrected_dates_gem], dtype=int)


            # Convert hour and month components back to datetime objects
            hour_objects = np.array([np.datetime64(f'1970-01-01T{hour:02d}:00:00') for hour in hours])
            month_objects = np.array([np.datetime64(f'1970-{month:02d}-01') for month in months])
            hour_objects_gem = np.array([np.datetime64(f'1970-01-01T{hour:02d}:00:00') for hour in hours_gem])
            month_objects_gem = np.array([np.datetime64(f'1970-{month:02d}-01') for month in months_gem])

            # Create var_comp array with wind speed as the first variable, hour_objects as the second variable, and month_objects as the third variable
            var_comp = np.array([WS, hour_objects, month_objects])
            var_comp_gem = np.array([WS_gem, hour_objects_gem, month_objects_gem])

            # Define bins for hours and months (assuming bins_all is a dictionary)
            bins_all = {'hours': np.arange(0, 25), 'months': np.arange(1, 14)}  
            
            # Define variables (assuming wind speed is the independent variable)
            variables = ['wind_speed', 'hours', 'months']

            # Define percentiles. Don't really need this
            percentiles = [25, 50, 75]

            # Call calculate_binning_mean function
            freq_bin, int_bin, _,_,_,_ = utils.calculate_binning_mean_time(var_comp, bins_all, variables, percentiles)
            freq_bin_gem, int_bin_gem,_,_,_,_ = utils.calculate_binning_mean_time(var_comp_gem, bins_all, variables, percentiles)
            
            amf_int_per_station.append(int_bin)
            amf_freq_per_station.append(freq_bin)
            
            if counter==0:
                amf_diu_per_station=ds.groupby("time.hour").mean()
                amf_ann_per_station=ds.groupby("time.month").mean()
                counter+=1
            else:
                amf_diu_per_station+=ds.groupby("time.hour").mean()
                amf_ann_per_station+=ds.groupby("time.month").mean()

            if sim in gem_int_per_station:
                gem_int_per_station[sim].append(int_bin_gem)
                gem_freq_per_station[sim].append(freq_bin_gem)
                gem_diu_per_station[sim] += ds_gem.groupby("time.hour").mean()
                gem_ann_per_station[sim] += ds_gem.groupby("time.month").mean()
            else:
                gem_int_per_station[sim] = [int_bin_gem]
                gem_freq_per_station[sim] = [freq_bin_gem]
                gem_diu_per_station[sim] = ds_gem.groupby("time.hour").mean()
                gem_ann_per_station[sim] = ds_gem.groupby("time.month").mean()
        


    # Compute Errors and make figures
    print("Saving Figures")
    aa=np.concatenate([arr.flatten() for arr in amf_ws])
    XY=float(aa.shape[0])
    X=float(len(amf_ws))
    Y=XY/X

    month_names=['J','F','M','A','M','J','J','A','S','O','N','D']
    hour_names=np.arange(1,25)

    # Diurnal cycle AMF
    diu_amf = amf_diu_per_station['WS']/len(amf_freq_per_station)
    plt.figure()
    plt.plot(hour_names,diu_amf,color='k',label='AMF')
    # Set labels and title                                                                                                                                 
    plt.xlabel('Local Hour')
    plt.xticks(hour_names[::3]+2, hour_names[::3]+2)
    plt.ylabel(r'$\overline{u_{h}}$'+' ('+varunit_dic['WS']+')')
    plt.ylim([3,4.5])
    #plt.legend()
    plt.savefig(path_out+"amf_diu",dpi=dpi)
    plt.close()

    # Annual cycle AMF
    #print(np.amax(amf_ann_per_station['WS']))
    ann_amf = amf_ann_per_station['WS']/len(amf_freq_per_station)
    plt.figure()
    plt.plot(np.arange(1,13),ann_amf,color='k',label='AMF')
    plt.xlabel('Month')
    plt.ylabel(r'$\overline{u_{m}}$'+' ('+varunit_dic['WS']+')')
    plt.ylim([3,4.5])
    #plt.legend()
    plt.xticks(np.arange(1,13), month_names)
    plt.savefig(path_out+"amf_ann",dpi=dpi)
    plt.close()

    # Intensity and freq for AMF
    int_bin_amf = np.mean(amf_int_per_station,axis=0)
    freq_bin_amf = np.sum(amf_freq_per_station,axis=0)/len(amf_freq_per_station)

    # Create the heatmap
    plt.figure(dpi=dpi)
    plt.imshow(int_bin_amf, cmap='viridis', aspect='auto',vmin=0, vmax=5.)

    # Add color bar
    plt.colorbar(label=r'$\overline{u_{h,m}}$'+' ('+varunit_dic['WS']+')')

    # Annotate each cell with the frequency
    for m in range(12):
        for h in range(24):
            plt.text(m,h,f'{int(freq_bin[h, m])}', ha='center', va='center', color='white',fontsize=11)
            
    # Set labels and title
    plt.xlabel('Month')
    plt.ylabel('Local Hour')
    plt.xticks(np.arange(0,12), month_names)
    plt.yticks(hour_names[::3]+2, hour_names[::3]+2)
    plt.title('AMF')
    plt.savefig(path_out+"AMF_WS_ann-diu",dpi=dpi)
    plt.close()

    bbox = dict(boxstyle='round', fc='grey', ec='white', alpha=0.5)

    bias_error = {}
    diu_error = {}
    diu_error_nb = {}
    ann_error = {}
    ann_error_s = {}
    ann_error_nb = {}
    maxWS=30.
    nbins=30
    error_names=['bias','diu','ann','diu-ann','fi','reg']
    error_names=[r'$\epsilon_{bias}$',r'$\epsilon_{diu}$',r'$\epsilon_{ann}$',r'$\epsilon_{ann,diu}$',r'$\epsilon_{fi}$',r'$\epsilon_{reg}$']
    errors_all=np.zeros((len(error_names),len(gem_int_per_station.keys()),2))
    histo_sims={}
    histo_sims_unbiased={}
    for c_key,key in enumerate(gem_int_per_station.keys()):
        ###############################################
        # Mean and distrubtion plotting.
        ###############################################
        # need to flatten
        amf_ws = np.concatenate([arr.flatten() for arr in amf_ws])
        sim_ws = np.concatenate([arr.flatten() for arr in gem_ws[key]])
        bins=np.histogram(amf_ws,range=(0,maxWS),bins=nbins)[1]
        bins_c=0.5*(bins[:-1]+bins[1:])
        
        # biased
        histo_amf=np.histogram(amf_ws,range=(0,maxWS),bins=nbins)[0]
        histo_sims[key]=np.histogram(sim_ws,range=(0,maxWS),bins=nbins)[0]
        error= bins_c*(histo_sims[key]-histo_amf)
        norm=1./XY
        errors_all[4,c_key,0]=norm*np.sum(np.abs(error))

        # unbiased
        histo_amf_unbiased=np.histogram(amf_ws-np.mean(amf_ws),range=(0,maxWS),bins=nbins)[0]
        histo_sims_unbiased[key]=np.histogram(sim_ws-np.mean(sim_ws),range=(0,maxWS),bins=nbins)[0]
        error= 0.5*(bins[:-1]+bins[1:])*(histo_sims_unbiased[key]-histo_amf_unbiased)
        norm=1./XY
        errors_all[4,c_key,1]=norm*np.sum(np.abs(error))

        plt.figure(dpi=dpi)
        plt.step(bins[:-1],histo_amf,where='post',label="AMF",color='k')
        plt.step(bins[:-1],histo_sims[key],where='post',label=p.sim_name[key],color='r')
        plt.axvline(x=np.mean(amf_ws),color='k',ls='--')
        plt.axvline(x=np.mean(sim_ws),color='r',ls='--')
        #plt.legend(loc='best')
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim([0,22])
        plt.savefig(path_out+p.sim_name[key]+"_dist",dpi=dpi)
        plt.close()

        #print("Mean error for sim : "+ p.sim_name[key])
        #print(np.abs(np.mean(sim_ws) - np.mean(amf_ws)))
        #print("\n")
        bias_error[key]=np.abs(np.mean(sim_ws) - np.mean(amf_ws))

        ###############################################
        # diu stuff here
        ###############################################
        gem_diu = gem_diu_per_station[key]['WS']/X
        hours = np.arange(1,25)
        plt.figure(dpi=dpi)
        plt.plot(hours,diu_amf,label="AMF",color='k',lw=lw)
        plt.plot(hours,gem_diu,label=p.sim_name[key],lw=lw,color='r')
        plt.xlabel("Local Hour")
        plt.xticks(hour_names[::3]+2, hour_names[::3]+2)
        plt.savefig(path_out+p.sim_name[key]+"_diu_error",dpi=dpi)
        plt.close()

        diu_error[key] = gem_diu-diu_amf
        diu_error_nb[key] = (gem_diu-np.mean(sim_ws))-(diu_amf-np.mean(amf_ws))
        
        ###############################################
        # ann stuff here
        ###############################################
        L = np.asarray(gem_int_per_station[key]).shape[0]
        gem_ann = gem_ann_per_station[key]['WS']/L
        months = np.arange(1,13)
        plt.figure(dpi=dpi)
        plt.plot(months,ann_amf,label="AMF",color='k',lw=lw)
        plt.plot(months,gem_ann,label=p.sim_name[key],lw=lw,color='r')
        plt.xlabel("Month")
        plt.savefig(path_out+p.sim_name[key]+"_ann_error",dpi=dpi)
        plt.close()

        ann_error[key] = gem_ann-ann_amf
        ann_error_nb[key] = (gem_ann-np.mean(sim_ws))-(ann_amf-np.mean(amf_ws))
        
        ###############################################
        # ann_diu stuff here
        ###############################################
        L = np.asarray(gem_int_per_station[key]).shape[0]
        int_bin_gem = np.sum(np.asarray(gem_int_per_station[key]),axis=0)/L
        freq_bin_gem = np.sum(np.asarray(gem_freq_per_station[key]),axis=0)

        # Compute absolute error
        error = freq_bin_amf*(int_bin_gem-int_bin_amf)
        abs_error = freq_bin_amf*np.abs(int_bin_gem-int_bin_amf)
        norm=1./XY
        errors_all[0,c_key,0]=norm*np.abs(np.sum(error))
        errors_all[1,c_key,0]=norm*np.sum(np.abs(np.sum(error,axis=0)))
        errors_all[2,c_key,0]=norm*np.sum(np.abs(np.sum(error,axis=1)))
        errors_all[3,c_key,0]=norm*np.sum(np.abs(error))
        
        #print("Diurnal error for sim: " + p.sim_name[key])
        #print("\n")

        # Create the heatmap
        plt.figure(dpi=dpi)
        lims=2.5
        plt.imshow(error/freq_bin_amf, cmap='seismic', aspect='auto',vmin=-lims,vmax=lims)
        plt.colorbar(label=r'$\epsilon_{ann,diu}$'+' ('+varunit_dic['WS']+')')
        plt.xlabel('Month')
        plt.ylabel('Local Hour')
        plt.xticks(np.arange(0,12), month_names)
        plt.yticks(hour_names[::3]+2, hour_names[::3]+2)
        plt.title(p.sim_name[key]+ ' - biased')
        plt.savefig(path_out+p.sim_name[key]+"_error_diu-ann",dpi=dpi)
        plt.close()

        # Compute all error
        error = freq_bin_amf*((int_bin_gem-np.mean(sim_ws))-(int_bin_amf-np.mean(amf_ws)))
        abs_error = freq_bin_amf*np.abs(int_bin_gem-np.mean(sim_ws))-(int_bin_amf-np.mean(amf_ws))
        errors_all[0,c_key,1]=norm*np.abs(np.sum(error))
        errors_all[1,c_key,1]=norm*np.sum(np.abs(np.sum(error,axis=0)))
        errors_all[2,c_key,1]=norm*np.sum(np.abs(np.sum(error,axis=1)))
        errors_all[3,c_key,1]=norm*np.sum(np.abs(error))
        #print("Diurnal error for sim: " + p.sim_name[key])
        #print("\n")

        # Create the heatmap
        plt.figure(dpi=dpi)
        plt.imshow(error/freq_bin_amf, cmap='seismic', aspect='auto',vmin=-lims,vmax=lims)
        plt.colorbar(label=r'$\epsilon^{\prime}_{ann,diu}$'+' ('+varunit_dic['WS']+')')
        plt.xlabel('Month')
        plt.ylabel('Local Hour')
        plt.xticks(np.arange(0,12), month_names)
        plt.yticks(hour_names[::3]+2, hour_names[::3]+2)
        plt.title(p.sim_name[key]+ ' - unbiased')
        plt.savefig(path_out+p.sim_name[key]+"_error_diu-ann_nb",dpi=dpi)
        plt.close()

    # Plot diu and ann error cycles
    plt.figure(dpi=dpi)
    plt.plot([0,25],[0,0],color='grey',linewidth=.75)
    for key in diu_error:
        plt.plot(np.arange(1,25),diu_error[key], label=p.sim_name[key], color=p.sim_color[key])
        print(np.mean(diu_error[key]))
    plt.xlabel('Local Hour')
    plt.xticks(hour_names[::3]+2, hour_names[::3]+2)
    plt.ylabel(r"$\epsilon_{diu}$"+' ('+varunit_dic['WS']+')')
    plt.ylim([-1.2,1.2])
    plt.xlim([0,25])
    plt.text(1, 1., 'biased', bbox=bbox)
    plt.savefig(path_out+"diu_error",dpi=dpi)
    
    plt.figure(dpi=dpi)
    plt.plot([0,13],[0,0],color='grey',linewidth=.75)
    for key in ann_error:
        plt.plot(np.arange(1,13),ann_error[key], label=p.sim_name[key], color=p.sim_color[key])
    plt.xlabel('Month')
    plt.ylabel(r"$\epsilon_{ann}$"+' ('+varunit_dic['WS']+')')
    plt.ylim([-1.2,1.2])
    plt.xticks(np.arange(1,13), month_names)
    plt.xlim([0,13])
    plt.text(1, 1., 'biased', bbox=bbox)
    plt.savefig(path_out+"ann_error",dpi=dpi)

    plt.figure(dpi=dpi)
    plt.plot([0,13],[0,0],color='grey',linewidth=.75)
    for key in ann_error_nb:
        plt.plot(np.arange(1,13),ann_error_nb[key], label=p.sim_name[key], color=p.sim_color[key])
    plt.xlabel("Month")
    plt.ylabel(r"$\epsilon^{\prime}_{ann}$"+' ('+varunit_dic['WS']+')')
    plt.ylim([-1.2,1.2])
    plt.xticks(np.arange(1,13), month_names)
    plt.xlim([0,13])
    plt.text(1, 1., 'unbiased', bbox=bbox)
    plt.savefig(path_out+"ann_error_nb",dpi=dpi)

    plt.figure(dpi=dpi)
    plt.plot([0,25],[0,0],color='grey',linewidth=.75)
    for key in diu_error:
        plt.plot(np.arange(1,25),diu_error_nb[key], label=p.sim_name[key], color=p.sim_color[key])
    plt.xlabel("Local Hour")
    plt.xticks(hour_names[::3]+2, hour_names[::3]+2)
    plt.ylabel(r"$\epsilon^{\prime}_{diu}$"+' ('+varunit_dic['WS']+')')
    plt.ylim([-1.2,1.2])
    plt.xlim([0,25])
    plt.text(1, 1., 'unbiased', bbox=bbox)
    plt.savefig(path_out+"diu_error_nb",dpi=dpi)

    # Plot dist in on figure
    plt.figure(dpi=dpi)
    amf_ws = np.concatenate([arr.flatten() for arr in amf_ws])
    plt.step(bins[:-1],histo_amf,where='post',label="AMF",color='black',linewidth=lw)
    plt.xlabel(r'$\overline{u^{o}_k}$'+' ('+varunit_dic['WS']+')')
    plt.ylabel(r'$N_{k}$')
    plt.yscale("log")
    plt.xlim([0,22])
    plt.savefig(path_out+"AMF_fi",dpi=dpi)
    plt.close()

    # Plot dist in on figure
    maxA=0
    plt.figure(dpi=dpi)
    amf_ws = np.concatenate([arr.flatten() for arr in amf_ws])
    for key in gem_int_per_station.keys():
        plt.step(bins[:-1],bins_c*(histo_sims_unbiased[key]-histo_amf_unbiased),where='post',label=p.sim_name[key], color=p.sim_color[key],linewidth=lw/2)
        maxA=np.max((maxA,np.max(np.abs(bins_c*(histo_sims_unbiased[key]-histo_amf_unbiased)))))
    plt.xlabel(r'$\overline{u^{o}_k}$'+' ('+varunit_dic['WS']+')')
    #plt.ylabel(r'$(N^{s}_{k}-N^{o}_{k}) \cdot u^{o}_{k}$'+' ('+varunit_dic['WS']+')')
    plt.ylabel('absolute error ('+varunit_dic['WS']+')')
    plt.xlim([0,22])
    limA=10**np.ceil(np.log10(maxA))
    plt.ylim([-30000,30000])
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.text(15, 24000.,"unbiased", bbox=bbox)
    plt.savefig(path_out+"error_fi_unbiased",dpi=dpi)
    plt.close()
    plt.close('all')
     
    plt.figure(dpi=dpi)
    for key in gem_int_per_station.keys():
        err=(histo_sims_unbiased[key]-histo_amf_unbiased)/(histo_amf_unbiased+histo_sims_unbiased[key])
        plt.step(bins[:-1],100.*err,where='post',label=p.sim_name[key], color=p.sim_color[key],linewidth=lw/2)
    plt.xlabel(r'$\overline{u^{o}_k}$'+' ('+varunit_dic['WS']+')')
    plt.ylabel('relative error (%)')
    plt.xlim([0,22])
    plt.ylim([-100,100])
    #plt.legend()
    plt.text(15, 80.,"unbiased", bbox=bbox)
    plt.savefig(path_out+"error_fi_unbiased-rel")
    plt.close()

    # Plot unbias dist in on figure
    plt.figure(dpi=dpi)
    maxA=0
    amf_ws = np.concatenate([arr.flatten() for arr in amf_ws])
    for key in gem_int_per_station.keys():
        plt.step(bins[:-1],bins_c*(histo_sims[key]-histo_amf),where='post',label=p.sim_name[key], color=p.sim_color[key],linewidth=lw/2)     
        maxA=np.max((maxA,np.max(np.abs(bins_c*(histo_sims[key]-histo_amf)))))
    plt.xlabel(r'$\overline{u^{o}_k}$'+' ('+varunit_dic['WS']+')')
    plt.ylabel(r'$(N^{s}_{k}-N^{o}_{k}) \cdot u^{o}_{k}$'+' ('+varunit_dic['WS']+')')
    plt.ylabel('absolute error ('+varunit_dic['WS']+')')
    plt.xlim([0,22])
    limA=10**np.ceil(np.log10(maxA))
    plt.ylim([-300000,300000])
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.text(17, 160000.,"biased", bbox=bbox)
    #plt.text(.99, .99, 'biased', ha='right', va='bottom', transform=ax.transAxes, bbox=bbox)
    plt.savefig(path_out+"error_fi_biased",dpi=dpi)
    plt.close()
    plt.close('all')

    plt.figure(dpi=dpi)
    for key in gem_int_per_station.keys():
        err=abs(histo_sims[key]-histo_amf)/(histo_amf+histo_sims[key])
        err=(histo_sims[key]-histo_amf)/(histo_amf+histo_sims[key])
        plt.step(bins[:-1],100.*err,where='post',label=p.sim_name[key], color=p.sim_color[key],linewidth=lw/2)
    #plt.yscale("log")
    plt.xlabel(r'$\overline{u^{o}_k}$'+' ('+varunit_dic['WS']+')')
    plt.ylabel('relative error (%)')
    #plt.yscale("log")
    plt.xlim([0,22])
    #plt.legend()
    plt.ylim([-100,100])
    plt.text(17, 80.,"biased", bbox=bbox)
    plt.savefig(path_out+"error_fi_biased-rel")
    plt.close()

    plt.close('all')

    nn=6
    errors_norm=cp.deepcopy(errors_all)
    for ii in range(2):
        if ii==0:
            title='biased'
        else:
            title='unbiased'
        errors_all[5,:,ii]=np.load('add_reg_'+title+'.npy')
        for jj,jl in enumerate(error_names):
            errors_norm[jj,:,ii]=errors_all[jj,:,ii]/np.max(errors_all[jj,:,ii])

    for ii in range(2):
        if ii==0:
            title='biased'
            nn1=1
        else:
            title='unbiased'
            nn1=2
        for c_key,key in enumerate(gem_int_per_station.keys()):
            plt.plot(np.arange(nn1,nn+1),errors_all[nn1-1:nn,c_key,ii],color=p.sim_color[key],label=p.sim_name[key],marker='o')
        plt.ylabel(r'$\epsilon$'+' ('+varunit_dic['WS']+')')
        plt.xticks(np.arange(1,nn+1), error_names[:nn], rotation=45)
        plt.yscale("log")
        plt.ylim([0.00005,2.])
        plt.xlim([0,nn+1])
        plt.text(.25, .8,title, bbox=bbox)
        plt.savefig(path_out+"abs_errors_"+title,dpi=dpi)
        plt.close()
        
        for c_key,key in enumerate(gem_int_per_station.keys()):
            plt.plot(np.arange(nn1,nn+1),errors_norm[nn1-1:nn,c_key,ii],color=p.sim_color[key],label=p.sim_name[key],marker='o')
        plt.ylabel(r'$\hat{\epsilon}$'+' ('+varunit_dic['WS']+')')
        plt.xticks(np.arange(1,nn+1), error_names[:nn], rotation=45)
        plt.ylim([0,1.1])
        plt.xlim([0,nn+1])
        plt.text(.25,1.,title, bbox=bbox)
        plt.savefig(path_out+"norm_errors_"+title,dpi=dpi)
        plt.close()

if __name__ == '__main__':
    main()
