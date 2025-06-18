import config
import utils
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator
import glob
from tqdm import tqdm
from datetime import datetime, timedelta
import xarray as xr
import warnings
import copy as cp
import sys
import pdb
import copy
from matplotlib.lines import Line2D

def main():
    # Setup configuration and simulation parameters
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
    path_out = p.path_out

    # Make output directories if they do not exist
    os.system('mkdir -p ' + path_out)
    os.system('mkdir -p ' + path_out + 'binning_mean_AMF/')
    os.system('mkdir -p ' + path_out + 'binning_mean_perc_error/')
    os.system('mkdir -p ' + path_out + 'binning_mean_perc_reg/')
    os.system('mkdir -p ' + path_out + 'binning_mean_error/')
    os.system('mkdir -p ' + path_out + 'binning_mean_reg-error/')
    os.system('mkdir -p ' + path_out + 'distribution/')
    os.system('mkdir -p ' + path_out + 'error_diu_ann/')
    os.system('mkdir -p ' + path_out + 'ann_error/')
    os.system('mkdir -p ' + path_out + 'diu_error/')
    os.system('mkdir -p ' + path_out + 'distribution_errors/')
    os.system('mkdir -p ' + path_out + 'combined_errors/')
    os.system('mkdir -p ' + path_out + 'combined_ann_diu_errors/')

    vars_ind={}
    for ivar,var in enumerate(variables):
        vars_ind[var]=ivar
    stations=np.load(path_stats+'stats-selected.npy')

    sim_name=p.sim_name
    sources=['AMFc-BEL_withz0','GEM']

    # Set up default plotting parameters
    MPL_DEFAULT_PARAMS = {
        'font.size': 19,
        # 'figure.figsize': (10, 8),
        # 'axes.titlesize': 20,
        # 'axes.labelsize': 18,
        # 'xtick.labelsize': 16,
        # 'ytick.labelsize': 16,
        'legend.fontsize': 15,
        'lines.linewidth': 2.5,
        'figure.dpi': 300,
        'figure.constrained_layout.use': True,
    }
    plt.rcParams.update(MPL_DEFAULT_PARAMS)

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
    flux_names=['LE','H']

    ##########################################
    # Load data for each simulation and regime
    ##########################################
    for var_flux in flux_names:
        print(var_flux," case")
        add=var_flux
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
        for c_var,var in enumerate(variables):
            all_dat_t=np.concatenate((data_gem['AMFc-BEL_withz0-all-'+'-'+z0][:,c_var],data_gem['GEM-all-NAM-11m_GEM511_ERA5_CLASS_USGS_Delage_SN8'+'-'+z0][:,c_var])).ravel()
            minv=np.floor(np.min(all_dat_t))
            bins=[]
            bins.append(minv)
            for c_per,per in enumerate(percentiles[1:]):
                bins.append(np.percentile(all_dat_t,per))
                cond= (all_dat_t>=bins[c_per]) & (all_dat_t<bins[c_per+1])
            bins_all[var]=np.asarray(bins)

        # Bins for the stability: Neutral when lower than 0.02 and extreme when higher than 0.4
        bins_all['ZL'][1]=-.4
        bins_all['ZL'][2]=-.02
        bins_all['ZL'][3]=.02
        bins_all['ZL'][4]=.4

        variables_plot = [var_flux,'ZL','USTAR'] # ['WS','ZL','USTAR','TA','H','LE','PA','RH']
        reg='all'
        add=z0+'_'+var_flux
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

                freq_bin,int_bin,max_bin,max5_bin,std_bin,per_bin=utils.calculate_binning_mean(var_gem,bins_all,variables_plot,percentiles,prob=True)
                freq_bin_a,int_bin_a,max_bin,max5_bin,std_bin,per_bin=utils.calculate_binning_mean(var_amf,bins_all,variables_plot,percentiles,prob=True)

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
                        cmap = copy.copy(plt.get_cmap('viridis'))
                    else:
                        cmap = copy.copy(plt.get_cmap('cubehelix_r'))
                    if term=='tot':
                        cmap = copy.copy(plt.get_cmap('cool'))
                    cmap.set_bad( "gray", alpha=0 )
                    if term=='int':
                        Z=int_bin_a.T
                        levels = MaxNLocator(nbins=11).tick_values(-100,250)
                        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                    if term=='freq':
                        Z=freq_bin_a.T
                        norm=mcolors.LogNorm(vmin=0.0005, vmax=.5, clip=True)
                    if term=='tot':
                        Z=int_bin_a.T*freq_bin_a.T
                        tt=Z.max()
                        bb=Z.min()
                        levels = MaxNLocator(nbins=11).tick_values(bb,tt)
                        norm=mcolors.SymLogNorm(linthresh=1, linscale=1, base=10, clip=True)

                    plt.ylabel(varname_dic[variables_plot[2]] + ' (' + varunit_dic[variables_plot[2]] + ')')
                    plt.xlabel('Stability') # variable_plot[1] is stability
                    labelsx = []
                    aa = bins_all[variables_plot[1]]
                    for ii in np.arange(aa.shape[0]-1) :
                        aa1=str("{:.0e}".format(aa[ii]))
                        aa2=str("{:.0e}".format(aa[ii+1]))
                        labelsx.append('[' + aa1 + ',' + aa2+ ')')
                    labelsx=['VU','U','N','S','VS'] # ['very unstable','unstable','neutral','stable','very stable']
                    labelsy = []
                    aa = bins_all[variables_plot[2]]
                    for ii in np.arange(aa.shape[0]-1):
                        labelsy.append(str(round(aa[ii+1],2)))

                    plt.xticks(np.arange(len(labelsx))+.5, labelsx, rotation=0)
                    plt.yticks(np.arange(len(labelsy))+.5, labelsy, rotation=0)
                    plt.pcolormesh(Z, cmap=cmap, norm=norm)
                    for m in range(5):
                        for h in range(5):
                            val = Z[h, m]
                            if np.ma.is_masked(val):
                                val = "--"
                            else:
                                val = f"{val:.2f}"
                            plt.text(m+0.5, h+0.5, val, ha='center', va='center', color='black',fontsize=11)
                    cbar = plt.colorbar(extend='max')
                    if term == 'int':
                        cbar.set_label(rf'$\overline{{{var_flux}}}$' + ' (' + varunit_dic[variables_plot[0]] + ')')
                    elif term == 'tot' :
                        cbar.set_label(rf'$p \overline{{{var_flux}}}$' + ' (' + varunit_dic[variables_plot[0]] + ')')
                    elif term == 'freq':
                        cbar.set_label(term_names[term] + ' ' )
                    fig_name = path_out+'binning_mean_AMF/'+'binning_mean_AMF_'+term+'_'+add+'.png'
                    plt.savefig(fig_name, bbox_inches='tight') # binning mean AMF
                    plt.close()
                ####################################
                for c_term,term in enumerate(terms):
                    Z = var_perc_errors[c_term,:].T*100
                    Z_v = var_errors[c_term,:].T
                    if term == 'tot' :
                        Zmax = np.ceil(np.abs(Z).max())
                    Zmax=100
                    levels = MaxNLocator(nbins=11, symmetric=True).tick_values(-1*Zmax, Zmax)
                    levels = levels[levels != 0]
                    cmap   = copy.copy(plt.get_cmap('PiYG'))
                    norm   = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                    cmap.set_bad("gray", alpha=0)
                    plt.title(sim_name[sim] + ' - ' + var_flux)
                    plt.ylabel(varname_dic[variables_plot[2]] + ' (' + varunit_dic[variables_plot[2]] + ')')
                    plt.xticks(np.arange(len(labelsx))+.5, labelsx)
                    plt.yticks(np.arange(len(labelsy))+.5, labelsy)
                    for m in range(5):
                        for h in range(5):
                            plt.text(m+0.5,h+0.5,f'{Z_v[h, m]:.3f}', ha='center', va='center', color='black',fontsize=11)
                    plt.pcolormesh(Z, cmap=cmap, norm=norm)
                    cbar = plt.colorbar(extend='both')
                    if term == 'int':
                        cbar.set_label(rf'$p^{{o}} \cdot \Delta \overline{{{var_flux}}}$ (' + varunit_dic[variables_plot[0]] + ') (%)')
                    elif term == 'tot' :
                        cbar.set_label(rf'$\Delta (p \cdot \overline{{{var_flux}}})$ (' + varunit_dic[variables_plot[0]] + ') (%)')
                    elif term == 'freq':
                        cbar.set_label(rf'$\overline{{{var_flux}^{{o}}}} \cdot \Delta p$ (' + varunit_dic[variables_plot[0]] + ') (%)')
                    elif term == 'residual':
                        cbar.set_label(rf'$\Delta \overline{{{var_flux}}} \cdot \Delta p$ (' + varunit_dic[variables_plot[0]] + ') (%)')
                    fig_name = path_out + 'binning_mean_perc_error/' + 'binning_mean_perc_error_' + source + '_' + term + '_' + sim_name[sim] +'_'+add+ '.png'
                    plt.savefig(fig_name, bbox_inches='tight')

                    plt.close()
                    if term == 'tot' :
                        Z = np.sum(abs(var_errors[1:,:,:]),axis=0).T
                        Zmax=100
                        levels = MaxNLocator(nbins=11).tick_values(0, Zmax)
                        cmap   = copy.copy(plt.get_cmap('cool'))
                        norm=mcolors.LogNorm(vmin=0.0005, vmax=Zmax, clip=True)
                        cmap.set_bad( "gray", alpha=0 )
                        plt.title(sim_name[sim] + ' - ' + var_flux)
                        plt.ylabel(varname_dic[variables_plot[2]] + ' (' + varunit_dic[variables_plot[2]] + ')')
                        plt.xticks(np.arange(len(labelsx))+.5, labelsx)
                        plt.yticks(np.arange(len(labelsy))+.5, labelsy)
                        plt.pcolormesh(Z, cmap=cmap, norm=norm)
                        for m in range(5):
                            for h in range(5):
                                plt.text(m+0.5,h+0.5,f'{Z[h, m]:.2f}', ha='center', va='center', color='black',fontsize=11)
                        cbar = plt.colorbar(extend='max')
                        cbar.set_label(r'$\epsilon_{reg}$'+ ' (' + varunit_dic[variables_plot[0]] + ')')
                        fig_name = path_out + 'binning_mean_perc_reg/' + 'binning_mean_perc_reg-error_' + source + '_' + term + '_' + sim_name[sim] +'_'+add+ '.png'
                        plt.savefig(fig_name, bbox_inches='tight')
                        plt.close()


                for c_term,term in enumerate(terms):
                    Z = var_errors[c_term,:].T
                    if term == 'tot' :
                        Zmax = np.ceil(np.abs(Z).max())
                    Zmax=.3
                    levels = MaxNLocator(nbins=11, symmetric=True).tick_values(-Zmax, Zmax)
                    levels = levels[levels != 0]
                    cmap   = copy.copy(plt.get_cmap('PiYG'))
                    norm   = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                    cmap.set_bad( "gray", alpha=0 )
                    plt.title(sim_name[sim])
                    plt.ylabel(varname_dic[variables_plot[2]] + ' (' + varunit_dic[variables_plot[2]] + ')')
                    plt.xticks(np.arange(len(labelsx))+.5, labelsx)
                    plt.yticks(np.arange(len(labelsy))+.5, labelsy)
                    for m in range(5):
                        for h in range(5):
                            plt.text(m+0.5,h+0.5,f'{Z[h, m]:.2f}', ha='center', va='center', color='black',fontsize=11)
                    plt.pcolormesh(Z, cmap=cmap, norm=norm)
                    cbar = plt.colorbar(extend='both')
                    if term == 'int':
                        cbar.set_label(rf'$p^{{o}} \cdot \Delta \overline{{{var_flux}}}$')
                    elif term == 'tot' :
                        cbar.set_label(rf'$\Delta (p \cdot \overline{{{var_flux}}})$')
                    elif term == 'freq':
                        cbar.set_label(rf'$\overline{{{var_flux}^{{o}}}} \cdot \Delta p$' + ' (' + varunit_dic[variables_plot[0]] + ')')
                    elif term == 'residual':
                        cbar.set_label(rf'$\Delta \overline{{{var_flux}}} \cdot \Delta p$' + ' (' + varunit_dic[variables_plot[0]] + ')')
                    fig_name = path_out+'binning_mean_error/'+ 'binning_mean_error_' + source + '_' + term + '_' + sim_name[sim] +'_'+add+ '.png'
                    plt.savefig(fig_name, bbox_inches='tight')

                    plt.close()
                    if term == 'tot' :
                        Z = np.sum(abs(var_errors[1:,:,:]),axis=0).T
                        Zmax=.5
                        levels = MaxNLocator(nbins=11).tick_values(0, Zmax)
                        cmap   = copy.copy(plt.get_cmap('cool'))
                        norm=mcolors.LogNorm(vmin=0.0005, vmax=Zmax, clip=True)
                        cmap.set_bad( "gray", alpha=0 )
                        plt.title(term_error_names[term] + r' ($\epsilon_{reg}=$' + str(np.round(np.ma.mean(abs(Z)),2)) + ' '+varunit_dic[variables_plot[0]] + ')')
                        plt.title(sim_name[sim])
                        plt.ylabel(varname_dic[variables_plot[2]] + ' (' + varunit_dic[variables_plot[2]] + ')')
                        plt.xticks(np.arange(len(labelsx))+.5, labelsx)
                        plt.yticks(np.arange(len(labelsy))+.5, labelsy)
                        plt.pcolormesh(Z, cmap=cmap, norm=norm)
                        for m in range(5):
                            for h in range(5):
                                plt.text(m+0.5,h+0.5,f'{Z[h, m]:.2f}', ha='center', va='center', color='black',fontsize=11)
                        cbar = plt.colorbar(extend='max')
                        cbar.set_label(r'$\epsilon_{reg}$'+ ' (' + varunit_dic[variables_plot[0]] + ')')
                        fig_name = path_out+'binning_mean_reg-error/'+'binning_mean_reg-error_' + source + '_' + term + '_' + sim_name[sim] +'_'+add+ '.png'
                        plt.savefig(fig_name, bbox_inches='tight')
                        plt.close()
                        add_met.append(np.mean(Z))

        np.save('add_reg_'+var_flux+'.npy',add_met)

    ######################################################
    # AMF data processing (Computing the mean, diurnal, and annual cycle errors)
    ######################################################
    # Get amf paths
    amf_data_path = glob.glob(path_in_r+"/AMFc-BEL_withz0_data_*.npy")
    amf_date_path = glob.glob(path_dates+"/AMF_dates_*.npy")

    # Sort paths
    amf_data_path.sort()
    amf_date_path.sort()
    N = len(amf_data_path)

    gem_flux = {} # Store LE and H to compute mean and show distribution
    gem_int_per_station = {} # to store LE_{h,m} and H_{h,m} intensity
    gem_freq_per_station = {} # to store LE_{h,m} and H_{h,m} freq
    gem_diu_per_station = {} # to store LE_{h} and H_{h} diurnal cycle
    gem_ann_per_station = {} # to store LE_{m} and H_{m} annual cycle

    for sim in tqdm(simulation_names, desc='Computing errors', unit='simulation'):
        amf_le = []
        amf_int_per_station_le = []
        amf_freq_per_station_le = []
        amf_diu_per_station_le = []
        amf_ann_per_station_le = []
        amf_hflux = [] # named hflux to avoid confusion with hours
        amf_int_per_station_hflux = []
        amf_freq_per_station_hflux = []
        amf_diu_per_station_hflux = []
        amf_ann_per_station_hflux = []
        counter=0
        for date_p,data_p in zip(amf_date_path, amf_data_path):
            # Load data
            data = np.load(data_p)
            dates = np.load(date_p)

            # Switch paths to gem
            data_p_gem = data_p.replace("AMFc-BEL_withz0_data","GEM_data").replace("1.npy",sim+"_1.npy")
            data_gem = np.load(data_p_gem)
            dates_gem =  dates # Use AMF dates

            # Correct dates
            reference_date=datetime(1971,1,1) # This is the reference date used in the GEM model
            corrected_dates = []
            corrected_dates_gem = []

            for i in range(len(dates)):
                first_date_timestamp = dates[i]
                first_date_timestamp_gem = dates_gem[i]
                timestamp_delta = timedelta(seconds=first_date_timestamp)
                timestamp_delta_gem = timedelta(seconds=first_date_timestamp_gem)
                corrected_dates.append(reference_date + timestamp_delta)
                corrected_dates_gem.append(reference_date + timestamp_delta_gem)

            HFLUX = data[:,4]
            LE = data[:,5]
            HFLUX_gem = data_gem[:,4]
            LE_gem = data_gem[:,5]

            # Save fluxes
            amf_hflux.append(HFLUX)
            amf_le.append(LE)
            if sim in gem_flux:
                gem_flux[sim]['LE'].append(LE_gem)
                gem_flux[sim]['H'].append(HFLUX_gem)
            else:
                gem_flux[sim] = {'LE': [LE_gem], 'H': [HFLUX_gem]}

            # Make dataframes to make quick diurnal and yearly cycles
            ds_le = xr.Dataset({
                'LE': (('time',), LE),
                }, coords={'time': corrected_dates})
            ds_hflux = xr.Dataset({
                'H': (('time',), HFLUX),
                }, coords={'time': corrected_dates})

            ds_gem_le = xr.Dataset({
                    'LE': (('time',), LE_gem),
                    }, coords={'time': corrected_dates_gem})
            ds_gem_hflux= xr.Dataset({
                    'H': (('time',), HFLUX_gem),
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
            var_comp_le = np.array([LE, hour_objects, month_objects])
            var_comp_hflux = np.array([HFLUX, hour_objects, month_objects])
            var_comp_gem_le = np.array([LE_gem, hour_objects_gem, month_objects_gem])
            var_comp_gem_hflux = np.array([HFLUX_gem, hour_objects_gem, month_objects_gem])

            # Define bins for hours and months (assuming bins_all is a dictionary)
            bins_all = {'hours': np.arange(0, 25), 'months': np.arange(1, 14)}

            # Define variables (assuming wind speed is the independent variable)
            variables = ['fluxes', 'hours', 'months'] # The first variable is not used in the function calculate_binning_mean_time

            # Define percentiles. Don't really need this
            percentiles = [25, 50, 75]

            # Call calculate_binning_mean function
            freq_bin_le, int_bin_le, _,_,_,_ = utils.calculate_binning_mean_time(var_comp_le, bins_all, variables, percentiles)
            freq_bin_hflux, int_bin_hflux, _,_,_,_ = utils.calculate_binning_mean_time(var_comp_hflux, bins_all, variables, percentiles)
            freq_bin_gem_le, int_bin_gem_le,_,_,_,_ = utils.calculate_binning_mean_time(var_comp_gem_le, bins_all, variables, percentiles)
            freq_bin_gem_hflux, int_bin_gem_hflux,_,_,_,_ = utils.calculate_binning_mean_time(var_comp_gem_hflux, bins_all, variables, percentiles)

            amf_int_per_station_le.append(int_bin_le)
            amf_int_per_station_hflux.append(int_bin_hflux)
            amf_freq_per_station_le.append(freq_bin_le)
            amf_freq_per_station_hflux.append(freq_bin_hflux)

            if counter==0:
                amf_diu_per_station_le=ds_le.groupby("time.hour").mean()
                amf_ann_per_station_le=ds_le.groupby("time.month").mean()

                amf_diu_per_station_hflux=ds_hflux.groupby("time.hour").mean()
                amf_ann_per_station_hflux=ds_hflux.groupby("time.month").mean()
                counter+=1
            else:
                amf_diu_per_station_le+=ds_le.groupby("time.hour").mean()
                amf_ann_per_station_le+=ds_le.groupby("time.month").mean()

                amf_diu_per_station_hflux+=ds_hflux.groupby("time.hour").mean()
                amf_ann_per_station_hflux+=ds_hflux.groupby("time.month").mean()

            if sim in gem_int_per_station:
                gem_int_per_station[sim]['LE'].append(int_bin_gem_le)
                gem_freq_per_station[sim]['LE'].append(freq_bin_gem_le)
                gem_diu_per_station[sim]['LE'] += ds_gem_le.groupby("time.hour").mean()
                gem_ann_per_station[sim]['LE'] += ds_gem_le.groupby("time.month").mean()

                gem_int_per_station[sim]['H'].append(int_bin_gem_hflux)
                gem_freq_per_station[sim]['H'].append(freq_bin_gem_hflux)
                gem_diu_per_station[sim]['H'] += ds_gem_hflux.groupby("time.hour").mean()
                gem_ann_per_station[sim]['H'] += ds_gem_hflux.groupby("time.month").mean()
            else:
                gem_int_per_station[sim] = {}
                gem_freq_per_station[sim] = {}
                gem_diu_per_station[sim] = {}
                gem_ann_per_station[sim] = {}

                gem_int_per_station[sim]['LE'] = [int_bin_gem_le]
                gem_freq_per_station[sim]['LE'] = [freq_bin_gem_le]
                gem_diu_per_station[sim]['LE'] = ds_gem_le.groupby("time.hour").mean()
                gem_ann_per_station[sim]['LE'] = ds_gem_le.groupby("time.month").mean()

                gem_int_per_station[sim]['H'] = [int_bin_gem_hflux]
                gem_freq_per_station[sim]['H'] = [freq_bin_gem_hflux]
                gem_diu_per_station[sim]['H'] = ds_gem_hflux.groupby("time.hour").mean()
                gem_ann_per_station[sim]['H'] = ds_gem_hflux.groupby("time.month").mean()


    # Compute Errors and make figures
    print("Saving Figures")
    aa_le=np.concatenate([arr.flatten() for arr in amf_le])
    XY_le=float(aa_le.shape[0])
    X_le=float(len(amf_le))
    Y_le=XY_le/X_le
    aa_hflux=np.concatenate([arr.flatten() for arr in amf_hflux])
    XY_hflux=float(aa_hflux.shape[0])
    X_hflux=float(len(amf_hflux))
    Y_hflux=XY_hflux/X_hflux

    month_names=['J','F','M','A','M','J','J','A','S','O','N','D']
    hour_names=np.arange(0,24) # Local hours from 0 to 23

    # -----------------------------------------------
    # Diurnal cycle AMF averaged across all stations
    # -----------------------------------------------
    diu_amf_le = amf_diu_per_station_le['LE']/len(amf_freq_per_station_le)
    diu_amf_hflux = amf_diu_per_station_hflux['H']/len(amf_freq_per_station_hflux)
    def plot_diu_amf_cycle(diu_data, flux_name, filename):
        """Plot the diurnal cycle of AMF data."""
        plt.figure()
        plt.plot(hour_names,diu_data,linestyle='-',color='k',label='Flux')
        plt.xlabel('Local Hour')
        plt.xticks(hour_names[::3]+2, hour_names[::3]+2)
        plt.ylabel(flux_name+' ('+varunit_dic[flux_name]+')')
        plt.ylim(-30, 155)
        plt.savefig(path_out+'combined_ann_diu_errors/'+filename)
        plt.close()
    plot_diu_amf_cycle(diu_amf_le, 'LE', 'amf_diu_LE')
    plot_diu_amf_cycle(diu_amf_hflux, 'H', 'amf_diu_H')

    # -----------------------------------------------
    # Annual cycle AMF averaged across all stations
    # -----------------------------------------------
    ann_amf_le = amf_ann_per_station_le['LE']/len(amf_freq_per_station_le)
    ann_amf_hflux = amf_ann_per_station_hflux['H']/len(amf_freq_per_station_hflux)
    def plot_ann_amf_cycle(ann_data, flux_name, filename):
        """Plot the annual cycle of AMF data."""
        plt.figure()
        plt.plot(np.arange(1,13),ann_data,linestyle='-',color='k',label='Flux')
        plt.xlabel('Month')
        plt.ylabel(flux_name+' ('+varunit_dic[flux_name]+')')
        plt.xticks(np.arange(1,13), month_names)
        plt.ylim(-5, 110)
        plt.savefig(path_out+'combined_ann_diu_errors/'+filename)
        plt.close()
    plot_ann_amf_cycle(ann_amf_le, 'LE', 'amf_ann_LE')
    plot_ann_amf_cycle(ann_amf_hflux, 'H', 'amf_ann_H')

    # -----------------------------------------------
    # Combined monthly and hourly latent fluxes from AMF
    # -----------------------------------------------
    # Intensity and freq for AMF
    int_bin_amf_le = np.mean(amf_int_per_station_le,axis=0)
    freq_bin_amf_le = np.sum(amf_freq_per_station_le,axis=0)/len(amf_freq_per_station_le)

    # Create the heatmap
    plt.figure()
    plt.imshow(int_bin_amf_le, cmap='viridis', aspect='auto')
    levels = MaxNLocator(nbins=11).tick_values(0, int_bin_amf_le.max())
    cmap = plt.get_cmap('viridis')
    norm_cmap = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    im = plt.imshow(int_bin_amf_le, cmap=cmap, norm=norm_cmap, aspect='auto')
    cbar = plt.colorbar(im, boundaries=levels, ticks=levels, spacing='proportional', label=r'$\overline{LE_{h,m}}$'+' ('+varunit_dic['LE']+')')

    # Annotate each cell with the frequency
    for m in range(12):
        for h in range(24):
            value = int_bin_amf_le[h, m]
            plt.text(m, h, f'{int(freq_bin_amf_le[h, m])}', ha='center', va='center', color="white", fontsize=11)

    # Set labels and title
    ticks = np.linspace(0, 250, 5)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{int(t)}" for t in ticks])
    plt.xlabel('Month')
    plt.ylabel('Local Hour')
    plt.xticks(np.arange(0,12), month_names)
    plt.yticks(hour_names[::3]+2, hour_names[::3]+2)
    plt.title('AMF')
    plt.savefig(path_out+'error_diu_ann/'+"AMF_LE_ann-diu")
    plt.close()

    # -----------------------------------------------
    # Combined monthly and hourly sensible fluxes from AMF
    # -----------------------------------------------
    # Intensity and freq for AMF
    int_bin_amf_hflux = np.mean(amf_int_per_station_hflux,axis=0)
    freq_bin_amf_hflux = np.sum(amf_freq_per_station_hflux,axis=0)/len(amf_freq_per_station_hflux)

    # Create the heatmap
    plt.figure()
    levels = np.linspace(-int_bin_amf_hflux.max(), int_bin_amf_hflux.max(), 12)
    cmap = plt.get_cmap('PuOr')
    norm_cmap = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    im = plt.imshow(int_bin_amf_hflux, cmap=cmap, norm=norm_cmap, aspect='auto')
    cbar = plt.colorbar(im, boundaries=levels, ticks=levels, spacing='proportional', label=r'$\overline{H_{h,m}}$'+' ('+varunit_dic['H']+')')

    # Annotate each cell with the frequency
    for m in range(12):
        for h in range(24):
            value = int_bin_hflux[h, m]
            text_color = "black"
            if value > 60:
                text_color = "white"
            plt.text(m, h, f'{int(freq_bin_hflux[h, m])}', ha='center', va='center', color=text_color, fontsize=11)

    # Set labels and title
    cbar.set_ticks(np.linspace(-150, 150, 5))
    plt.xlabel('Month')
    plt.ylabel('Local Hour')
    plt.xticks(np.arange(0,12), month_names)
    plt.yticks(hour_names[::3]+2, hour_names[::3]+2)
    plt.title('AMF')
    plt.savefig(path_out+'error_diu_ann/'+"AMF_H_ann-diu")
    plt.close()

    ###################################################
    # Comparing simulations with AMF
    ###################################################

    bbox = dict(boxstyle='round', fc='grey', ec='white', alpha=0.5)

    bias_error = {}
    diu_error = {}
    ann_error = {}
    ann_error_s = {}
    maxLE=1100
    maxH=950
    nbins=30
    error_names=['bias','diu','ann','diu-ann','fi','reg']
    error_names=[r'$\epsilon_{bias}$',r'$\epsilon_{diu}$',r'$\epsilon_{ann}$',r'$\epsilon_{ann,diu}$',r'$\epsilon_{fi}$',r'$\epsilon_{reg}$']
    errors_all=np.zeros((len(error_names),len(gem_int_per_station.keys()),2)) # Shape is (number of errors, number of simulations, LE and H)
    histo_sims={}
    histo_sims_unbiased={}
    for c_key,key in enumerate(gem_int_per_station.keys()):
        ###############################################
        # Mean and distribution plotting.
        ###############################################
        # need to flatten
        amf_le = np.concatenate([arr.flatten() for arr in amf_le])
        sim_le = np.concatenate([arr.flatten() for arr in gem_flux[key]['LE']])
        bins_le=np.histogram(amf_le,range=(0,maxLE),bins=nbins)[1]
        bins_c_le=0.5*(bins_le[:-1]+bins_le[1:])
        bins2_le=np.histogram(amf_le,range=(-maxLE,maxLE),bins=2*nbins)[1]
        bins_c2_le=0.5*(bins2_le[:-1]+bins2_le[1:])

        amf_hflux = np.concatenate([arr.flatten() for arr in amf_hflux])
        sim_hflux = np.concatenate([arr.flatten() for arr in gem_flux[key]['H']])
        bins_hflux=np.histogram(amf_hflux,range=(0,maxH),bins=nbins)[1]
        bins_c_hflux=0.5*(bins_hflux[:-1]+bins_hflux[1:])
        bins2_hflux=np.histogram(amf_hflux,range=(-maxH,maxH),bins=2*nbins)[1]
        bins_c2_hflux=0.5*(bins2_hflux[:-1]+bins2_hflux[1:])

        # biased
        if key not in histo_sims:
            histo_sims[key] = {}
        histo_amf_le=np.histogram(amf_le,range=(0,maxLE),bins=nbins)[0]
        histo_sims[key]['LE']=np.histogram(sim_le,range=(0,maxLE),bins=nbins)[0]
        error_le= bins_c_le*(histo_sims[key]['LE']-histo_amf_le)
        norm_le=1./XY_le
        errors_all[4,c_key,0]=norm_le*np.sum(np.abs(error_le)) # 0 is for LE

        histo_amf_hflux=np.histogram(amf_hflux,range=(0,maxH),bins=nbins)[0]
        histo_sims[key]['H']=np.histogram(sim_hflux,range=(0,maxH),bins=nbins)[0]
        error_hflux= bins_c_hflux*(histo_sims[key]['H']-histo_amf_hflux)
        norm_hflux=1./XY_hflux
        errors_all[4,c_key,1]=norm_hflux*np.sum(np.abs(error_hflux)) # 1 is for H

        plt.figure()
        plt.step(bins_le[:-1],histo_amf_le,where='post',label="AMF LE",color='#1f77b4')
        plt.step(bins_le[:-1],histo_sims[key]['LE'],where='post',label=f"{p.sim_name[key]} LE",color='#ff7f0e')
        plt.axvline(x=np.mean(amf_le),color='#1f77b4',ls='--')
        plt.axvline(x=np.mean(sim_le),color='#ff7f0e',ls='--')
        plt.step(bins_hflux[:-1],histo_amf_hflux,where='post',label="AMF H",color='#2ca02c')
        plt.step(bins_hflux[:-1],histo_sims[key]['H'],where='post',label=f"{p.sim_name[key]} H",color='#d62728')
        plt.axvline(x=np.mean(amf_hflux),color='#2ca02c',ls='--')
        plt.axvline(x=np.mean(sim_hflux),color='#d62728',ls='--')
        plt.legend(loc='best')
        plt.xscale("log")
        plt.yscale("log")
        plt.savefig(path_out+"distribution/"+p.sim_name[key]+"_dist") # TODO: Change colours
        plt.close()

        bias_error[key] = {}
        bias_error[key]['LE']=np.abs(np.mean(sim_le) - np.mean(amf_le))
        bias_error[key]['H']=np.abs(np.mean(sim_hflux) - np.mean(amf_hflux))

        ###############################################
        # diu stuff here
        ###############################################
        gem_diu_le = gem_diu_per_station[key]['LE']/X_le
        gem_diu_hflux = gem_diu_per_station[key]['H']/X_hflux
        hours = np.arange(1,25)
        plt.figure()
        plt.plot(hours,diu_amf_le,label="AMF LE",color='b')
        plt.plot(hours,gem_diu_le['LE'],label=f"{p.sim_name[key]} LE",color='b',ls='--')
        plt.plot(hours,diu_amf_hflux,label="AMF H",color='r')
        plt.plot(hours,gem_diu_hflux['H'],label=f"{p.sim_name[key]} H",color='r',ls='--') # TODO: Verify why ['H'] is needed
        plt.xlabel("Local Hour")
        plt.xticks(hour_names[::3]+2, hour_names[::3]+2)
        plt.legend()
        plt.savefig(path_out+"diu_error/"+p.sim_name[key]+"_diu_error")
        plt.close()

        diu_error[key] = {}
        diu_error[key]['LE'] = gem_diu_le['LE'] - diu_amf_le
        diu_error[key]['H'] = gem_diu_hflux['H'] - diu_amf_hflux

        ###############################################
        # ann stuff here
        ###############################################
        L_le = np.asarray(gem_int_per_station[key]['LE']).shape[0]
        gem_ann_le = gem_ann_per_station[key]['LE']/L_le
        L_hflux = np.asarray(gem_int_per_station[key]['H']).shape[0]
        gem_ann_hflux = gem_ann_per_station[key]['H']/L_hflux
        months = np.arange(1,13)
        plt.figure()
        plt.plot(months,ann_amf_le,label="AMF LE",color='b')
        plt.plot(months,gem_ann_le['LE'],label=f"{p.sim_name[key]} LE",color='b',ls='--')
        plt.plot(months,ann_amf_hflux,label="AMF H",color='r')
        plt.plot(months,gem_ann_hflux['H'],label=f"{p.sim_name[key]} H",color='r',ls='--')
        plt.xlabel('Month')
        plt.ylabel('Fluxes'+' ('+varunit_dic['LE']+')')
        plt.xticks(np.arange(1,13), month_names)
        plt.legend()
        plt.savefig(path_out+"ann_error/"+p.sim_name[key]+"_ann_error")
        plt.close()

        ann_error[key] = {}
        ann_error[key]['LE'] = gem_ann_le['LE']-ann_amf_le
        ann_error[key]['H'] = gem_ann_hflux['H']-ann_amf_hflux

        ###############################################
        # ann_diu stuff here
        ###############################################
        L_le = np.asarray(gem_int_per_station[key]['LE']).shape[0]
        int_bin_gem_le = np.sum(np.asarray(gem_int_per_station[key]['LE']),axis=0)/L_le
        freq_bin_gem_le = np.sum(np.asarray(gem_freq_per_station[key]['LE']),axis=0)
        L_hflux = np.asarray(gem_int_per_station[key]['H']).shape[0]
        int_bin_gem_hflux = np.sum(np.asarray(gem_int_per_station[key]['H']),axis=0)/L_hflux
        freq_bin_gem_hflux = np.sum(np.asarray(gem_freq_per_station[key]['H']),axis=0)

        # Compute absolute error
        error_le = freq_bin_amf_le*(int_bin_gem_le-int_bin_amf_le)
        abs_error_le = freq_bin_amf_le*np.abs(int_bin_gem_le-int_bin_amf_le)
        norm_le=1./XY_le
        errors_all[0,c_key,0]=norm_le*np.abs(np.sum(error_le))
        errors_all[1,c_key,0]=norm_le*np.sum(np.abs(np.sum(error_le,axis=0)))
        errors_all[2,c_key,0]=norm_le*np.sum(np.abs(np.sum(error_le,axis=1)))
        errors_all[3,c_key,0]=norm_le*np.sum(np.abs(error_le))

        error_hflux = freq_bin_amf_hflux*(int_bin_gem_hflux-int_bin_amf_hflux)
        abs_error_hflux = freq_bin_amf_hflux*np.abs(int_bin_gem_hflux-int_bin_amf_hflux)
        norm_hflux=1./XY_hflux
        errors_all[0,c_key,1]=norm_hflux*np.abs(np.sum(error_hflux))
        errors_all[1,c_key,1]=norm_hflux*np.sum(np.abs(np.sum(error_hflux,axis=0)))
        errors_all[2,c_key,1]=norm_hflux*np.sum(np.abs(np.sum(error_hflux,axis=1)))
        errors_all[3,c_key,1]=norm_hflux*np.sum(np.abs(error_hflux))

        ###############################################
        # Heatmaps of errors for LE and H
        ###############################################
        # Create the heatmap for LE
        plt.figure()
        LIMIT = 170
        levels = np.linspace(-LIMIT, LIMIT, 12)
        cmap = copy.copy(plt.get_cmap('seismic'))
        cmap.set_bad("gray", alpha=0)  # Set bad values to gray
        norm_cmap = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        im = plt.imshow(error_le/freq_bin_amf_le, cmap=cmap, norm=norm_cmap, aspect='auto')
        cbar = plt.colorbar(im, boundaries=levels, ticks=levels, spacing='proportional', label=r'$\epsilon_{ann,diu}$'+' ('+varunit_dic['LE']+')')
        cbar.set_ticks(np.linspace(-150, 150, 5))
        plt.xlabel('Month')
        plt.ylabel('Local Hour')
        plt.xticks(np.arange(0,12), month_names)
        plt.yticks(hour_names[::3]+2, hour_names[::3]+2)
        plt.title(p.sim_name[key] + ' - LE')
        plt.savefig(path_out+"error_diu_ann/"+p.sim_name[key]+"_error_diu-ann_LE")
        plt.close()

        # Create the heatmap for H
        plt.figure()
        norm_cmap = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        im = plt.imshow(error_hflux/freq_bin_amf_hflux, cmap=cmap, norm=norm_cmap, aspect='auto')
        cbar = plt.colorbar(im, boundaries=levels, ticks=levels, spacing='proportional', label=r'$\epsilon_{ann,diu}$'+' ('+varunit_dic['H']+')')
        cbar.set_ticks(np.linspace(-150, 150, 5))
        plt.xlabel('Month')
        plt.ylabel('Local Hour')
        plt.xticks(np.arange(0,12), month_names)
        plt.yticks(hour_names[::3]+2, hour_names[::3]+2)
        plt.title(p.sim_name[key] + ' - H')
        plt.savefig(path_out+"error_diu_ann/"+p.sim_name[key]+"_error_diu-ann_H")
        plt.close()

    # Plot diu and ann error cycles
    def plot_diu_error_cycle(diu_data, flux_name):
        """Plot the diurnal error cycle of AMF data."""
        plt.figure()
        plt.plot([-1,24],[0,0],color='grey',linewidth=.75)
        for key in diu_error:
            plt.plot(hour_names,diu_error[key][flux_name], label=f"{p.sim_name[key]} FLUX", color=p.sim_color[key], linewidth=1.5)
        plt.xlabel('Local Hour')
        plt.xticks(hour_names[::3]+2, hour_names[::3]+2)
        plt.ylabel(r"$\epsilon_{diu}$"+' ('+varunit_dic['LE']+')')
        plt.ylim([-100,100])
        plt.xlim([-0.5, 23.5])
        plt.text(0.5, 82, flux_name, bbox=bbox)
        plt.savefig(path_out + 'combined_ann_diu_errors/' + "diu_error_" + flux_name)
        plt.close()
    plot_diu_error_cycle(diu_error, 'LE')
    plot_diu_error_cycle(diu_error, 'H')

    def plot_ann_error_cycle(ann_data, flux_name):
        """Plot the annual error cycle of AMF data."""
        plt.figure()
        plt.plot([0,13],[0,0],color='grey',linewidth=.75)
        for key in ann_data:
            plt.plot(np.arange(1,13),ann_data[key][flux_name], label=f"{p.sim_name[key]} FLUX", color=p.sim_color[key], linewidth=1.5)
        plt.xlabel('Month')
        plt.ylabel(r"$\epsilon_{ann}$"+' ('+varunit_dic['LE']+')')
        plt.xticks(np.arange(1,13), month_names)
        plt.ylim([-64, 64])
        plt.xlim([0.5, 12.5])
        plt.text(1, 52, flux_name, bbox=bbox)
        plt.savefig(path_out + 'combined_ann_diu_errors/' + "ann_error_" + flux_name)
        plt.close()
    plot_ann_error_cycle(ann_error, 'LE')
    plot_ann_error_cycle(ann_error, 'H')

    # Create a legend for combined_ann_diu_errors
    colors = [p.sim_color[key] for key in gem_int_per_station.keys()]
    labels = [p.sim_name[key] for key in gem_int_per_station.keys()]
    handles = [Line2D([0], [0], color='k', linestyle='-', label='FLUXES')]
    handles += [Line2D([0], [0], color=color, linewidth=1.5, label=label) for color, label in zip(colors, labels)]
    labels = ['AMF'] + labels
    fig_legend = plt.figure(figsize=(3, 2))
    fig_legend.legend(handles, labels, loc='center', ncol=2)
    fig_legend.savefig(path_out + 'combined_ann_diu_errors/legend_only.png', bbox_inches='tight')
    plt.close(fig_legend)

    # Plot dist in on figure
    plt.figure()
    amf_le = np.concatenate([arr.flatten() for arr in amf_le])
    plt.step(bins_le[:-1],histo_amf_le,where='post',label="AMF",color='black')
    plt.xlabel(r'$\overline{LE^{o}_k}$'+' ('+varunit_dic['LE']+')')
    plt.ylabel(r'$N_{k}$')
    plt.yscale("log")
    # plt.xlim([0,22])
    plt.ylim()
    plt.savefig(path_out+'distribution_errors/'+"AMF_fi_LE")
    plt.close()

    # Plot dist in on figure
    plt.figure()
    amf_hflux = np.concatenate([arr.flatten() for arr in amf_hflux])
    plt.step(bins_hflux[:-1],histo_amf_hflux,where='post',label="AMF",color='black') # Check for negative values
    plt.xlabel(r'$\overline{H^{o}_k}$'+' ('+varunit_dic['H']+')')
    plt.ylabel(r'$N_{k}$')
    plt.yscale("log")
    # plt.xlim([0,22])
    plt.savefig(path_out+'distribution_errors/'+"AMF_fi_H")
    plt.close()

    # Plot bias dist in on figure
    plt.figure()
    maxA_le=0
    amf_le = np.concatenate([arr.flatten() for arr in amf_le])
    for key in gem_int_per_station.keys():
        plt.step(bins_le[:-1],bins_c_le*(histo_sims[key]['LE']-histo_amf_le),where='post',label=p.sim_name[key], color=p.sim_color[key])     
        maxA_le=np.max((maxA_le,np.max(np.abs(bins_c_le*(histo_sims[key]['LE']-histo_amf_le)))))
    plt.xlabel(r'$\overline{LE^{o}_k}$'+' ('+varunit_dic['LE']+')')
    plt.ylabel('Absolute error ('+varunit_dic['LE']+')') # r'$(N^{s}_{k}-N^{o}_{k}) \cdot LE^{o}_{k}$'+' ('+varunit_dic['LE']+')'
    # plt.xlim([0,22])
    plt.ylim([-4000000,4000000])
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # plt.text(130, 270000.,"LE", bbox=bbox)
    # plt.text(.99, .99, 'biased', ha='right', va='bottom', transform=ax.transAxes, bbox=bbox)
    plt.savefig(path_out+'distribution_errors/'+"error_fi_LE")
    plt.close()
    plt.close('all')

    plt.figure()
    maxA_hflux=0
    amf_hflux = np.concatenate([arr.flatten() for arr in amf_hflux])
    for key in gem_int_per_station.keys():
        plt.step(bins_hflux[:-1],bins_c_hflux*(histo_sims[key]['H']-histo_amf_hflux),where='post',label=p.sim_name[key], color=p.sim_color[key])
        maxA_hflux=np.max((maxA_hflux,np.max(np.abs(bins_c_hflux*(histo_sims[key]['H']-histo_amf_hflux)))))
    plt.xlabel(r'$\overline{H^{o}_k}$'+' ('+varunit_dic['H']+')')
    plt.ylabel('Absolute error ('+varunit_dic['H']+')') # r'$(N^{s}_{k}-N^{o}_{k}) \cdot H^{o}_{k}$'+' ('+varunit_dic['H']+')'
    # plt.xlim([0,22])
    plt.ylim([-4000000,4000000])
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # plt.text(140, 120000.,"H", bbox=bbox)
    # plt.text(.99, .99, 'biased', ha='right', va='bottom', transform=ax.transAxes, bbox=bbox)
    plt.savefig(path_out+'distribution_errors/'+"error_fi_H")
    plt.close()
    plt.close('all')

    plt.figure()
    for key in gem_int_per_station.keys():
        denom_le = histo_amf_le + histo_sims[key]['LE']
        with np.errstate(divide='ignore', invalid='ignore'):
            err_le = (histo_sims[key]['LE'] - histo_amf_le) / denom_le
            err_le = np.where(denom_le == 0, np.nan, err_le)
        plt.step(bins_le[:-1],100.*err_le,where='post',label=p.sim_name[key], color=p.sim_color[key])
    plt.xlabel(r'$\overline{LE^{o}_k}$'+' ('+varunit_dic['LE']+')')
    plt.ylabel('Relative error (%)')
    plt.ylim([-100,100])
    plt.savefig(path_out+'distribution_errors/'+"error_fi_LE-rel")
    plt.close()
    plt.close('all')

    plt.figure()
    for key in gem_int_per_station.keys():
        denom_hflux = histo_amf_hflux + histo_sims[key]['H']
        with np.errstate(divide='ignore', invalid='ignore'):
            err_hflux = (histo_sims[key]['H'] - histo_amf_hflux) / denom_hflux
            err_hflux = np.where(denom_hflux == 0, np.nan, err_hflux)
        plt.step(bins_hflux[:-1],100.*err_hflux,where='post',label=p.sim_name[key], color=p.sim_color[key])
    plt.xlabel(r'$\overline{H^{o}_k}$'+' ('+varunit_dic['H']+')')
    plt.ylabel('Relative error (%)')
    plt.ylim([-100,100])
    plt.savefig(path_out+'distribution_errors/'+"error_fi_H-rel")
    plt.close()
    plt.close('all')

    # Create a legend for the distribution errors
    colors = [p.sim_color[key] for key in gem_int_per_station.keys()]
    labels = [p.sim_name[key] for key in gem_int_per_station.keys()]
    handles = [Line2D([0], [0], color='k', linestyle='-', label='AMF')] + [Line2D([0], [0], color=color, linestyle='-', label=label) for color, label in zip(colors, labels)]
    labels = ['AMF'] + labels
    fig_legend = plt.figure(figsize=(3, 2))
    fig_legend.legend(handles, labels, loc='center', ncol=2)
    fig_legend.savefig(path_out + 'distribution_errors/legend_only.png', bbox_inches='tight')
    plt.close(fig_legend)

    nn=6
    errors_norm=cp.deepcopy(errors_all)
    for ii in ['LE','H']:
        if ii=='LE':
            errors_all[5,:,0]=np.load('add_reg_LE.npy')
            for jj,jl in enumerate(error_names):
                errors_norm[jj,:,0]=errors_all[jj,:,0]/np.max(errors_all[jj,:,0])
        elif ii=='H':
            errors_all[5,:,1]=np.load('add_reg_H.npy')
            for jj,jl in enumerate(error_names):
                errors_norm[jj,:,1]=errors_all[jj,:,1]/np.max(errors_all[jj,:,1])

    for ii in range(2):
        nn1=1
        if ii==0:
            title='LE'
            for c_key,key in enumerate(gem_int_per_station.keys()):
                plt.plot(np.arange(nn1,nn+1),errors_all[nn1-1:nn,c_key,0],color=p.sim_color[key],label=p.sim_name[key],marker='o')
        else:
            title='H'
            for c_key,key in enumerate(gem_int_per_station.keys()):
                plt.plot(np.arange(nn1,nn+1),errors_all[nn1-1:nn,c_key,1],color=p.sim_color[key],label=p.sim_name[key],marker='o')
        plt.ylabel(r'$\epsilon$'+' ('+varunit_dic['LE']+')')
        plt.xticks(np.arange(1,nn+1), error_names[:nn], rotation=45)
        plt.yscale("log")
        plt.ylim([0.2,40])
        yticks = [0.2, 1, 10]
        yticklabels = [r'$2 \times 10^{-1}$', r'$10^0$', r'$10^1$']
        plt.yticks(yticks, yticklabels)
        plt.xlim([0,nn+1])
        plt.text(.25, 22, title, bbox=bbox)
        plt.savefig(path_out+'combined_errors/'+"abs_errors_"+title)
        plt.close()

        if ii==0:
            for c_key,key in enumerate(gem_int_per_station.keys()):
                plt.plot(np.arange(nn1,nn+1),errors_norm[nn1-1:nn,c_key,0],color=p.sim_color[key],label=p.sim_name[key],marker='o')
        else:
            for c_key,key in enumerate(gem_int_per_station.keys()):
                plt.plot(np.arange(nn1,nn+1),errors_norm[nn1-1:nn,c_key,1],color=p.sim_color[key],label=p.sim_name[key],marker='o')
        plt.ylabel(r'$\hat{\epsilon}$'+' ('+varunit_dic['LE']+')')
        plt.xticks(np.arange(1,nn+1), error_names[:nn], rotation=45)
        plt.ylim([0,1.1])
        plt.xlim([0,nn+1])
        plt.text(.25, 1, title, bbox=bbox)
        plt.savefig(path_out+'combined_errors/'+"norm_errors_"+title)
        plt.close()

    # Create a legend for the combined errors
    colors = [p.sim_color[key] for key in gem_int_per_station.keys()]
    labels = [p.sim_name[key] for key in gem_int_per_station.keys()]
    handles = [Line2D([0], [0], color=color, marker='o', linestyle='-', label=label) for color, label in zip(colors, labels)]
    fig_legend = plt.figure(figsize=(3, 2))
    fig_legend.legend(handles, labels, loc='center', ncol=2)
    fig_legend.savefig(path_out + 'combined_errors/legend_only.png', bbox_inches='tight')
    plt.close(fig_legend)

if __name__ == '__main__':
    main()
