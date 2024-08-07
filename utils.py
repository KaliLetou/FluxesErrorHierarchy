import numpy as np

class variables:
  """
  This class has information about variables including long names and units
  """
  names={'WS':'Wind Speed','LE':'Latent heat flux','H':'Sensible heat flux','SWC':'Soil water content','P':'Precipitation','USTAR':'Friction velocity','WD':'Wind direction',
         'SW_IN':'Shortwave radiation, incoming',         'SW_OUT':'Shortwave radiation, outgoing',         'LW_IN':'Longwave radiation, incoming',         'LW_OUT':'Longwave radiation, outgoing','G':'Soil heat flux', 'ZL':'Stability parameter','PA':'Atmospheric pressure','RH':'Relative humidity','TS':'Skin temperature','TA':'2-m temperature','VPD':'Vapor Pressure Deficit','Rn':'Net radiation','SWn':'Net shortwave radiation','LWn':'Net longwave radiation'}
  units={'bilan':r'W m$^{-2}$','WS':'m s$^{-1}$', 'LE':r'W m$^{-2}$','H':r'W m$^{-2}$','SWC':'%','P':'mm','USTAR':'m s$^{-1}$','WD':'degrees','SW_IN': r'W m$^{-2}$','SW_OUT':r'W m$^{-2}$','LW_IN':r'W m$^{-2}$','LW_OUT':r'W m$^{-2}$','G':r'W m$^{-2}$','ZL':'-','PA':'kPa','RH':'%','TA':r'$K$','TS':r'$^{\circ}$C','VPD':'hPa','Rn':r'W m$^{-2}$','SWn':r'W m$^{-2}$','LWn':r'W m$^{-2}$'}
  symbols={'WS':r'$u$','USTAR':r'$u_{*}$','ZL':r'$zL^{-1}$','TA':r'$T_a$','TA':r'$T_a$','PA':r'$p$','RH':r'$RH$','LE':r'$LE$','H':r'$H$','P':r'$P$'}
  min_value={'WS':0, 'LE':-1000.,'H':-1000.,'SWC':0,'P':0,'USTAR':0.0001,'WD':0,'SW_IN':-10.,'SW_OUT':-10.,'LW_IN':-10,'LW_OUT':-10,'G':-400.,'ZL':-5000000.,'PA':70.,'RH':0,'TA':-50.,'TS':-50.,'VPD':0}
  max_value={'WS':50,'LE':1000.,'H':1000.,'SWC':100,'P':200,'USTAR':20.,'WD':3600,'SW_IN':1500,'SW_OUT':1500,'LW_IN':1500,'LW_OUT':900,'G':400,'ZL':5000000.,'PA':110.,'RH':100.1,'TA':50.,'TS':50.,'VPD':100}
  hm_limits={'WS':[0,8],'USTAR':[0,.6],'ZL':[-.75,.75],'TA':[0,25],'PA':[99,102],'RH':[0,100],'LE':[0,125],'H':[-50,250],'P':[0,.2]}
  mm_limits={'WS':[0,5],'USTAR':[0,.5],'ZL':[-.1,.5],'TA':[0,25],'PA':[99,102],'RH':[0,100],'LE':[0,125],'H':[-50,250],'P':[0,.2]}
  pdf_limits={'WS':[0,20],'USTAR':[0,1.6],'ZL':[-.75,.75],'TA':[-25,45],'PA':[85,110],'RH':[0,100],'LE':[-100,600],'H':[-200,600],'P':[0,20]}


class constants:
  """ ADD HERE ONLY UNIVERSAL CONSTANTS                                                                                                                      
  """
  import datetime as dt
  celcius_to_kelvin_conversion_constant = 273.15
  stefan_boltzmann_constant = 5.67 * 10**(-8)     # units : W * m**(-2) * K**(-4)
  reference_date=dt.datetime(1971,1,1)
  cutoff_date=dt.datetime(2021,12,31)             # last date for the GEM simulations
  min_samples=365.25*24
  missing_value=1.e+20
  dpi=200
  extension='.png'
  min_data_percentage=50.
  min_nbr_of_yrs     =1.
  central_time_output=1.5
  budget_close_value=30.
  bins_vars = {'WS':np.arange(0,50.,5.),'LE':np.arange(-500,1100.,100.),'H':np.arange(-500.,1100.,100.),'SWC':np.arange(0,1.0,.04),'P':np.arange(0,1.0,.04),'USTAR':np.arange(0,5.0,.5),'WD':np.arange(0,1.0,.04),
               'SW_IN':np.arange(0,1200.0,100.),'SW_OUT':np.arange(0,400.0,40.),'LW_IN':np.arange(100,600,50.), 'LW_OUT':np.arange(100,600.,50.),'G':np.arange(-100,100.,20.), 'ZL':np.arange(0,1.0,.04),'PA':np.arange(0,1.0,.04),'RH':np.arange(0,1.0,.04),'TS':np.arange(-15.,45.,5.),'TA':np.arange(-15.,45.,5.),'VPD':np.arange(0,1.0,.04)}
  data_names={'AMF':'AmeriFlux','GEM':'CRCM6/GEM5'}
  term_names={'int':r'$\overline{u}$','freq':r'$p$','tot':r'$p \overline{u}$'}
  term_error_names={'int':r'$N^{o} \cdot \Delta \overline{u}$','freq':r'$\overline{u^{o}} \cdot \Delta N$','tot':r'$\Delta (N \cdot \overline{u})$','residual':r'$\Delta \overline{u} \cdot \Delta N$'}
  sea_months={'ANN':[1,2,3,4,5,6,7,8,9,10,11,12],'DJF':[12,1,2],'MAM':[3,4,5],'JJA':[6,7,8],'SON':[9,10,11]}

def get_mask(array,mask=None):
    missing_m = np.ma.masked_array(array,array==constants.missing_value).mask
    nan_m= np.ma.masked_array(array, np.isnan(array)).mask
    if mask==None:
        array_out = np.ma.masked_array(array,np.logical_or(nan_m,missing_m))
    else:
        array_out = np.ma.masked_array(array,np.logical_or(mask,nan_m,missing_m))
    return array_out

def calculate_binning_mean(var_comp,bins_all,variables,percentiles,prob=None):
    """                                                                                                                                                                                                              
    This function separates a given variable (var_comp[0,:]) into different regimes that are defined by binning two other variables: var_comp[1,:] and var_comp[2,:].
    
    INPUT: 
  	1) var_comp: is an array with two dimensions. The first dimension has a shape 3 corresponding to the number of variables. The second dimension has values of each variables and can contain different times/grid points.
	2) bins_all: is a dictionary with arrays describing the binning of each variable. It is calculated using the calculate_bins function.
	3) variables: include the names of the variables to use in the binning process. Variables[0] is the independent variable and the other two are the dependent variables.
	4) percentiles: for each regime (i,j) we calculate various statistics include the mean, std, max and variaous percentiles.

    OUTPUT:
	- 2-dimensional array with the number (freq_bin), mean intensity (int_bin),maximum (max_bin), etc. max5_bin,std_bin,per_bin

    Created: 02/08/2021                                                                                                                                                                                                                       
    Last Modification:                                                                                                                                                                                                                    
    - 02/08/2021:
                                                                                                                                                                                                                              
    """
    min_n=5 # minimum number of values in a bin so the statistics are calculated.

    #print('\n','\n','+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    #print('  --> CALCULATING BINNING ')
    #print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    #print('Independent variable: ', variables[0])
    var_comp=get_mask(var_comp)
    prm=var_comp[0,:][var_comp[0,:].mask==False].flatten()
    #print('Large-scale variables:', variables[1:])
    #for ivar,var in enumerate(variables):
    #  print('  -> ',var, '(min, max): ', np.min(var_comp[ivar,:]),np.max(var_comp[ivar,:]))

    #print('\n','  --> Digitizing large-scale variables')
    ind1=np.digitize(var_comp[1,:][var_comp[1,:].mask==False].flatten(), bins_all[variables[1]])
    ind2=np.digitize(var_comp[2,:][var_comp[2,:].mask==False].flatten(), bins_all[variables[2]])
    #pdb.set_trace()
    # Define arrays for the output
    freq_bin=np.zeros((bins_all[variables[1]].shape[0]-1,bins_all[variables[2]].shape[0]-1))
    int_bin=np.zeros((bins_all[variables[1]].shape[0]-1,bins_all[variables[2]].shape[0]-1))
    max_bin=np.zeros((bins_all[variables[1]].shape[0]-1,bins_all[variables[2]].shape[0]-1))
    max5_bin=np.zeros((bins_all[variables[1]].shape[0]-1,bins_all[variables[2]].shape[0]-1))
    std_bin=np.zeros((bins_all[variables[1]].shape[0]-1,bins_all[variables[2]].shape[0]-1))
    per_bin=np.zeros((len(percentiles),bins_all[variables[1]].shape[0]-1,bins_all[variables[2]].shape[0]-1))

    #print('\n','  --> Calculating intensity, frequency and other statistiques at each regime')
    for cb1,bin1 in enumerate(bins_all[variables[1]][:-1]):
        for cb2,bin2 in enumerate(bins_all[variables[2]][:-1]):
          #print(cb1,cb2)
          freq_bin[cb1,cb2]=np.sum((ind1==cb1+1) & (ind2==cb2+1))
          #pdb.set_trace()
          if freq_bin[cb1,cb2]<=min_n:
            int_bin[cb1,cb2]=constants.missing_value
            max5_bin[cb1,cb2]=constants.missing_value
            std_bin[cb1,cb2]=constants.missing_value
            for c_per,per in enumerate(percentiles):
              per_bin[c_per,cb1,cb2]=constants.missing_value
          else:
            int_bin[cb1,cb2]=np.ma.mean(prm[(ind1==cb1+1) & (ind2==cb2+1)])
            max5_bin[cb1,cb2]=np.ma.mean(np.sort(prm[(ind1==cb1+1) & (ind2==cb2+1)])[-5:])
            std_bin[cb1,cb2]=np.ma.std(prm[(ind1==cb1+1) & (ind2==cb2+1)])
            for c_per,per in enumerate(percentiles):
              per_bin[c_per,cb1,cb2]=np.percentile(prm[(ind1==cb1+1) & (ind2==cb2+1)],float(per))

    int_bin=get_mask(int_bin)
    std_bin=get_mask(std_bin)
    if prob:
        freq_bin=freq_bin/(np.sum(freq_bin))
    return freq_bin,int_bin,max_bin,max5_bin,std_bin,per_bin

def error_decomposition_taylor(var_tot,terms,Im,Nm,Io,No):
    """
    This functions calculate errors for the difference between two functions of the form f=I*N is used.

    Created: 01/08/2020
    Last Modification:
    - 01/08/2020:
       I
    """
    factor=1.
    for term in terms:
      if term=='tot':
        var_tot[0,:]=factor*(Im*Nm-Io*No)
      if term=='int':
        var_tot[2,:]=factor*(Im-Io)*No
      if term=='freq':
        var_tot[1,:]=factor*(Nm-No)*Io
      if term=='residual':
        var_tot[3,:]=(Im-Io)*(Nm-No)
    #pdb.set_trace()
    return var_tot
    
def percent_error_decomposition_taylor(var_tot,terms,Im,Nm,Io,No):
    """
    This functions calculate errors for the difference between two functions of the form f=I*N is used.

    Created: 01/08/2020
    Last Modification:
    - 01/08/2020:
       I
    """
    factor=1.
    for term in terms:
      if term=='tot':
        var_tot[0,:]=factor*(Im*Nm-Io*No) / Io*No 
      if term=='int':
        var_tot[2,:]=factor*(Im-Io)*No / (Io*No)  
      if term=='freq':
        var_tot[1,:]=factor*(Nm-No)*Io / (Io*No)
      if term=='residual':
        var_tot[3,:]=(Im-Io)*(Nm-No)
    #pdb.set_trace()
    return var_tot

#############
def error_illdefined(var_tot,terms,Im,Nm,Io,No):
    var_tot=get_mask(var_tot)
    for bb in [True,False]:
      if bb==False:
        bb1=True
      else:
        bb1=False
      condition=((Io*No).mask==bb) & ((Im*Nm).mask==bb1)
      #print(bb,np.sum(condition))
      if np.sum(condition)>0:
        if bb==False:
          for n_stat,stat in enumerate(['tot','freq']):
            ind=np.where(np.asarray(terms)==stat)[0][0]
            var_tot[ind,:][condition]=-(Io*No)[condition].data
            var_tot.mask[ind,:][condition]=False
          for n_stat,stat in enumerate(['int','residual']):
            ind=np.where(np.asarray(terms)==stat)[0][0]
            var_tot[ind,:][condition]=-(Io*No)[condition].data*0
            var_tot[ind,:][condition].mask=False
        else:
          for stat in ['tot','freq']:
            ind=np.where(np.asarray(terms)==stat)[0][0]
            #pdb.set_trace()
            var_tot[ind,:][condition]=(Im*Nm)[condition].data
            var_tot[ind,:][condition].mask=False
          for stat in ['int','residual']:
            ind=np.where(np.asarray(terms)==stat)[0][0]
            var_tot[ind,:][condition]=(Im*Nm)[condition].data*0
            var_tot[ind,:][condition].mask=False
    return var_tot

# Modified binnin function
def calculate_binning_mean_time(var_comp, bins_all, variables, percentiles, prob=None):
    """                                                                                                                                                                                                              
    This function separates a given variable (var_comp[0,:]) into different regimes that are defined by binning two other variables: var_comp[1,:] and var_comp[2,:].
    XY
    INPUT: 
        1) var_comp: is an array with two dimensions. The first dimension has a shape 3 corresponding to the number of variables. The second dimension has values of each variables and can contain different times/grid points.
        2) bins_all: is a dictionary with arrays describing the binning of each variable. It is calculated using the calculate_bins function.
        3) variables: include the names of the variables to use in the binning process. Variables[0] is the independent variable and the other two are the dependent variables.
        4) percentiles: for each regime (i,j) we calculate various statistics include the mean, std, max and variaous percentiles.

    OUTPUT:
        - 2-dimensional array with the number (freq_bin), mean intensity (int_bin), maximum (max_bin), etc. max5_bin,std_bin,per_bin

    Created: 02/08/2021                                                                                                                                                                                                                       
    Last Modification:                                                                                                                                                                                                                    
    - 02/08/2021:
    - 12/04/2023: TW. Modified to work to datetime binning.
                                                                                                                                                                                                                              
    """
    min_n = 5 # minimum number of values in a bin so the statistics are calculated.

    #print('\n', '\n', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    #print('  --> CALCULATING BINNING ')
    #print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    #print('Independent variable: ', variables[0])
    prm = var_comp[0, :]
    #print('Large-scale variables:', variables[1:])
    #for ivar, var in enumerate(variables):
    #    print('  -> ', var, '(min, max): ', np.min(var_comp[ivar, :]), np.max(var_comp[ivar, :]))

    #print('\n', '  --> Digitizing large-scale variables')
    # Extract hour and month components from var_comp[1,:] and var_comp[2,:]
    hours = np.array([date.hour for date in var_comp[1, :]])
    months = np.array([date.month for date in var_comp[2, :]])
    # Bin the hour and month components separately
    hour_indices = np.digitize(hours, bins_all[variables[1]])
    month_indices = np.digitize(months, bins_all[variables[2]])

    # Define arrays for the output
    freq_bin = np.zeros((bins_all[variables[1]].shape[0] - 1, bins_all[variables[2]].shape[0] - 1))
    int_bin = np.zeros((bins_all[variables[1]].shape[0] - 1, bins_all[variables[2]].shape[0] - 1))
    max_bin = np.zeros((bins_all[variables[1]].shape[0] - 1, bins_all[variables[2]].shape[0] - 1))
    max5_bin = np.zeros((bins_all[variables[1]].shape[0] - 1, bins_all[variables[2]].shape[0] - 1))
    std_bin = np.zeros((bins_all[variables[1]].shape[0] - 1, bins_all[variables[2]].shape[0] - 1))
    per_bin = np.zeros((len(percentiles), bins_all[variables[1]].shape[0] - 1, bins_all[variables[2]].shape[0] - 1))

    #print('\n', '  --> Calculating intensity, frequency and other statistics at each regime')
    for cb1, bin1 in enumerate(bins_all[variables[1]][:-1]):
        for cb2, bin2 in enumerate(bins_all[variables[2]][:-1]):
            # Calculate frequency in the bin
            freq_bin[cb1, cb2] = np.sum((hour_indices == cb1 + 1) & (month_indices == cb2 + 1))
            # If the frequency is less than the minimum number, set statistics to missing values
            if freq_bin[cb1, cb2] <= min_n:
                int_bin[cb1, cb2] = np.nan
                max_bin[cb1, cb2] = np.nan
                max5_bin[cb1, cb2] = np.nan
                std_bin[cb1, cb2] = np.nan
                for c_per, per in enumerate(percentiles):
                    per_bin[c_per, cb1, cb2] = np.nan
            else:
                # Calculate mean, max, max5, and std in the bin
                int_bin[cb1, cb2] = np.mean(prm[(hour_indices == cb1 + 1) & (month_indices == cb2 + 1)])
                max_bin[cb1, cb2] = np.max(prm[(hour_indices == cb1 + 1) & (month_indices == cb2 + 1)])
                max5_bin[cb1, cb2] = np.mean(np.sort(prm[(hour_indices == cb1 + 1) & (month_indices == cb2 + 1)])[-5:])
                std_bin[cb1, cb2] = np.std(prm[(hour_indices == cb1 + 1) & (month_indices == cb2 + 1)])
                # Calculate percentiles in the bin
                for c_per, per in enumerate(percentiles):
                    per_bin[c_per, cb1, cb2] = np.percentile(prm[(hour_indices == cb1 + 1) & (month_indices == cb2 + 1)], float(per))

    int_bin = np.ma.masked_invalid(int_bin)
    std_binhist = np.ma.masked_invalid(std_bin)
    if prob:
        freq_bin = freq_bin / (np.sum(freq_bin))
    return freq_bin, int_bin, max_bin, max5_bin, std_bin, per_bin

"""                                                                                                                                                                                                 
This function uses the Obukhov length scale L (l_array) to determine the stability in the surface later based on the paper of                                                                       
Gryning et al (https://link.springer.com/article/10.1007/s10546-007-9166-9)                                                                                                                         
"""
def get_stability(zl_array,reg,threshold=0.02):
    if reg=='all':
        select=zl_array>-10**(10)
    if reg=='neutral':
        select=np.abs(zl_array)<=threshold
    if reg=='stable':
        #select=(zl_array>=0) & (zl_array<500.)                                                                           
        select=zl_array>threshold
    if reg=='unstable':
        #select=(zl_array<0) & (zl_array>-500.)                                                                       
        select=zl_array<-threshold
    return select
