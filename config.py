import datetime as dt

class params:
    path_in_r = '/home/whittaker/Documents/validateBL/SAVE-WINDS/data/'
    path_in = '/home/whittaker/Documents/validateBL/SAVE-WINDS/AMF_corrected/'
    path_out = 'png/'
    variables_dic=[[ 'WS', 'USTAR', 'ZL', 'TA', 'H', 'LE','PA','RH']]
    simulations=['NAM-11m_ERA5_GEM5_CLASS_NEWVEG_newP3-SCPF_SN8hrs_Lakes','NAM-11m_GEM511_ERA5_CLASS_USGS_Delage_SN8','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_CLASSusesISBA-Z0','NAM-11m_GEM511_ERA5_ISBA_USGS_SN8','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD']
    
    i_date=dt.datetime(2015,1,1)
    f_date=dt.datetime(2020,12,31)
    
    sim_color={"NAM-11m_GEM511_ERA5_CLASS_USGS_Delage_SN8":'r',"NAM-11m_GEM511_ERA5_ISBA_USGS_SN8":'b',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8":'orange',"NAM-11m_ERA5_GEM5_CLASS_NEWVEG_newP3-SCPF_SN8hrs_Lakes":'green',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD":'magenta','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_newDPTH_FVAP_5Ksat_5yrs':'pink','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD_ALPHA45':'pink','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_CLASSusesISBA-Z0':'purple'}
    sim_name={"NAM-11m_GEM511_ERA5_CLASS_USGS_Delage_SN8":'GEM51-C-DEL',"NAM-11m_GEM511_ERA5_ISBA_USGS_SN8":'GEM51-I-BEL',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8":'GEM51-C-BEL',"NAM-11m_ERA5_GEM5_CLASS_NEWVEG_newP3-SCPF_SN8hrs_Lakes":'GEM50-C-DEL',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD":'GEM51-C-DEL-TOFD','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_newDPTH_FVAP_5Ksat_5yrs':'GEM51-C-DEL-TOFD-LE','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD_ALPHA45':'GEM51-C-DEL-TOFD-45','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_CLASSusesISBA-Z0':'GEM51-C-BEL-Iz0'}
    sim_line={"NAM-11m_GEM511_ERA5_CLASS_USGS_Delage_SN8":'solid',"NAM-11m_GEM511_ERA5_ISBA_USGS_SN8":'dashed',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8":'solid',"NAM-11m_ERA5_GEM5_CLASS_NEWVEG_newP3-SCPF_SN8hrs_Lakes":'solid',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD":'solid','GEM51C-DEL-TOFD':'solid','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD_ALPHA45':'solid','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_CLASSusesISBA-Z0':'solid'}
    sources_col={'AMF':'grey','AMFc-BEL_withoutz0':'k','AMFc-BEL_withz0':'k'}
    
    sea_months={'DJF':[12,1,2],'MAM':[3,4,5],'JJA':[6,7,8],'SON':[9,10,11]}
    hm_limits={'WS':[0,8],'USTAR':[0,.6],'ZL':[-.75,.75],'TA':[0,25],'PA':[99,102],'RH':[0,100],'LE':[0,125],'H':[-50,250],'P':[0,.2]}
    mm_limits={'WS':[0,5],'USTAR':[0,.5],'ZL':[-.1,.5],'TA':[0,25],'PA':[99,102],'RH':[0,100],'LE':[0,125],'H':[-50,250],'P':[0,.2]}
    pdf_limits={'WS':[0,20],'USTAR':[0,1.6],'ZL':[-.75,.75],'TA':[-25,45],'PA':[85,110],'RH':[0,100],'LE':[-100,600],'H':[-200,600],'P':[0,20]}
    pdf_res=100
    #var_names={'WS':r'$u$','USTAR':r'$u_{*}$','ZL':r'$zL^{-1}$','TA':r'$T_a$','TA':r'$T_a$','PA':r'$p$','RH':r'$RH$','LE':r'$LE$','H':r'$H$','P':r'$P$'}
    #var_units={'WS':r'm s$^{-1}$','USTAR':r'm s$^{-1}$','ZL':r'-','TA':r'$K$','H':r'W m$^{-2}$','LE':r'W m$^{-2}$','PA':r'$kPa$','RH':r'$\%$','P':r'mm h$^{-1}$'}
