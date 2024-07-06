import datetime as dt

class params:
    # Paths
    path_data = '/pampa/diluca/SAVE-WINDS/data/'
    path_dates = '/pampa/diluca/SAVE-WINDS/dates/'
    path_z0 = '/pampa/diluca/SAVE-WINDS/z0_GEM/'
    path_stats = '/pampa/diluca/SAVE-WINDS/stats/'
    path_amf_corrected = '/pampa/diluca/SAVE-WINDS/AMF_corrected/'
    path_out = '/pampa/diluca/SAVE-WINDS/png/'
    variables_dic=[[ 'WS', 'USTAR', 'ZL', 'TA', 'H', 'LE','PA','RH']]
    simulations=['NAM-11m_ERA5_GEM5_CLASS_NEWVEG_newP3-SCPF_SN8hrs_Lakes','NAM-11m_GEM511_ERA5_CLASS_USGS_Delage_SN8','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_CLASSusesISBA-Z0','NAM-11m_GEM511_ERA5_ISBA_USGS_SN8','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD']
    
    i_date=dt.datetime(2015,1,1)
    f_date=dt.datetime(2020,12,31)

    # Simulation parameters
    sim_color={"NAM-11m_GEM511_ERA5_CLASS_USGS_Delage_SN8":'r',"NAM-11m_GEM511_ERA5_ISBA_USGS_SN8":'b',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8":'orange',"NAM-11m_ERA5_GEM5_CLASS_NEWVEG_newP3-SCPF_SN8hrs_Lakes":'green',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD":'magenta','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_newDPTH_FVAP_5Ksat_5yrs':'pink','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD_ALPHA45':'pink','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_CLASSusesISBA-Z0':'purple'}
    sim_name={"NAM-11m_GEM511_ERA5_CLASS_USGS_Delage_SN8":'GEM51-C-DEL',"NAM-11m_GEM511_ERA5_ISBA_USGS_SN8":'GEM51-I-BEL',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8":'GEM51-C-BEL',"NAM-11m_ERA5_GEM5_CLASS_NEWVEG_newP3-SCPF_SN8hrs_Lakes":'GEM50-C-DEL',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD":'GEM51-C-DEL-TOFD','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_newDPTH_FVAP_5Ksat_5yrs':'GEM51-C-DEL-TOFD-LE','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD_ALPHA45':'GEM51-C-DEL-TOFD-45','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_CLASSusesISBA-Z0':'GEM51-C-BEL-Iz0'}
    sim_line={"NAM-11m_GEM511_ERA5_CLASS_USGS_Delage_SN8":'solid',"NAM-11m_GEM511_ERA5_ISBA_USGS_SN8":'dashed',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8":'solid',"NAM-11m_ERA5_GEM5_CLASS_NEWVEG_newP3-SCPF_SN8hrs_Lakes":'solid',"NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD":'solid','GEM51C-DEL-TOFD':'solid','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_TOFD_ALPHA45':'solid','NAM-11m_GEM511_ERA5_CLASS_USGS_SN8_CLASSusesISBA-Z0':'solid'}
    sources_col={'AMF':'grey','AMFc-BEL_withoutz0':'k','AMFc-BEL_withz0':'k'}

    #Other parameters
    #regimes=['all','neutral','stable','unstable']
    #sources=['AMF','AMFc-BEL_withoutz0','AMFc-BEL_withz0','GEM']
    #z0s=['all','low','mid','hig','p01-p1']
    #seasons=['DJF','MAM','JJA','SON']
    regimes=['all']
    sources=['AMF','AMFc-BEL_withoutz0','AMFc-BEL_withz0','GEM']
    z0s=['all']
    seasons=['DJF','MAM','JJA','SON']

