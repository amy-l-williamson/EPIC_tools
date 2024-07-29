#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 13:29:31 2023

@author: amy
"""

import pandas as pd
import numpy as np
import warnings
import glob


def process_log(log_dir,out_dir,fname,handle,target_id = 'all', save_non_events=False):
    def add_to_EI_dataframe(parent_df,cols,l):
        df = pd.DataFrame(columns=cols)
        df.loc[len(df)]=[int(l[1]),int(l[2]),float(l[3]),float(l[4]),float(l[5]) ,float(l[6]),l[7],
                          float(l[8]),float(l[9]),float(l[10]),float(l[11]),float(l[12]),float(l[13]),
                          int(l[14]),int(l[15]),int(l[16]),int(l[17]),float(l[18]) ,float(l[19]),int(l[20]),
                          int(l[21]),int(l[22]),int(l[23]),int(l[24]),int(l[25]),int(l[26]),float(l[27]),
                          int(l[28]),float(l[29]),float(l[30]),int(l[31]),float(l[32]),int(l[33]),
                          int(l[34]),int(l[35]),float(l[36]),float(l[37]),float(l[38]),float(l[39]),
                          int(l[40]),int(l[41]),int(l[42]),l[43],l[44].split('\n')[0]]
        
        with warnings.catch_warnings(): 
            warnings.filterwarnings("ignore", category=FutureWarning)
            parent_df= pd.concat([parent_df,df])
            
        return(parent_df)

    def add_to_EIT_dataframe(parent_df,cols,l):
        df = pd.DataFrame(columns=cols)
        df.loc[len(df)]= [int(l[1]),int(l[2]),int(l[3]),int(l[4]),l[5],l[6],l[7],l[8],float(l[9]),float(l[10]),
                          l[11],int(l[12]),int(l[13]),float(l[14]),float(l[15]),int(l[16]),float(l[17]),float(l[18]),
                          l[19] ,float(l[20]),int(l[21]),float(l[22]),int(l[23]),int(l[24]),int(l[25]),
                          int(l[26]),int(l[27]),int(l[28]),int(l[29]),int(l[30]),float(l[31]),float(l[32]),float(l[33]),int(l[34]),
                          float(l[35]),float(l[36]),float(l[37]),int(l[38]),float(l[39]),float(l[40]), float(l[41]), float(l[42]),
                          float(l[43]), float(l[44]), float(l[45]),  float(l[46]), float(l[47]), float(l[48]), float(l[49]),
                          float(l[50]),  float(l[51])]
        with warnings.catch_warnings(): 
            warnings.filterwarnings("ignore", category=FutureWarning)
            parent_df= pd.concat([parent_df,df])
        return(parent_df)

    def add_to_LE_dataframe(parent_df,cols,l):
        df = pd.DataFrame(columns=cols)
        df.loc[len(df)]= [int(l[1]),int(l[2]), int(l[3]),float(l[4]),float(l[5]),float(l[6]),l[7],float(l[8]),
                          float(l[9]),float(l[10]),l[11],float(l[12]),float(l[13]),float(l[14]),int(l[15]),int(l[16])]

        with warnings.catch_warnings(): 
            warnings.filterwarnings("ignore", category=FutureWarning)
            parent_df= pd.concat([parent_df,df])
        return(parent_df)

    
    # OPEN THE LOG FILE
    try:   file1 = open(log_dir + fname, 'r');  Lines = file1.readlines()
    except: print('ERROR: cannot open log file. Exiting.'); 
   
   
    EI_cols = ['eventid', 'version', 'event lat', 'event lon', 'depth',
               'mag','time','latu','lonu','depu', 'magu','timeu','lk','nTb','nSb', 'nT', 'nS', 'ave',
               'rms','fitOK','splitOK','near','statrig','active','inact', 'nsta','percent','prcntOK',
               'mindist','maxdist','distOK','azspan','MOK','nSOK','LOK','Tdif','tpave','pdave', 'TF',
               'TOK','AZOK','AOK','Ast','alert_time']
    EI_df = pd.DataFrame(columns=EI_cols)

    EIT_cols = ['eventid', 'version', 'update', 'order', 'sta', 'chan',
               'net','loc','lat','lon','trigger time','rsmp', 'tsmp', 'log taup', 'taup snr', 'dsmp',
               'log pd','pd snr','assoc','tpmag','utpm','pdmag','updm','uch','ukm','upd','ups','utp',
               'uts','tel','tsec','distkm','azimuth','TF','tterr', 'azerror', 'incid', 'plen', 'sps',
               'toffset', 'arrtime', 'protime', 'fndtime', 'quetime', 'sndtime', 'e2time', 'buftime', 
               'alert', 'zc', 'ne_to_z',  'acc_range']
    EIT_df = pd.DataFrame(columns=EIT_cols)
    
    LE_cols = ['eventid','ver', 's','lat0','lon0','dep0','time0',
                'lat','lon','dep','time', 'ddist',  'avefit',  'rmsfit',  'nT',  'nS']
    LE_df = pd.DataFrame(columns = LE_cols)
    


    #-----------------------------------------------------------------------------#    
    for k in range(len(Lines)):
        l0 = list(filter(None, Lines[k].split('|')))[-1];    l = list(filter(None, l0.split(' ')))
        if ('E:I:' in l) or ('E:I:F:' in l):
            
            if (target_id == 'all') & (save_non_events==False) & (int(l[1])>=0): EI_df= add_to_EI_dataframe(EI_df,EI_cols,l)
            elif (target_id == 'all') & (save_non_events==True):                 EI_df= add_to_EI_dataframe(EI_df,EI_cols,l)
            elif l[1]==target_id:                                                EI_df= add_to_EI_dataframe(EI_df,EI_cols,l)
           
        if ('E:I:T:' in l):
            if (target_id == 'all') & (save_non_events==False) & (int(l[1])>=0):  EIT_df= add_to_EIT_dataframe(EIT_df,EIT_cols,l)
            elif (target_id == 'all') & (save_non_events==True):                  EIT_df= add_to_EIT_dataframe(EIT_df,EIT_cols,l)
            elif l[1]==target_id:                                                 EIT_df= add_to_EIT_dataframe(EIT_df,EIT_cols,l)
        
        if ('L:E:' in l):
            if (target_id == 'all') & (save_non_events==False) & (int(l[1])>=0):  LE_df= add_to_LE_dataframe(LE_df,LE_cols,l)
            elif (target_id == 'all') & (save_non_events==True):                  LE_df= add_to_LE_dataframe(LE_df,LE_cols,l)
            elif l[1]==target_id:                                                 LE_df= add_to_LE_dataframe(LE_df,LE_cols,l)
        
       
    # #----------------------------------------------------------------------------
    
    if save_non_events==True:  
        unique_events = np.unique(EI_df['eventid'])
    else:
        unique_events = np.unique(EI_df['eventid'])
        idx = np.where(unique_events >=0)[0]
        unique_events = unique_events[idx]
    
    for i in range(len(unique_events)):
       
        EI_subset_df = EI_df.iloc[ np.where(EI_df['eventid'] == (unique_events[i]))].reset_index(drop=True)
        EI_subset_df.to_csv(out_dir+handle+'_EI_eid'+str(unique_events[i])+'.txt',index=False,sep='\t')
        
        EIT_subset_df = EIT_df.iloc[np.where(EIT_df['eventid'] ==(unique_events[i]))].reset_index(drop=True)
        EIT_subset_df.to_csv(out_dir+handle+'_EIT_eid'+str(unique_events[i])+'.txt',index=False,sep='\t')
        
        LE_subset_df = LE_df.iloc[np.where(LE_df['eventid'] == (unique_events[i]))].reset_index(drop=True)
        LE_subset_df.to_csv(out_dir+handle+'_LE_eid'+str(unique_events[i])+'.txt',index=False,sep='\t')
        

    return()



def make_run_file(log_dir,trigger_fname,epic_id):
    #  required modules
    import pandas as pd       # will work to remove
    import numpy as np
    from datetime import datetime

    station_trigger_log_file = log_dir+'/'+trigger_fname
    
    output_filename = log_dir+'/'+epic_id+'.run'
    #-----------------------------------------------------------------------------#
    v =  np.genfromtxt(station_trigger_log_file,usecols=[1],skip_header=1)
    o =  np.genfromtxt(station_trigger_log_file,usecols=[3],skip_header=1)
    
    station_sta   =  np.genfromtxt(station_trigger_log_file,usecols=[4], dtype='str',skip_header=1)
    station_chan  =  np.genfromtxt(station_trigger_log_file,usecols=[5], dtype='str',skip_header=1)
    station_net   =  np.genfromtxt(station_trigger_log_file,usecols=[6], dtype='str',skip_header=1)
    station_loc   =  np.genfromtxt(station_trigger_log_file,usecols=[7], dtype='str',skip_header=1)
    station_t_str =  np.genfromtxt(station_trigger_log_file,usecols=[10],dtype='str',skip_header=1)
    station_lon   =  np.genfromtxt(station_trigger_log_file,usecols=[9],skip_header=1)
    station_lat   =  np.genfromtxt(station_trigger_log_file,usecols=[8],skip_header=1)
    station_t_str =  np.genfromtxt(station_trigger_log_file,usecols=[10],dtype='str',skip_header=1)
    tdata =          np.genfromtxt(station_trigger_log_file,usecols=[11],skip_header=1)
    station_pd_snr = np.genfromtxt(station_trigger_log_file,usecols=[17],dtype='str',skip_header=1)
    station_updm = np.genfromtxt(station_trigger_log_file,usecols=[22],skip_header=1)
    station_t     = []
    version       = []
    order         = []
    
    
    for i in range(len(station_t_str)):
        station_t = np.append(station_t,datetime.strptime(station_t_str[i], '%Y-%m-%dT%H:%M:%S.%f').timestamp())
        version   = np.append(version,int(v[i]))
        order     = np.append(order,int(o[i]))
    station_logpd = np.genfromtxt(station_trigger_log_file,usecols=[16],skip_header=1)
    station_tterr = np.genfromtxt(station_trigger_log_file,usecols=[34],skip_header=1)
    
    #--------------------------------------------------#
    df = pd.DataFrame({'version':version,
                       'order':order,
                       'station':station_sta,
                       'channel':station_chan,
                       'network':station_net,
                       'location':station_loc,
                       'longitude':station_lon,
                       'latitude':station_lat,
                       'trigger time':station_t,
                       'tterr':station_tterr,
                       'logPd':station_logpd,
                       'Pd SNR':station_pd_snr,
                       'updm':station_updm,
                       'tdata':tdata})
    df = df.astype({'version':'int'})
    df = df.astype({'order':'int'})
    df.to_csv(output_filename,index=False,sep='\t')
    #--------------------------------------------------#    
