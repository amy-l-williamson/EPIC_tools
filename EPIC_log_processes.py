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


def grab_empty_log_dataframes():
    
    '''
    N: New trigger information
    
        sta:     station,        chan:    channel,     net:    network, 
        loc:     location,       lat:     latitude,    lon:    longitude, 
        time:    trigger time,   up:      update index,
        c:       Component code (Z,E,N),
        rsp:     The sample offset from the trigger sample of the most recent data available to the waveform processor
        tsp:     The sample offset from the trigger sample of the maximum taup measurement on the c-component
        tp:      The maximum taup measurement on the c-component
        tpsnr:   The signal to noise for the maximum taup measurement on the c-component
        dsp:     The sample offset from the trigger sample of the maximum displacement value on the c-component
        pd:      The maximum displacement on the c-component
        pdsnr:   The signal to noise for the maximum displacement on the c-component
        vsp:     The sample offset from the trigger sample of the maximum velocity value on the c-component
        pv:      The maximum velocity on the c-component
        pvsnr:   The signal to noise for the maximum velocity on the c-component
        asp:     The sample offset from the trigger sample of the maximum acceleration value on the c-component
        pa:      The maximum acceleration on the c-component
        pasnr:   The signal to noise for the maximum accelration on the c-component
        assoc:   one of:
                 UNASSOC: the trigger is unassociated
                 NO_ZCOMP: no Z-component measurements are available yet
                 PD_SMALL: trigger cannot be used because log(Pd) < MinPD
                 PD_LARGE: trigger cannot be used because log(Pd) > MaxPD
                 PV_SMALL: trigger cannot be used because log(Pv) < MinPV
                 PV_LARGE: trigger cannot be used because log(Pv) > MaxPV
                 PA_SMALL: trigger cannot be used because log(Pa) < MinPA
                 TP_SMALL: trigger cannot be used because log(Taup) < MinTP
                 TP_LARGE: trigger cannot be used because log(Taup) > MaxTP
                 NEAR_SRC: trigger is associated to an event using the near-source association criterion
                 TTWINDOW: trigger is associated to an event using the TT-window association criterion
                 MULTISTA: trigger is associated to an event using the multiStationEvent method
                 TWOSTATN: trigger is associated to an event using the twoStationEvent method
                 ONESTATN: trigger is associated to an event using the singleStationEvent method
       tel:      1 - trigger time is within a teleseismic event window, (window_start-teleseismic_secs <= trigger_time < window_end+teleseismic_secs)
                 2 - trigger time is within +/- teleseismic_sec of the predicted phase arrival time, computed on distance to the trigger station.
                 0 - trigger time has not been identified as teleseismic.
       tsec:     If tel is 1 or 2, this is the trigger time minus predicted teleseismic arrival time.
       plen:     the number of samples in the packet that contains the trigger sample
       sps:      the sample rate for this channel
       offset:   the offset of the trigger (seconds) from the first sample in the packet
       arrtime:  the arrival time of the trigger-packet at the waveform processor (packet-arrival-time) minus the trigger time.
       protime:  the time that the trigger-packet was passed to the WP processing thread minus the packet-arrival-time
       fndtime:  the time that the trigger was detected minus the packet-arrival-time
       quetime:  the time that the trigger entered the WP message-send-queue minus the packet-arrival-time
       sndtime:  the time that the trigger message was passed to the activemq send function minus the packet-arrival-time
       e2time:   the time that trigger message arrived at the Event Associator and was placed in the trigger-buffer minus the packet-arrival-time
       buftime:  the time that the trigger was removed from the trigger-buffer and inserted into the trigger pool minus the packet-arrival-time
       zc:       Minimum number of zero-crossings
    '''
    N_cols = ['sta', 'chan', 'net', 'loc', 'lat', 'lon', 'time', 'up', 'c',  'rsp', 
                'tsp', 'tp', 'tpsnr', 'dsp', 'pd', 'pdsnr',  'vsp','pv','pvsnr',  'asp', 'pa', 'pasnr', 'assoc', 'tel', 'tsec', 
                'plen', 'sps', 'toffset', 'arrtime', 'protime', 'fndtime','quetime', 'sndtime',  'e2time', 'buftime', 'zc',
                'ne_to_z',  'acc_range']
    N_df = pd.DataFrame(columns =  N_cols)
    #################################################################################################
    
    '''
    U: Updated trigger information
    sta:     station,        chan:    channel,     net:    network, 
    loc:     location,       lat:     latitude,    lon:    longitude, 
    time:    trigger time,   up:      update index,
    c:       Component code (Z,E,N),
    rsp:     The sample offset from the trigger sample of the most recent data available to the waveform processor
    tsp:     The sample offset from the trigger sample of the maximum taup measurement on the c-component
    tp:      The maximum taup measurement on the c-component
    tpsnr:   The signal to noise for the maximum taup measurement on the c-component
    dsp:     The sample offset from the trigger sample of the maximum displacement value on the c-component
    pd:      The maximum displacement on the c-component
    pdsnr:   The signal to noise for the maximum displacement on the c-component
    vsp:     The sample offset from the trigger sample of the maximum velocity value on the c-component
    pv:      The maximum velocity on the c-component
    pvsnr:   The signal to noise for the maximum velocity on the c-component
    asp:     The sample offset from the trigger sample of the maximum acceleration value on the c-component
    pa:      The maximum acceleration on the c-component
    pasnr:   The signal to noise for the maximum accelration on the c-component
    assoc:   see above (N dataframe)
    zc:      Minimum number of zero-crossings
    '''
    U_cols =  ['sta', 'chan', 'net', 'loc', 'lat', 'lon', 'time', 'up', 'c',  'rsp',  
                'tsp', 'tp', 'tpsnr', 'dsp', 'pd', 'pdsnr',  'vsp','pv', 'pvsnr',  'asp', 'pa', 'pasnr', 'assoc','zc']
    U_df = pd.DataFrame(columns = U_cols)
    #################################################################################################
    # R: lines in log file
    R_cols=['unkn0','unkn1','unkn2','unkn2','sta','chan','net','loc','lat','lon','time',
             'unkn3','unkn4','unkn5','unkn6','unkn7','unkn8','unkn9','unkn10','assoc','unkn11','unkn12',
             'unkn13','unkn14','unkn15','unkn16','unkn17','unkn18','unkn19','unkn20','unkn21','unkn22']
    R_df = pd.DataFrame(columns=R_cols)
    #################################################################################################
    # format all the different summary files
    '''
    E:I: EVENT INFORMATION
       eventid:   Event id: < 0 for events that are never alerted; > 0 for event that are alerted
       ver:       Event version (update) number
       evlat:     latitude,       evlon:     longitude,       dep:       depth, 
       mag:       magnitude,      time:      origin time, 
       latu:      latitude uncertainty,       lonu:      longitude uncertainty
       depu:      depth uncertainty,          magu:      magnitude uncertainty, 
       timeu:     origin time uncertainty, 
       lk:        event likelihood
       nTb:       number of triggers in the build,          nSb: number of stations in the build
       nT:        number of triggers in the current alert,   nS:  number of stations in the current alert
       ave:       average station travel time error, 
       rms:       root-mean-square of station travel time error,
       fitok:     (ave < MaxMisfit)
       splitok:   1: this is not a split event, 
                  0: this event was rejected because it is a split event
       near:      number of near stations, counting clusters as one station
       statrig:   number of near stations that triggered, counting clusters as one station
       active:    number of near stations that are active, counting individual stations
       inact:     number of near stations that are inactive, counting individual stations
       nsta:      number of near active plus inactive stations
       percnt:    statrig / near,  
       prcntok:   (percnt >= TrigStaPercent)
       mindist:   minimum associated station distance, 
       maxdist:   maximum associated station distance, 
       distok:    (min_dist < MaxMagDistkm)
       azspan:    maximum station azimuth minus minimum station azimuth
       Mok:       (AlertMinMag < magnitude < AlertMaxMag)
       nSok:      (nT >= AlertMinTrigs && nS >= AlertMinStats && (if within a AlertMinStaRegion: nS >= region_min_sta))
       Lok:       (current location is not on the edge of the search grid)
       Tdif:      log10(pdaverage) - (-3.728 + 2.817 * log10(tpaverage))
       tpave:     log10(tpaverage)
       pdave:     log10(pdaverage)
       TF:        trigger teleseismic code from the multi-window/filter-band teleseismic discriminator
                  0 - trigger is NOT teleseismic
                  1 - trigger is teleseismic
                  2 - waiting for window length of TFMinWindowIndex
                  3 - insufficient data before the trigger for the filter window.
       Tok:      (log10(pdaverage) < (-3.728 + 2.817 * log10(tpaverage))
       Azok:     Event azimuth span >= minAzSpan
       Aok:      Alert passed all criteria and can be sent
       Ast:      Alert actually sent (alerts are not sent more that once per second)
       alerttime:Time that the alert criteria were checked
    '''
    EI_cols = ['eventid', 'version', 'event lat', 'event lon', 'depth',
               'mag','time','latu','lonu','depu', 'magu','timeu','lk','nTb','nSb', 'nT', 'nS', 'ave',
               'rms','fitOK','splitOK','near','statrig','active','inact', 'nsta','percent','prcntOK',
               'mindist','maxdist','distOK','azspan','MOK','nSOK','LOK','Tdif','tpave','pdave', 'TF',
               'TOK','AZOK','AOK','Ast','alert_time']
    EI_df = pd.DataFrame(columns=EI_cols)
    #################################################################################################
    # E:I:T: lines in log file
    '''
    eventid:      Event id
    ver:          Event version 
    (update)      number
    order:        The order that the triggers were associated with the event
    sta:          station name,      chan:    channel name, 
    net:          network name,      loc:     location name, 
    lat:          station latitude,  lon:     station longitude
    trigger_time: trigger time
    rsmp:         the number of data samples more recent than the trigger sample 
                  that were available to the waveform processor at detection time
    tsmp:         the sample offset from the trigger of the maximum value of Taup.
    log_taup:     log(maximum Taup),
    taup_snr:     SNR at the time of the maximum Taup
    dsmp:         the sample offset from the trigger of the maximum value of displacement (Pd)
    log_pd:       log(abs(maximum displacement)), pd_snr: SNR at the time of the maximum displacement
    pd_snr:       Pd signal-to-noise at dsmp
    assoc:        association code 
    tpmag:        Magnitude computed from log_taup measurements only
    utpm:         tpmag is used
    pdmag:        Magnitude computed from log_pd measurements only
    updm:         pdmag is used
    uch:          channel is okay for magnitude computation
                  (net=="CI" || net=="AZ" || chan_2nd_char='L' || chan_2nd_char='N' || chan_2nd_char=='H')
    ukm:          trigger distance is okay for magnitude computation (see E2Magnitude.cc)
    upd:          (Pdmag > 0 && Pdmag < 9)
    ups:          (Pd-SNR >= MinPdSNR)
    utp:          (Tpmag > 0 && Tpmag < 9)
    uts:          (Tp-SNR >= MinTaupSNR)
    tel:          1 - trigger time is within a teleseismic event window, (window_start-teleseismic_secs <= trigger_time < window_end+teleseismic_secs)
                  2 - trigger time is within +/- teleseismic_sec of the predicted phase arrival time, computed using distance to the trigger station.
                  0 - trigger time has not been identified as teleseismic.
    tsec:         If tel is 1 or 2, this is the trigger time minus predicted teleseismic arrival time.
    distkm:       trigger distance, azimuth: trigger azimuth, tterr: travel-time error
    azimuth:      the event to station azimuth
    TF:           trigger teleseismic code from multi-window/filter-band discriminator
                  0 - trigger is NOT teleseismic
                  1 - trigger is teleseismic
                  2 - waiting for window length of TFMinWindowIndex
                  3 - insufficient data before the trigger for the filter window.
    tterr:        trigger time minus predicted trigger time using the current event location
    plen:         trigger packet length, 
    sps:          channel sample rate
    toffset:      the offset of the trigger (seconds) from the first sample in the packet
    arrtime:      the arrival time of the trigger-packet at the waveform processor (packet-arrival-time) minus the trigger time.
    protime:      the time that the trigger-packet was passed to the WP processing thread minus the packet-arrival-time
    fndtime:      the time that the trigger was detected minus the packet-arrival-time
    quetime:      the time that the trigger entered the WP message-send-queue minus the packet-arrival-time
    sndtime:      the time that the trigger message was passed to the activemq send function minus the packet-arrival-time
    e2time:       the time that trigger message arrived at the Event Associator and was placed in the trigger-buffer minus the packet-arrival-time
    buftime:      the time that the trigger was removed from the trigger-buffer and inserted into the trigger pool minus the packet-arrival-time
    alert:        the time of the alert minus the packet-arrival-time
    zc
    ne_to_z
    acc_range
    '''
    EIT_cols = ['eventid', 'version', 'update', 'order', 'sta', 'chan',
               'net','loc','lat','lon','trigger time','rsmp', 'tsmp', 'log taup', 'taup snr', 'dsmp',
               'log pd','pd snr','assoc','tpmag','utpm','pdmag','updm','uch','ukm','upd','ups','utp',
               'uts','tel','tsec','distkm','azimuth','TF','tterr', 'azerror', 'incid', 'plen', 'sps',
               'toffset', 'arrtime', 'protime', 'fndtime', 'quetime', 'sndtime', 'e2time', 'buftime', 
               'alert', 'zc', 'ne_to_z',  'acc_range']
    EIT_df = pd.DataFrame(columns=EIT_cols)
    
    #################################################################################################
    '''
    L:E: Location Algorithm
    eventid:      event located, 
    ver:          event version, 
    nT:           number of triggers, 
    s:            where location routine was called in the code (1,2,or 3)
    lat0, lon0, dep0, time0:    initial values; 
    lat, lon, dep, time:        new location, 
    ddist:        distance(km) from old to new location
    avefit:       average trigger travel-time error, 
    rmsfit:       route-mean-square of the travel time error, 
    nT:           number of triggers, 
    nS:           number of stations

    '''
    LE_cols = ['eventid','ver', 's','lat0','lon0','dep0','time0',
                'lat','lon','dep','time', 'ddist',  'avefit',  'rmsfit',  'nT',  'nS']
    LE_df = pd.DataFrame(columns = LE_cols)
    
    #################################################################################################
    '''
    L:T: Location Triggers
    U:      trigger was used in the location, 
    dist:   trigger distance to new location, 
    tt:     travel-time to new location, 
    tterr:  travel-time error

    '''
    LT_cols = ['eventid', 'ver','nT', 'index', 'sta', 'chan', 'net', 
                                                   'loc', 'lat','lon', 'U',  'dist', 'tt',   'tterr']
    LT_df = pd.DataFrame(columns = LT_cols)

    #################################################################################################
    '''
    B: Event Build Summary
    B:N: starting station
    loop:      search algorithm loop (1 or 2), 
    ord:       station-added index; 
    net, sta, ch, loc: station tested
    evlat, evlon, dep, evtime: current event location, 
    ttok:       tt accepted, 
    tt:         station travel time (must be <= MultiStaMaxDist/PVelocity)
    dok:        distance accepted, dist: station-to-event distance (must be <= MultiStaMaxDist)
    msta:       station with max station-to-station tt-residual; 
    maxd:       distance to this sta; 
    mttres:     tt-residual to this sta; 
    mttok:      mttres accepted
    trig:       add trigger, 
    ns:         number of stations, 
    nns:        num stations with trigger, 
    nt:         num triggers, 
    lo:         located,
    avefit:     location average error, 
    rmsfit:     rms error
    '''
    BN_cols = ['eventid', 'loop', 'ord','net', 'sta', 'ch', 'loc',
                'evlat', 'evlon',  'dep', 'evtime', 'ttok', 'tt', 'dok', 'dist','msta', 'maxd',
                'mttres', 'mttok', 'trig', 'ns', 'nns', 'nt', 'lo', 'avefit',  'rmsfit']
    BN_df = pd.DataFrame(columns = BN_cols )
    #################################################################################################
    '''
    E:S Station Count Summary
        eventid, ver, evlat, evlon, time: counting triggered stations for this event
        mindist:       minimum event-to-station distance, 
        maxdist:       maximum station distance
        percnt:        the percentage of active stations that triggered (sta_trig_cnt/near_sta_cnt)
        near_sta_cnt:  number of active stations, counting a cluster as one station
        sta_trig_cnt:  number of triggered stations, counting a triggered cluster as one station
        active:        number of activem stations, 
        inactive:      number of inactive stations, 
        nsta:          total number of stations
    '''
    # E:S: lines in log file
    ES_cols = ['eventid', 'ver', 'evlat',  'evlon',    'time', 
                'mindist', 'maxdist', 'percnt', 'near_sta_cnt','sta_trig_cnt', 'active',   'inactive', 
                'nsta']
    ES_df = pd.DataFrame(columns = ES_cols)
    #################################################################################################
    
    '''
    E:C:    Station Count Details
    eventid, ver:         counting triggered stations for this event
    sta, net, lat, lon:   a station whose distance to the event is less than the distance to the fartherest triggered station
    cluster:              station is a member of a cluster, 
    dist:                 station to event distance, 
    tt:                   predicted travel time to station
    time:                 last_data_time = time of most recent gmpeak packet from the station
    time_check:           last_data_time - (event_time + tt - MaxStationDelay). If this is positive, the station is counted as active
    active:               counted as active, 
    trig:                 triggered, 
    clu:                  part of a cluster, 
    ctrig:                the cluster triggered (40% of cluster stations triggered)
    '''
    EC_cols = ['eventid','ver','sta','net', 'lat', 'lon', 
                'cluster', 'dist', 'tt', 'time',  'time_check', 'active', 'trig', 'clu', 'ctrig']
    EC_df = pd.DataFrame( columns = EC_cols)
    #################################################################################################
    # A:S: lines in log file
    AS_cols = ['eventid', 'ver','evlat','evlon', 'time', 'mindist', 'maxdist',
                                    'percnt', 'near_sta_cnt', 'sta_trig_cnt', 'active', 'inactive', 'nsta']
    AS_df = pd.DataFrame(columns = AS_cols)
    #################################################################################################
    # A:C: lines in log file
    AC_cols = ['eventid', 'ver', 'sta', 'net','lat', 'lon', 'cluster', 'dist', 
                                             'tt', 'time', 'time_check','active', 'trig', 'clu', 'ctrig']
    AC_df = pd.DataFrame(columns = AC_cols)
    #################################################################################################
   
    # F: lines in log file
    F_cols = ['sta', 'chan','net', 'loc',  'lat','lon','trigger-time',  'U', 'C',  'lead',
                'lag', 'npts',  'Pgv1',  'T1',  'Pgv2', 'T2',  'Pgv3', 'T3',  'Pgv4', 'T4', 'Pgv5', 'T5',  'Pgv6', 'T6',  'Pgv7', 'T7',
                'Pgv8', 'T8',  'Pgv9', 'T9', 'TTest', 'TF']
    F_df = pd.DataFrame(columns = F_cols)


    
    return(N_df,U_df,R_df,EI_df,EIT_df,LE_df,LT_df,BN_df,ES_df,EC_df,AS_df,AC_df,F_df,
           N_cols,U_cols,R_cols,EI_cols,EIT_cols,LE_cols,LT_cols,BN_cols,ES_cols,
           EC_cols,AS_cols,AC_cols,F_cols,R_cols)


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

def add_to_LT_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [int(l[1]),int(l[2]), int(l[3]), int(l[4]), l[5], l[6], l[7], l[8], float(l[9]), 
                      float(l[10]), int(l[11]),  float(l[12]), float(l[13]),   float(l[14])]

    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)

def add_to_BN_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [int(l[1]),int(l[2]),int(l[3]), l[4], l[5], l[6], l[7],float(l[8]), float(l[9]),  
                      float(l[10]), l[11], int(l[12]), float(l[13]), int(l[14]), float(l[15]),l[16], 
                      float(l[17]),  float(l[18]), int(l[19]), float(l[20]),float(l[21]), float(l[22]), 
                      float(l[23]), float(l[24]),  float(l[25]),  (l[26]).split('\n')[0]]

    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)

def add_to_ES_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [int(l[1]),int(l[2]),float(l[3]),float(l[4]),l[5],float(l[6]), float(l[7]),
                      float(l[8]) , int(l[9]),int(l[10]), int(l[11]), int(l[12]), int(l[13])]

    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)

def add_to_EC_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [int(l[1]),int(l[2]),l[3],l[4], float(l[5]), float(l[6]), l[7], float(l[8]),
                      float(l[9]), l[10],int(l[11]), int(l[12]), int(l[13]), int(l[14]), int(l[15])]

    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)

def add_to_AS_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [int(l[1]),int(l[2]),float(l[3]),float(l[4]), l[5], float(l[6]), float(l[7]),
                      float(l[8]), int(l[9]),int(l[10]), int(l[11]), int(l[12]), int(l[13])]

    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)

def add_to_AC_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [int(l[1]),int(l[2]),l[3],l[4], float(l[5]), float(l[6]), l[7], float(l[8]), 
                      float(l[9]), l[10], int(l[11]), int(l[12]), int(l[13]), int(l[14]), int(l[15])]

    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)

def add_to_N_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [l[1], l[2], l[3], l[4], float(l[5]), float(l[6]), l[7], int(l[8]), l[9],  
                      int(l[10]),  int(l[11]), float(l[12]), float(l[13]), int(l[14]), float(l[15]), 
                      float(l[16]), int(l[17]),float(l[18]), float(l[19]),  int(l[20]), float(l[21]), 
                      float(l[22]), l[23], int(l[24]), float(l[25]), int(l[26]), float(l[27]), 
                      float(l[28]), float(l[29]), float(l[30]), float(l[31]), float(l[32]),float(l[33]),  
                      float(l[34]), float(l[35]), float(l[36]), float(l[37]),float(l[38])]
    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)


def add_to_U_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [l[1], l[2], l[3], l[4], float(l[5]), float(l[6]), l[7], int(l[8]), l[9],  int(l[10]),
                      int(l[11]), float(l[12]), float(l[13]), int(l[14]), float(l[15]), float(l[16]), 
                      int(l[17]),float(l[18]), float(l[19]),  int(l[20]), float(l[21]), float(l[22]), l[23], 
                      float(l[24])]
    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)

def add_to_F_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [l[1], l[2], l[3], l[4], float(l[5]), float(l[6]), l[7], int(l[8]), l[9],  float(l[10]),
                     float(l[11]), int(l[12]),  float(l[13]),  l[14],  float(l[15]), l[16],  float(l[17]), 
                     l[18],  float(l[19]), l[20], float(l[21]), l[22],  float(l[23]), l[24],  float(l[25]), 
                     l[26],  float(l[27]), l[28],  float(l[29]), l[30], float(l[31]), l[32]]
    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)

def add_to_R_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [int(l[1]), int(l[2]), int(l[3]), int(l[4]), str(l[5]), str(l[6]), str(l[7]), str(l[8]), 
                      float(l[9]),float(l[10]),str(l[11]), int(l[12]), int(l[13]),float(l[14]), float(l[15]), 
                      int(l[16]), float(l[17]),float(l[18]), str(l[19]),l[20],l[21],float(l[22]), int(l[23]), 
                      int(l[24]), int(l[25]),int(l[26]), int(l[27]), int(l[28]), int(l[29]), int(l[30]),
                      (l[30]),l[31].split('\n')[0]]

    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)

def process_log(log_dir,out_dir,fname,handle,target_id = 'all', save_non_events=False, basic_process = True):
    
    
    # OPEN THE LOG FILE
    try:   file1 = open(log_dir+fname, 'r');  Lines = file1.readlines()
    except: print('ERROR: cannot open log file. Exiting.'); 
    
    # get all the summary dataframes
    (N_df,U_df,R_df,EI_df,EIT_df,LE_df,LT_df,BN_df,ES_df,EC_df,AS_df,AC_df,F_df,
           N_cols,U_cols,R_cols,EI_cols,EIT_cols,LE_cols,LT_cols,BN_cols,ES_cols,
           EC_cols,AS_cols,AC_cols,F_cols,R_cols)  = grab_empty_log_dataframes()
    
    
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
        
        if basic_process == False:
        
            if ('L:E:' in l):
                if (target_id == 'all') & (save_non_events==False) & (int(l[1])>=0):  LE_df= add_to_LE_dataframe(LE_df,LE_cols,l)
                elif (target_id == 'all') & (save_non_events==True):                  LE_df= add_to_LE_dataframe(LE_df,LE_cols,l)
                elif l[1]==target_id:                                                 LE_df= add_to_LE_dataframe(LE_df,LE_cols,l)
            
            
            
            
            if ('N:' in l[0]) & (len(l[0])==2) :  N_df= add_to_N_dataframe(N_df,N_cols,l)
            if ('U:' in l[0]) & (len(l[0])==2) :  U_df= add_to_U_dataframe(U_df,U_cols,l)
            if ('F:' in l[0]) & (len(l[0])==2): F_df= add_to_F_dataframe(F_df,F_cols,l)
            if ('R:' in l) & (len(l)>5):        R_df= add_to_R_dataframe(R_df,R_cols,l)
            
            
            
            
            if ('L:T:' in l):
                if (target_id == 'all') & (save_non_events==False) & (int(l[1])>=0):  LT_df= add_to_LT_dataframe(LT_df,LT_cols,l)
                elif (target_id == 'all') & (save_non_events==True):                  LT_df= add_to_LT_dataframe(LT_df,LT_cols,l)
                elif l[1]==target_id:                                                 LT_df= add_to_LT_dataframe(LT_df,LT_cols,l)
            
            if ('B:' in l[0]) and (len(l[0]) ==2) and (len(l[1])>1):
               if (target_id == 'all') & (save_non_events==False) & (int(l[1])>=0):  BN_df= add_to_BN_dataframe(BN_df,BN_cols,l)
               elif (target_id == 'all') & (save_non_events==True):                  BN_df= add_to_BN_dataframe(BN_df,BN_cols,l)
               elif l[1]==target_id:                                                 BN_df= add_to_BN_dataframe(BN_df,BN_cols,l)
            
            if ('E:S:' in l):
                if (target_id == 'all') & (save_non_events==False) & (int(l[1])>=0):  ES_df= add_to_ES_dataframe(ES_df,ES_cols,l)
                elif (target_id == 'all') & (save_non_events==True):                  ES_df= add_to_ES_dataframe(ES_df,ES_cols,l)
                elif l[1]==target_id:                                                 ES_df= add_to_ES_dataframe(ES_df,ES_cols,l)
            
            if ('E:C:' in l):
                if (target_id == 'all') & (save_non_events==False) & (int(l[1])>=0):  EC_df= add_to_EC_dataframe(EC_df,EC_cols,l)
                elif (target_id == 'all') & (save_non_events==True):                  EC_df= add_to_EC_dataframe(EC_df,EC_cols,l)
                elif l[1]==target_id:                                                 EC_df= add_to_EC_dataframe(EC_df,EC_cols,l)
            
            if ('A:S:' in l):
                if (target_id == 'all') & (save_non_events==False) & (int(l[1])>=0):  AS_df= add_to_AS_dataframe(AS_df,AS_cols,l)
                elif (target_id == 'all') & (save_non_events==True):                  AS_df= add_to_AS_dataframe(AS_df,AS_cols,l)
                elif l[1]==target_id:                                                 AS_df= add_to_AS_dataframe(AS_df,AS_cols,l)
            
            if ('A:C:' in l):
                if (target_id == 'all') & (save_non_events==False) & (int(l[1])>=0):  AC_df= add_to_AC_dataframe(AC_df,AC_cols,l)
                elif (target_id == 'all') & (save_non_events==True):                  AC_df= add_to_AC_dataframe(AC_df,AC_cols,l)
                elif l[1]==target_id:                                                 AC_df= add_to_AC_dataframe(AC_df,AC_cols,l)
                
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
        if basic_process == False:
            
            LT_subset_df = LT_df.iloc[np.where(LT_df['eventid'] == (unique_events[i]))].reset_index(drop=True)
            LT_subset_df.to_csv(out_dir+handle+'_LT_eid'+str(unique_events[i])+'.txt',index=False,sep='\t')
            
            BN_subset_df = BN_df.iloc[np.where(BN_df['eventid'] == (unique_events[i]))].reset_index(drop=True)
            BN_subset_df.to_csv(out_dir+handle+'_BN_eid'+str(unique_events[i])+'.txt',index=False,sep='\t')
            
            ES_subset_df = ES_df.iloc[np.where(ES_df['eventid'] == (unique_events[i]))].reset_index(drop=True)
            ES_subset_df.to_csv(out_dir+handle+'_ES_eid'+str(unique_events[i])+'.txt',index=False,sep='\t')
            
            EC_subset_df = EC_df.iloc[np.where(EC_df['eventid'] == (unique_events[i]))].reset_index(drop=True)
            EC_subset_df.to_csv(out_dir+handle+'_EC_eid'+str(unique_events[i])+'.txt',index=False,sep='\t')
            
            AS_subset_df = AS_df.iloc[np.where(AS_df['eventid'] == (unique_events[i]))].reset_index(drop=True)
            AS_subset_df.to_csv(out_dir+handle+'_AS_eid'+str(unique_events[i])+'.txt',index=False,sep='\t')
            
            AC_subset_df = AC_df.iloc[np.where(AC_df['eventid'] == (unique_events[i]))].reset_index(drop=True)
            AC_subset_df.to_csv(out_dir+handle+'_AC_eid'+str(unique_events[i])+'.txt',index=False,sep='\t')
    
    if basic_process == False:
        U_df.to_csv(out_dir+handle+'_U.txt',index=False,sep='\t')
        N_df.to_csv(out_dir+handle+'_N.txt',index=False,sep='\t')
        F_df.to_csv(out_dir+handle+'_F.txt',index=False,sep='\t')
    
        R_df.to_csv(out_dir+handle+'_R.txt',index=False,sep='\t')

    return(EI_df, EIT_df,LE_df,LT_df,BN_df,ES_df,EC_df,AS_df,AC_df,U_df,N_df,F_df,R_df)

def add_to_M_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [l[3],l[4].split('\n')[0]]

    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)


def add_to_S_dataframe(parent_df,cols,l):
    df = pd.DataFrame(columns=cols)
    df.loc[len(df)]= [l[1],l[2],l[3],l[4],l[5],l[6],-1]

    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", category=FutureWarning)
        parent_df= pd.concat([parent_df,df])
    return(parent_df)


def get_active_stations_log(log_dir,out_dir,fname,handle):
    # OPEN THE LOG FILE
    try:   file1 = open(log_dir+fname, 'r');  Lines = file1.readlines()
    except: print('ERROR: cannot open log file. Exiting.'); 
    
    M_df_cols = ['network','station']
    M_df = pd.DataFrame(columns = M_df_cols)
    
    S_df_cols = ['network','station','channel','loc','sta lat','sta lon','install timestamp']
    S_df = pd.DataFrame(columns = S_df_cols)
    
    
    for k in range(len(Lines)):
        l0 = list(filter(None, Lines[k].split('|')))[-1];    l = list(filter(None, l0.split(' ')))
        
        # if ('M:' in l[0]) & (len(l[0])==2) : 
        #     if (l[1]=='Adding') & (l[2]=='station:'):
        #         M_df = add_to_M_dataframe(M_df,M_df_cols,l)
                
        if ('S:' in l[0]) & (len(l[0])==2) : 
            S_df = add_to_S_dataframe(S_df,S_df_cols,l)
            
    
    # remove duplicates
    S_df = S_df.drop_duplicates()
    
    S_df.to_csv(out_dir+'/'+handle+'_S.txt',index=False,sep='\t')
    #M_df.to_csv(out_dir+handle+'_M.txt',index=False,sep='\t')
    

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
