#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: amy
"""

import  pandas            as pd       
import  numpy             as np
from    scipy             import interpolate


def get_two_station_location(sta_df):
 
    station_01_lon = sta_df['longitude'].iloc[0]    ;  station_02_lon = sta_df['longitude'].iloc[1]
    station_01_lat = sta_df['latitude'].iloc[0]     ;  station_02_lat = sta_df['latitude'].iloc[1]
    station_01_OT  = sta_df['trigger time'].iloc[0] ;  station_02_OT = sta_df['trigger time'].iloc[1]

    if   station_01_OT <= station_02_OT:   eq_lat = (station_01_lat*2 + station_02_lat)/3;    eq_lon = (station_01_lon*2 + station_02_lon)/3
    elif station_02_OT <= station_01_OT:   eq_lat = (station_02_lat*2 + station_01_lat)/3;    eq_lon = (station_02_lon*2 + station_01_lon)/3
    return([ eq_lon, eq_lat])



def EPIC_locate(trigger_lons,trigger_lats,trigger_times,evlon,evlat):
   
    # open travel time model
    tt_file               = '/Users/amy/projects/EPICdb/EPICoffline/misc/h2p+ak135.080'
    tt_mod                = np.genfromtxt(tt_file,skip_header=1)
    ttf                    = interpolate.interp1d(tt_mod[:,0], tt_mod[:,1])
   
    
    # various semi-permanent variables
    GridSize      = 200              ;  GridSpacing = 1                                 # grid size and spacing (km) on either side of epicenter
    R             = 6378.137         ;  ff = 1./298.257                                 # // flattening factor
    evlatr        = evlat*np.pi/180. ;  r  = R*(1 - ff*np.power(np.sin(evlatr), 2));      # // radius - radius at lat [m]
    mpd           = r*np.pi/180.     ;  f  = mpd*np.cos(evlatr)                           # // mpd - meters per degree

    # construct the search grid about the initial earthquake location (lon0,lat0)
    grid_y   = grid_x = np.arange(-1*GridSize,GridSize+GridSpacing,GridSpacing)
    grid_lon = evlon + grid_x/f;     grid_lat = evlat + grid_y/mpd 
    xx,yy=np.meshgrid(grid_x ,grid_x )
    
    
    n = len(trigger_lons)              # number of active stations
    p = len(grid_x)                    # number of grid nodes in x (or y) direction
    m = p*p                            # total number of grid nodes  (p*p)
    
    #-------------------------------------------------------------------------#
    
    
    # # // get the station coordinates on the local grid centered at lon0, lat0
    stay = np.zeros(n);    stax = np.zeros(n)
    for i in range(n):     stay[i] = (trigger_lats[i]-evlat)*mpd;    stax[i]=  (trigger_lons[i]-evlon)*f
    
    
    # station_distance is a   (m x n) array with the each station's distance (km) from each grid node in the grid search
    station_distance = np.sqrt(  np.square(np.subtract(np.tile(  np.ravel(xx),(n,1)).T,   np.tile(stax,(m,1)) )) 
                               + np.square(np.subtract(np.tile(  np.ravel(yy),(n,1)).T,   np.tile(stay,(m,1)) ))  )  
    
    
    # travel_time is a   (m x n) array with the each station's travel time (s) from grid node to station location
    travel_time       = ttf(np.ravel(station_distance)).reshape(m,n)
    
    
    # average_OT is a (m , ) array with the origin time at each grid node
    #         basically this says that if any grid node is the epicenter, then this is the expected origin time
    average_OT        = np.mean( np.subtract(np.tile(trigger_times,(m,1)) , travel_time),axis=1) 
    
    # trigger_time_calc is a   (m x n) array  and is the (forward) modeled trigger time for each station from each grid node
    trigger_time_calc = np.add(np.tile(average_OT ,(n,1)).T,travel_time )         

    # tt_error is a (m x n) aray and is the misfit  between the observed and forward modeled travel times for each station at each grid node       
    tt_error          = np.abs(np.subtract(trigger_time_calc,np.tile(trigger_times,(m,1))))**2   
    
    # get the likelihood (L(m|d)) L for the model m given the observation of the data d
    # this can actually use a sigma value- I ignore right now, (sigma =1) but could have a step variable or station variable sigma....
    
    # rho is a (m x n) aray 
    rho               = np.exp(-0.5*(tt_error))  
    
    # like is a (m,1) array with the likelihood value (between 0 - 1) at each grid node
    like              = np.prod(rho,axis=1)
    
    
    # misfit grid
    trig_ot =np.subtract(np.tile(trigger_times,(m,1)) , travel_time)
    rms = np.zeros(m)
    for kk in range(trig_ot.shape[1]):    rms+= np.power(trig_ot[:,kk] - average_OT,2)
    rms = rms/n 
    
    
    xx,yy=np.meshgrid(grid_lon ,grid_lat )
    grid_lon_ravel = np.ravel(xx);    grid_lat_ravel = np.ravel(yy)
    
    best_location  =   np.where(like == np.nanmax(like))[0]
    likelihood_lon =   grid_lon_ravel[best_location][0]  
    likelihood_lat =   grid_lat_ravel[best_location][0]
    
    best_location  =   np.where(rms == np.nanmin(rms))[0]
    rms_lon        =   grid_lon_ravel[best_location][0]  
    rms_lat        =   grid_lat_ravel[best_location][0]
    
    bestOT         = average_OT[best_location[0]]
    bestrms        = np.nanmin(rms)
    
    
    return(likelihood_lon,likelihood_lat, rms_lon, rms_lat, bestOT, bestrms)
    

def run_locate(run_file):

    #
    run_df        =   pd.read_csv(run_file,sep='\t')   # read the .run file
    
    # -------     output data frames
    cols = ['version','evlat0','evlon0','event lat','event lon','origin time','misfit']
    output_event_summary = pd.DataFrame(columns = cols)
    
    
    # allow for simulation to stop early
    max_version = run_df['version'].iloc[-1]
    version = 0
    while version < max_version+1:
        
        sta_df         = run_df.iloc[np.where((run_df['version']==version)  & 
                                              (run_df['tterr'] > -999))[0]].reset_index(drop=True).sort_values(by=['order']) 
        
        # location and trigger time of associated stations
        trigger_lons   = np.array(sta_df['longitude'])  
        trigger_lats   = np.array(sta_df['latitude'])
        trigger_times  = np.array(sta_df['trigger time']) 
        
        
        if version ==0:
            (evlon,evlat)  =  get_two_station_location(sta_df)
            
        else:
            evlon = output_event_summary['event lon'].iloc[version-1]
            evlat = output_event_summary['event lat'].iloc[version-1]
            
        
        # run the location algorithm for this step
        (likelihood_lon,likelihood_lat, rms_lon, rms_lat, bestOT, bestrms) = EPIC_locate(trigger_lons,trigger_lats,trigger_times,
                                                                   evlon,evlat)
    
        output_event_summary.loc[len(output_event_summary)] = [version, evlat, evlon, rms_lat, rms_lon, bestOT, bestrms]
        print(str(version)+' of '+str(max_version)+': start lon,lat '+str(np.round(evlon,3))+', '+str(np.round(evlat,3))+' | end lon,lat '+str(np.round(rms_lon,3))+', '+str(np.round(rms_lat,3)))
        version+=1
    return(output_event_summary)
