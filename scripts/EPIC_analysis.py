#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 12:46:19 2024

@author: amy
"""

def get_loc_error(catalog_lon,catalog_lat,event_lons,event_lats):
    from    obspy.geodetics   import gps2dist_azimuth
    from numpy import zeros, nan
    
    loc_error = zeros(len(event_lons ))*nan
    for i in range(len(event_lons)):    
        loc_error[i],_,_ = gps2dist_azimuth(catalog_lat,  catalog_lon,
                                                                            event_lats[i],event_lons[i])
    loc_error=loc_error/1000
    return(loc_error)


def get_alert_seconds(catalog_time,alert_times):
    from numpy     import zeros
    from datetime  import datetime
    
    alert_seconds = zeros(len(alert_times))
    for i in range(len(alert_times)):
        try:      alert_seconds[i] = datetime.strptime(alert_times[i], '%Y-%m-%dT%H:%M:%S.%f').timestamp() - catalog_time
        except:   alert_seconds[i] = datetime.strptime(alert_times[i], '%Y-%m-%dT%H:%M:%S').timestamp() - catalog_time
    return(alert_seconds)   