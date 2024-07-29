#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: amy
"""
import run_EPIC_offline
import EPIC_log_processes
import EPIC_analysis
import pandas as pd

'''----------------------------------------------------------------------------#
# Step 01: get the log file and put it in  {path}/EPICoffline/log
# on eew1:   QUEREY_EPIC_db.sh   EPIC_id
#
#
#   Here, I will use the Oct 2022 M5.1 Alum Rock earthquake as an example.
#   I put the log file  in {path}/EPICoffline/log
#
#-----------------------------------------------------------------------------'''
log_dir  = '/Users/amy/projects/EPICdb/EPICoffline/log/' 
fname    =  'eew-nc-prod2_20221025.log'



'''-----------------------------------------------------------------------------#
#  the log file for this event contains everything that happened on that day. 
#  This is a ton of extra info that we do not need right now. So we will run a 
#  script that pulls out the relevant information for our target event. What we
#  need is the log file and the eventID that we care about.  
#
#  
#-----------------------------------------------------------------------------'''
eventid = '52918'
tag     = 'example'   # any name
out_dir = '/Users/amy/projects/EPICdb/EPICoffline/log/' # where to save the parsed log file


#EPIC_log_processes.process_log(log_dir,out_dir,fname,tag,target_id = str(eventid), save_non_events=False)



'''----------------------------------------------------------------------------#
# the script EPIC_log_processes.process_log() usually runs within a few seconds. But
# for big files, it can take minutes.
#
# The output are a few text files (that we will read with pandas) containing the 
# event and trigger info for that run.
#-----------------------------------------------------------------------------'''


event_summary_log =pd.read_csv(out_dir+tag+'_EI_eid'+str(eventid)+'.txt',sep='\t')


trigger_fname = tag+'_EIT_eid'+str(eventid)+'.txt'
event_trigger_log =pd.read_csv(out_dir+trigger_fname,sep='\t')


'''-----------------------------------------------------------------------------#
# make a run file
#-----------------------------------------------------------------------------'''


EPIC_log_processes.make_run_file(log_dir,trigger_fname,eventid)




run_file = '/Users/amy/projects/EPICdb/EPICoffline/log/52918.run'
output_event_summary = run_EPIC_offline.run_locate(run_file)










#-----------------------------------------------------------------------------#
# we can plot the results of this event if we want
#  USGS event page:https://earthquake.usgs.gov/earthquakes/eventpage/nc73799091/executive






lon = -121.672;  lat = 37.312; 
origin_time = '2022-10-25T18:42:02'
from datetime import datetime
origin_timestamp = datetime.strptime(origin_time, '%Y-%m-%dT%H:%M:%S').timestamp() 
event_summary_log['d'] = EPIC_analysis.get_loc_error(lon,lat, event_summary_log['event lon'], event_summary_log['event lat'])
event_summary_log['t'] = EPIC_analysis.get_alert_seconds(origin_timestamp,event_summary_log['alert_time'])

output_event_summary['d']= EPIC_analysis.get_loc_error(lon,lat,output_event_summary['event lon'],output_event_summary['event lat'])


def plot_location_error(event_summary_log,lon,lat,origin_time):
    import matplotlib.pyplot as plt
    event_summary_log['d'] = EPIC_analysis.get_loc_error(lon,lat, event_summary_log['event lon'], event_summary_log['event lat'])
    event_summary_log['t'] = EPIC_analysis.get_alert_seconds(origin_timestamp,event_summary_log['alert_time'])
    plot_color = '#0099cc'
    fig,ax = plt.subplots(figsize=(8,4))
    ax.set_ylabel('location error (km)',       fontsize=14)
    ax.set_xlabel('seconds since origin time', fontsize=14)
    ax.grid(linestyle='dashed')
    ax.plot(event_summary_log['t'],event_summary_log['d'],c=plot_color)
    ax.scatter(event_summary_log['t'],event_summary_log['d'],c=plot_color,edgecolor='k',s=50)
    
plot_location_error(event_summary_log,lon,lat,origin_time)
    
    
