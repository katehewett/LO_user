"""
Code to extract corrosive volume from LO domain 

running testing code with two history files on kh personal mac, looks like:
run extract_corrosive_volume_rev1 -gtx cas7_t0_x4b -ro 1 -0 2017.12.12 -1 2017.12.13 -lt lowpass -test True

apogee:
python extract_corrosive_volume_rev1.py -gtx cas7_t0_x4b -ro 3 -0 2013.01.01 -1 2013.12.31 -lt lowpass > a_corr.log & 

First created: January 2023
July 2024, rev1: 
* using argparse so can set gtx and lt 
* added thresholds (0.5, 1, 1.7 Oag)

2 days is ~38 seconds (long(er) than extract_hypoxic_volume_rev1.py because of CO2sys)

"""

# imports
from lo_tools import Lfun, zfun, zrfun
from lo_tools import extract_argfun as exfun
#Ldir = exfun.intro() # this handles the argument passing

from subprocess import Popen as Po
from subprocess import PIPE as Pi

from time import time
import sys 
import os
import argparse
import pandas as pd
import xarray as xr
import numpy as np
import datetime 

pid = os.getpid()
print('Calculating corrosive volumes '.center(60,'='))
print('PID for this job = ' + str(pid))
print('started process at: ' + str(datetime.datetime.now()))

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_lo8b
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
parser.add_argument('-lt', '--list_type', type=str) # list type: hourly, daily, weekly, lowpass
# Optional args for testing and proc:
parser.add_argument('-fb','--false_bottom', default = False, type = Lfun.boolean_string) # places 200m false bottom 
parser.add_argument('-Nproc', type=int, default=10)
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided
argsd = args.__dict__
for a in ['gtagex']:
    if argsd[a] == None:
        print('*** Missing required argument: ' + a)
        sys.exit()
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
        
def messages(mess_str, stdout, stderr):
    # utility function to help with subprocess errors
    try:
        if len(stdout) > 0:
            print(mess_str)
            print(stdout.decode())
    except TypeError:
        pass
    try:
        if len(stderr) > 0:
            print(mess_str)
            print(stderr.decode())
    except TypeError:
        pass
        
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]
    
# output location and filenames + temp dir to accumulate individual extractions
dd_str = '_' + Ldir['ds0'] + '_' + Ldir['ds1']
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'corrosive_volume' / ('LO' + dd_str)
temp_vol_dir = out_dir / 'temp_vol_dir'

Lfun.make_dir(out_dir, clean=False)
Lfun.make_dir(temp_vol_dir, clean=True)

if args.false_bottom == True:
    vol_fn_final = 'LO_corrosive_volume_' + args.list_type + dd_str + '_false_bottom.nc'
elif args.false_bottom == False:
    vol_fn_final = 'LO_corrosive_volume_' + args.list_type + dd_str + '.nc'
        
# get list of files to work on and check before entering loop over all files 
fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])

# put the fields that do not change with time to a dictionary 
tt00 = time()

FF = dict()
G, S, T = zrfun.get_basic_info(fn_list[0])
ds = xr.open_dataset(fn_list[0], decode_times=False) 
FF['Lon'] = G['lon_rho']
FF['Lat'] = G['lat_rho']
FF['DA'] = G['DX'] * G['DY']      # cell horizontal area 
h = ds.h.values  
FF['h'] = h
mask_rho = ds.mask_rho.values
FF['mask_rho']=mask_rho
del ds, h, mask_rho

print('Time to get initial fields = %0.2f sec' % (time()-tt0))
print('Time check: ' + str(datetime.datetime.now()))
sys.stdout.flush()

# loop over all jobs
tt0 = time()
N = len(fn_list)
proc_list = []
for ii in range(N):
    # Launch a job and add its process to a list.
    fn = fn_list[ii]
    ii_str = ('0000' + str(ii))[-5:]
    out_fn = temp_vol_dir / ('CC_' + ii_str + '.nc')
    cmd_list = ['python3', 'get_one_corrosive_volume_rev1.py',
            '-lt',args.list_type,
            #'-fb',str(args.false_bottom),
            '-in_fn',str(fn),
            '-out_fn', str(out_fn)]
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    proc_list.append(proc)
    # If we have accumulated Nproc jobs, or are at the end of the
    # total number of jobs, then stop and make sure all the jobs
    # in proc_list have finished, using the communicate method.
    if ((np.mod(ii,Ldir['Nproc']) == 0) and (ii > 0)) or (ii == N-1):
        for proc in proc_list:
            stdout, stderr = proc.communicate()
            if len(stdout) > 0:
                print('\nSTDOUT:')
                print(stdout.decode())
                sys.stdout.flush()
            if len(stderr) > 0:
                print('\nSTDERR:')
                print(stderr.decode())
                sys.stdout.flush()
        # Then initialize a new list.
        proc_list = []
    # Print screen output about progress.
    if (np.mod(ii,10) == 0) and ii>0:
        print(str(ii), end=', ')
        print('Time spent calculating corrosive vol = %0.2f sec' % (time()-tt0))
        print('Time check: ' + str(datetime.datetime.now()))
        sys.stdout.flush()
    if (np.mod(ii,50) == 0) and (ii > 0):
        print('') # line feed
        print('Time spent calculating corrosive vol = %0.2f sec' % (time()-tt0))
        print('Time check: ' + str(datetime.datetime.now()))
        sys.stdout.flush()
    if (ii == N-1):
        print(str(ii))
        print('Time spent calculating corrosive vol = %0.2f sec' % (time()-tt0))
        print('Time check: ' + str(datetime.datetime.now()))
        sys.stdout.flush()

print('Total processing vol calculations = %0.2f sec' % (time()-tt0))

tt0 = time()
# concatenate the records into one file
# This bit of code is a nice example of how to replicate a bash pipe
pp1 = Po(['ls', str(temp_vol_dir)], stdout=Pi)
pp2 = Po(['grep','CC'], stdin=pp1.stdout, stdout=Pi)
temp_fn = str(temp_vol_dir)+'/'+vol_fn_final # this is all the maps put to one
cmd_list = ['ncrcat','-p', str(temp_vol_dir), '-O', temp_fn]
proc = Po(cmd_list, stdin=pp2.stdout, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if len(stdout) > 0:
    print('\nSTDOUT:')
    print(stdout.decode())
    sys.stdout.flush()
if len(stderr) > 0:
    print('\nSTDERR:')
    print(stderr.decode())
    sys.stdout.flush()
 
print('Time to concat = %0.2f sec' % (time()-tt0))     
sys.stdout.flush()
  
"""
Next we want to repackage these results into one NetCDF file per section, with all times.

Variables in the NetCDF files:
- hyp_dz (mild, severe, anoxic) is depth of the corrosive layer(s) in each cell (t, x, y) [same for all other variables]
- corrosive_dz is the undersaturated layer (t, x, y)

To save space, we save the variables that don't change over time once 
- DA is the area of each cell (x,y)
- h is bathymetric depth 
- ocean_time is a vector of time in seconds since (typically) 1/1/1970.
- Lat and Lon on rho points 
    
A useful tool is isel():
a = ds.isel(p=np.arange(10,15))
"""
print('now saving vars')

tt0 = time()    

ds1 = xr.open_dataset(temp_fn)

ds1['lat_rho'] = (('eta_rho', 'xi_rho'),FF['Lat'],{'units':'degree_north'})
ds1['lat_rho'].attrs['standard_name'] = 'grid_latitude_at_cell_center'
ds1['lat_rho'].attrs['long_name'] = 'latitude of RHO-points'
ds1['lat_rho'].attrs['field'] = 'lat_rho'

ds1['lon_rho'] = (('eta_rho', 'xi_rho'),FF['Lon'],{'units': 'degree_east'})
ds1['lon_rho'].attrs['standard_name'] = 'grid_longitude_at_cell_center'
ds1['lon_rho'].attrs['long_name'] = 'longitude of RHO-points'
ds1['lon_rho'].attrs['field'] = 'lon_rho'

ds1['h'] = (('eta_rho', 'xi_rho'),FF['h'],{'units': 'm'})
ds1['h'].attrs['standard_name'] = 'sea_floor_depth'
ds1['h'].attrs['long_name'] = 'time_independent bathymetry'
ds1['h'].attrs['field'] = 'bathymetry'
ds1['h'].attrs['grid'] =  args.gtagex

ds1['mask_rho'] = (('eta_rho', 'xi_rho'),FF['mask_rho'],{'units': 'mask'})
ds1['mask_rho'].attrs['standard_name'] = 'land_sea_mask_at_cell_center'
ds1['mask_rho'].attrs['long_name'] = 'mask on RHO-points'
ds1['mask_rho'].attrs['flag_values'] = np.array([0.,1.])
ds1['mask_rho'].attrs['flag_meanings'] = 'land water'
ds1['mask_rho'].attrs['grid'] =  args.gtagex

ds1['DA'] = (('eta_rho', 'xi_rho'),FF['DA'],{'units': 'm^2'})
ds1['DA'].attrs['standard_name'] = 'cell area'
ds1['DA'].attrs['long_name'] = 'cell horizontal area'
ds1['DA'].attrs['grid'] =  args.gtagex
    
this_fn = out_dir / (vol_fn_final)
ds1.to_netcdf(this_fn)

# clean up 
#Lfun.make_dir(temp_vol_dir, clean=True)
#temp_vol_dir.rmdir()

print('Time to save and clean-up = %0.2f sec' % (time()-tt0)) 
sys.stdout.flush()

print('Time to run extract_corrosive_volume_rev1.py = %0.2f sec' % (time()-tt00)) 
print('Time check finish: ' + str(datetime.datetime.now()))
sys.stdout.flush()


