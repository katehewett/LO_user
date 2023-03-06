"""
Code to extract grid info and zetas for plotting 

Future: encorporate to box - 
takes about 1-2 seconds to run 2 history fiels on mac 
~4-5seconds on perigee 

Testing Feb. 2023 - present
"""

# imports
from lo_tools import Lfun, zfun, zrfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

from subprocess import Popen as Po
from subprocess import PIPE as Pi

from time import time
import sys 
import pandas as pd
import xarray as xr
import numpy as np
import pickle

gctag = Ldir['gridname']

fn_list = Lfun.get_fn_list('daily', Ldir, Ldir['ds0'], Ldir['ds1'])

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'bottom_layer'
out_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
temp_dir = out_dir0 / ('temp_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)
Lfun.make_dir(temp_dir, clean=True)

# loop over all jobs
tt0 = time()
N = len(fn_list)
proc_list = []
for ii in range(N):
    # Launch a job and add its process to a list.
    fn = fn_list[ii]
    ii_str = ('0000' + str(ii))[-5:]
    out_fn = temp_dir / ('GZ_' + ii_str + '.nc')
    # use subprocesses
    cmd_list = ['python3', 'get_one_grid_layer.py',
            '-lt','daily',
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
        sys.stdout.flush()
    if (np.mod(ii,50) == 0) and (ii > 0):
        print('') # line feed
        sys.stdout.flush()
    if (ii == N-1):
        print(str(ii))
        sys.stdout.flush()

print('Total processing time = %0.2f sec' % (time()-tt0))

# concatenate the records into one file
# This bit of code is a nice example of how to replicate a bash pipe
pp1 = Po(['ls', str(temp_dir)], stdout=Pi)
pp2 = Po(['grep','GZ'], stdin=pp1.stdout, stdout=Pi)
fn_p = 'GZinfo_'+str(Ldir['ds0'])+'_'+str(Ldir['ds1']+'.nc')
temp_fn = str(temp_dir)+'/'+fn_p # this is all the maps put to one
cmd_list = ['ncrcat','-p', str(temp_dir), '-O', temp_fn]
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
        
"""
Next we want to repackage these results into one NetCDF file per section, with all times.

Variables in the NetCDF files:
(x,y)
- DA is the area of each cell  
- h is bathymetric depth 
- Lat and Lon on rho points 

(t,x,y)
- zr_bot bottom depth (with SSH = 0, for depth calcs)
- zeta surface (SSH ~= 0)

(t)
- ocean_time is a vector of time in seconds since (typically) 1/1/1970.    
"""

tt0 = time()
# add the variables that dont change w/ time and save
fn = fn_list[0]
dsg = xr.open_dataset(fn, decode_times=False)
G, S, T = zrfun.get_basic_info(fn)
DA = G['DX'] * G['DY']
#open dataset and save with time varying variables extracted in get_one_grid_layer.py
ds1 = xr.open_dataset(temp_fn)

ds1['lon_rho'] = dsg.lon_rho
ds1['lat_rho'] = dsg.lat_rho
#ds1['lon_psi'] = dsg.lon_psi
#ds1['lat_psi'] = dsg.lat_psi
ds1['mask_rho'] = dsg.mask_rho
ds1['h'] = dsg.h
ds1['DA'] = (('eta_rho', 'xi_rho'), DA, {'units':'m^2', 'long_name': 'cell horizontal area '})

this_fn = out_dir / (fn_p)
ds1.to_netcdf(this_fn)

print('Time to add grid fields = %0.2f sec' % (time()-tt0))


