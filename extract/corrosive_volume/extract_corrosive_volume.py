"""
Code to extract corrosive volumes from LO domain.

Was taking a long time to run on perigree and then run would terminate after processing 100-ish history files. 
Ran hypoxic volume calcs seperate while figuring out Oag calcs... 
Need to combine hypoxic and corrosive vol calcs

2 history files ran in 43 seconds (~36 seconds to do whole LO domain all layers loop in get_one_corrosive_volume.py)

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
hv_dir = Ldir['LOo'] / 'extract' / 'hypoxic_volume'

fn_list = Lfun.get_fn_list('daily', Ldir, Ldir['ds0'], Ldir['ds1'])

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'corrosive_volume'
out_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
temp_dir = out_dir0 / ('temp_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)
Lfun.make_dir(temp_dir, clean=True)

# still working on this if section
#if Ldir['testing']:
#    fn_list = fn_list[:3]

# loop over all jobs
tt0 = time()
N = len(fn_list)
proc_list = []
for ii in range(N):
    # Launch a job and add its process to a list.
    fn = fn_list[ii]
    ii_str = ('0000' + str(ii))[-5:]
    out_fn = temp_dir / ('CC_' + ii_str + '.nc')
    # use subprocesses
    cmd_list = ['python3', 'get_one_corrosive_volume.py',
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
pp2 = Po(['grep','CC'], stdin=pp1.stdout, stdout=Pi)
fn_p = 'Volumes_Oag_'+str(Ldir['ds0'])+'_'+str(Ldir['ds1']+'.nc')
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
- corrosive_dz is depth of the the undersaturated layer (t, x, y) in each cell 
    
Need to save a seperate file with the variables that don't change over time - this is to save space
- DA is the area of each cell (x,y) < doesn't DA stay the same thru ocean_time? 
- h is bathymetric depth 
- ocean_time is a vector of time in seconds since (typically) 1/1/1970.
- Lat and Lon on rho points 
    
A useful tool is isel():
a = ds.isel(p=np.arange(10,15))
"""

ds1 = xr.open_dataset(temp_fn)
this_fn = out_dir / (fn_p)
ds1.to_netcdf(this_fn)



