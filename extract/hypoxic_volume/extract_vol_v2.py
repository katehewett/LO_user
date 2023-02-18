"""
Code to extract hypoxic volumes from LO domain.

Baesd code on extract/box and /tef2 code: extract_sections; extract_one_section
updated code to make the volume calculations run faster 

To test on mac:
run extract_vol_v2 -gtx cas6_v0_live -ro 1 -0 2022.08.08 -1 2022.08.09 -lt daily 

Testing January 2023 - present
"""
# 
# % clear 
# % reset
#import os 
#os.chdir('/Users/katehewett/Documents/LO_user/extract/hypoxic_volume')
#run extract_vol_v2 -gtx cas6_v0_live -ro 1 -0 2022.08.08 -1 2022.08.09 -lt daily 
#-job LO_oxygen_WA 

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

# sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p') # buncha pickled sections // still working on this
# sect_df = pd.read_pickle(sect_df_fn)

fn_list = Lfun.get_fn_list('daily', Ldir, Ldir['ds0'], Ldir['ds1'])

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'hypoxic_volume'
out_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
temp_dir = out_dir0 / ('temp_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)
Lfun.make_dir(temp_dir, clean=True)

# still working on this if section
if Ldir['testing']:
    fn_list = fn_list[:3]

#if Ldir['get_bio']:
#    vn_type = 'bio'
#else:
#    vn_type = 'salt'

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
    cmd_list = ['python3', 'get_one_volume.py',
            #'-hv_df_fn', str(hv_df_fn),
            #'#-job','LO_oxygen_WA'
            '-in_fn',str(fn),
            '-out_fn', str(out_fn)]
            #'-vn_type', vn_type]
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
temp_fn = str(temp_dir)+'/all.nc'
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

We will follow the structure of the output of LO/tef/extract_sections.py so that we can mostly
recycle the subsequent processing code:

Variables in the NetCDF files:
- salt is daily salinity in each cell (t, z, x-or-y) [same for all other variables]
- q is hourly transport in each cell (t, z, x-or-y)
- vel is velocity in each cell (t, z, x-or-y) positive to East or North
- DA is the area of each cell (t, z, x-or-y) hence: q = vel * DA
- z0 is the average z-position of cell centers (assumes SSH=0), useful for plotting
- DA0 is the average cross-sectional area of each cell (assumes SSH=0)
- h is depth on the section (x-or-y) positive down
- zeta is SSH on the section (t, x-or-y) positive up
- ocean_time is a vector of time in seconds since (typically) 1/1/1970.
    
A useful tool is isel():
a = ds.isel(p=np.arange(10,15))
"""

ds1 = xr.open_dataset(temp_fn)
S = zrfun.get_basic_info(fn_list[0], only_S=True)
eta = ds1.zeta.values.squeeze() # packed (t, p)
NT, NP = eta.shape
hh = ds1.h.values.squeeze().reshape(1,NP) * np.ones((NT,1))
zw = zrfun.get_z(hh, eta, S, only_w=True)
dz = np.diff(zw, axis=0) # NOTE: this is packed (z,t,p)
DZ = np.transpose(dz, (1,0,2)) # packed (t,z,p)

## Help on this too b/c not sure how to pack back in -- confusion continued from lines 35ish this coe
#hv_list = list(hv_df.sn.unique())
#hv_list.sort()
#for sn in hv_list:
#    ii = np.where(hv_df.sn == sn)[0]
#    this_ds = ds1.isel(p=ii)
#    # add DZ
#    this_DZ = DZ[:,:,ii]
#    this_ds['DZ'] = (('time','z','p'), this_DZ)
#    this_fn = out_dir / (sn + '.nc')
#    this_ds.to_netcdf(this_fn)



