"""
Function to do the calculation for one hypoxic volume for a single history file.
This does the whole domain. 

"""

import numpy as np
import sys

from argparse import ArgumentParser
from lo_tools import Lfun, zfun, zrfun
import xarray as xr
from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
from time import time
from lo_tools import zrfun

parser = ArgumentParser()
parser.add_argument('-in_fn', type=str)             # path to history file
parser.add_argument('-out_fn', type=str)            # path to outfile (temp directory)
parser.add_argument('-lt', '--list_type', type=str) # list type: hourly, daily, weekly, lowpass
parser.add_argument('-fb','--false_bottom', default = False, type = Lfun.boolean_string) # places 200m false bottom 
args = parser.parse_args()

tt0 = time()
ds = xr.open_dataset(args.in_fn, decode_times=False)
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat

G, S, T = zrfun.get_basic_info(args.in_fn)   # grid information

# Fields that do not change with time
lon = ds.lon_rho.values                
lat = ds.lat_rho.values
DA = G['DX'] * G['DY']                       # cell horizontal area 
CC = dict()                                  # this is for holding fields extracted on sections
CC['lon'] = lon
CC['lat'] = lat
CC['DA'] = DA
h = ds.h.values       
CC['h'] = h
mask_rho = ds.mask_rho.values
CC['mask_rho']=mask_rho
z_rho, z_w = zrfun.get_z(h, 0*h, S)
dzr = np.diff(z_w, axis=0)
print('Time to get initial fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

tt0 = time()

if args.false_bottom == True:
    oxy = ds.oxygen.values.squeeze()
    oxy[z_rho<-200] = 9999        # put big O2 value so that does not include in depths below 200m in the hyp vol calculation 
elif args.false_bottom == False:
    oxy = ds.oxygen.values.squeeze()
    
# hypoxia thresholds 
dzrm = np.ma.masked_where(oxy>106.6,dzr) 
mild_dz = dzrm.sum(axis=0)

dzrm = np.ma.masked_where(oxy>60.9,dzr) 
hyp_dz = dzrm.sum(axis=0)

dzrm = np.ma.masked_where(oxy>21.6,dzr) 
severe_dz = dzrm.sum(axis=0)

dzrm = np.ma.masked_where(oxy>0,dzr) 
anoxic_dz = dzrm.sum(axis=0)

print('Time to get hypoxia dz = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

CC['mild_dz'] = mild_dz
CC['hyp_dz'] = hyp_dz
CC['severe_dz'] = severe_dz
CC['anoxic_dz'] = anoxic_dz

# put them in a dataset, ds1
NR, NC = CC['hyp_dz'].shape        
ot = ds.ocean_time.values          # an array with dtype='datetime64[ns]'

ds1 = Dataset()
ds1['ocean_time'] = (('ocean_time'), ot, {'calendar':ds.ocean_time.calendar,'units':ds.ocean_time.units})

ds1['mild_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['mild_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of mild hypoxic layer'})
ds1['hyp_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['hyp_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of hypoxic layer'})
ds1['severe_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['severe_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of severe hypoxic layer'})
ds1['anoxic_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['anoxic_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of anoxic layer'})

ds1.to_netcdf(args.out_fn, unlimited_dims='ocean_time')

