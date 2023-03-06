"""
Function to do the calculation to grab grid info and zeta z_rho's for each history file 
to be used in plotting - and needs to get cleaned up!
"""
import xarray as xr
import numpy as np
import sys


from argparse import ArgumentParser
from lo_tools import zfun, zrfun

from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
from time import time

parser = ArgumentParser()
parser.add_argument('-in_fn', type=str) # path to history file
parser.add_argument('-out_fn', type=str) # path to outfile (temp directory)
parser.add_argument('-lt', '--list_type', default = 'daily', type=str) # list type: hourly, daily, weekly
args = parser.parse_args()
  
tt0 = time()            
# 1. Get some grid information
ds = xr.open_dataset(args.in_fn, decode_times=False)
G, S, T = zrfun.get_basic_info(args.in_fn)
h = ds.h.values
mask_rho = ds.mask_rho.values # 0 = land
z_rho, z_w = zrfun.get_z(h, 0*h, S) # use 0 for SSH
zbot_rho = z_rho[0,:,:].squeeze()      # just the bottom layer 
zbot_rho2 = np.ma.masked_where(mask_rho==0,zbot_rho) # mask land
zbot_rho2 = zfun.fillit(zbot_rho2) # nan's where land

dzr = np.diff(z_w, axis=0) # vertical thickness of all cells [m]
dzr = np.sum(dzr,axis=0)   # the total dzr 
dzrm = np.ma.masked_where(mask_rho==0,dzr) # mask land
dzrm = zfun.fillit(dzrm) # nan's where land

# 2. put them in a dataset, ds1
NT, NR, NC = np.shape(ds.zeta) 
ot = ds.ocean_time.values          # an array with dtype='datetime64[ns]'

ds1 = Dataset()
ds1['ocean_time'] = (('ocean_time'), ot, {'long_name':ds.ocean_time.long_name,'units':ds.ocean_time.units})
ds1['zbot_rho'] = (('ocean_time','eta_rho', 'xi_rho'), zbot_rho2.reshape(NT,NR,NC), {'units':'m', 'long_name': 'depth of bottom with SSH=0'})
ds1['dzr'] = (('ocean_time','eta_rho', 'xi_rho'), dzrm.reshape(NT,NR,NC), {'units':'m', 'long_name': 'vertical thickness of water on rho_point'})
ds1['zeta'] = ds.zeta

ds1.to_netcdf(args.out_fn, unlimited_dims='ocean_time')

print('Time to get time varying fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()
