"""
Function to do the calculation for one hypoxic volume for a single history file.
"""

from argparse import ArgumentParser
from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
from pickle import dump
from pandas import read_pickle, to_datetime
from lo_tools import zrfun
from time import time

import numpy as np

parser = ArgumentParser()
#parser.add_argument('-sect_df_fn', type=str) # path to sect_df
parser.add_argument('-in_fn', type=str) # path to history file
parser.add_argument('-out_fn', type=str) # path to outfile (temp directory)
#parser.add_argument('-vn_type', type=str) # 'salt' or 'bio'
args = parser.parse_args()

#hv_df = read_pickle(args.hv_df_fn)

ds = open_dataset(args.in_fn, decode_times=False)
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat

vn_list = ['oxygen']
# aa = [-126, -122.5, 46, 49] # WA Shelf, future - want to pass this in the arg passer job_defn?

G, S, T = zrfun.get_basic_info(args.in_fn)    # grid information

lon_rho = ds.lon_rho.values                
lat_rho = ds.lat_rho.values

DX = 1/ds.pm.values
DY = 1/ds.pn.values

DA = DX * DY          # cell horizontal area 

CC = dict()           # this is for holding fields extracted on sections

CC['lon_rho'] = lon_rho
CC['lat_rho'] = lat_rho

CC['DA'] = DA

# Fields that do not change with time
# 
h = ds.h.values       
CC['h'] = h

mask_rho = ds.mask_rho.values
CC['mask_rho']=mask_rho

# Fields that do change with time
#
# First: tracers and zeta
#for vn in vn_list:
#    aa = ds[vn].values.squeeze()
#    CC[vn] = aa; # whole shebang // not cropped to job yet - improve once get working

zeta = ds.zeta.values.squeeze()
z_rho, z_w = zrfun.get_z(h, zeta, S)
#CC['zeta_rho'] = z_rho
#CC['zeta_w'] = z_w

dzr = np.diff(z_w, axis = 0)
oxy = ds.oxygen.values.squeeze()
dzrm = np.ma.masked_where(oxy>61,dzr) 
hyp_dz = dzrm.sum(axis=0)

CC['hyp_dz'] = hyp_dz

#dump(CC,open(args.out_fn),'wb')

# put them in a dataset, ds1
NR, NC = CC['hyp_dz'].shape
ot = ds.ocean_time.values # an array with dtype='datetime64[ns]'
#attrs = {'units':ds.ocean_time.units} # ds time object has no attribute 'units' for me
ds1 = Dataset()
ds1['ocean_time'] = (('ocean_time'), ot, {'long_name':ds.ocean_time.long_name,'units':ds.ocean_time.units})

ds1['hyp_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['hyp_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of hypoxic layer'})

ds1['DA'] = (('eta_rho', 'xi_rho'), CC['DA'], {'units':'m^2', 'long_name': 'cell horizontal area '})
ds1['mask_rho'] = (('eta_rho', 'xi_rho'), CC['DA'], {'flag_values':[0., 1.],'flag_meanings':'land water','long_name': 'mask on RHO-points'})
ds1['h'] = (('eta_rho', 'xi_rho'), CC['h'], {'units':ds.h.units, 'long_name': ds.h.long_name})

ds1['Lat'] = (('eta_rho', 'xi_rho'), CC['lat_rho'], {'units':'degree_north','long_name': 'latitude of RHO-points'})
ds1['Lon'] = (('eta_rho', 'xi_rho'), CC['lon_rho'], {'units':'degree_east','long_name': 'longitude of RHO-points'})

ds1.to_netcdf(args.out_fn, unlimited_dims='ocean_time')

