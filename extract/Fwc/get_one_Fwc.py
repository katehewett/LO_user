"""
Function to do the calculation for a single file sent
by extract_Fwc.py 
This script: 
calculates (1) Rig (2) Fwc using a threshold method 

"""

import numpy as np
import sys

from argparse import ArgumentParser
from lo_tools import zrfun
from xarray import open_dataset, Dataset
from time import time
from lo_tools import zrfun

import gsw 

parser2 = ArgumentParser()
parser2.add_argument('-in_fn', '--input_box_files', type=str) # path to clipped box files (temp directory)
parser2.add_argument('-out_fn', '--output_Fwc_files', type=str) # path to save Fwc files (temp directory)
parser2.add_argument('-his_fn', '--input_his_file', type=str) # path to one history file for grid info > z's
parser2.add_argument('-lt', '--list_type', type=str) # list type: hourly, daily, weekly, lowpass      
args2 = parser2.parse_args()

#############################make this an argument
threshold = 32.5
#############################

# 1. Get some grid information
tt0 = time()
dsb = open_dataset(args2.input_box_files, decode_times=False) # TODO UPDATE
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat

# add z variables
# NOTE: this only works if you have salt as a saved field!
G, S, T = zrfun.get_basic_info(args2.input_his_file)   # grid information # TODO UPDATE

lat = np.nanmean(G['lat_rho'])
lon = np.nanmean(G['lon_rho']) # for gsw calcs

NT, NZ, NR, NC = np.shape(dsb.salt)
h = dsb.h.values
zeta = dsb.zeta.values
mask_rho = dsb.mask_rho

z_rho, z_w = zrfun.get_z(h, zeta, S) # (30, 1111, 356) (31, 1111, 356)
dzr = np.diff(z_w, axis=0) # vertical thickness of all cell [m]

z_rho = np.expand_dims(z_rho,axis=0) # a little unnecessary step but helps make it easier below 
z_w = np.expand_dims(z_w,axis=0)

tempC = dsb.temp.values # (1, 30, 1111, 356)
SP = dsb.salt.values

print('Time to get initial fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# 2. Calc SA and SIG0: there's a holdup with sending gsw LARGE files - 
# it can accomodate 4D variables, but has an upper limit 
# and can be slow. Not a problem with NT = 1 tho if change to chunks 2xcheck
tt0 = time()

P = gsw.p_from_z(z_rho,lat)
SA = gsw.SA_from_SP(SP, P, lon, lat)
CT = gsw.CT_from_pt(SA, tempC)
SIG0 = gsw.sigma0(SA,CT)   

print('Time to apply GSW = %0.2f sec' % (time()-tt0))
sys.stdout.flush()
     
# 3. Form the map of the thickness of Freshwater and Sigma26.5
tt0 = time()
CC = dict()

dzrm = dzr.copy()
dzrm[SP>32.5] = 0
Fwc_dz = dzrm.sum(axis=0)
CC['Fwc_dz'] = Fwc_dz

dzrm1 = dzr.copy()
dzrm1[SIG0<26.5] = 0
SIG26h_dz = dzrm1.sum(axis=0)
CC['SIG26h_dz'] = SIG26h_dz

print('Time to calc Fwc and sig 26.5 = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

#save some values 
ds1 = Dataset()
ds1['ocean_time'] = dsb.ocean_time

ds1['zeta'] = (('ocean_time', 'eta_rho', 'xi_rho'), dsb.zeta.values, {'units':'m', 'long_name':'free-surface'})
ds1['Fwc_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['Fwc_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of layer SA<=32.5'})
ds1['SIG26h_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['SIG26h_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of layer SIG0>=26.5'})

ds1.to_netcdf(args2.output_cline_files, unlimited_dims='ocean_time')

print('saved!')