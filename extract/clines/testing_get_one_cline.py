"""
Function to do the calculation for a single clipped file
calculates GSW vars > drho/dz dT/dz
finds max drho(dT)/dz + returns the max depth + values

THIS IS A TESTING ONE OFF FILE for get_one_cline.py 
update original when working and delete this file

Time to get initial fields = 0.72 sec
Time to apply GSW = 0.82 sec

"""

import numpy as np
import sys

from argparse import ArgumentParser
from lo_tools import zfun, zrfun
from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
from time import time
from lo_tools import Lfun, zrfun, zfun

import gsw 

######## ######## ######## 
# remove this testing after get the pycnocline to work here. and replace with the 
# arg pass tags in get_one_cline.py 
box_fn_in = '/Users/katehewett/Documents/LO_output/extract/cas6_v0_live/clines/shelf_box_2022.08.08_2022.08.09/temp_box_dir/box_000001.nc'
his_in = '/Users/katehewett/Documents/LO_roms/cas6_v0_live/f2022.08.08/ocean_his_0021.nc'
######## ######## ######## 

tt0 = time()
dsb = open_dataset(box_fn_in, decode_times=False) # TODO UPDATE
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat

# add z variables
# NOTE: this only works if you have salt as a saved field!
G, S, T = zrfun.get_basic_info(his_in)   # grid information # TODO UPDATE
NT, NZ, NETA, NXI = np.shape(dsb.salt)
h = dsb.h.values
zeta = dsb.zeta.values
lat = dsb.lat_rho.values
lon = dsb.lon_rho.values
mask_rho = dsb.mask_rho

zeta = dsb.zeta.values
z_rho, z_w = zrfun.get_z(h, zeta, S) # (30, 1111, 356) (31, 1111, 356)

NW = z_w.shape[0]

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
       
dp = np.diff(SIG0,axis=1)
dT = np.diff(CT,axis=1) 
dz = np.diff(z_rho,axis=1)









  