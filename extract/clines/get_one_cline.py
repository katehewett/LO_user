"""
Function to do the calculation for a single clipped file
calculates GSW vars > drho/dz dT/dz
finds max drho(dT)/dz + returns the max depth + values

"""

import numpy as np
import sys

from argparse import ArgumentParser
from lo_tools import zfun, zrfun
from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
from time import time
from lo_tools import zrfun
import gsw 

parser2 = ArgumentParser()
parser2.add_argument('-in_fn', '--input_box_files', type=str) # path to clipped box files (temp directory)
parser2.add_argument('-out_fn', '--output_cline_files', type=str) # path to save cline files (temp directory)
parser2.add_argument('-his_fn', '--input_his_file', type=str) # path to one history file for grid info > z's
parser2.add_argument('-lt', '--list_type', type=str) # list type: hourly, daily, weekly, lowpass         
args2 = parser2.parse_args()

# 1. Get some grid information 
tt0 = time()
dsb = open_dataset(args2.input_box_files, decode_times=False)
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat

# add z variables
# NOTE: this only works if you have salt as a saved field!
G, S, T = zrfun.get_basic_info(args2.input_his_file)   # grid information
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
    
g = 9.8     
po = 1025 
dpdz = dp/dz
N2 = -(g/po)*dpdz     # radians / second
dTdz = dp/dz

C = np.argmax(N2,axis=1,keepdims=True)
z_N2max = np.take_along_axis(z_w.values,C+1,axis=1)
val_N2max = np.take_along_axis(N2,C,axis=1)

B = np.argmax(dTdz,axis=1,keepdims=True)
z_dTdzmax = np.take_along_axis(z_w.values,B+1,axis=1)
val_dTdzmax = np.take_along_axis(dTdz,B,axis=1)

# put to 
#ds1 = Dataset()

#ds1.to_netcdf(args2.output_cline_files, unlimited_dims='ocean_time')





