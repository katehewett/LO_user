"""
Function to do the calculation for one corrosive volume (ARAG<0.5, 1, 1.7) for a single 
file created using extract_corrosive_box.py
This does the whole domain. Need to add portion that allows a lat/lon box 
"""

import xarray as xr
import numpy as np
import gsw
import sys
import PyCO2SYS as pyco2

from argparse import ArgumentParser
from lo_tools import Lfun, zfun, zrfun

from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
from time import time

parser2 = ArgumentParser()
parser2.add_argument('-in_fn', '--input_box_files', type=str) # path to clipped box files (temp directory)
parser2.add_argument('-out_fn', '--output_vol_files', type=str) # path to save cline files (temp directory)
parser2.add_argument('-his_fn', '--input_his_file', type=str) # path to one history file for grid info > z's
parser2.add_argument('-lt', '--list_type', type=str) # list type: hourly, daily, weekly, lowpass     
args2 = parser2.parse_args()

testing = False # True => process fewer layers for ARAG

# 1. Get some grid information
tt0 = time()
ds = xr.open_dataset(args2.input_box_files, decode_times=False)
h = ds.h.values
mask_rho = ds.mask_rho.values

G, S, T = zrfun.get_basic_info(args2.input_his_file)
#z_rho, z_w = zrfun.get_z(h, zeta, S)  # we saved zeta 
z_rho, z_w = zrfun.get_z(h, 0*h, S) # use 0 for SSH
dzr = np.diff(z_w, axis=0) # vertical thickness of all cell [m]

lon = np.nanmean(G['lon_rho']) # it's just for the GSW calcs and faster than passing args and clipping 
lat = np.nanmean(G['lat_rho'])

print('Time to get initial fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# 2. Prepare for carbon variable calculations
tt0 = time()
p = gsw.p_from_z(z_rho, lat) # pressure [dbar]
SP = ds.salt.values.squeeze() # practical salinity
# Note: SP and other fields loaded this way will be packed z,y,x
# and have nan's where masked. The singleton time dimension is removed
# by the squeeze() call.
# Make sure salinity is not negative. Could be a problem for pyco2.
SP[SP<0] = 0
PT = ds.temp.values.squeeze() # potential temperature [degC]
ALK = ds.alkalinity.values.squeeze() # alkalinity [milli equivalents m-3 = micro equivalents L-1]
TIC = ds.TIC.values.squeeze() # # TIC [millimol C m-3 = micromol C L-1]
# We then use gsw routines to calculate in situ density and temperature.
SA = gsw.SA_from_SP(SP, p, lon, lat) # absolute salinity [g kg-1]
CT = gsw.CT_from_pt(SA, PT) # conservative temperature [degC]
rho = gsw.rho(SA, CT, p) # in situ density [kg m-3]
ti = gsw.t_from_CT(SA, CT, p) # in situ temperature [degC]
# Convert from micromol/L to micromol/kg using in situ dentity because these are the
# units expected by pyco2.
ALK1 = 1000 * ALK / rho
TIC1 = 1000 * TIC / rho
# I'm not sure if this is needed. In the past a few small values of these variables had
# caused big slowdowns in the MATLAB version of CO2SYS.
ALK1[ALK1 < 100] = 100
TIC1[TIC1 < 100] = 100
print('Time to load and prepare carbon fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# 3. Calculate aragonite saturation and a map of corrosive thickness.
# We do this one layer at a time because pyco2 got bogged down if we fed it
# the whole 3-D arrays at once.
tt0 = time()
ARAG = np.nan * np.ones(SP.shape) # Initialize array (z,y,z) to hold results.
nz, nr, nc = SP.shape # handy dimension sizes
amat = np.nan * np.ones((nr,nc)) # Initialize array (y,x) for single layer.
if testing:
    nlay = 1 # number of layers to process
else:
    nlay = nz
for ii in range(nlay):
    # ii = nz-1 # hand override to test surface layer
    tt00 = time()
    # Note that by using the [mask_rho==1] operator on each of the layers
    # we go from a 2-D array with nan's to a 1-D vector with no nan's. We
    # do this to speed up the calculation. Note that the cas6 grid is about half
    # landmask.
    aALK = ALK1[ii,:,:].squeeze()[mask_rho==1]
    aTIC = TIC1[ii,:,:].squeeze()[mask_rho==1]
    aTemp = ti[ii,:,:].squeeze()[mask_rho==1]
    aPres = p[ii,:,:].squeeze()[mask_rho==1]
    aSalt = SP[ii,:,:].squeeze()[mask_rho==1]
    # Note: here is where to get info on the inputs, outputs, and units:
    # https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/
    CO2dict = pyco2.sys(par1=aALK, par1_type=1, par2=aTIC, par2_type=2,
        salinity=aSalt, temperature=aTemp, pressure=aPres,
        total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)
    aARAG = CO2dict['saturation_aragonite']
    if False:
        # Check by comparing to the old way Parker was calling the function.
        from PyCO2SYS import CO2SYS
        CO2dict_alt = CO2SYS(aALK, aTIC, 1, 2, aSalt, aTemp, aTemp,
            aPres, aPres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
        # Accounting for the defaults, this is the same as the pyco2.sys() call above,
        # except for the opt_buffers_mode, which I do not think we need.
        aARAG_alt = CO2dict_alt['OmegaARout']
        # Then look at aARAG - aARAG_alt by hand.
        # RESULT: the results for the deepest layer are only different by -2e-6
        # for a field with mean = .55, so basically identical. The difference is
        # considerably smaller for the surface layer.
    # Then write the arag field into the ARAG array, indexing with [mask_rho==1].
    aamat = amat.copy()
    aamat[mask_rho==1] = aARAG
    ARAG[ii,:,:] = aamat 
    #print('  ii = %d' % (ii))
    #print('  Time to get one slice = %0.2f sec' % (time()-tt00))
    sys.stdout.flush()
print('Time to calculate ARAG for all layers = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# 4. Form the map of the thickness of corrosive water [m].
tt0 = time()
CC = dict()

dzrm = dzr.copy()
dzrm[ARAG>0.5] = 0
corrosive_severe_dz = dzrm.sum(axis=0)
CC['corrosive_severe_dz'] = corrosive_severe_dz

dzrm1 = dzr.copy()
dzrm1[ARAG>1] = 0
corrosive_int_dz = dzrm1.sum(axis=0)
CC['corrosive_int_dz'] = corrosive_int_dz

dzrm2 = dzr.copy()
dzrm2[ARAG>1.7] = 0
corrosive_mild_dz = dzrm2.sum(axis=0)
CC['corrosive_mild_dz'] =corrosive_mild_dz

# put them in a dataset, ds1
NR, NC = CC['corrosive_int_dz'].shape        
ot = ds.ocean_time.values          # an array with dtype='datetime64[ns]'

ds1 = Dataset()
ds1['ocean_time'] = (('ocean_time'), ot, {'calendar':ds.ocean_time.calendar,'units':ds.ocean_time.units})

ds1['corrosive_severe_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['corrosive_severe_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of layer Oag<0.5'})
ds1['corrosive_int_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['corrosive_int_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of layer Oag<1'})
ds1['corrosive_mild_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['corrosive_mild_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of layer Oag<1.7'})

ds1.to_netcdf(args2.output_vol_files, unlimited_dims='ocean_time')

print('Time to get corrosive_dz and save = %0.2f sec' % (time()-tt0))
sys.stdout.flush()
