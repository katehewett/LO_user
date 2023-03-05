"""
Function to do the calculation for one layer T S O2 TIC Alkalinity 
and then calculate and pass back ARAG and pH total for a single history file.
This does the whole domain. Need to add portion that allows a lat/lon box 
encorporate w/ box  ?

on mac ~4-5seconds for 2 history files whole domain 
Each history file runs ~12-15 seconds on perigee

"""
import xarray as xr
import numpy as np
import gsw
import sys
import PyCO2SYS as pyco2

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

testing = False # True => process fewer layers for ARAG

vn_list = ['salt', 'temp', 'oxygen', 'TIC', 'alkalinity']
            
# 1. Get some grid information
tt0 = time()
ds = xr.open_dataset(args.in_fn, decode_times=False)
G, S, T = zrfun.get_basic_info(args.in_fn)
lon = ds.lon_rho.values
lat = ds.lat_rho.values
DA = G['DX'] * G['DY']
# we pack some things in the CC dict that will eventually be written out to a
# NetCDF file.
CC = dict()
CC['lon'] = lon
CC['lat'] = lat
CC['DA'] = DA
h = ds.h.values
CC['h'] = h
mask_rho = ds.mask_rho.values
CC['mask_rho']=mask_rho
z_rho, z_w = zrfun.get_z(h, 0*h, S) # use 0 for SSH
z_rho = z_rho[0,:,:].squeeze()      # just the bottom layer 
#dzr = np.diff(z_w, axis=0) # vertical thickness of all celle [m]
print('Time to get initial fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

tt0 = time()
# 2. Tracers and zeta 
for vn in vn_list:
    aa = ds[vn].values[:,0,:,:].squeeze() # just the bottom layer 
    CC[vn] = aa
aa = ds.zeta.values.squeeze()
CC['zeta'] = aa
print('Time to grab tracers and zeta = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# 2. Prepare for carbon variable calculations
tt0 = time()
p = gsw.p_from_z(z_rho, lat) # pressure [dbar]
SP = CC['salt'].squeeze() # ds.salt.values.squeeze() # practical salinity
# Note: SP and other fields loaded this way will be packed z,y,x
# and have nan's where masked. The singleton time dimension is removed
# by the squeeze() call.
# Make sure salinity is not negative. Could be a problem for pyco2.
SP[SP<0] = 0
PT = CC['temp'].squeeze() # ds.temp.values.squeeze() # potential temperature [degC]
ALK = CC['alkalinity'].squeeze() #ds.alkalinity.values.squeeze() # alkalinity [milli equivalents m-3 = micro equivalents L-1]
TIC = CC['TIC'].squeeze() #ds.TIC.values.squeeze() # # TIC [millimol C m-3 = micromol C L-1]
# We then use gsw routines to calculate in situ density and temperature.
SA = gsw.SA_from_SP(SP, p, lon, lat) # absolute salinity [g kg-1]
CC['SA'] = SA
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

# 3. Calculate aragonite saturation and pH a map of bottom values.
tt0 = time()
ARAG = np.nan * np.ones(SP.shape) # Initialize array (y,z) to hold results.
pH = np.nan * np.ones(SP.shape)
nr, nc = SP.shape                 # handy dimension sizes
amat = np.nan * np.ones((nr,nc))  # Initialize array (y,x) for single layer.
amat2 = np.nan * np.ones((nr,nc)) 

# Note that by using the [mask_rho==1] operator on each of the layers
# we go from a 2-D array with nan's to a 1-D vector with no nan's. We
# do this to speed up the calculation. Note that the cas6 grid is about half
# landmask.
aALK = ALK1.squeeze()[mask_rho==1]
aTIC = TIC1.squeeze()[mask_rho==1]
aTemp = ti.squeeze()[mask_rho==1]
aPres = p.squeeze()[mask_rho==1]
aSalt = SP.squeeze()[mask_rho==1]
# Note: here is where to get info on the inputs, outputs, and units:
# https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/
CO2dict = pyco2.sys(par1=aALK, par1_type=1, par2=aTIC, par2_type=2,
salinity=aSalt, temperature=aTemp, pressure=aPres,
total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)
aARAG = CO2dict['saturation_aragonite']
apH = CO2dict['pH_total']

# Then write the arag and pH fields into an array, indexing with [mask_rho==1].
aamat = amat.copy()
aamat2 = amat2.copy()
aamat[mask_rho==1] = aARAG
aamat2[mask_rho==1] = apH
ARAG = aamat 
pH = aamat2
#print('  ii = %d' % (ii))
#print('  Time to get one slice = %0.2f sec' % (time()-tt00))
sys.stdout.flush()
print('Time to calculate ARAG for all layers = %0.2f sec' % (time()-tt0))
#print('saving:', args.in_fn)
sys.stdout.flush()

CC['pH'] = pH
CC['ARAG'] = ARAG

# 4. put them in a dataset, ds1
NR, NC = CC['salt'].shape        
ot = ds.ocean_time.values          # an array with dtype='datetime64[ns]'

ds1 = Dataset()
ds1['ocean_time'] = (('ocean_time'), ot, {'long_name':ds.ocean_time.long_name,'units':ds.ocean_time.units})

ds1['salt'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['salt'].reshape(1,NR,NC), 
{'units':'none', 'long_name': 'practical salinity'})

ds1['SA'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['SA'].reshape(1,NR,NC), 
{'units':'g kg^-1', 'long_name': 'absolute salinity'})

ds1['temp'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['temp'].reshape(1,NR,NC), 
{'units':ds.temp.units, 'long_name': ds.temp.long_name})

ds1['TIC'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['TIC'].reshape(1,NR,NC), 
{'units':ds.TIC.units, 'long_name': ds.TIC.long_name})

ds1['Alkalinity'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['alkalinity'].reshape(1,NR,NC), 
{'units':ds.alkalinity.units, 'long_name': ds.alkalinity.long_name})

ds1['oxygen'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['oxygen'].reshape(1,NR,NC), 
{'units':ds.oxygen.units, 'long_name': ds.oxygen.long_name})

ds1['pH'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['pH'].reshape(1,NR,NC), 
{'units':'total scale', 'long_name':'pH total scale'})

ds1['ARAG'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['ARAG'].reshape(1,NR,NC), 
{'long_name':'aragonite saturation state'})

ds1['DA'] = (('eta_rho', 'xi_rho'), CC['DA'], {'units':'m^2', 'long_name': 'cell horizontal area '})

ds1['mask_rho'] = (('eta_rho', 'xi_rho'), CC['DA'], 
{'flag_values':[0., 1.],'flag_meanings':'land water','long_name': 'mask on RHO-points'})

ds1['h'] = (('eta_rho', 'xi_rho'), CC['h'], {'units':ds.h.units, 'long_name': ds.h.long_name})

ds1['Lat'] = (('eta_rho', 'xi_rho'), CC['lat'], {'units':'degree_north','long_name': 'latitude of RHO-points'})
ds1['Lon'] = (('eta_rho', 'xi_rho'), CC['lon'], {'units':'degree_east','long_name': 'longitude of RHO-points'})

ds1.to_netcdf(args.out_fn, unlimited_dims='ocean_time')

