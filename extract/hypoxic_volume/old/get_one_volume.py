"""
Function to do the calculation for one hypoxic volume for a single history file.
This does the whole domain. Need to add portion that allows a lat/lon box 
"""

import numpy as np
import gsw
import sys
import PyCO2SYS as pyco2

from argparse import ArgumentParser
from lo_tools import zfun, zrfun
from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
from time import time
from lo_tools import zrfun

parser = ArgumentParser()
parser.add_argument('-in_fn', type=str) # path to history file
parser.add_argument('-out_fn', type=str) # path to outfile (temp directory)
parser.add_argument('-lt', '--list_type', default = 'daily', type=str) # list type: hourly, daily, weekly
args = parser.parse_args()

tt0 = time()
ds = open_dataset(args.in_fn, decode_times=False)
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
# hypoxia thresholds 
oxy = ds.oxygen.values.squeeze()
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

tt0 = time()
# carbon claculations
p = gsw.p_from_z(z_rho, lat)           # pressure [dbar]; adjusted in zrfun.get_z(h, 0*h, S)
SP = ds.salt.values.squeeze()
PT = ds.temp.values.squeeze()
ALK = ds.alkalinity.values.squeeze()
TIC = ds.TIC.values.squeeze()

SA = gsw.SA_from_SP(SP, p, lon, lat)
CT = gsw.CT_from_pt(SA, PT)
rho = gsw.rho(SA, CT, p)              # in situ density
ti = gsw.t_from_CT(SA, CT, p)         # in situ temperature

# convert from umol/L to umol/kg using in situ dentity
ALK1 = 1000 * ALK / rho
TIC1 = 1000 * TIC / rho

## I'm not sure if this is needed
#ALK1[ALK1 < 100] = np.nan                 # Q from dm_pfun.py: why? 
#TIC1[TIC1 < 100] = np.nan                 # Q from dm_pfun.py: why? 

# do we need to check for neg salinities / temps

# calculate aragonite saturation:
# For CO2SYS: All temperatures are in °C, 
#             all salinities are on the PSS, 
#             and all pressures are in dbar. 
ARAG = np.nan * np.ones(SP.shape)
nz, nr, nc = SP.shape
amat = np.nan * np.ones((nr,nc))
for ii in range(nz):
# for ii in range(2):
    tt00 = time()
    aALK = ALK1[ii,:,:].squeeze()[mask_rho==1]
    aTIC = TIC1[ii,:,:].squeeze()[mask_rho==1]
    aTemp = ti[ii,:,:].squeeze()[mask_rho==1]
    aPres = p[ii,:,:].squeeze()[mask_rho==1]
    aSalt = SP[ii,:,:].squeeze()[mask_rho==1]
    # note: still need to consult:
    # https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/
    # to make sure the other inputs are handled correctly (e.g. what does 50 mean?)
    CO2dict = pyco2.sys(par1=aALK, par1_type=1, par2=aTIC, par2_type=2,
        salinity=aSalt, temperature=aTemp, pressure=aPres, opt_buffers_mode=0)
    # CO2dict = CO2SYS(aALK, aTIC, 1, 2, aSalt, aTemp, aTemp,
    #     aPres, aPres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)             # assumptions from dm_pfun.py
    aARAG = CO2dict['saturation_aragonite']
    aamat = amat.copy()
    aamat[mask_rho==1] = aARAG
    ARAG[ii,:,:] = aamat 
    print('  ii = %d' % (ii))
    print('  Time to get one slice = %0.2f sec' % (time()-tt00))
    sys.stdout.flush()

dzrm = dzr.copy()
dzrm[ARAG>1] = 0
corrosive_dz = dzrm.sum(axis=0)
CC['corrosive_dz'] = corrosive_dz
print('Time to get corrosive_dz = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# put them in a dataset, ds1
NR, NC = CC['hyp_dz'].shape        
ot = ds.ocean_time.values          # an array with dtype='datetime64[ns]'

ds1 = Dataset()
ds1['ocean_time'] = (('ocean_time'), ot, {'long_name':ds.ocean_time.long_name,'units':ds.ocean_time.units})

ds1['mild_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['mild_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of mild hypoxic layer'})
ds1['hyp_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['hyp_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of hypoxic layer'})
ds1['severe_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['severe_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of severe hypoxic layer'})
ds1['anoxic_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['anoxic_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of anoxic layer'})

ds1['corrosive_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['corrosive_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of undersaturated layer'})

#ds1['DA'] = (('eta_rho', 'xi_rho'), CC['DA'], {'units':'m^2', 'long_name': 'cell horizontal area '})
#ds1['mask_rho'] = (('eta_rho', 'xi_rho'), CC['DA'], {'flag_values':[0., 1.],'flag_meanings':'land water','long_name': 'mask on RHO-points'})
#ds1['h'] = (('eta_rho', 'xi_rho'), CC['h'], {'units':ds.h.units, 'long_name': ds.h.long_name})

#ds1['Lat'] = (('eta_rho', 'xi_rho'), CC['lat'], {'units':'degree_north','long_name': 'latitude of RHO-points'})
#ds1['Lon'] = (('eta_rho', 'xi_rho'), CC['lon'], {'units':'degree_east','long_name': 'longitude of RHO-points'})

ds1.to_netcdf(args.out_fn, unlimited_dims='ocean_time')

