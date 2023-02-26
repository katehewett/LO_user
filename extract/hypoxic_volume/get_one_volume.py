"""
Function to do the calculation for one hypoxic volume for a single history file.
This does the whole domain. Need to add portion that allows a lat/lon box 
"""

from argparse import ArgumentParser
from lo_tools import zfun, zrfun
from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
#from time import time
from lo_tools import zrfun
import numpy as np

parser = ArgumentParser()
parser.add_argument('-in_fn', type=str) # path to history file
parser.add_argument('-out_fn', type=str) # path to outfile (temp directory)
parser.add_argument('-lt', '--list_type', default = 'daily', type=str) # list type: hourly, daily, weekly
args = parser.parse_args()

ds = open_dataset(args.in_fn, decode_times=False)
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat

G, S, T = zrfun.get_basic_info(args.in_fn)    # grid information

lon = ds.lon_rho.values                
lat = ds.lat_rho.values

DX = 1/ds.pm.values
DY = 1/ds.pn.values

DA = DX * DY                                 # cell horizontal area 

CC = dict()                                  # this is for holding fields extracted on sections

CC['lon'] = lon
CC['lat'] = lat

CC['DA'] = DA

# Fields that do not change with time
# 
h = ds.h.values       
CC['h'] = h

mask_rho = ds.mask_rho.values
CC['mask_rho']=mask_rho

# Fields that do change with time
#

zeta = ds.zeta.values.squeeze()
z_rho, z_w = zrfun.get_z(h, zeta, S)

# hypoxia thresholds 
dzr = np.diff(z_w, axis = 0)
oxy = ds.oxygen.values.squeeze()

dzrm = np.ma.masked_where(oxy>106.6,dzr) 
mild_dz = dzrm.sum(axis=0)
del dzrm 

dzrm = np.ma.masked_where(oxy>60.9,dzr) 
hyp_dz = dzrm.sum(axis=0)
del dzrm 

dzrm = np.ma.masked_where(oxy>21.6,dzr) 
severe_dz = dzrm.sum(axis=0)

dzrm = np.ma.masked_where(oxy>0,dzr) 
anoxic_dz = dzrm.sum(axis=0)

CC['mild_dz'] = mild_dz
CC['hyp_dz'] = hyp_dz
CC['severe_dz'] = severe_dz
CC['anoxic_dz'] = anoxic_dz

# add in Oag, see toy_Oag.py
import gsw
# Pressure calcs: Lpres
ZZ = z_rho-zeta                         # zeta adjusted, see p_ref comment below
Lpres = gsw.p_from_z(ZZ, lat)           # pressure [dbar]

# Note gsw.p_from_z uses p_ref = 0 and requires z's to be neg. So need to adjust zeta to 'zero' for pressure calcs. There may be a better gsw function for this w/ adjustable p_ref, but prob takes longer than subtracting the two arrays (?) 

# grab and convert physical variables + alkalinity and TIC from LO history files 
SP = ds.salt.values.squeeze()
TI = ds.temp.values.squeeze()
ALK = ds.alkalinity.values.squeeze()
TIC = ds.TIC.values.squeeze()

SA = gsw.SA_from_SP(SP, Lpres, lon, lat)  # Q from dm_pfun.py: isn't LO output SA? 
CT = gsw.CT_from_pt(SA, TI)
rho = gsw.rho(SA, CT, Lpres)              # in situ density
Ltemp = gsw.t_from_CT(SA, CT, Lpres)      # in situ temperature

# convert from umol/L to umol/kg using in situ dentity
Lalkalinity = 1000 * ALK / rho
Lalkalinity[Lalkalinity < 100] = np.nan   # Q from dm_pfun.py: why? 
 
LTIC = 1000 * TIC / rho
LTIC[LTIC < 100] = np.nan                 # Q from dm_pfun.py: why? 

Lpres = zfun.fillit(Lpres)               
Ltemp = zfun.fillit(Ltemp)
# zfun.fillit ensures a is an array with nan's for masked values
# instead of a masked array  

# calculate aragonite saturation:
# For CO2SYS: All temperatures are in °C, 
#             all salinities are on the PSS, 
#             and all pressures are in dbar. 
from PyCO2SYS import CO2SYS   

ARAG = np.full(np.shape(SP),np.nan)
A = np.shape(Lalkalinity)
for ii in range(A[0]): 
    aALK = Lalkalinity[ii,:,:].squeeze()
    aTIC = LTIC[ii,:,:].squeeze()
    aTemp = Ltemp[ii,:,:].squeeze()
    aPres = Lpres[ii,:,:].squeeze()
    aSalt = SP[ii,:,:].squeeze()
    
    CO2dict = CO2SYS(aALK, aTIC, 1, 2, aSalt, aTemp, aTemp,
    aPres, aPres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)             # assumptions from dm_pfun.py
    
    aARAG = CO2dict['OmegaARout']
    aARAG = aARAG.reshape((aSalt.shape))                         # reshape 
    ARAG[ii,:,:] = np.expand_dims(aARAG, axis=0)

dzr = np.diff(z_w, axis = 0)
dzrm = np.ma.masked_where(ARAG>1,dzr) 
corrosive_dz = dzrm.sum(axis=0)

CC['corrosive_dz'] = corrosive_dz

# put them in a dataset, ds1
NR, NC = CC['hyp_dz'].shape        # this feels sloppy 
ot = ds.ocean_time.values          # an array with dtype='datetime64[ns]'

ds1 = Dataset()
ds1['ocean_time'] = (('ocean_time'), ot, {'long_name':ds.ocean_time.long_name,'units':ds.ocean_time.units})

ds1['mild_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['mild_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of mild hypoxic layer'})
ds1['hyp_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['hyp_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of hypoxic layer'})
ds1['severe_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['severe_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of severe hypoxic layer'})
ds1['anoxic_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['anoxic_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of anoxic layer'})

ds1['corrosive_dz'] = (('ocean_time', 'eta_rho', 'xi_rho'), CC['corrosive_dz'].reshape(1,NR,NC), {'units':'m', 'long_name': 'Thickness of undersaturated layer'})

ds1['DA'] = (('eta_rho', 'xi_rho'), CC['DA'], {'units':'m^2', 'long_name': 'cell horizontal area '})
ds1['mask_rho'] = (('eta_rho', 'xi_rho'), CC['DA'], {'flag_values':[0., 1.],'flag_meanings':'land water','long_name': 'mask on RHO-points'})
ds1['h'] = (('eta_rho', 'xi_rho'), CC['h'], {'units':ds.h.units, 'long_name': ds.h.long_name})

ds1['Lat'] = (('eta_rho', 'xi_rho'), CC['lat'], {'units':'degree_north','long_name': 'latitude of RHO-points'})
ds1['Lon'] = (('eta_rho', 'xi_rho'), CC['lon'], {'units':'degree_east','long_name': 'longitude of RHO-points'})

ds1.to_netcdf(args.out_fn, unlimited_dims='ocean_time')

