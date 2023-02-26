"""
Testing Oag calc's for volume extractions 

Whole LO domain all levels: 
Entering CO2sys calcs = 4.92 sec
Total processing time = 11.23 sec < 1 depth layer 

Entering CO2sys calcs = 4.92 sec
Total processing time = 175.14 sec sec < all 30 layers 

"""

# imports
import os
import sys 
import xarray as xr
import numpy as np

from lo_tools import Lfun, zfun, zrfun
from xarray import open_dataset, Dataset
from time import time

os.chdir('/Users/katehewett/Documents/LO_roms/cas6_v0_live/f2022.08.08')
ds = open_dataset('ocean_his_0018.nc', decode_times=False)

# grid information 
G, S, T = zrfun.get_basic_info('/Users/katehewett/Documents/LO_roms/cas6_v0_live/f2022.08.08/ocean_his_0018.nc')                 # grid information

lon = G['lon_rho']
lat = G['lat_rho']
    
h = ds.h.values  
zeta = ds.zeta.values.squeeze()          # surface elevation 
z_rho, z_w = zrfun.get_z(h, zeta, S)     

tt0 = time()
import gsw

# Pressure calcs: Lpres
ZZ = z_rho-zeta                         # this is the faux zeta
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

print('Entering CO2sys calcs = %0.2f sec' % (time()-tt0))

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
    aPres, aPres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)             # from dm_pfun.py
    
    aARAG = CO2dict['OmegaARout']
    aARAG = aARAG.reshape((aSalt.shape))
    #y = np.expand_dims(aARAG, axis=0)
    ARAG[ii,:,:] = np.expand_dims(aARAG, axis=0)

# to test later, can you input 3D to CO2SYS? 
#CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, SP, Ltemp, Ltemp,
#Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)             # from dm_pfun.py 

dzr = np.diff(z_w, axis = 0)
oxy = ds.oxygen.values.squeeze()
dzrm = np.ma.masked_where(ARAG>1,dzr) 
corrosive_dz = dzrm.sum(axis=0)

print('Total processing time = %0.2f sec' % (time()-tt0))
