"""
Testing Oag calc's for volume extractions 

Whole LO domain all levels: 
Entering CO2sys calcs = 4.92 sec
Total processing time = 11.23 sec < 1 depth layer 

Entering CO2sys calcs = 4.92 sec
Total processing time = 175.14 sec sec < all 30 layers 
[173.59 sec < with corrosive calc]
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

#print('Entering CO2sys calcs = %0.2f sec' % (time()-tt0))

# calculate aragonite saturation:
# For CO2SYS: All temperatures are in °C, 
#             all salinities are on the PSS, 
#             and all pressures are in dbar. 

from PyCO2SYS import CO2SYS               


# take a subset of vars and see if we can put 3d to co2sys and compare times
sALK = Lalkalinity[:,567:600,220:280]
sTIC = LTIC[:,567:600,220:280]
sTemp = Ltemp[:,567:600,220:280]
sPres = Lpres[:,567:600,220:280]
sSalt = SP[:,567:600,220:280].squeeze()
    
tt0 = time()
ARAG = np.full(np.shape(sSalt),np.nan)
A = np.shape(Lalkalinity)
for ii in range(A[0]): 
    aALK = sALK[ii,:,:].squeeze()
    aTIC = sTIC[ii,:,:].squeeze()
    aTemp = sTemp[ii,:,:].squeeze()
    aPres = sPres[ii,:,:].squeeze()
    aSalt = sSalt[ii,:,:].squeeze()
    
    CO2dict = CO2SYS(aALK, aTIC, 1, 2, aSalt, aTemp, aTemp,
    aPres, aPres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)             # from dm_pfun.py
    
    aARAG = CO2dict['OmegaARout']
    aARAG = aARAG.reshape((aSalt.shape))
    #y = np.expand_dims(aARAG, axis=0)
    ARAG[ii,:,:] = np.expand_dims(aARAG, axis=0)

# to test later, can you input 3D to CO2SYS? 
#CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, SP, Ltemp, Ltemp,
#Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)             # from dm_pfun.py 
print('Loop, total processing time = %0.2f sec' % (time()-tt0))
    

falkalinity = np.matrix.flatten(Lalkalinity)
fSP = np.matrix.flatten(SP)
fTIC = np.matrix.flatten(LTIC)
ftemp = np.matrix.flatten(Ltemp)
fpres = np.matrix.flatten(Lpres)

tt0 = time()
CO2dict = CO2SYS(falkalinity, fTIC, 1, 2, fSP, ftemp, ftemp,
    fpres, fpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)             # from dm_pfun.py
    
ARAG2 = CO2dict['OmegaARout']
ARAG2 = ARAG2.reshape((SP.shape))

print('3D, total processing time = %0.2f sec' % (time()-tt0))





