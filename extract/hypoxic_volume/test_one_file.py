"""
Code to test the speed and meory use of a single hypoxic-corrosive volume extraction

"""
import xarray as xr
import numpy as np
import gsw
import PyCO2SYS as pyco2
from pathlib import Path
from lo_tools import zrfun
from time import time
import sys

#fn = Path('/Users/pm8/Documents/LO_roms/cas6_v0_live/f2021.07.03/ocean_his_0021.nc')
fn = Path('/Users/katehewett/Documents/LO_roms/cas6_v0_live/f2022.08.08/ocean_his_0021.nc')

tt0 = time()
ds = xr.open_dataset(fn, decode_times=False)
G, S, T = zrfun.get_basic_info(fn)    # grid information
lon = ds.lon_rho.values
lat = ds.lat_rho.values
DA = G['DX'] * G['DY']
CC = dict()
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
dzrm = dzr.copy()
dzrm[oxy>106.6] = 0
mild_dz = dzrm.sum(axis=0)
CC['mild_dz'] = mild_dz
print('Time to get mild_dz = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

tt0 = time()
# carbon claculations
p = gsw.p_from_z(z_rho, lat)           # pressure [dbar]
SP = ds.salt.values.squeeze()
PT = ds.temp.values.squeeze()
ALK = ds.alkalinity.values.squeeze()
TIC = ds.TIC.values.squeeze()

SA = gsw.SA_from_SP(SP, p, lon, lat)
CT = gsw.CT_from_pt(SA, PT)
rho = gsw.rho(SA, CT, p)           # in situ density
ti = gsw.t_from_CT(SA, CT, p)      # in situ temperature

# convert from umol/L to umol/kg using in situ dentity
ALK1 = 1000 * ALK / rho
TIC1 = 1000 * TIC / rho

# I'm not sure if this is needed
#ALK1[ALK1 < 100] = np.nan                 # Q from dm_pfun.py: why? 
#TIC1[TIC1 < 100] = np.nan                 # Q from dm_pfun.py: why? 

# calculate aragonite saturation:
# For CO2SYS: All temperatures are in °C, 
#             all salinities are on the PSS, 
#             and all pressures are in dbar. 
ARAG = np.nan * np.ones(SP.shape)
nz, nr, nc = SP.shape
amat = np.nan * np.ones((nr,nc))

mask_arag = mask_rho                       # kh added h to reduce input to co2sys
mask_arag[h<400] = 1

for ii in range(nz):
# for ii in range(2):
    tt00 = time()
    #aALK = ALK1[ii,:,:].squeeze()[mask_rho==1]    
    #aTIC = TIC1[ii,:,:].squeeze()[mask_rho==1]
    #aTemp = ti[ii,:,:].squeeze()[mask_rho==1]
    #aPres = p[ii,:,:].squeeze()[mask_rho==1]
    #aSalt = SP[ii,:,:].squeeze()[mask_rho==1]        
    
    aALK = ALK1[ii,:,:].squeeze()[mask_arag==1]    
    aTIC = TIC1[ii,:,:].squeeze()[mask_arag==1]
    aTemp = ti[ii,:,:].squeeze()[mask_arag==1]
    aPres = p[ii,:,:].squeeze()[mask_arag==1]
    aSalt = SP[ii,:,:].squeeze()[mask_arag==1]
    
    if any(aSalt<0):
        ALK1[aSalt<0] = -9        #co2sys nan flag
        aTIC[aSalt<0] = -9        
        aTemp[aSalt<0] = -9        
        aPres[aSalt<0] = -9        
        aSalt[aSalt<0] = -9  
    if any(aTemp<0):
        ALK1[aSalt<0] = -9        
        aTIC[aSalt<0] = -9        
        aTemp[aSalt<0] = -9        
        aPres[aSalt<0] = -9        
        aSalt[aSalt<0] = -9   
    if any(aPres<0):
        ALK1[aSalt<0] = -9        
        aTIC[aSalt<0] = -9        
        aTemp[aSalt<0] = -9        
        aPres[aSalt<0] = -9        
        aSalt[aSalt<0] = -9   
        
    # note: still need to consult:
    # https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/
    # to make sure the other inputs are handled correctly (e.g. what does 50 mean?)
    CO2dict = pyco2.sys(par1=aALK, par1_type=1, par2=aTIC, par2_type=2,
        salinity=aSalt, temperature=aTemp, pressure=aPres, opt_buffers_mode=0)
    # CO2dict = CO2SYS(aALK, aTIC, 1, 2, aSalt, aTemp, aTemp,
    #     aPres, aPres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)             # assumptions from dm_pfun.py
    aARAG = CO2dict['saturation_aragonite']
    aamat = amat.copy()
    #aamat[mask_rho==1] = aARAG
    aamat[mask_arag==1] = aARAG
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