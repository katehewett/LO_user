"""
Function to do the calculation for a single clipped file
calculates GSW vars > drho/dz dT/dz
finds max drho(dT)/dz + returns the max depth + values

THIS IS A TESTING 1 OFF FILE for
get_one_cline.py 

Time to get initial fields = 0.65 sec
Time to apply GSW = 0.81 sec
Time to calc dp/dz and dt/dz = 0.09 sec

# TODO make a function for the masked data flag when finding mean 

Time to get initial fields = 2.32 sec
Time to apply GSW = 4.79 sec
Time to calc SML + assign T/S= 1.33 sec
Time to calc BML + assign T/S= 1.29 sec
Time to calc max dT/dz(drho/dz) and assign SA(CT)= 1.26 sec
Time to calc interior gradients and T/S= 2.11 sec

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

#plotting things
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import cmcrameri.cm as cmc

#box_fn_in = '/Users/katehewett/Documents/LO_output/extract/cas6_v0_live/clines/shelf_box_2022.08.08_2022.08.09/temp_box_dir/box_000001.nc'
#his_in = '/Users/katehewett/Documents/LO_roms/cas6_v0_live/f2022.08.08/ocean_his_0021.nc'

box_fn_in = '/Users/katehewett/Documents/LO_output/extract/cas7_t0_x4b/clines/shelf_box_2017.12.12_2017.12.13/temp_box_dir/box_000000.nc'
his_in = '/Users/katehewett/Documents/LO_roms/cas7_t0_x4b/f2017.12.12/ocean_his_0021.nc'

threshold = 0.03
######## ######## ######## 

tt0 = time()
dsb = open_dataset(box_fn_in, decode_times=False) # TODO UPDATE
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat

# add z variables
# NOTE: this only works if you have salt as a saved field!
G, S, T = zrfun.get_basic_info(his_in)   # grid information # TODO UPDATE

lat = np.nanmean(G['lat_rho'])
lon = np.nanmean(G['lon_rho']) # for gsw calcs

NT, NZ, NR, NC = np.shape(dsb.salt)
h = dsb.h.values
zeta = dsb.zeta.values
mask_rho = dsb.mask_rho

z_rho, z_w = zrfun.get_z(h, zeta, S) # (30, 1111, 356) (31, 1111, 356)
NW = z_w.shape[0]

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
  
# grabbing SIG0 for all depths above -250m to see po over time for the box just curious

print('Time to apply GSW = %0.2f sec' % (time()-tt0))
sys.stdout.flush()
     
# 1. calc threshold SML 
tt0 = time()       

Dsurf = SIG0[:,-1,:,:]
sthresh = np.abs(SIG0 - Dsurf)

# SML calc 
# flip along z axis so that ocean packed surface:bottom
# find where the first value passes the threshold ("B")
fSIG0 = np.flip(SIG0, axis=1)
fCT = np.flip(CT, axis=1)
fSA = np.flip(SA, axis=1)
fsthresh = np.flip(sthresh, axis=1)
fz_rho = np.flip(z_rho, axis=1)
fz_w = np.flip(z_w, axis=1)

# the flag returned by argmin acts on SIG0, which is on 
# z_rho. The base of the sml, with our thresh, occurs somewhere 
# between those two z_rho points at B and B-1. 
# Assign the midpoint fz_w[:,B,:,:] as the base. 
# Note: if sml is found to be zero, the depth is saved as 
# the first z_w point in the flipped "fz_w" array. 
# So when we calc the thickness, zeta - SML base, 
# we get zero for the SML thickness.
B = np.argmin(fsthresh<threshold,axis=1,keepdims=True)
SML_z = np.take_along_axis(fz_w,B,axis=1)

# if SML is zero and in water keep SIG0 value just as the the surface 
# value else if in water and non zero we calc the avg SIG0 for the 
# two z_rhos at B and B-1:
Bval = np.take_along_axis(fSIG0,B,axis=1) 
bmask = (B.squeeze()!=0) & (mask_rho.values == 1) 
b1 = np.take_along_axis(fSIG0,B-1,axis=1)
b2 = np.take_along_axis(fSIG0,B,axis=1)
bb = (b1+b2)/2
Bval[:,:,bmask] = bb[:,:,bmask]

# save mean value of T and S in the SML; 1st initialize: 
# note: when SML = 0; we still save the surface CT(SA); and save the land as NaNs
SML_CT = np.expand_dims(fCT[:,0,:,:],axis=1)    
SML_SA = np.expand_dims(fSA[:,0,:,:],axis=1)

# below the SML set values to NaN  & find the mean CT(SA) in the SML                         
fCTsurf = np.copy(fCT)
fCTsurf[fz_rho<SML_z] = np.NaN          
# included this step s/t don't get error msg when sum an empty slice
masked_data = np.ma.masked_array(fCTsurf, np.isnan(fCTsurf)) 
mCTsurf = np.ma.average(masked_data, axis=1,keepdims=True)
#mCTsurf = np.nanmean(fCTsurf,axis=1,keepdims=True)
del masked_data

fSAsurf = np.copy(fSA)
fSAsurf[fz_rho<SML_z] = np.NaN  
# included this step s/t don't get error msg when sum an empty slice
masked_data = np.ma.masked_array(fSAsurf, np.isnan(fSAsurf)) 
mSAsurf = np.ma.average(masked_data, axis=1,keepdims=True)        
#mSAsurf = np.nanmean(fSAsurf,axis=1,keepdims=True)
del masked_data

SML_CT[:,:,bmask] = mCTsurf[:,:,bmask]
SML_SA[:,:,bmask] = mSAsurf[:,:,bmask]

print('Time to calc SML + assign T/S= %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# BML calc :: similar to SML but not using flipped data 
tt0 = time()  

Dbot = SIG0[:,0,:,:]
bthresh = np.abs(SIG0 - Dbot)

C = np.argmin(bthresh<threshold,axis=1,keepdims=True)
BML_z = np.take_along_axis(z_w,C,axis=1)

Cval = np.take_along_axis(SIG0,C,axis=1) 
cmask = (C.squeeze()!=0) & (mask_rho.values == 1) 
c1 = np.take_along_axis(SIG0,C-1,axis=1)
c2 = np.take_along_axis(SIG0,C,axis=1)
cc = (c1+c2)/2
Cval[:,:,cmask] = cc[:,:,cmask]

BML_CT = np.expand_dims(CT[:,0,:,:],axis=1)    
BML_SA = np.expand_dims(SA[:,0,:,:],axis=1)
                        
CTbot = np.copy(CT)
CTbot[z_rho>BML_z] = np.NaN      
# included this step s/t don't get error msg when sum an empty slice
masked_data = np.ma.masked_array(CTbot, np.isnan(CTbot)) 
mCTbot = np.ma.average(masked_data, axis=1,keepdims=True)      
#mCTbot = np.nanmean(CTbot,axis=1,keepdims=True)
del masked_data

SAbot = np.copy(SA)
SAbot[z_rho>BML_z] = np.NaN          
# included this step s/t don't get error msg when sum an empty slice
masked_data = np.ma.masked_array(SAbot, np.isnan(SAbot)) 
mSAbot = np.ma.average(masked_data, axis=1,keepdims=True)      
#mSAbot = np.nanmean(SAbot,axis=1,keepdims=True)

BML_CT[:,:,cmask] = mCTbot[:,:,cmask]
BML_SA[:,:,cmask] = mSAbot[:,:,cmask]

print('Time to calc BML + assign T/S= %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# find loco of max drho/dz and dT/dz and grab value
tt0 = time()

dp = np.diff(SIG0,axis=1)
dT = np.diff(CT,axis=1)
dz = np.diff(z_rho,axis=1)
dpdz = dp/dz
dTdz = dT/dz

D = np.argmax(np.abs(dpdz),axis=1,keepdims=True)
Dz = np.take_along_axis(z_w,D+1,axis=1)
Dval = np.take_along_axis(dpdz,D,axis=1)

d1 = np.take_along_axis(SA,D,axis=1)         # grab SA(CT) at dpdz max
d2 = np.take_along_axis(SA,D-1,axis=1)
Dsa = (d1+d2)/2

d1 = np.take_along_axis(CT,D,axis=1)
d2 = np.take_along_axis(CT,D-1,axis=1)
Dct = (d1+d2)/2

F = np.argmax(np.abs(dTdz),axis=1,keepdims=True)     
Fz = np.take_along_axis(z_w,F+1,axis=1)
Fval = np.take_along_axis(dTdz,F,axis=1)

f1 = np.take_along_axis(SA,F,axis=1)         # grab SA(CT) at dTdz max
f2 = np.take_along_axis(SA,F-1,axis=1)
Fsa = (f1+f2)/2

f1 = np.take_along_axis(CT,F,axis=1)
f2 = np.take_along_axis(CT,F-1,axis=1)
Fct = (d1+d2)/2

print('Time to calc max dT/dz(drho/dz) and assign SA(CT)= %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# grab mean CT(SA)(dpdz)(dTdz) in interior layer 
tt0 = time()

iCT = np.copy(CT)
iSA = np.copy(SA)
idpdz = np.copy(dpdz)
idTdz = np.copy(dTdz)

iCT[(z_rho[:,:,:,:]>SML_z[:,:,:,:])] = np.nan
iCT[(z_rho[:,:,:,:]<BML_z[:,:,:,:])] = np.nan
masked_data = np.ma.masked_array(iCT, np.isnan(iCT)) 
INT_CT = np.ma.average(masked_data, axis=1,keepdims=True) 
#INT_CT = np.nanmean(iCT,axis=1,keepdims=True)

iSA[(z_rho[:,:,:,:]>SML_z[:,:,:,:])] = np.nan
iSA[(z_rho[:,:,:,:]<BML_z[:,:,:,:])] = np.nan
masked_data = np.ma.masked_array(iSA, np.isnan(iSA)) 
INT_SA = np.ma.average(masked_data, axis=1,keepdims=True) 
#INT_SA = np.nanmean(iSA,axis=1,keepdims=True)

idpdz[(z_w[:,1:NW-1,:,:]>SML_z[:,:,:,:])] = np.nan
idpdz[(z_w[:,1:NW-1,:,:]<BML_z[:,:,:,:])] = np.nan
masked_data = np.ma.masked_array(idpdz, np.isnan(idpdz)) 
INT_dpdz = np.ma.average(masked_data, axis=1,keepdims=True) 
#INT_dpdz = np.nanmean(idpdz,axis=1,keepdims=True)

idTdz[(z_w[:,1:NW-1,:,:]>SML_z[:,:,:,:])] = np.nan
idTdz[(z_w[:,1:NW-1,:,:]<BML_z[:,:,:,:])] = np.nan
masked_data = np.ma.masked_array(idTdz, np.isnan(idTdz)) 
INT_dTdz = np.ma.average(masked_data, axis=1,keepdims=True)
#INT_dTdz = np.nanmean(idTdz,axis=1,keepdims=True)

print('Time to calc interior gradients and T/S= %0.2f sec' % (time()-tt0))
sys.stdout.flush()

#save some values 
tt0 = time()

po = np.nanmean(np.where(z_rho>-250,SIG0,np.nan))
print('po= %0.2f kg/m3' % (po))

SAstack = np.hstack((SML_SA,INT_SA,Dsa,BML_SA))
CTstack = np.hstack((SML_CT,INT_CT,Dct,BML_CT))

ds1 = Dataset()
ds1['ocean_time'] = dsb.ocean_time

ds1['zeta'] = (('ocean_time', 'lon_rho', 'lat_rho'), dsb.zeta.values, {'units':'m', 'long_name':'free-surface'})
ds1['zSML'] = (('ocean_time', 'lon_rho', 'lat_rho'), np.squeeze(SML_z,axis=1), {'units':'m', 'long_name': 'depth of bottom of SML threshold'})
ds1['zBML'] = (('ocean_time', 'lon_rho', 'lat_rho'), np.squeeze(BML_z,axis=1), {'units':'m', 'long_name': 'depth of bottom of SML threshold'})
ds1['zDGmax'] = (('ocean_time', 'lon_rho', 'lat_rho'), np.squeeze(Dz,axis=1), {'units':'m', 'long_name': 'depth of max density gradient'})
ds1['DGmax'] = (('ocean_time', 'lon_rho', 'lat_rho'), np.squeeze(Dval,axis=1), {'units':'kg/m3/m', 'long_name': 'value of max density gradient'})
ds1['DGint'] = (('ocean_time', 'lon_rho', 'lat_rho'), np.squeeze(INT_dpdz,axis=1), {'units':'kg/m3/m', 'long_name': 'mean density gradient between SML and BML'})

ds1.coords['zone'] = [0,1,2,3]
ds1['SA'] = (('ocean_time', 'zone', 'lon_rho', 'lat_rho'), SAstack, {'units':'g/kg', 'long_name': 'mean Absolute Salinity by zone'})
ds1['SA'].attrs['zone_values'] = [0,1,2,3]
ds1['SA'].attrs['zone_meanings'] = ['SML','INT(SML:BML)','DGmax','BML'] 

ds1['CT'] = (('ocean_time', 'zone', 'lon_rho', 'lat_rho'), CTstack, {'units':'degrees celcius', 'long_name': 'mean Conservative Temperature by zone'})
ds1['CT'].attrs['zone_values'] = [0,1,2,3]
ds1['CT'].attrs['zone_meanings'] = ['SML','INT(SML:BML)','DGmax','BML']

sys.exit()

######## ######## ######## 
# remove this testing after get the pycnocline to work here. 
##### PRINTING THINGS 

# get indices for plotting // ilon and ilat will be the index of the mlon mlat inputs::
mlon = -124.06
mlat = 45.135

Lon = lon[0,:]
Lat = lat[:,0]
# error checking
if (mlon < Lon[0]) or (mlon > Lon[-1]):
    print('ERROR: lon out of bounds ' + moor_fn.name)
    sys.exit()
if (mlat < Lat[0]) or (mlat > Lat[-1]):
    print('ERROR: lat out of bounds ' + moor_fn.name)
    sys.exit()
# get indices
ilon = zfun.find_nearest_ind(Lon, mlon)
ilat = zfun.find_nearest_ind(Lat, mlat)


# initialize plot 
plt.close('all')
fs=12
plt.rc('font', size=fs)
fig = plt.figure(figsize=(16,8))

ax1 = plt.subplot2grid((1,4),(0,0),colspan=1,rowspan=1)
ax2 = plt.subplot2grid((1,4),(0,1),colspan=1,rowspan=1)
ax3 = plt.subplot2grid((1,4),(0,2),colspan=1,rowspan=1)
ax4 = plt.subplot2grid((1,4),(0,3),colspan=1,rowspan=1)

# plot T and SIG0
Tmin = np.floor(np.nanmin(CT[:,:,ilat,ilon])) - 0.5
Tmax = np.floor(np.nanmax(CT[:,:,ilat,ilon])) + 0.5
tticks = np.arange(Tmin,Tmax+1,2)

Tplot = ax1.plot(CT[:,:,ilat,ilon].squeeze(),z_rho[:,:,ilat,ilon].squeeze(),color='blue',marker='.',label = 'CT') 
#ax1.plot(SML_CT[:,:,ilat,ilon].squeeze(),SML_z[:,:,ilat,ilon].squeeze(),color='blue',marker='o',label = 'CT') 
ax1.set_ylabel('depth m')
ax1.set_xlabel('CT',color='blue')
ax1.tick_params(axis='x',color='blue')
ax1.xaxis.label.set_color('blue')
ax1.set_xlim([Tmin, Tmax])
ax1.set_xticks(tticks)
ax1.tick_params(labelcolor='blue')
#ax1.axhline(y = z_dtmax[:,:,ilon,ilat], color = 'k', linestyle = '--') 

axsig = ax1.twiny()
SIGmin = np.floor(np.min(SIG0[:,:,ilat,ilon])) - 0.5
SIGmax = np.floor(np.max(SIG0[:,:,ilat,ilon])) + 0.5
sigticks = np.arange(SIGmin,SIGmax+1,1)

Splot = axsig.plot(SIG0[:,:,ilat,ilon].squeeze(),z_rho[:,:,ilat,ilon].squeeze(),color='red',marker='.',label = 'SIG0')
axsig.set_xlabel('SIG0',color='red')
axsig.tick_params(axis='x',color='red')
axsig.xaxis.label.set_color('red')
axsig.set_xlim([SIGmin, SIGmax])
axsig.set_xticks(sigticks)
axsig.tick_params(labelcolor='red')
#axsig.axhline(y = z_dpmax[:,:,ir,ic], color = 'k', linestyle = '-') 
axsig.plot(Bval[:,:,ilat,ilon].squeeze(),SML_z[:,:,ilat,ilon].squeeze(),'x')
axsig.plot(Cval[:,:,ilat,ilon].squeeze(),BML_z[:,:,ilat,ilon].squeeze(),'x')

# plot thresh
pp = ax2.plot(np.abs(sthresh[:,:,ilat,ilon].squeeze()),z_rho[:,:,ilat,ilon].squeeze(),color='green',marker='.',label = 'surf') 
ax2.plot(np.abs(bthresh[:,:,ilat,ilon].squeeze()),z_rho[:,:,ilat,ilon].squeeze(),color='orange',marker='.',label = 'bot') 
ax2.set_ylabel('depth m')
ax2.set_xlabel('surf:bot',color='blue')
ax2.set_title('thresholds (orange B; green S)')
ax2.tick_params(axis='x',color='blue')
ax2.xaxis.label.set_color('blue')
#ax2.set_xlim([Tmin, Tmax])
#ax2.set_xticks(tticks)
ax2.tick_params(labelcolor='blue')
ax2.axvline(x=threshold)

# plot dp/dz 
ax3.plot(np.abs(dpdz[:,:,ilat,ilon].squeeze()),z_w[:,1:NW-1,ilat,ilon].squeeze(),color='green',marker='.',label = 'dpdz') 
ax3.plot(np.abs(Dval[:,:,ilat,ilon].squeeze()),Dz[:,:,ilat,ilon].squeeze(),color='black',marker='o',label = 'max') 
ax3.set_ylabel('depth m')
ax3.set_xlabel('dp/dz')
ax3.tick_params(axis='x')

# plot dp/dz 
ax4.plot(np.abs(dTdz[:,:,ilat,ilon].squeeze()),z_w[:,1:NW-1,ilat,ilon].squeeze(),color='green',marker='.',label = 'dpdz') 
ax4.plot(np.abs(Fval[:,:,ilat,ilon].squeeze()),Fz[:,:,ilat,ilon].squeeze(),color='black',marker='o',label = 'max') 
ax4.set_ylabel('depth m')
ax4.set_xlabel('dT/dz')
ax4.tick_params(axis='x')

# po = np.where(z_rho>-250,SIG0,np.nan)  

