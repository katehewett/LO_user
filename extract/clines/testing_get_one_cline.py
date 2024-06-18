"""
Function to do the calculation for a single clipped file
calculates GSW vars > drho/dz dT/dz
finds max drho(dT)/dz + returns the max depth + values

THIS IS A TESTING 1 OFF FILE for
get_one_cline.py 

Time to get initial fields = 0.65 sec
Time to apply GSW = 0.81 sec
Time to calc dp/dz and dt/dz = 0.09 sec

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
NT, NZ, NR, NC = np.shape(dsb.salt)
h = dsb.h.values
zeta = dsb.zeta.values
lat = dsb.lat_rho.values
lon = dsb.lon_rho.values
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
 
po = np.where(z_rho>-250,SIG0,np.nan)    
# grabbing SIG0 for all depths above -250m to see po over time for the box just curious

print('Time to apply GSW = %0.2f sec' % (time()-tt0))
sys.stdout.flush()
     
# 1. calc dp/dz and dT/dz 
tt0 = time()       

dp = np.diff(SIG0,axis=1)
dT = np.diff(CT,axis=1) 
dz = np.diff(z_rho,axis=1)

dpdz = dp/dz
dtdz = dT/dz

C = np.argmin(dpdz,axis=1,keepdims=True)
z_dpmax = np.take_along_axis(z_w,C+1,axis=1)
val_dpmax = np.take_along_axis(dpdz,C,axis=1)

D = np.argmax(dtdz,axis=1,keepdims=True)
z_dtmax = np.take_along_axis(z_w,D+1,axis=1)
val_dtmax = np.take_along_axis(dtdz,D,axis=1)

print('Time to calc dp/dz and dt/dz = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

######## ######## ######## 
# remove this testing after get the pycnocline to work here. 
##### PRINTING THINGS 

# get indices for plotting // ilon and ilat will be the index of the mlon mlat inputs::
mlon = -124.4
mlat = 44.65
Lon = G['lon_rho'][0,:]
Lat = G['lat_rho'][:,0]
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

ax1 = plt.subplot2grid((1,3),(0,0),colspan=1,rowspan=1)
ax2 = plt.subplot2grid((1,3),(0,1),colspan=1,rowspan=1)
ax3 = plt.subplot2grid((1,3),(0,2),colspan=1,rowspan=1)

# plot T and SIG0
Tmin = np.floor(np.min(CT[:,:,ilon,ilat]))-1
Tmax = np.floor(np.max(CT[:,:,ilon,ilat]))+1
tticks = np.arange(Tmin,Tmax+1,2)

Tplot = ax1.plot(CT[:,:,ilon,ilat].squeeze(),z_rho[:,:,ilon,ilat].squeeze(),color='blue',marker='.',label = 'CT') 
ax1.set_ylabel('depth m')
ax1.set_xlabel('CT',color='blue')
ax1.tick_params(axis='x',color='blue')
ax1.xaxis.label.set_color('blue')
ax1.set_xlim([Tmin, Tmax])
ax1.set_xticks(tticks)
ax1.tick_params(labelcolor='blue')
#ax1.axhline(y = z_dtmax[:,:,ilon,ilat], color = 'k', linestyle = '--') 

axsig = ax1.twiny()
SIGmin = np.floor(np.min(SIG0[:,:,ilon,ilat]))-0.5
SIGmax = np.floor(np.max(SIG0[:,:,ilon,ilat]))+0.5
sigticks = np.arange(SIGmin,SIGmax+1,1)

Splot = axsig.plot(SIG0[:,:,ilon,ilat].squeeze(),z_rho[:,:,ilon,ilat].squeeze(),color='red',marker='.',label = 'SIG0')
axsig.set_xlabel('SIG0',color='red')
axsig.tick_params(axis='x',color='red')
axsig.xaxis.label.set_color('red')
axsig.set_xlim([SIGmin, SIGmax])
axsig.set_xticks(sigticks)
axsig.tick_params(labelcolor='red')
#axsig.axhline(y = z_dpmax[:,:,ir,ic], color = 'k', linestyle = '-') 

# plot dp/dz
pp = ax2.plot(dpdz[:,:,ilon,ilat].squeeze(),z_w[:,1:NW-1,ilon,ilat].squeeze(),color='red',marker='.',label = 'dpdz') 
ax2.set_ylabel('depth m')
ax2.set_xlabel('drho/dz',color='red')
ax2.tick_params(axis='x',color='red')
ax2.xaxis.label.set_color('red')
#ax2.set_xlim([Tmin, Tmax])
#ax2.set_xticks(tticks)
ax2.tick_params(labelcolor='red')
ax2.plot(val_dpmax[:,:,ilon,ilat],z_dpmax[:,:,ilon,ilat],color='black',marker='x')

# plot dt/dz
pt = ax3.plot(dtdz[:,:,ilon,ilat].squeeze(),z_w[:,1:NW-1,ilon,ilat].squeeze(),color='blue',marker='.',label = 'dpdz') 
ax3.set_ylabel('depth m')
ax3.set_xlabel('dT/dz',color='blue')
ax3.tick_params(axis='x',color='blue')
ax3.xaxis.label.set_color('blue')
#ax2.set_xlim([Tmin, Tmax])
#ax2.set_xticks(tticks)
ax3.tick_params(labelcolor='blue')
ax3.plot(val_dtmax[:,:,ilon,ilat],z_dtmax[:,:,ilon,ilat],color='black',marker='x')