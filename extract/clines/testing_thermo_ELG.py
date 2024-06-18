"""
Function to do the calculation for a single clipped file
calculates GSW vars > drho/dz dT/dz
finds max drho(dT)/dz + returns the max depth + values

THIS IS A TESTING FILE for ELG method to find thermocline 
at first doing one location // enter lat and lon
then modify to fit into get one cline - 
do winter first and then summer 

If you want this code to work then you need to have run test_extract_clines.py 
for the shelf box job & have access to one of these history files, and 
will need to adjust the path on your computer 
Summer Test 
box_fn_in = '/Users/katehewett/Documents/LO_output/extract/cas6_v0_live/clines/shelf_box_2022.08.08_2022.08.09/temp_box_dir/box_000001.nc'
his_in = '/Users/katehewett/Documents/LO_roms/cas6_v0_live/f2022.08.08/ocean_his_0021.nc'

Winter Test 
box_fn_in = '/Users/katehewett/Documents/LO_output/extract/cas7_t0_x4b/clines/shelf_box_2017.12.12_2017.12.12/temp_box_dir/box_000000.nc'
his_in = '/Users/katehewett/Documents/LO_roms/cas7_t0_x4b/f2017.12.12/ocean_his_0021.nc'

The history files are ~ random in that I had 4 on my personal computer and wanted one summer and one winter 

TODO: look at this in an estuary and in shallow water 

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
from math import log2
import pandas as pd 

######## ######## ######## 
# remove this testing after get the pycnocline to work here. 
# and replace with the arg pass tags in get_one_cline.py 

#plotting things
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import cmcrameri.cm as cmc

#box_fn_in = '/Users/katehewett/Documents/LO_output/extract/cas6_v0_live/clines/shelf_box_2022.08.08_2022.08.09/temp_box_dir/box_000001.nc'
#his_in = '/Users/katehewett/Documents/LO_roms/cas6_v0_live/f2022.08.08/ocean_his_0021.nc'

box_fn_in = '/Users/katehewett/Documents/LO_output/extract/cas7_t0_x4b/clines/shelf_box_2017.12.12_2017.12.12/temp_box_dir/box_000000.nc'
his_in = '/Users/katehewett/Documents/LO_roms/cas7_t0_x4b/f2017.12.12/ocean_his_0021.nc'
######## ######## ######## 

# 1. Get initial fields 
tt0 = time()
dsb = open_dataset(box_fn_in, decode_times=False) # TODO UPDATE
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat

# add z variables
# NOTE: this only works if you have salt as a saved field!
G, S, T = zrfun.get_basic_info(his_in)   # grid information to get S 
                                         # G is the whole domain. not the clipped box 

lat = dsb.lat_rho.values    
lon = dsb.lon_rho.values

# get indices for plotting // ilon and ilat will be the index of the mlon mlat inputs::
mlon = -124.1#-124.502
mlat = 45.135
#mlon = -124.1#.502
#mlat = 45.135
Lon = lon[0,:]
Lat = lat[:,0]
# error checking
if (mlon < Lon[0]) or (mlon > Lon[-1]):
    print('ERROR: lon out of bounds ')
    sys.exit()
if (mlat < Lat[0]) or (mlat > Lat[-1]):
    print('ERROR: lat out of bounds ')
    sys.exit()
# get indices
ilon = zfun.find_nearest_ind(Lon, mlon)
ilat = zfun.find_nearest_ind(Lat, mlat)

NT, NZ, NR, NC = np.shape(dsb.salt)
zeta = dsb.zeta.values
h = dsb.h.values

z_rho, z_w = zrfun.get_z(h, zeta, S) # (30, 1111, 356) (31, 1111, 356)
NW = z_w.shape[0]

#### isolate the one profile to test method #### 
z_w = z_w[:,ilat,ilon]
z_rho = z_rho[:,ilat,ilon]

tempC = dsb.temp.values[:,:,ilat,ilon]
SP = dsb.salt.values[:,:,ilat,ilon]

lat = lat[ilat,ilon]
lon = lon[ilat,ilon]

print('Time to get initial fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# 2. Calc SA and SIG0
tt0 = time()

P = gsw.p_from_z(z_rho,lat).squeeze()
SA = gsw.SA_from_SP(SP, P, lon, lat).squeeze()
CT = gsw.CT_from_pt(SA, tempC).squeeze()
SIG0 = gsw.sigma0(SA,CT).squeeze()   

fSIG0 = np.flipud(SIG0)     # flip so orientated surface to bottom 
fCT = np.flipud(CT)
fz_w = np.flipud(z_w)
fz_rho = np.flipud(z_rho)

print('Time to apply GSW = %0.2f sec' % (time()-tt0))
sys.stdout.flush()
     
# 3. Make some assumptions about the WC profiles 
tt0 = time()

# id min gradient above -3m (i think this step will be more important in obs data but not lp model output data)
dt = diff(fCT)
dz = diff(fz_rho)
dtdz = dt/dz

ztrim = -3
Gsurf_min = np.nanmin(dtdz[np.where(fz_w>=ztrim)])

iz1 = zfun.find_nearest_ind(dtdz, Gsurf_min)
z1 = fz_rho[iz1]  

# Temp variablity across the profile rep by 
Tmin = np.nanmin(CT)
Tmax = np.nanmax(CT)
T1 = fCT[iz1]

dTwc  = Tmax - Tmin      # this is the temp variabilty across the whole wc

T10 = Tmax - 0.1*dTwc 
T70 = Tmax - 0.7*dTwc

iT10 = zfun.find_nearest_ind(fCT, T10)   # grab the locations of the 10% and 70% del T whole WC
iT70 = zfun.find_nearest_ind(fCT, T70)

z01 = fz_rho[iT10]
z07 = fz_rho[iT70]
T01 = fCT[iT10]
T07 = fCT[iT70]

if z01 == z07:
    print('z10 equals z70 exit')
    sys.exit()
if z01<z07:
    print('weird depths')
    sys.exit()
    
# the vertical gradient will be strongest in the thermocline 
# and lowest in the surface mixed layer (and deep layer) 

# Main part of thermocline (next do pycnocline) 
# the data between (z0.1 and z0.7) is rearranged into [pi,i=0,1,2,...,I] 
# with po = p(z0.1) = D0.1 and pI = p(z0.7) = D0.7 
# And we calculate vertical gradients between pi(i=1,2,...,I) and po
zi = fz_rho[iT10:iT70+1]       # We flipped it! the data is packed surface to bottom
ti = fCT[iT10:iT70+1]
NGD = np.shape(zi)[0]          # number of data points in the zone bounded by z0.1 : z0.7

# if it is an easy profile, could calculate gradient like this:
GT17 = ((T01-T07)/(z01-z07))                 
# but if it's jumpy then may not be representative, and instead caculate like this:
# And recall, only doing the calc between i=1:I (we exclude i=0=po + remember we flipped data is packed surface on top)
GTi = ((T01-ti[1:NGD])/(z01-zi[1:NGD]))  

# I think we might need to caclulate the gradients like this (was it a typo?)
dtdz17 = diff(ti)/diff(zi)  

# The median of those gradients is used to represent a 
# characteristic gradient; pycnocline gradient                               
# If small, TD <= 0.001degC/m (GD <= ?? need threshold) then MLD extends to base of profile 
GT = np.nanmedian(GTi)

print('Time to calc characteristic gradient = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# 4. Exponential Leap-forward Gradient("ELG") method for MLD: 
# Let NG = the number of data points between z1 and z0.7 
# Let N = floor(log2(Ng)) which provides an N that is 
# much smaller than Ng which is req. for the ELG method
# A note on z1: NOt sure what depth is best to set this to
# 0-, 3-, or 10m?? 
# Let's set it, ztrim = 0 = z1, for now and explore a little.
# might have to treat seperate - 10m = z1 and then look at surface
# shallower than 10m ?? 

# initialize plot 
plt.close('all')
fs=12
plt.rc('font', size=fs)
fig = plt.figure(figsize=(16,8))

ax1 = plt.subplot2grid((1,3),(0,0),colspan=1,rowspan=1)
ax2 = plt.subplot2grid((1,3),(0,1),colspan=1,rowspan=1)
ax3 = plt.subplot2grid((1,3),(0,2),colspan=1,rowspan=1)

# plot T and SIG0
Tmin = np.floor(np.min(CT))-1
Tmax = np.floor(np.max(CT))+1
tticks = np.arange(Tmin,Tmax+1,2)

Tplot = ax1.plot(fCT,fz_rho,color='blue',marker='.',label = 'CT') 
ax1.set_ylabel('depth m')
ax1.set_xlabel('CT',color='blue')
ax1.tick_params(axis='x',color='blue')
ax1.xaxis.label.set_color('blue')
ax1.set_xlim([Tmin, Tmax])
ax1.set_xticks(tticks)
ax1.tick_params(labelcolor='blue')
#ax1.axhline(y = z_dtmax[:,:,ilon,ilat], color = 'k', linestyle = '--') 
ax1.axhline(z01)
ax1.axhline(z07)

axsig = ax1.twiny()

Splot = axsig.plot(fSIG0,fz_rho,color='red',marker='.',label = 'SIG0')
axsig.set_xlabel('SIG0',color='red')
axsig.tick_params(axis='x',color='red')
axsig.xaxis.label.set_color('red')
axsig.tick_params(labelcolor='red')
#axsig.axhline(y = z_dpmax[:,:,ir,ic], color = 'k', linestyle = '-') 

# plot dt/dz
pp = ax2.plot(dtdz,fz_w[1:NW-1],color='blue',marker='.',label = 'dpdz') 
ax2.set_ylabel('depth m')
ax2.set_xlabel('dT/dz',color='blue')
ax2.tick_params(axis='x',color='blue')
ax2.xaxis.label.set_color('blue')
#ax2.set_xlim([Tmin, Tmax])
#ax2.set_xticks(tticks)
ax2.tick_params(labelcolor='blue')
#ax2.plot(val_dpmax[:,:,ilat,ilon],z_dpmax[:,:,ilat,ilon],color='black',marker='x')
ax2.axhline(z01)
ax2.axhline(z07)

#gradient filter
ax3.axhline(z01)
ax3.axhline(z07)
#ax3.plot(GTi,zi[1:NGD],color='pink',marker='x',label = 'dpdz') 
ax3.plot(dtdz/GT,fz_w[1:NW-1],color='black',marker='x')
#ax3.plot(dtdz/np.nanmedian(dtdz17),fz_w[1:NW-1],color='green',marker='x')
ax3.axvline(0.5)

# plot dt/dz
#pt = ax3.plot(dtdz,z_w[1:NW-1],color='blue',marker='.',label = 'dpdz') 
#ax3.set_ylabel('depth m')
#ax3.set_xlabel('dT/dz',color='blue')
#ax3.tick_params(axis='x',color='blue')
#ax3.xaxis.label.set_color('blue')
#ax2.set_xlim([Tmin, Tmax])
#ax2.set_xticks(tticks)
ax3.tick_params(labelcolor='blue')
#ax3.plot(val_dtmax[:,:,ilat,ilon],z_dtmax[:,:,ilat,ilon],color='black',marker='x')



    



