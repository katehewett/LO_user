"""
Function to do the calculation for a single clipped file
calculates GSW vars > drho/dz dT/dz
finds max drho(dT)/dz + returns the max depth + values

THIS IS A TESTING FILE for threshold method to find SML 
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

box_fn_in = '/Users/katehewett/Documents/LO_output/extract/cas6_v0_live/clines/shelf_box_2022.08.08_2022.08.09/temp_box_dir/box_000001.nc'
his_in = '/Users/katehewett/Documents/LO_roms/cas6_v0_live/f2022.08.08/ocean_his_0021.nc'

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
#mlon = -124.502
#mlat = 45.135

mlon = -124.06
mlat = 45.135

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

SIGd = gsw.sigma0(SA,CT+0.8)

fSIG0 = np.flipud(SIG0)     # flip so orientated surface to bottom 
fCT = np.flipud(CT)
fz_w = np.flipud(z_w)
fz_rho = np.flipud(z_rho)

sS = fSIG0[0] #surface 
sTp = gsw.sigma0(SA[-1],CT[-1]+0.8)
sTm = gsw.sigma0(SA[-1],CT[-1]-0.8)

sm8 = fSIG0[0]-fSIG0

print('Time to apply GSW = %0.2f sec' % (time()-tt0))
sys.stdout.flush()
     
# 3. id pycnocline(thermocline) gradient
tt0 = time()   

# Vertical density difference rep total variabilty in profile , using D and lil p's for SIG
Dmin = np.nanmin(fSIG0)
Dmax = np.nanmax(fSIG0)
dDwc = Dmax - Dmin        
D10 = Dmin + 0.1*dDwc 
D70 = Dmin + 0.7*dDwc

iD10 = zfun.find_nearest_ind(fSIG0, D10)   # grab the locations of the 10% and 70% del SIG0
iD70 = zfun.find_nearest_ind(fSIG0, D70)

z01 = fz_rho[iD10]
z07 = fz_rho[iD70]
D01 = fSIG0[iD10]
D07 = fSIG0[iD70]

if z01 == z07:
    print('z10 equals z70 exit')
    sys.exit()
if z01<z07:
    print('weird depths')
    sys.exit()

# Main part of pycnocline(thermocline) 
# the data between (z0.1 and z0.7) is rearranged into [pi,i=0,1,2,...,I] 
# with po = p(z0.1) = D0.1 and pI = p(z0.7) = D0.7 
# And we calculate vertical gradients between pi(i=1,2,...,I) and po
zi = fz_rho[iD10:iD70+1]       # We flipped it! the data is packed surface to bottom
pi = fSIG0[iD10:iD70+1]
NGD = np.shape(zi)[0]          # number of data points in the zone bounded by z0.1 : z0.7

# if it is an easy profile, could calculate gradient like this:
GD17 = -((D01-D07)/(z01-z07))                 
# but if it's jumpy then may not be representative, and instead caculate like this:
# And recall, only doing the calc between i=1:I (we exclude i=0=po + remember we flipped data is packed surface on top)
GDi = -((D01-pi[1:NGD])/(z01-zi[1:NGD]))  
  
# I think we should actually caclulate the gradients like this:
dpdz17 = -diff(pi)/diff(zi)  

# The median of those gradients is used to represent a 
# characteristic gradient; pycnocline gradient                               
# If small, TD <= 0.001degC/m (GD <= ?? need threshold) then MLD extends to base of profile 
GD = np.nanmedian(GDi)


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

# DATA IS FLIPPED : packed surface to bottom 

tt0 = time()  

if max(z_rho)>0:
    print('z_rho[-1]<0')
    
ztrim = 0
# need an if statement to flag if water depth < ztrim (when 10m)
# need an if statement to flag if ztrim is above the surface 

iz = zfun.find_nearest_ind(fz_rho, ztrim)
z1 = fz_rho[iz]

zk = fz_rho[iz:iD70+1] 
pk = fSIG0[iz:iD70+1]

NG = np.shape(zk)[0]
N = np.floor(log2(NZ-1-NG-1)).astype('int64') # if NG = 20; N = 4

DnD = np.nan * np.ones((NZ))   # initialize 
for k in range(NG-1): 
    print(k)
    for ni in range(N-1): 
        print(ni)
        n = np.arange(0,N,1)
        tn = 2**n
        k2 = k+tn
    
        DnD[k]=np.nanmean((fSIG0[k] - fSIG0[k+k2])/(fz_rho[k]-fz_rho[k+k2]))

# 
dp = diff(fSIG0)
dz = diff(fz_rho)
dT = diff(fCT)
dpdz = dp/dz
dtdz = dT/dz

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

axsig = ax1.twiny()
SIGmin = np.floor(np.min(SIG0))-0.5
SIGmax = np.floor(np.max(SIG0))+0.5
sigticks = np.arange(SIGmin,SIGmax+1,1)

Splot = axsig.plot(fSIG0,fz_rho,color='red',marker='.',label = 'SIG0')
axsig.set_xlabel('SIG0',color='red')
axsig.tick_params(axis='x',color='red')
axsig.xaxis.label.set_color('red')
axsig.set_xlim([SIGmin, SIGmax])
axsig.set_xticks(sigticks)
axsig.tick_params(labelcolor='red')
#axsig.axhline(y = z_dpmax[:,:,ir,ic], color = 'k', linestyle = '-') 
plt.plot(D01,z01,'o')
plt.plot(D07,z07,'o')
plt.axhline(z01)
plt.axhline(z07)
axsig.set_ylim([np.nanmin(fz_rho)-1, np.nanmax(fz_rho)+1])

# plot dp/dz
#pp = ax2.plot(-dpdz,fz_w[1:NW-1],color='red',marker='.',label = 'dpdz') 
ax2.plot(-sm8,fz_rho,color='orange',marker='.')
ax2.set_ylabel('depth m')
ax2.set_xlabel('drho/dz',color='red')
ax2.tick_params(axis='x',color='red')
ax2.xaxis.label.set_color('red')
#ax2.set_xlim([Tmin, Tmax])
#ax2.set_xticks(tticks)
ax2.tick_params(labelcolor='red')
#ax2.plot(val_dpmax[:,:,ilat,ilon],z_dpmax[:,:,ilat,ilon],color='black',marker='x')
#ax2.axhline(z01)
#ax2.axhline(z07)
ax2.set_ylim([np.nanmin(fz_rho)-1, np.nanmax(fz_rho)+1])
#ax2.axvline(x = 0.001*1.1, color = 'k', linestyle = '-') 
#ax2.axvline(x = 0.002*1.1, color = 'k', linestyle = '-') 
ax2.axvline(x = 0.03)
ax2.axvline(x = 0.05)
ax2.axvline(x = 0.08)

#gradient filter
ax3.axhline(z01)
ax3.axhline(z07)
ax3.plot(-dpdz/GD,fz_w[1:NW-1],color='black',marker='x',label = 'dpdz') 
ax3.plot(-DnD,fz_rho,color='green',marker='x',label = 'dpdz') 
ax3.axvline(0.1)
ax3.set_ylim([np.nanmin(fz_rho)-1, np.nanmax(fz_rho)+1])
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



    



