"""
Function to do the calculation for a single file sent
by extract_box_test.py 
This script: 
(1) finds the layers within some threshold of the bottom. 
Made for some DO work with Jack Barth; 
beta for surface / bot extractions for indicators 
(2) calculates:
(1) GSW vars (2) SML using a threshold method 
(3) the average values of t s alk oxygen n03 nh4 tic in the sml
(4) values are saved in a temp file and then concat'd in extract_clines

Time to get initial fields = 3.32 sec
Time to apply GSW = 1.26 sec
Time to calc SML + assign T/S +CHEM = 1.09 sec

TODO make a function for the masked data flag when finding mean
make flags for pycnocline and thermocline and pass threshold too 

Right now using 0.125 kg/m3 = thresh (previously 0.03 kg/m3)

"""

import numpy as np
import sys

from argparse import ArgumentParser
from lo_tools import zrfun
from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
from time import time
from lo_tools import Lfun, zrfun
import warnings

import gsw 

import matplotlib.pyplot as plt
from lo_tools import plotting_functions as pfun
from matplotlib.patches import Polygon
from cmcrameri import cm
from matplotlib.colors import BoundaryNorm

'''parser2 = ArgumentParser()
parser2.add_argument('-in_fn', '--input_box_files', type=str) # path to clipped box files (temp directory)
parser2.add_argument('-out_fn', '--output_cline_files', type=str) # path to save cline files (temp directory)
parser2.add_argument('-his_fn', '--input_his_file', type=str) # path to one history file for grid info > z's
parser2.add_argument('-lt', '--list_type', type=str) # list type: hourly, daily, weekly, lowpass   
parser2.add_argument('thresh','--threshold', type=int, default = 15) # threshold (distance in meters): default 15m    
args2 = parser2.parse_args()


##########################################################
threshold = args.threshold
##########################################################

tt0 = time()
##########################################################
# STEP 1: setup and get initial fields
##########################################################
dsb = open_dataset(args2.input_box_files, decode_times=False)
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat'''

## TESTING START 
tt0 = time() 

threshold = 15
Ldir = Lfun.Lstart()

# load regions and no data list for that year and make a plot
in_dir = Ldir['roms_out'] / 'cas7_t1_x11ab' / ('f2025.08.11')
fn_in = in_dir / 'ocean_avg_0001.nc'

out_dir = Ldir['parent'] / 'LKH_output' / 'pnw_oxygen_maps' / 'early_eval' 
Lfun.make_dir(out_dir, clean=False)

dsb = open_dataset(fn_in, decode_times=False)

# add z variables
# NOTE: this only works if you have salt as a saved field!
#G, S, T = zrfun.get_basic_info(args2.input_his_file)   # grid information 
G, S, T = zrfun.get_basic_info(fn_in)   # grid information 

lat = np.nanmean(G['lat_rho'])
lon = np.nanmean(G['lon_rho']) # for gsw calcs

NT, NZ, NR, NC = np.shape(dsb.salt)
h = dsb.h.values
hthresh = dsb.h.values + threshold
zeta = dsb.zeta.values
mask_rho = dsb.mask_rho

z_rho, z_w = zrfun.get_z(h, zeta, S) # (30, 1111, 356) (31, 1111, 356)
NW = z_w.shape[0]

z_rho = np.expand_dims(z_rho,axis=0) # a little unnecessary step but helps make it easier below 
z_w = np.expand_dims(z_w,axis=0)

tempC = dsb.temp.values # (1, 30, 1111, 356)
SP = dsb.salt.values
OXY = dsb.oxygen.values

print('Time to get initial fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

'''##########################################################
# STEP 2. Calc SA and SIG0: there's a holdup with sending gsw LARGE files - 
# it can accomodate 4D variables, but has an upper limit 
# and can be slow. Not a problem with NT = 1 tho if change to chunks 2xcheck
##########################################################
tt0 = time()

P = gsw.p_from_z(z_rho,lat)
SA = gsw.SA_from_SP(SP, P, lon, lat)
CT = gsw.CT_from_pt(SA, tempC)
SIG0 = gsw.sigma0(SA,CT)   

print('Time to apply GSW = %0.2f sec' % (time()-tt0))
sys.stdout.flush()'''

##########################################################     
# STEP 2: calc threshold layer 
##########################################################
tt0 = time()       
plt.close('all')
fs=11
plt.rc('font', size=fs)
height_of_image = 15 
width_of_image = 8 
fig1 = plt.figure(figsize=(height_of_image,width_of_image))
fig1.set_size_inches(height_of_image,width_of_image, forward=False)
frac=0.047 * (height_of_image / (width_of_image/13))
fig1.tight_layout()

xrho = dsb.lon_rho
yrho = dsb.lat_rho
mask = dsb.mask_rho #land 0 water 1

cmap1=cm.roma
levels = np.arange(0, 15.5, 0.5)
norm1 = BoundaryNorm(levels, ncolors=cmap1.N, clip=False)
cmap1.set_extremes(over = 'LightPink',under='black')

dbot = np.squeeze(z_w[:,1,:,:]-z_w[:,0,:,:])

# shelf box
pwest = -126
peast = -123
pnorth = 48.75
psouth = 42
vertices = [(pwest, psouth), (pwest, pnorth), (peast, pnorth), (peast, psouth)]
polygon = Polygon(vertices, closed=True, facecolor='none', edgecolor='peru', linewidth=2, linestyle = '-')

ax = plt.subplot2grid((2,2), (0,0), colspan=1, rowspan=2)
pfun.add_coast(ax)
pfun.dar(ax)

pcm = ax.pcolormesh(xrho,yrho,dbot,cmap=cmap1, norm = norm1) 

ax.axis([-130, -122, 42, 52])
ax.add_patch(polygon)
#ax.plot(x, y, marker='o', ms=5, mec='crimson', mfc='lightpink', linestyle = 'none', label = 'NDBC location')
ax.set_title('Main LiveOcean domain')

ax2 = plt.subplot2grid((2,2), (0,1), colspan=2, rowspan=2)
pfun.add_coast(ax2)
pfun.dar(ax2)
ax2.axis([pwest, peast, psouth, pnorth])
#ax2.plot(x, y, marker='o', ms=5, mec='crimson', mfc='lightpink', linestyle = 'none', label = 'NDBC location')
ax2.set_title('Shelf Box', color = 'peru', weight = 'bold')

pcm2 = ax2.pcolormesh(xrho,yrho,dbot,cmap=cmap1, norm = norm1) 

ax2.contour(xrho,yrho,h, [50],
colors=['white'], linewidths=1, linestyles='solid',alpha=0.4)
ax2.contour(xrho,yrho,h, [200],
colors=['black'], linewidths=1, linestyles='solid',alpha=0.4)

cbar2 = fig1.colorbar(pcm2, extend = 'both')
cbar2 = fig1.colorbar(pcm, extend = 'both')

sys.exit() 

##########################################################     
# STEP 2: calc threshold layer 
##########################################################
tt0 = time()       

Dsurf = SIG0[:,-1,:,:]
dthresh = np.abs(h - h15)

# flip along z axis so that ocean packed surface:bottom
# find where the first value passes the threshold ("B")
fSIG0 = np.flip(SIG0, axis=1)
fCT = np.flip(CT, axis=1)
fSA = np.flip(SA, axis=1)
fALK = np.flip(ALK, axis=1)
fTIC = np.flip(TIC, axis=1)
fOXY = np.flip(OXY, axis=1)
fNO3 = np.flip(NO3, axis=1)
fNH4 = np.flip(NH4, axis=1)
fsthresh = np.flip(sthresh, axis=1)
fz_rho = np.flip(z_rho, axis=1)
fz_w = np.flip(z_w, axis=1)
fdz = -diff(fz_w,axis=1) # 1, 30, 1302, 663 LO domain 

# the flag returned by argmin acts on SIG0, which is on 
# z_rho. The base of the sml, with our thresh, occurs somewhere 
# between those two z_rho points at B and B-1. 
# Assign the midpoint fz_w[:,B,:,:] as the base. 
# Note: if sml is found to be zero, the depth is saved as 
# the first z_w point in the flipped "fz_w" array. 
# So when we calc the thickness, zeta - SML base, 
# we get zero for the SML thickness.
B = np.argmin(fsthresh<threshold,axis=1,keepdims=True)

# Add on after looking at data: there are some cases where 
# the whole w/c is w/in and under the threshold (0.125 kg/m3). 
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    fmax = np.nanmax(fsthresh,axis=1,keepdims=True)
B[fmax<threshold] = 30

SML_z = np.take_along_axis(fz_w,B,axis=1)

##########################################################     
# STEP 4: Find weighted means of vars in SML
# if SML is zero and in water keep SIG0 value just as the the surface 
# value else if in water and non zero we calc the avg SIG0 for the SML
# note: when SML = 0; we still save the surface CT(SA); and save the land as NaNs
##########################################################     
# save mean value of T and S in the SML; 1st initialize: 
# TODO there is a better way to do this, being lazy and writing everything out
SML_CT = np.expand_dims(fCT[:,0,:,:],axis=1)    
SML_SA = np.expand_dims(fSA[:,0,:,:],axis=1)

SML_ALK = np.expand_dims(fALK[:,0,:,:],axis=1)
SML_TIC = np.expand_dims(fTIC[:,0,:,:],axis=1)
SML_NH4 = np.expand_dims(fNH4[:,0,:,:],axis=1)
SML_NO3 = np.expand_dims(fNO3[:,0,:,:],axis=1)
SML_OXY = np.expand_dims(fOXY[:,0,:,:],axis=1)

bmask = (B.squeeze()!=0) & (mask_rho.values == 1)  # will use at end of this step
broadcast_mask_rho = (mask_rho.values == 1)[np.newaxis, np.newaxis, :, :]

# below the SML set values to NaN  & find the weighted mean CT(SA) in the SML      
# calc weights first for averaging 
condition = fz_rho<SML_z
dz_masked = np.ma.masked_array(fdz,condition)
SML_thickness = np.ma.sum(dz_masked, axis=1, keepdims=True)
weights_masked = dz_masked/SML_thickness

# retain values above the SML; values below will be nans 
T_masked = np.ma.masked_array(fCT,condition)
S_masked = np.ma.masked_array(fSA,condition)

ALK_masked = np.ma.masked_array(fALK,condition) # PLUS
TIC_masked = np.ma.masked_array(fTIC,condition)
NO3_masked = np.ma.masked_array(fNO3,condition)
NH4_masked = np.ma.masked_array(fNH4,condition)
OXY_masked = np.ma.masked_array(fOXY,condition)

# inside takes a weighted average of each var, outside extracts 
mCT_SML = np.ma.getdata(np.ma.average(T_masked, weights = weights_masked, axis=1, keepdims=True))
mSA_SML = np.ma.getdata(np.ma.average(S_masked, weights = weights_masked, axis=1, keepdims=True))

mALK_SML = np.ma.getdata(np.ma.average(ALK_masked, weights = weights_masked, axis=1, keepdims=True)) # PLUS
mTIC_SML = np.ma.getdata(np.ma.average(TIC_masked, weights = weights_masked, axis=1, keepdims=True))
mNO3_SML = np.ma.getdata(np.ma.average(NO3_masked, weights = weights_masked, axis=1, keepdims=True))
mNH4_SML = np.ma.getdata(np.ma.average(NH4_masked, weights = weights_masked, axis=1, keepdims=True))
mOXY_SML = np.ma.getdata(np.ma.average(OXY_masked, weights = weights_masked, axis=1, keepdims=True))

# replace non zeros with data we just calc'd 
SML_CT[broadcast_mask_rho] = mCT_SML[broadcast_mask_rho]
SML_SA[broadcast_mask_rho] = mSA_SML[broadcast_mask_rho]
SML_ALK[broadcast_mask_rho] = mALK_SML[broadcast_mask_rho]
SML_TIC[broadcast_mask_rho] = mTIC_SML[broadcast_mask_rho]
SML_NO3[broadcast_mask_rho] = mNO3_SML[broadcast_mask_rho]
SML_NH4[broadcast_mask_rho] = mNH4_SML[broadcast_mask_rho]
SML_OXY[broadcast_mask_rho] = mOXY_SML[broadcast_mask_rho]

SMLT = np.ma.getdata(SML_thickness)

print('Time to calc SML + assign T/S +CHEM = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

#save some values 
tt0 = time()

ds1 = Dataset()
ds1['ocean_time'] = dsb.ocean_time

ds1['zeta'] = (('ocean_time', 'eta_rho', 'xi_rho'), dsb.zeta.values, {'units':'m', 'long_name':'free-surface'})
ds1['zSML'] = (('ocean_time', 'eta_rho', 'xi_rho'), np.squeeze(SML_z,axis=1), {'units':'m', 'long_name': str('depth of bottom with SML threshold = ' + str(threshold))})
ds1['SML_thickness'] = (('ocean_time', 'eta_rho', 'xi_rho'), np.squeeze(SMLT,axis=1), {'units':'m', 'long_name': str('SML thickness with threshold = ' + str(threshold))})

ds1['SA'] = (('ocean_time', 'eta_rho', 'xi_rho'), np.squeeze(SML_SA,axis=1), {'units':'g/kg', 'long_name': 'mean Absolute Salinity'})
ds1['CT'] = (('ocean_time', 'eta_rho', 'xi_rho'), np.squeeze(SML_CT,axis=1), {'units':'degrees celcius', 'long_name': 'mean Conservative Temperature'})

ds1['ALK'] = (('ocean_time', 'eta_rho', 'xi_rho'), np.squeeze(SML_ALK,axis=1), {'standard_name':'sea_water_alkalinity_expressed_as_mole_equivalent','long_name': 'total alkalinity','units':'milliequivalents meter-3'})
ds1['TIC'] = (('ocean_time', 'eta_rho', 'xi_rho'), np.squeeze(SML_TIC,axis=1), {'standard_name':'mole_concentration_of_total_inorganic_carbon_in_sea_water','long_name': 'total inorganic carbon','units':'millimole_carbon meter-3'})
ds1['NH4'] = (('ocean_time', 'eta_rho', 'xi_rho'), np.squeeze(SML_NH4,axis=1), {'standard_name':'mole_concentration_of_ammonium_expressed_as_nitrogen_in_sea_water','long_name': 'ammonium concentration','units':'millimole_nitrogen meter-3'})
ds1['NO3'] = (('ocean_time', 'eta_rho', 'xi_rho'), np.squeeze(SML_NO3,axis=1), {'standard_name':'mole_concentration_of_nitrate_expressed_as_nitrogen_in_seawater','long_name': 'nitrate concentration','units':'millimole_nitrogen meter-3'})
ds1['oxygen'] = (('ocean_time', 'eta_rho', 'xi_rho'), np.squeeze(SML_OXY,axis=1), {'standard_name':'mole_concentration_of_dissolved_oxygen_in_sea_water','long_name': 'dissolved oxygen concentration','units':'millimole_oxygen meter-3'})

ds1.to_netcdf(args2.output_cline_files, unlimited_dims='ocean_time')

print('saved!')