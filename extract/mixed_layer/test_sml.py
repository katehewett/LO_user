"""
Code to test the speed and memory use of one sml calculation for the LO domain

Quesiton: this version is using 0 for SSH, and I don't know if that is best; 
see line in step #1:
z_rho, z_w = zrfun.get_z(h, 0*h, S) # [can we use 0 for SSH??]

"""
import xarray as xr
import numpy as np
import gsw
from lo_tools import Lfun, zrfun
from time import time
import sys
import gsw

#plotting things
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import cmcrameri.cm as cmc

Ldir = Lfun.Lstart()
fn = Ldir['roms_out'] / 'cas6_v0_live' / 'f2022.08.08' / 'ocean_his_0021.nc'

yplotting = True 

# 1. Get some grid information
tt0 = time()
ds = xr.open_dataset(fn, decode_times=False)
G, S, T = zrfun.get_basic_info(fn)
lon = ds.lon_rho.values
lat = ds.lat_rho.values
DA = G['DX'] * G['DY']
# we pack some things in the CC dict that will eventually be written out to a
# NetCDF file.
CC = dict()
CC['lon'] = lon
CC['lat'] = lat
CC['DA'] = DA
h = ds.h.values
CC['h'] = h
mask_rho = ds.mask_rho.values
CC['mask_rho']=mask_rho
z_rho, z_w = zrfun.get_z(h, 0*h, S) # [can we use 0 for SSH??]

dzr = np.diff(z_w, axis=0) # vertical thickness of all celle [m]
print('Time to get initial fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# 2. naming variables and calc SA; CT; density 
tt0 = time()

temp = ds.temp.values.squeeze() # just one history file; squeezing time dimen
SP = ds.salt.values.squeeze()
P = gsw.p_from_z(z_rho,lat)
SA = gsw.SA_from_SP(SP, P, lon, lat)
CT = gsw.CT_from_pt(SA, temp)
SIG0 = gsw.sigma0(SA,CT)

# 3. Calculate a map of the thickness of SML using
# notes:
# Temp threshold N.CA/OR/Peru/NW Africa used 0.02(0.05) degC (lentz 92)
# threshold difference = 0.01 kg/m3 (defacto standard in oceans?)
# threshold gradient = drho / dz <= 0.01 kg/m4 (?)
# split merge methods et al. 
# for now we are going to use a threshold difference of 0.02 degC
# because there are more temperature data near surface to check this with obs
# although most casts also have S so can check using density too 

#dz = np.diff(z_rho,axis=1)
#drho = np.diff(SIG0,axis=1)
#dT = np.diff(CT,axis=1)
#dS = np.diff(SA,axis=1)
#dpdz = drho/dz


# grabbing a random location on shelf  (300 300 deeper ; 360 360 shallower)
ts = CT[:,300,300].squeeze()
sigs = SIG0[:,300,300].squeeze()
zs = z_rho[:,300,300].squeeze()

NZ = np.shape(ts)[0]
t0s = ts[NZ-1]
t0b = ts[0]
sig0s = sigs[NZ-1]
sig0b = sigs[0]

dT_s = np.round(ts-t0s,2)
dT_b = np.round(ts-t0b,2)

dSIG_s = np.round(sigs-sig0s,2)
dSIG_b = np.round(sigs-sig0b,2)

# define function to find consecutive true statements in a list 
# put this in a seperate file that we can call 

def find_consecutive_true(list1):
    """
    Finds consecutive true statements in a list.

    Args:
        list1 (list): A list of boolean values.

    Returns:
        list: A list of lists, where each sublist contains the indices of the
        consecutive true statements in the original list.
    """

    result = []
    start = 0
    for i in range(len(list1)):
        if list1[i]:
            if i == start:
                result.append([i])
            else:
                result[-1].append(i)
        else:
            start = i + 1
    return result

if dT_s[NZ-1]<-0.05:
    print('sml 0')
else:
    A = dT_s>=-0.05 #boolean for temp threshold 
    if np.any(A): # there are true values 
        # list1 = np.where(A==True) 
        AC = find_consecutive_true(A==True)
        zidx=np.min(AC[-1])
        zSML = zs[zidx]
        
              
if yplotting==True:
    # initialize plot 
    plt.close('all')
    fs=12
    plt.rc('font', size=fs)
    fig = plt.figure(figsize=(16,8))

    ax1 = plt.subplot2grid((1,4),(0,0),colspan=1,rowspan=1)
    ax2 = plt.subplot2grid((1,4),(0,1),colspan=1,rowspan=1)
    ax3 = plt.subplot2grid((1,4),(0,2),colspan=1,rowspan=1)
    ax4 = plt.subplot2grid((1,4),(0,3),colspan=1,rowspan=1)
    
    axnum = 0 
    tplot = plt.gcf().axes[axnum].plot(ts,zs,color='blue',marker='x') 
    plt.gcf().axes[axnum].set_title('CT')
    plt.gcf().axes[axnum].set_ylabel('depth m')

    axnum = 1
    rhoplot = plt.gcf().axes[axnum].plot(sigs,zs,color='pink',marker='x') 
    plt.gcf().axes[axnum].plot(sigs[zidx],zs[zidx],color='red',marker='o') 
    plt.gcf().axes[axnum].set_title('SIG0')
    
    axnum = 2 
    plt.gcf().axes[axnum].plot(dT_s,zs,color='black',marker='x')
    plt.gcf().axes[axnum].set_title('diff T surf > bot')
    plt.gcf().axes[axnum].axvline(x = -0.02, color = 'r', label = 'axvline - full height')
    plt.gcf().axes[axnum].axvline(x = -0.05, color = 'tab:pink', label = 'axvline - full height')
    
    axnum = 3 
    plt.gcf().axes[axnum].plot(dSIG_s,zs,color='black',marker='x')
    plt.gcf().axes[axnum].set_title('diff sig0 surf > bot')
    plt.gcf().axes[axnum].axvline(x = 0.01, color = 'r', label = 'axvline - full height')
    



