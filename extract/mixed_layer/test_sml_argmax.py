"""
Code to test the speed and memory use of one sml calculation for the LO domain
using argmax to find indicies 

running testing code with two history files on kh personal mac
run test_sml_whole -gtx cas6_v0_live -ro 1 -0 2022.08.08 -1 2022.08.09

TOdos: (1) this version is using 0 for SSH, and I don't know if that is best; 
see line in step #1:
z_rho, z_w = zrfun.get_z(h, 0*h, S) # [can we use 0 for SSH??]
(2) explore thresh and grad thresh with observations woib
(3) add testing line 
(4) place function outside this code consecutive. and is there a better option??

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

gridyah = True
if gridyah == True:
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
    dzr = np.diff(z_w, axis=0)          # vertical thickness of all cell [m]

    # 2. set vars; just one history file; squeezing time dimen
    temp = ds.temp.values.squeeze()  
    SP = ds.salt.values.squeeze()
    P = gsw.p_from_z(z_rho,lat)
    SA = gsw.SA_from_SP(SP, P, lon, lat)
    CT = gsw.CT_from_pt(SA, temp)
    SIG0 = gsw.sigma0(SA,CT)
    mask_rho = ds.mask_rho 

    NZ, NETA, NXI = CT.shape
    # NZ = np.shape(CT)[0]
    # NETA = np.shape(CT)[1]
    # NXI = np.shape(CT)[2] 

    print('Time to get initial fields = %0.2f sec' % (time()-tt0))
    sys.stdout.flush()

# 3. Calculate a map of the thickness of SML using temp and sig0 thresholds
# notes: Using 0.05 deg C (lentz 92) and 0.01 kg/m3 (standard?) as thresholds 
# for SML calc. There are other methods too. Using these for now and then will 
# compare w/ obs (to-do) Using gradient diff or split/merge might be excessive? 

tt0 = time()

# flipud so we find the first instance for the SML false 
SIGf = np.flipud(SIG0)
Tf = np.flipud(CT)
dzrf = np.flipud(dzr)
z_wf = np.flipud(z_w)
z_rhof = np.flipud(z_rho)

sig0s = SIGf[0,:,:]                      
t0s = Tf[0,:,:]

dSIG_s = np.round(SIGf-sig0s,2)
dT_s = np.round(Tf-t0s,2)

Tbool = dT_s>-0.05
Sbool = dSIG_s<0.02

A = np.argmax(Tbool==False,axis=0)    # grab the first index where false 
B = np.argmax(Sbool==False,axis=0)

print('Time to get initial SML = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

yplotting = True 
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
    Tplot = plt.gcf().axes[axnum].plot(CT[:,300,300].squeeze(),z_rho[:,300,300].squeeze(),color='pink',marker='x') 
    Tplot = plt.gcf().axes[axnum].plot(Tf[:,300,300].squeeze(),z_rhof[:,300,300].squeeze(),color='blue',marker='o')
    plt.gcf().axes[axnum].plot(Tf[A[300,300],300,300].squeeze(),z_rhof[A[300,300],300,300],color='black',marker='*')
    plt.gcf().axes[axnum].set_title('CT')
    plt.gcf().axes[axnum].set_ylabel('depth m')

    axnum = 1
    rhoplot = plt.gcf().axes[axnum].plot(SIG0[:,300,300].squeeze(),z_rho[:,300,300].squeeze(),color='pink',marker='x') 
    rhoplot = plt.gcf().axes[axnum].plot(SIGf[:,300,300].squeeze(),z_rhof[:,300,300].squeeze(),color='blue',marker='o')
    plt.gcf().axes[axnum].plot(SIGf[B[300,300],300,300].squeeze(),z_rhof[B[300,300],300,300],color='black',marker='*') 
    plt.gcf().axes[axnum].set_title('SIG0')
    



 