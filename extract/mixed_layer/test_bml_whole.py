"""
Code to test the speed and memory use of one sml  calculation for the LO domain

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

sig0s = SIG0[NZ-1,:,:]                      
t0s = CT[NZ-1,:,:]

dSIG_s = np.round(SIG0-sig0s,2)
dT_s = np.round(CT-t0s,2)

Tdzrs = np.ma.masked_where(dT_s<-0.05,dzr)
Tdz = Tdzrs.sum(axis=0)                       # SML thickness temp threshold 

Sdzrs = np.ma.masked_where(dSIG_s>0.01,dzr)
SIGdz = Sdzrs.sum(axis=0)                     # SML thickness sigma threshold 

Sdzrs = dSIG_s>0.01

print('Time to get initial SML = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# 4. function(s) for mixed layer calcs 

tt0 = time()

def find_consecutive_true(list1):
    """
    Finds consecutive true statements in a list.
    Args: A list of boolean values.
    Returns: A list of lists, where each sublist contains the indices of the
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
 
print('time to build function = %0.2f sec' % (time()-tt0))
sys.stdout.flush()
    
# 5. Check for locations that threshold was met then skipped a cell and 
#    then was met again (e.g. true true false true sums to 3 and not 2 (or 4) 
#    comes up more with sig0 

tt0 = time()

#x = np.zeros_like(CT[1,300:350,300:350]) #initialize and below smaller testing size 30 x 50 x 50 
#bSdzrm = Sdzrs[:,300:350,300:350]
#NETA = np.shape(x)[0]
#NXI = np.shape(x)[1] 

#initialize flag arrays 
xs = np.ones_like(mask_rho.values)*mask_rho.values 
xt = np.ones_like(mask_rho.values)*mask_rho.values 
xf = np.ones_like(mask_rho.values)*mask_rho.values

# Initialize array (y,x) for SML calcs 
zSML_sig = np.nan * np.ones((NETA,NXI))   
dSML_sig = np.nan * np.ones((NETA,NXI)) 

zSML_T = np.nan * np.ones((NETA,NXI)) 
dSML_T = np.nan * np.ones((NETA,NXI)) 

SML_sig = np.nan * np.ones((NETA,NXI)) 
SML_T = np.nan * np.ones((NETA,NXI)) 


exit
for i in range(NETA-1):
    for j in range(NXI-1):
        if mask_rho[i,j]==1:
            print(i)
            print(j)
            AS = find_consecutive_true(Sdzrs.mask[:,i,j]==False)
            xs[i,j] = len(AS) # if x > 1 then need to look at sml calc again   
            if xs[i,j]>1:
                dSML_sig[i,j] = sum(Sdzrs[:,i,j].data[AS[-1]])     
                SML_sig[i,j] = SIG0[np.min(AS[-1]),i,j]
                #zSML_sig[i,j] = Sdzrs[:,i,j].data[AS[-1]]
        
                AT = find_consecutive_true(Tdzrs.mask[:,i,j]==False)
                xt[i,j] = len(AT) # if x > 1 then need to look closer 
                if xt[i,j]>1:
                    dSML_T[i,j] = sum(Tdzrs[:,i,j].data[AT[-1]])  
                    SML_T[i,j] = CT[np.min(AT[-1]),i,j]
         
print('time to id flags = %0.2f sec' % (time()-tt0))
sys.stdout.flush()    

itesting = False
if itesting == True:
    ts = CT[:,300,300].squeeze()
    sigs = SIG0[:,300,300].squeeze()
    zs = z_rho[:,300,300].squeeze()
    
    dSIG_s = dSIG_s[:,300,300].squeeze()
    dT_s = dT_s[:,300,300].squeeze()

    Tdz = Tdz[300,300] 

    SIGdz = SIGdz[300,300] 
    
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
    #plt.gcf().axes[axnum].plot(SML_T[300,300],zSML_T[300,300],color='red',marker='o') 

    axnum = 1
    rhoplot = plt.gcf().axes[axnum].plot(sigs,zs,color='pink',marker='x') 
    #plt.gcf().axes[axnum].plot(SML_sig[300,300],zSML_sig[300,300],color='red',marker='o') 
    plt.gcf().axes[axnum].set_title('SIG0')
    
    axnum = 2 
    plt.gcf().axes[axnum].plot(dT_s,zs,color='black',marker='x')
    plt.gcf().axes[axnum].set_title('diff T surf')
    plt.gcf().axes[axnum].axvline(x = -0.05, color = 'black', label = 'axvline - full height')
    plt.gcf().axes[axnum].axvline(y = -dSML_T[300,300], color = 'black', label = 'axvline - full height')
        
    axnum = 3 
    plt.gcf().axes[axnum].plot(dSIG_s,zs,color='black',marker='x')
    plt.gcf().axes[axnum].set_title('diff sig0 surf')
    plt.gcf().axes[axnum].axvline(x = 0.01, color = 'black', label = 'axvline - full height')
    plt.gcf().axes[axnum].axvline(y = -dSML_sig[300,300], color = 'black', label = 'axvline - full height')
     

 

        
