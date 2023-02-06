"""
This is code for doing geting hypoxic volume.
Can I make it more generic in the next step to 
have flags for different goals: e.g. anoxic vol, 
arag saturation state? 

This is based on the old code from fjord that calculates 
hypoxic volume. 

Need to have run extract/box code first
Testing: Jan 2023

"""

# imports
# took these from extract_moor.py, extract_sections.py has a different way to deal with arg passing, but not sure which is best yet?
import sys
import argparse
from lo_tools import Lfun, zfun, zrfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import os
from time import time
import numpy as np
import xarray as xr

from datetime import datetime, timedelta
start_time = datetime.now()

pid = os.getpid()
print(' extract_box '.center(60,'='))
print('PID for this job = ' + str(pid))

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_l08b
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
parser.add_argument('-lt', '--list_type', type=str) # list type: hourly, daily, weekly

# Later can we add this to make the box smaller and add in arag sat state etc?
# select job name
#parser.add_argument('-job', type=str) # job name 

# Optional: set max number of subprocesses to run at any time
parser.add_argument('-Nproc', type=int, default=10)
# Optional: for testing
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided
argsd = args.__dict__
for a in ['gtagex']:
    if argsd[a] == None:
        print('*** Missing required argument: ' + a)
        sys.exit()
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]

# testing
if Ldir['testing']:
    Ldir['roms_out_num'] = 2
    Ldir['ds0'] = '2019.07.04'
    Ldir['ds1'] = '2019.07.06'
    Ldir['list_type'] = 'daily'
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]

# output location
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'hypoxic_volume'
Lfun.make_dir(out_dir)
# if add job names, then add this:
# bname = Ldir['job'] + '_' + Ldir['ds0'] + '_' + Ldir['ds1']
bname = 'HypoxicVolume_' + Ldir['ds0'] + '_' + Ldir['ds1']
out_fn = out_dir / (bname + '.nc')
out_fn.unlink(missing_ok=True)

# name the temp dir to accumulate individual extractions
temp_dir = out_dir / ('temp_' + bname)
Lfun.make_dir(temp_dir, clean=True)
    
# get list of history files to plot
fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])
#if Ldir['testing']:
#    fn_list = fn_list[:5]
#G, S, T = zrfun.get_basic_info(fn_list[0])
#Lon = G['lon_rho'][0,:]
#Lat = G['lat_rho'][:,0]

# lists of variables to process
dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
vn_list_2d = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
vn_list_2d_custom = ['DA']
vn_list_3d_t_custom = ['hyp_dz']

# make some things
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
DA = G['DX'] * G['DY']
ny, nx = DA.shape
h = G['h']
S = zrfun.get_basic_info(fn, only_S=True)
zr, zw = zrfun.get_z(h, 0*h, S)
dzr = np.diff(zw, axis=0)


dataset1 = xr.DataArray(fn)
dataset2 = xr.DataArray(out_fn)
# was ds1 = nc.Dataset(fn)
#     ds2 = nc.Dataset(out_fn, 'w')

# Create dimensions
for dname, the_dim in dataset1.dimensions.items():
    if dname in dlist:
        dataset2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
        
# Create variables and their attributes
# - first time
vn = 'ocean_time'
varin = dataset1[vn]
vv = dataset2.createVariable(vn, varin.dtype, varin.dimensions)
vv.long_name = varin.long_name
vv.units = varin.units
# - then static 2d fields
for vn in vn_list_2d:
    varin = dataset1[vn]
    vv = dataset1.createVariable(vn, varin.dtype, varin.dimensions)
    vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    vv[:] = dataset1[vn][:]
# - then static custom fields
for vn in vn_list_2d_custom:
    if vn == 'DA':
        vv = dataset2.createVariable(vn, float, ('eta_rho', 'xi_rho'))
        vv.long_name = 'Cell horizontal area'
        vv.units = 'm2'
        vv[:] = DA
#
# - then time-dependent custom 3d fields (processed into 2d)
for vn in vn_list_3d_t_custom:
    if vn == 'hyp_dz':
        vv = dataset2.createVariable(vn, float, ('ocean_time', 'eta_rho', 'xi_rho'))
        vv.long_name = 'Thickess of hypoxic water'
        vv.units = 'm'
        vv.time='ocean_time'
# Copy time dependent data
omat = np.nan * np.ones(h.shape)
omat = np.ma.masked_where(G['mask_rho']==0, omat)

tt = 0
NF = len(fn_list)
for fn in fn_list:
    if np.mod(tt,24)==0:
        print(' working on %d of %d' % (tt, NF))
        sys.stdout.flush()
        
    ds = xr.DataArray(fn)
    # was ds = nc.Dataset(fn)
    
    dataset2['ocean_time'][tt] = ds['ocean_time'][0]
        
    for vn in vn_list_3d_t_custom:
        if vn == 'hyp_dz':
            oxy = ds['oxygen'][0,:,:,:]
            dzrm = np.ma.masked_where(oxy > 61, dzr)
            hyp_dz = dzrm.sum(axis=0)
            dataset2[vn][tt,:,:] = hyp_dz
        
    tt += 1
    ds.close()

if args.testing:
    sys.path.append(os.path.abspath('../plotting'))
    import pfun
    import matplotlib.pyplot as plt
    plt.close('all')
    fs=16
    plt.rc('font', size=fs)
    fig = plt.figure(figsize=(10,12))
    ax = fig.add_subplot(111)
    cs = ax.pcolormesh(G['lon_psi'],G['lat_psi'],
        hyp_dz[1:-1,1:-1]/h[1:-1,1:-1], vmin=0, vmax=1)
    fig.colorbar(cs)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.contour(G['lon_rho'],G['lat_rho'],h, [50, 100, 150, 200, 2000],
        colors=['darkorange','plum','darkorchid'], linewidths=2, linestyles='solid')
    ax.text(.95,.16,'100 m',color='darkorange',weight='bold',transform=ax.transAxes,ha='right')
    ax.text(.95,.13,'200 m',color='plum',weight='bold',transform=ax.transAxes,ha='right')
    ax.text(.95,.1,'2000 m',color='darkorchid',weight='bold',transform=ax.transAxes,ha='right')
    ax.set_title('Hypoxic Depth Fraction')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xticks([-130, -128, -126, -124, -122])
    plt.show()
    plt.rcdefaults()
    
    print(' Saved Variables '.center(50,'='))
    for vn in dataset2.variables:
        print(vn)

dataset1.close()
dataset2.close()

# finale
import collections
result_dict = collections.OrderedDict()
result_dict['out_fn'] = out_fn
result_dict['date_string0'] = args.date_string0
result_dict['date_string1'] = args.date_string1
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'
print('')
for k in result_dict.keys():
    print('%s: %s' % (k, result_dict[k]))
    


