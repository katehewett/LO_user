"""
This code uses a job_list.py to clip history files to box locations
in the LO domain
From that temp file 'box_###.nc' we calculate a corrosive volume
in get_one_corrosive_volume_rev1.py. 

This code was created because extract_corrosive_volume was taking a 
very long time to run on apogee (and would stop on day 100 - 189). 
 
running testing code with two history files on kh personal mac, looks like:
run extract_clines -gtx cas7_t0_x4b -ro 1 -0 2017.12.12 -1 2017.12.13 -lt lowpass -job shelf_box -test True

scripts dev based on extract_box; extract_box_chunks.py and extract_hypoxic_volume.py

----
processing time and file size: 

get_one_cline runs for each file (e.g. per day or per hour): lowpass run takes
Time to get initial fields = 2.81 sec
Time to apply GSW = 2.21 sec
Time to calc SML + assign T/S= 0.71 sec
Time to calc BML + assign T/S= 0.66 sec
Time to calc max dT/dz(drho/dz) and assign SA(CT)= 0.69 sec
Time to calc interior gradients and T/S= 1.06 sec

And then the total time for extract clines (1 year of lowpass data):
Total processing time = 1299.55 sec 

~21 - 22  minutes for a year, which makes a 6.14 GB file 

---- 
Notes: * after the cat step... we open and reassign time independent vars and delete the temp folders. 
This step takes an additional 6-7 seconds. But we put that on the outside of the loops because we 
just need to do it 1x / job, which makes this run faster
* TODO nearshore/estuaries need to look closer at the N2-max 

"""
import sys
import argparse
from lo_tools import Lfun, zrfun, zfun
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import os
import xarray as xr
import numpy as np
from time import time

import gsw

pid = os.getpid()
print(' Finding pycnocline '.center(60,'='))
print('PID for this job = ' + str(pid))

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas6_v3_lo8b
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str) # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) # e.g. 2019.07.06
parser.add_argument('-lt', '--list_type', type=str) # list type: hourly, daily, weekly, lowpass
# select job name
parser.add_argument('-job', type=str) # job name
# for grabbing pycnocline
#parser.add_argument('-pycno', '--pycnocline', default=True, type=Lfun.boolean_string)
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
#if Ldir['testing']:
#    Ldir['roms_out_num'] = 1
#    Ldir['ds0'] = '2022.08.08'
#    Ldir['ds1'] = '2022.08.09'
#    Ldir['list_type'] = 'daily'
#    Ldir['job'] = 'shelf_box'
    
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]
    
# output location and filenames + temp dir to accumulate individual extractions
dd_str = Ldir['ds0'] + '_' + Ldir['ds1']
bb_str = '_'     # leaving this so can update to have other types of clines later in gtags
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'clines' / (Ldir['job'] + bb_str + dd_str)

temp_box_dir = out_dir / 'temp_box_dir'
temp_cline_dir = out_dir / 'temp_cline_dir'

Lfun.make_dir(out_dir, clean=True)
Lfun.make_dir(temp_box_dir, clean=True)
Lfun.make_dir(temp_cline_dir, clean=True)

cline_fn_final = out_dir / (Ldir['job'] + '_pycnocline' + bb_str + dd_str + '.nc')

# get list of files to work on and check before entering loop over all files 
fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])
G, S, T = zrfun.get_basic_info(fn_list[0])
Lon = G['lon_rho'][0,:]
Lat = G['lat_rho'][:,0]

def check_bounds(lon, lat):
    # error checking
    if (lon < Lon[0]) or (lon > Lon[-1]):
        print('ERROR: lon out of bounds ')
        sys.exit()
    if (lat < Lat[0]) or (lat > Lat[-1]):
        print('ERROR: lat out of bounds ')
        sys.exit()
    # get indices
    ilon = zfun.find_nearest_ind(Lon, lon)
    ilat = zfun.find_nearest_ind(Lat, lat)
    return ilon, ilat
    
# get the job_definitions module, looking first in LO_user (rn this is only in KH LO_user)
pth = Ldir['LO'] / 'extract' / 'clines'
upth = Ldir['LOu'] / 'extract' / 'clines'
if (upth / 'job_list.py').is_file():
    print('Importing job_list from LO_user')
    job_list = Lfun.module_from_file('job_list', upth / 'job_list.py')
else:
    print('Error: code only in LO_user atm')
    sys.exit()
    #print('Importing job_list from LO')
    #job_list = Lfun.module_from_file('job_list', pth / 'job_list.py')
aa, vn_list = job_list.get_cbox(Ldir['job'], Lon, Lat)
lon0, lon1, lat0, lat1 = aa
ilon0, ilat0 = check_bounds(lon0, lat0)
ilon1, ilat1 = check_bounds(lon1, lat1)

# do a final check to drop missing variables from the list
ds = xr.open_dataset(fn_list[0])
print('Original vn_list:')
print(' ' + vn_list)
vn_list = (',').join([item for item in vn_list.split(',') if item in ds.data_vars])
print('Trimmed vn_list:')
print(' ' + vn_list)
ds.close()

# NOTES from box chunks code:
# NOTE: ncks indexing is zero-based but is INCLUSIVE of the last point.
# NOTE: ncks extractions retain singleton dimensions

def messages(mess_str, stdout, stderr):
    # utility function to help with subprocess errors
    try:
        if len(stdout) > 0:
            print(mess_str)
            print(stdout.decode())
    except TypeError:
        pass
    try:
        if len(stderr) > 0:
            print(mess_str)
            print(stderr.decode())
    except TypeError:
        pass
        
# do the extractions - clip the box then feed the temp file to get_one_cline.py
# TODO is there a clever way that i don't have to loop thru NFN 2x: 1st time to clip
# second to run thru clipped files to calc clines?
NFN = len(fn_list)
proc_list = []
proc_list2 = []
tt0 = time()
print('Working on ' + cline_fn_final.name + ' (' + str(NFN) + ' times)')
for ii in range(NFN):
    fn = fn_list[ii]
    sys.stdout.flush()
    # extract one day at a time using ncks
    count_str = ('000000' + str(ii))[-6:]
    box_out_fn = temp_box_dir / ('box_' + count_str + '.nc')
    cmd_list1 = ['ncks',
        '-v', vn_list,
        '-d', 'xi_rho,'+str(ilon0)+','+str(ilon1), '-d', 'eta_rho,'+str(ilat0)+','+str(ilat1),
        '-d', 'xi_u,'+str(ilon0)+','+str(ilon1-1), '-d', 'eta_u,'+str(ilat0)+','+str(ilat1),
        '-d', 'xi_v,'+str(ilon0)+','+str(ilon1), '-d', 'eta_v,'+str(ilat0)+','+str(ilat1-1),
        '--mk_rec_dim', 'ocean_time']
    cmd_list1 += ['-O', str(fn), str(box_out_fn)]
            
    proc = Po(cmd_list1, stdout=Pi, stderr=Pi)
    proc_list.append(proc)

    # screen output about progress
    if (np.mod(ii,10) == 0) and ii>0:
        print(str(ii), end=', ')
        sys.stdout.flush()
    if (np.mod(ii,50) == 0) and (ii > 0):
        print('') # line feed
        sys.stdout.flush()
    if (ii == NFN-1):
        print(str(ii))
        sys.stdout.flush()
        
    # Nproc controls how many ncks subprocesses we allow to stack up
    # before we require them all to finish.
    if ((np.mod(ii,Ldir['Nproc']) == 0) and (ii > 0)) or (ii == NFN-1):
        for proc in proc_list:
            stdout, stderr = proc.communicate()
            messages('ncks messages:', stdout, stderr)
        # make sure everyone is finished before continuing
        proc_list = []
    
    ii += 1

print(' - Time for clipping box; calculating layers = %0.2f sec' % (time()- tt0))   

# Originally, I had get_one_cline_rev.py placed ~here but w/in the ii/NFN loop. 
# But it was grabbing the last box####.nc for get_one_cline_rev.py and not starting at box00000
# so breaking to it's own loop. Need to figure how to make this part faster and not have a 2x loop. 
# Maybe with grep and then put them in order
    
print('Working on extracting gradients...')
for ii in range(NFN):
    fn = fn_list[ii]
    sys.stdout.flush()
    # extract one day at a time using ncks
    count_str = ('000000' + str(ii))[-6:]
    box_out_fn = temp_box_dir / ('box_' + count_str + '.nc')
    cline_out_fn = temp_cline_dir / ('cline_' + count_str + '.nc')  
    print(str(box_out_fn))
    cmd_list2 = ['python3', 'get_one_cline_rev.py',
            '-lt',args.list_type,
            '-his_fn',str(fn_list[0]), 
            '-in_fn',str(box_out_fn),
            '-out_fn', str(cline_out_fn)]
    
    proc2 = Po(cmd_list2, stdout=Pi, stderr=Pi)
    proc_list2.append(proc2)

    # screen output about progress
    if (np.mod(ii,10) == 0) and ii>0:
        print(str(ii), end=', ')
        sys.stdout.flush()
    if (np.mod(ii,50) == 0) and (ii > 0):
        print('') # line feed
        sys.stdout.flush()
    if (ii == NFN-1):
        print(str(ii))
        sys.stdout.flush()
        
    # Nproc controls how many ncks subprocesses we allow to stack up
    # before we require them all to finish.
    if ((np.mod(ii,Ldir['Nproc']) == 0) and (ii > 0)) or (ii == NFN-1):
        for proc2 in proc_list2:
            stdout, stderr = proc2.communicate()
            messages('get_one_cline messages:', stdout, stderr)
        # make sure everyone is finished before continuing
        proc_list2 = []
            
    ii += 1

# don't think we need this step b/c running cas7, but check 
## Ensure that all days have the same fill value.  This was required for cas6_v3_lo8b
## when passing from 2021.10.31 to 2021.11.01 because they had inconsistent fill values,
## which leaks through the ncrcat call below.
## can comment out if using long hindcast
#tt1 = time()
#enc_dict = {'_FillValue':1e20}
#vn_List = vn_list.split(',')
#Enc_dict = {vn:enc_dict for vn in vn_List}
#for out_fn in list(temp_box_dir.glob('box_*.nc')):
#    ds = xr.load_dataset(box_out_fn) # need to load, not open, for overwrite
#    ds.to_netcdf(box_out_fn, encoding=Enc_dict)
#    ds.close()
#print(' - Time for adding fill value = %0.2f sec' % (time()- tt1))
#sys.stdout.flush()

# we now have box.nc's and cline.nc's in a temp dir and want to stitch together all the clines
# At the very end we will delete the temp dirs to save space. 
# The temp file PSB "pycnocline SML BML" will be saved in the temp dir too, because once concat'd 
# we are going to add vars that don't change w/ time: lat lon mask_rho h to the dataset
# and save the final outside the temp dir's

# concatenate the records into one file
# This bit of code is a nice example of how to replicate a bash pipe
pp1 = Po(['ls', str(temp_cline_dir)], stdout=Pi)
pp2 = Po(['grep','cline'], stdin=pp1.stdout, stdout=Pi)
fn_p = 'PSB_'+str(Ldir['ds0'])+'_'+str(Ldir['ds1']+'.nc')
temp_fn = str(temp_cline_dir)+'/'+fn_p # this is all the maps put to one
cmd_list = ['ncrcat','-p', str(temp_cline_dir), '-O', temp_fn]
proc = Po(cmd_list, stdin=pp2.stdout, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if len(stdout) > 0:
    print('\nSTDOUT:')
    print(stdout.decode())
    sys.stdout.flush()
if len(stderr) > 0:
    print('\nSTDERR:')
    print(stderr.decode())
    sys.stdout.flush()

    
print('saving vars')
"""
Saving variables that don't change over time - this is here to make the cat step faster
- h is bathymetric depth 
- ocean_time is a vector of time; [lowpass is converted to UTC and centered at 1200 UTC; daily is diff] 
- Lat and Lon on rho points 
"""

lon_rho = G['lon_rho'][ilat0:ilat1+1,ilon0:ilon1+1]
lat_rho = G['lat_rho'][ilat0:ilat1+1,ilon0:ilon1+1]
mask_rho = G['mask_rho'][ilat0:ilat1+1,ilon0:ilon1+1]
h = G['h'][ilat0:ilat1+1,ilon0:ilon1+1]

ds1 = xr.open_dataset(temp_fn)

ds1['lat_rho'] = (('eta_rho', 'xi_rho'),lat_rho,{'units':'degree_north'})
ds1['lat_rho'].attrs['standard_name'] = 'grid_latitude_at_cell_center'
ds1['lat_rho'].attrs['long_name'] = 'latitude of RHO-points'
ds1['lat_rho'].attrs['field'] = 'lat_rho'

ds1['lon_rho'] = (('eta_rho', 'xi_rho'),lon_rho,{'units': 'degree_east'})
ds1['lon_rho'].attrs['standard_name'] = 'grid_longitude_at_cell_center'
ds1['lon_rho'].attrs['long_name'] = 'longitude of RHO-points'
ds1['lon_rho'].attrs['field'] = 'lon_rho'
 
ds1['h'] = (('eta_rho', 'xi_rho'),h,{'units': 'm'})
ds1['h'].attrs['standard_name'] = 'sea_floor_depth'
ds1['h'].attrs['long_name'] = 'time_independent bathymetry'
ds1['h'].attrs['field'] = 'bathymetry'
ds1['h'].attrs['grid'] =  args.gtagex

ds1['mask_rho'] = (('eta_rho', 'xi_rho'),mask_rho,{'units': 'm'})
ds1['mask_rho'].attrs['standard_name'] = 'land_sea_mask_at_cell_center'
ds1['mask_rho'].attrs['long_name'] = 'mask on RHO-points'
ds1['mask_rho'].attrs['flag_values'] = np.array([0.,1.])
ds1['mask_rho'].attrs['flag_meanings'] = 'land water'
ds1['mask_rho'].attrs['grid'] =  args.gtagex

ds1.to_netcdf(cline_fn_final, unlimited_dims='ocean_time')

## clean up
Lfun.make_dir(temp_box_dir, clean=True)
temp_box_dir.rmdir()

Lfun.make_dir(temp_cline_dir, clean=True)
temp_cline_dir.rmdir()

print('Total processing time = %0.2f sec' % (time()-tt0))

