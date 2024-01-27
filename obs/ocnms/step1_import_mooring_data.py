"""
Code to process the OCNMS mooring data.

Takes processed OCNMS mooring data saved as .mat files and 
1) calculates conservative temperatutre; practical salinity; sigma0
2) converts oxygen values from mg/L to uM 
3) places in xarray Datasets in LO standard form 

Questions: 
(1) would saving with Z x time dimensions be better than time x Z? 
becasue that matches the LO mooring extraction format. Reason that this is time x Z 
is because it matched the orca formatting. 
(2) The depths are estimated here. Need to extract the surface elevation from 
each mooring and then estimate Z (fix after talking with Parker on Jan 23rd)

Comment: 2011 thru 2023 is being imported. data prior to 2011 hasn't been QC'd and 
needs extra O2 (and maybe sal) qc.


"""
#import mat73 # is mat73 better/worse than scipy.io?(all OCNMS are provided as .mat files)
import scipy.io as sio
import pandas as pd
import xarray as xr
import numpy as np
import gsw
from time import time

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-test','--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
testing = args.testing

# input location
source = 'ocnms'
otype = 'moor' # introducing a new "otype" beyond ctd and bottle

in_dir = Ldir['data'] / 'obs' / 'ocnms' / 'ocnms_mooring' / 'datasets'
#in_dir = '/Users/katehewett/Documents/LKH_output/OCNMS_processed/step2_sitestructs/daily' 
## need to update to above line

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
if not testing:
    Lfun.make_dir(out_dir, clean=True)

sn_name_dict = {
    'MB015':'Makah Bay 15m',
    'MB042':'Makah Bay 42m',
    'CA015':'Cape Alava 15m',
    'CA042':'Cape Alava 42m',
    'TH015':'Teahwhit Head 15m',
    'TH042':'Teahwhit Head 42m',
    'KL015':'Kalaloch 15m',
    'KL042':'Kalaloch 42m',
    'CE015':'Cape Elizabeth 15m',
    'CE042':'Cape Elizabeth 42m'
}

if testing:
    sn_list = ['CE042']
else:
    sn_list = list(sn_name_dict.keys())
    
sn_loc_dict = {
    'MB015': [48.32538333,-124.6768333],
    'MB042': [48.32396667,-124.7353833],
    'CA015': [48.16630000,-124.7568333],
    'CA042': [48.16601667,-124.8233667],
    'TH015': [47.87550000,-124.6194667],
    'TH042': [47.87615000,-124.7334167],
    'KL015': [47.60083333,-124.4284000],
    'KL027': [47.59456667,-124.4970667],
    'CE015': [47.35678333,-124.3481333],
    'CE042': [47.35313333,-124.4887333]
}


for sn in sn_list:
    print(sn)
    #in_fn = in_dir / (sn + '_ds_daily.nc')
    in_fn = in_dir / (sn + '_2011_2023_daily.mat')
    out_fn = out_dir / (sn + '_daily.nc')
    
    #ds0 = xr.open_dataset(in_fn)
    mat_data = sio.loadmat(in_fn)
    #matlab_datetime_values = mat_data['timestamp_UTC'][0]
    #times aren't loading from matlab, but the sampling time vector is known
    #cheat and need to fix:
    # Generate datetime array using pandas
    start_date = pd.to_datetime('2011-01-01 20:00:00',utc=True)
    end_date = pd.to_datetime('2023-12-31 20:00:00',utc=True)
    datetime_array_pandas = pd.date_range(start=start_date, end=end_date, freq='D')
    
    end_date = pd.to_datetime('2022-01-10')
    datetime_array_pandas = pd.date_range(start_date, end_date, freq='D')
    
    SP = mat_data['SAL']
    IT = mat_data['IT']
    OXY = mat_data['OXY']
    pressure = mat_data['P']
    Z_est = mat_data['Z_est']
    
    z = Z_est;
    
    # convert to LO standard units
    z = -ds0.depth.values # packed shallow to deep
    zz = z[::-1] # packed deep to shallow
    NZ = len(zz)
    t = ds0.time.values
    tt = pd.DatetimeIndex(t)
    NT = len(tt)
    Z = zz.reshape((1,NZ)) * np.ones((NT,1))
    # salt and temp
    SP = np.transpose(ds0.sal.values[::-1,:])
    IT = np.transpose(ds0.temp.values[::-1,:])
    lon = sn_loc_dict[sn][0]
    lat = sn_loc_dict[sn][1]
    P = gsw.p_from_z(Z, lat)
    # - do the conversions
    SA = gsw.SA_from_SP(SP, P, lon, lat)
    CT = gsw.CT_from_t(SA, IT, P)
    # oxygen converted mg/L to uM
    DO = (1000/32) * np.transpose(ds0.oxy.values[::-1,:])
    # compute sigma0
    # potential density relative to 0 dbar, minus 1000 kg m-3
    SIG0 = gsw.sigma0(SA,CT)
    
    # initialize new Dataset and fill
    coords = {'time':(('time'),tt),'z':(('z'),zz)}
    ds = xr.Dataset(coords=coords, attrs={'Station Name':sn_name_dict[sn],'lon':lon,'lat':lat})
    ds['SA'] = xr.DataArray(SA, dims=('time','z'),
        attrs={'units':'g kg-1', 'long_name':'Absolute Salinity'})
    ds['CT'] = xr.DataArray(CT, dims=('time','z'),
        attrs={'units':'degC', 'long_name':'Conservative Temperature'})
    ds['DO (uM)'] = xr.DataArray(DO, dims=('time','z'),
        attrs={'units':'uM', 'long_name':'Dissolved Oxygen'})
        
    # add sigma0
    ds['SIG0'] = xr.DataArray(SIG0, dims=('time','z'),
        attrs={'units':'kg m-3', 'long_name':'Sigma0'})
    
    if not testing:
        ds.to_netcdf(out_fn)
