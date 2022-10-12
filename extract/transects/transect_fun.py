"""
Functions used the by transect extractor and processing code.
"""
import pandas as pd
import netCDF4 as nc
import numpy as np
import pickle
import sys

from lo_tools import zfun, zrfun

# long list of variables to extract
vn_list = ['salt', 'temp', 'oxygen',
    'NO3', 'phytoplankton', 'zooplankton', 'detritus', 'Ldetritus',
    'TIC', 'alkalinity']

def get_sect_df(gridname):
    # section definitions
    # * x and y are latitude and longitude and we require sections to be NS or EW so
    # either x0=x1 or y0=y1
    sect_df = pd.DataFrame(columns=['x0', 'x1', 'y0', 'y1'])
    
    if gridname == 'cas6':
        # NH Line
        sect_df.loc['NHL1_NHL45',:] = [-124.0999985, -125.1166992,   44.65169907,   44.65169907]
        
    else:
        print('** sect_df not supported for this gridname **')
        sys.exit()

    return sect_df
    
def get_inds(x0, x1, y0, y1, G, verbose=False):
    # get grid indices and lon,lat for the TEF sections
    
    # determine the direction of the section
    # and make sure indices are *increasing*
    if (x0==x1) and (y0!=y1):
        sdir = 'NS'
        a = [y0, y1]; a.sort()
        y0 = a[0]; y1 = a[1]
    elif (x0!=x1) and (y0==y1):
        sdir = 'EW'
        a = [x0, x1]; a.sort()
        x0 = a[0]; x1 = a[1]
    else:
        print('Input points do not form a proper section')
        sdir='bad'
        sys.exit()
    
    # we assume a plaid grid, as usual
    if sdir == 'NS':
        lon = G['lon_u'][0,:].squeeze()
        lat = G['lat_u'][:,0].squeeze()
    elif sdir == 'EW':
        lon = G['lon_v'][0,:].squeeze()
        lat = G['lat_v'][:,0].squeeze()
        
    # we get all 4 i's or j's but only 3 are used
    i0, i1, fr = zfun.get_interpolant(np.array([x0]), lon)
    if np.isnan(fr):
        print('Bad x point')
        sys.exit()
    else:
        ii0 = int(i0)
    i0, i1, fr = zfun.get_interpolant(np.array([x1]), lon)
    if np.isnan(fr):
        print('Bad x point')
        sys.exit()
    else:
        ii1 = int(i1)
    j0, j1, fr = zfun.get_interpolant(np.array([y0]), lat)
    if np.isnan(fr):
        print('Bad y0 point')
        sys.exit()
    else:
        jj0 = int(j0)
    j0, j1, fr = zfun.get_interpolant(np.array([y1]), lat)
    if np.isnan(fr):
        print('Bad y1 point')
        sys.exit()
    else:
        jj1 = int(j1)

    # get mask and trim indices
    # Note: the mask in G = 1 on water points
    if sdir == 'NS':
        mask = G['mask_u'][jj0:jj1+1, ii0] == 1
        # Note: argmax finds the index of the first True in this case
        igood0 = np.argmax(mask)
        igood1 = np.argmax(mask[::-1])
        # check to see if section is "closed"
        if (igood0==0) | (igood1==0):
            print('Warning: not closed one or both ends')
        # keep only to end of water points, to allow for ocean sections
        Mask = mask[igood0:-igood1]
        # and change the indices to match.  These will be the indices
        # of the start and end points.
        jj0 = jj0 + igood0
        jj1 = jj1 - igood1
        if verbose:
            print('  sdir=%2s: jj0=%4d, jj1=%4d, ii0=%4d' % (sdir, jj0, jj1, ii0))
        Lat = lat[jj0:jj1+1]
        Lon = lon[ii0] * np.ones_like(Mask)
    elif sdir == 'EW':
        mask = G['mask_v'][jj0, ii0:ii1+1] == 1
        igood0 = np.argmax(mask)
        igood1 = np.argmax(mask[::-1])
        # check to see if section is "closed"
        if (igood0==0) | (igood1==0):
            print('Warning: not closed one or both ends')
        Mask = mask[igood0:-igood1]
        ii0 = ii0 + igood0
        ii1 = ii1 - igood1
        if verbose:
            print('  sdir=%2s: jj0=%4d, ii0=%4d, ii1=%4d' % (sdir, jj0, ii0, ii1))
        Lon = lon[ii0:ii1+1]
        Lat = lat[jj0] * np.ones_like(Mask)
    
    Lon = zfun.fillit(Lon)
    Lat = zfun.fillit(Lat)
    
    return ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask
    
def start_netcdf(fn, out_fn, NT, NX, NZ, Lon, Lat, Ldir, vn_list):
    # initialize the NetCDF files the extract_sections.py creates
    out_fn.unlink(missing_ok=True)
    ds = nc.Dataset(out_fn, 'w')
    # and some dicts of long names and units
    long_name_dict = dict()
    units_dict = dict()
    for vn in vn_list + ['ocean_time']:
        try:
            long_name_dict[vn] = ds.variables[vn].long_name
        except:
            long_name_dict[vn] = ''
        try:
            units_dict[vn] = ds.variables[vn].units
        except:
            units_dict[vn] = ''
    ds.close()
    # add custom dict fields
    long_name_dict['q'] = 'transport'
    units_dict['q'] = 'm3 s-1'
    long_name_dict['lon'] = 'longitude'
    units_dict['lon'] = 'degrees'
    long_name_dict['lat'] = 'latitude'
    units_dict['lat'] = 'degrees'
    long_name_dict['h'] = 'depth'
    units_dict['h'] = 'm'
    long_name_dict['z0'] = 'z on rho-grid with zeta=0'
    units_dict['z0'] = 'm'
    long_name_dict['DA0'] = 'cell area on rho-grid with zeta=0'
    units_dict['DA0'] = 'm2'
    long_name_dict['DA'] = 'cell area on rho-grid'
    units_dict['DA'] = 'm2'

    # initialize netcdf output file
    foo = nc.Dataset(out_fn, 'w')
    foo.createDimension('xi_sect', NX)
    foo.createDimension('s_rho', NZ)
    foo.createDimension('ocean_time', NT)
    foo.createDimension('sdir_str', 2)
    for vv in ['ocean_time']:
        v_var = foo.createVariable(vv, float, ('ocean_time',))
        v_var.long_name = long_name_dict[vv]
        v_var.units = units_dict[vv]
    for vv in vn_list + ['q', 'DA']:
        v_var = foo.createVariable(vv, float, ('ocean_time', 's_rho', 'xi_sect'))
        v_var.long_name = long_name_dict[vv]
        v_var.units = units_dict[vv]
    for vv in ['z0', 'DA0']:
        v_var = foo.createVariable(vv, float, ('s_rho', 'xi_sect'))
        v_var.long_name = long_name_dict[vv]
        v_var.units = units_dict[vv]
    for vv in ['lon', 'lat', 'h']:
        v_var = foo.createVariable(vv, float, ('xi_sect'))
        v_var.long_name = long_name_dict[vv]
        v_var.units = units_dict[vv]
    for vv in ['zeta']:
        v_var = foo.createVariable(vv, float, ('ocean_time', 'xi_sect'))
        v_var.long_name = 'Free Surface Height'
        v_var.units = 'm'

    # add static variables
    foo['lon'][:] = Lon
    foo['lat'][:] = Lat

    # add global attributes
    foo.gtagex = Ldir['gtagex']
    foo.date_string0 = Ldir['ds0']
    foo.date_string1 = Ldir['ds1']

    foo.close()

def add_fields(out_fn, temp_dir, sect_name, vn_list, S, NT):
    # unpack the data made by extract_one_time.py and load into NetCDF
    foo = nc.Dataset(out_fn, 'a')
    A_list = list(temp_dir.glob('A*.p'))
    A_list.sort()
    count = 0
    for A_fn in A_list:
        A = pickle.load(open(A_fn, 'rb'))
        C = A[sect_name]
        
        if count == 0:
            d = C['d']
            NX = len(d)
            NZ = S['N']
            h = C['h']
            foo['h'][:] = h
            z0 = zrfun.get_z(h, 0*h, S, only_rho=True)
            foo['z0'][:] = z0
            zw0 = zrfun.get_z(h, 0*h, S, only_w=True)
            DZ0 = np.diff(zw0, axis=0)
            DA0 = d.reshape((1, NX)) * DZ0
            foo['DA0'][:] = DA0
            zeta_arr = np.nan * np.ones((NT, NX))
            h_arr = np.nan * np.ones((NT, NX))
            vel_arr = np.nan * np.ones((NT, NZ, NX))
        zeta_arr[count,:] = C['zeta']
        h_arr[count,:] = h
        vel_arr[count,:,:] = C['vel']
        for vn in vn_list:
            foo[vn][count,:,:] = C[vn]
        foo['ocean_time'][count] = C['ot']
        count += 1
    z = zrfun.get_z(h_arr, zeta_arr, S, only_w=True)
    # initially z is packed (z,t,x)
    z = np.transpose(z, (1,0,2))
    # now it should be (t,z,x)
    DZ = np.diff(z, axis=1)
    DA = d.reshape((1, 1, NX)) * DZ
    q = vel_arr * DA
    foo['zeta'][:] = zeta_arr
    foo['q'][:] = q
    foo['DA'][:] = DA
    foo.close()
    
