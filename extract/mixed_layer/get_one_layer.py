"""
Function to do the calculation for depth of sml for a single history file.
This does the whole domain. 

Want to add portion that allows a lat/lon box 

"""

import numpy as np
import sys

from argparse import ArgumentParser
from lo_tools import zfun, zrfun
from xarray import open_dataset, Dataset
from numpy import nan, ones, diff
from time import time
from lo_tools import zrfun

parser = ArgumentParser()
parser.add_argument('-in_fn', type=str) # path to history file
parser.add_argument('-out_fn', type=str) # path to outfile (temp directory)
parser.add_argument('-lt', '--list_type', default = 'daily', type=str) # list type: hourly, daily, weekly
args = parser.parse_args()

tt0 = time()
ds = open_dataset(args.in_fn, decode_times=False)
# the decode_times=False part is important for correct treatment
# of the time axis later when we concatenate things in the calling function
# using ncrcat

G, S, T = zrfun.get_basic_info(args.in_fn)   # grid information

# Fields that do not change with time
lon = ds.lon_rho.values                
lat = ds.lat_rho.values
DA = G['DX'] * G['DY']                       # cell horizontal area 
CC = dict()                                  # this is for holding fields extracted on sections
CC['lon'] = lon
CC['lat'] = lat
CC['DA'] = DA
h = ds.h.values       
CC['h'] = h
mask_rho = ds.mask_rho.values
CC['mask_rho']=mask_rho
z_rho, z_w = zrfun.get_z(h, 0*h, S)
dzr = np.diff(z_w, axis=0)
print('Time to get initial fields = %0.2f sec' % (time()-tt0))
sys.stdout.flush()



