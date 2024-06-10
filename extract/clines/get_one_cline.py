"""
Function to do the calculation for a single clipped file
calculates GSW vars > drho/dz dT/dz
finds max drho(dT)/dz + returns the max depth + values

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
parser.add_argument('-in_fn', type=str) # path to clipped box files (temp directory)
parser.add_argument('-out_fn', type=str) # path to save cline files (temp directory)
parser.add_argument('-lt', '--list_type', type=str) # list type: hourly, daily, weekly, lowpass
            
args = parser.parse_args()

box_name = args.out_fn.name # name of temporary input file (a clipped box_######.nc) 

print('hello')
