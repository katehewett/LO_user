"""
Code to extract tef2 sections.

To test on mac:
run extract_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.04 -test True
run extract_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.04 -Nproc 10 -get_bio True

Doing this with subprocesses (Nproc = 10, on my mac) was about 2x as fast as doing
it sequentially within this program. The less-than-expected speedup may be because
each job only takes about a second, and there is overhead to spinning up new python
jobs because of imports.

Also, this is a memory-intensive calculation, so be careful about using Nproc > 10
(10 is the default in extract_argfun).

Performance: took about 1-2 sec per history file (Nproc = 10, on my mac).
- 58 sec per day with get_bio True (11 3-D variables)
- 24 sec per day with get_bio False (only salt)

"""

from lo_tools import Lfun, zrfun, zfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time
import sys
import pandas as pd
import xarray as xr
import numpy as np
import pickle

# gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
gctag = Ldir['gridname'] + '_' + Ldir['job']
hv_dir = Ldir['LOo'] / 'extract' / 'hypoxic_volume'

hv_df_fn = hv_dir / ('hv_df_' + gctag + '.p')
hv_df = pd.read_pickle(hv_df_fn)

