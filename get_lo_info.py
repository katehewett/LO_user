"""
This is the one place where you set the path structure of the LO code.
The info is stored in the dict Ldir.

All paths are pathlib.Path objects.

This program is meant to be loaded as a module by Lfun which then adds more
entries to the Ldir dict based on which model run you are working on.

Users should copy this to LO_user/get_lo_info.py, edit as needed, and make it into
their own GitHub repo.

"""
import os
from pathlib import Path

# defaults that should work on all machines
parent = Path(__file__).absolute().parent.parent
LO = parent / 'LO'
LOo = parent / 'LO_output'
LOu = parent / 'LO_user'
data = parent / 'LO_data'

# This is a new piece of information, to help with integration of
# Aurora Leeson's new LO_traps repo, 2023.11.03.
traps_name = 'traps00'
# In order for this to be more useful it would have to be integrated
# into Aurora's code.
# I'm not sure this is the best way to solve this problem.

# These are places where the ROMS history files are kept
roms_out = parent / 'LO_roms'
roms_out1 = parent / 'BLANK' # placeholder
roms_out2 = parent / 'BLANK'
roms_out3 = parent / 'BLANK'
roms_out4 = parent / 'BLANK'

# these are for klone
remote_user = 'BLANK'
remote_machine = 'BLANK'
remote_dir0 = 'BLANK'
local_user = 'BLANK'

# default for linux machines
which_matlab = '/usr/local/bin/matlab'

HOME = Path.home()
try:
    HOSTNAME = os.environ['HOSTNAME']
except KeyError:
    HOSTNAME = 'BLANK'
    
# debugging
# print('** from get_lo_info.py **')
# print('HOME = ' + str(HOME))
# print('HOSTNAME = ' + HOSTNAME)

if (str(HOME) == '/Users/katehewett'):
    lo_env = 'kh_mac'
    which_matlab = '/Applications/MATLAB_R2023b.app/bin/matlab'
    roms_out1 = Path('/Users/katehewett/Documents/LO_roms') # beeecareful!
    roms_out2 = Path('/Users/katehewett/Documents/apogee_parker/LO_roms') # special case
    roms_out3 = Path('/Users/katehewett/Documents/apogee_auroral/LO_roms') # special case

# parker
elif (str(HOME) == '/home/parker') & ('perigee' in HOSTNAME):
    lo_env = 'pm_perigee'
    roms_out1 = Path('/agdat1/parker/LO_roms')
    roms_out2 = Path('/agdat2/parker/LO_roms')
    roms_out3 = Path('/data1/auroral/LO_roms')
    roms_out4 = Path('/data2/parker/LiveOcean_roms/output')

elif (str(HOME) == '/home/parker') & ('apogee' in HOSTNAME):
    lo_env = 'pm_apogee'
    roms_out1 = Path('/dat1/parker/LO_roms')
    roms_out2 = Path('/dat2/parker/LO_roms')
    roms_out3 = Path('/dat2/jxiong/LO_roms')
    roms_out4 = Path('/pgdat2/parker/LO_roms')

elif ((str(HOME) == '/mmfs1/home/pmacc') or (str(HOME) == '/mmfs1/home/darrd')):
    lo_env = 'pm_klone'
    remote_user = 'parker'
    remote_machine = 'apogee.ocean.washington.edu'
    remote_dir0 = '/dat1/parker'
    local_user = 'pmacc'

# kate
elif (str(HOME) == '/home/kmhewett') & ('perigee' in HOSTNAME):
    lo_env = 'kh_perigee'
    roms_out1 = Path('/agdat1/parker/LO_roms')
    roms_out2 = Path('/agdat2/parker/LO_roms')
    roms_out3 = Path('/data1/parker/LO_roms')
    roms_out4 = Path('/data2/parker/LiveOcean_roms/output') # old liveocean
    #roms_out5 = Path('/data1/kmhewett/LO_roms')             

elif (str(HOME) == '/home/kmhewett') & ('apogee' in HOSTNAME):
    lo_env = 'kh_apogee'
    roms_out1 = Path('/dat1/parker/LO_roms')
    roms_out2 = Path('/dat2/parker/LO_roms')
    #roms_out3 = Path('/pgdat1/parker/LO_roms')
    roms_out3 = Path('/dat1/auroral/LO_roms')
    roms_out4 = Path('/pgdat2/parker/LO_roms')
    roms_out5 = Path('/dat1/kmhewett/LO_roms')

elif (str(HOME) == '/mmfs1/home/kmhewett'):
    lo_env = 'kh_klone'
    remote_user = 'kmhewett'
    remote_machine = 'apogee.ocean.washington.edu'
    remote_dir0 = '/dat1/kmhewett'
    local_user = 'kmhewett'

Ldir0 = dict()
Ldir0['lo_env'] = lo_env
Ldir0['parent'] = parent
Ldir0['LO'] = LO
Ldir0['LOo'] = LOo
Ldir0['LOu'] = LOu
Ldir0['data'] = data
Ldir0['roms_out'] = roms_out
Ldir0['roms_out1'] = roms_out1
Ldir0['roms_out2'] = roms_out2
Ldir0['roms_out3'] = roms_out3
Ldir0['roms_out4'] = roms_out4
#
Ldir0['remote_user'] = remote_user
Ldir0['remote_machine'] = remote_machine
Ldir0['remote_dir0'] = remote_dir0
Ldir0['local_user'] = local_user
#
Ldir0['traps_name'] = traps_name