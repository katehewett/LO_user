"""
Testing kuser flag 

run arg_test.py -g cas7 -t t2 -x x11b -r forecast -s continuation -grp macc -np 192 -vip True -done c11 

"""

import sys, os
import shutil
import argparse
from datetime import datetime, timedelta
from pathlib import Path
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time, sleep
import random
import string
from math import ceil

# Add the path to lo_tools by hand so that it we can import Lfun on klone
# without loenv. In general we write code to run on klone using only the
# default python3 installation.
pth = '/Users/katehewett/Documents/LO/lo_tools/lo_tools'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun

# >>> START Command Line Arguments <<<

# Arguments without defaults are required.

import sys, os
import shutil
import argparse
from datetime import datetime, timedelta
from pathlib import Path
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from time import time, sleep
import random
import string
from math import ceil

# Add the path to lo_tools by hand so that it we can import Lfun on klone
# without loenv. In general we write code to run on klone using only the
# default python3 installation.
pth = Path(__file__).absolute().parent.parent / 'lo_tools' / 'lo_tools'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun

# >>> START Command Line Arguments <<<

# Arguments without defaults are required.

parser = argparse.ArgumentParser()
# Typically when you use these at the command line you can use the short version,
# like "-g" but the long version works as well "--gridname". In the code the long name
# is what is used.

# Basic info to specify which model configuration to run
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas7
parser.add_argument('-t', '--tag', type=str)        # e.g. t2
parser.add_argument('-x', '--ex_name', type=str)    # e.g. x11b

# Set the run_type
parser.add_argument('-r', '--run_type', type=str, default='backfill')
# Choices: forecast or backfill

parser.add_argument('-s', '--start_type', type=str, default='continuation')
# Choices
# - new: only run 1 day, starting from ocean_ini.nc
# - perfect: start from ocean_rst.nc of previous day
# - continuation: start from ocean_his_[last one].nc of previous day
# - newperfect: first day = new, then perfect thereafter
# - newcontinuation: first day = new, then continuation thereafter

# If the run type is backfill you need to set the time range
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str) #  this is set to ds0 if omitted

# The next two/three flags are needed for compute resource allocation.
#
# Set how many cpu's (cores) to use.
parser.add_argument('-np', '--np_num', type=int) # e.g. 192, number of cores
#
# Choose which type of computer resource
parser.add_argument('-cpu','--cpu_choice', default='cpu-g2', type=str) # used in the sbatch script
# Choices: cpu-g2, compute, or ckpt-g2
#
# Choose which group (not needed when using --cpu_choice ckpt_g2)
parser.add_argument('-grp','--group_choice', type=str) # used in the sbatch script
# Choices: macc, coenv
# NOTE: Check with Kate or Parker before using macc or coenv.

# Optional flag used only for forecasts, so that they leave a clue if they are done for a given day.
# This allows a backup job in the crontab to exist but only run if needed. By using different
# done_tags we can have several different forecasts running, with backups, at the same time.
parser.add_argument('-done','--done_tag', type=str, default='00')

# Optional flags to facilitate testing.
parser.add_argument('-v', '--verbose', default=False, type=Lfun.boolean_string)
parser.add_argument('--short_roms', default=False, type=Lfun.boolean_string)
parser.add_argument('--run_dot_in', default=True, type=Lfun.boolean_string)
parser.add_argument('--run_roms', default=True, type=Lfun.boolean_string)

# You should use this for larger runs that use a full node (192 cores). 
parser.add_argument('-vip','--exclusive', default=False, type=Lfun.boolean_string)

# Flags to send output to kopah.
parser.add_argument('-k','--to_kopah', default=True, type=Lfun.boolean_string)
parser.add_argument('-ktest','--test_to_kopah', default=False, type=Lfun.boolean_string)
# This kopah flag sets the destination bucket in the macc group's kopah storage 
# that will recieve history files. If you have access to macc group kopah storage, then 
# Enter your username that is used on klone. 
# Please do not use -kuser pmacc unless you are Kate or Parker. 
parser.add_argument('-kuser','--kopah_user', type=str, default = None) 

# >>> END Command Line Arguments <<<

args = parser.parse_args()

# check for required arguments and other input errors
argsd = args.__dict__
for a in ['gridname', 'tag', 'ex_name', 'run_type', 'start_type', 'np_num', 'cpu_choice']:
    if argsd[a] == None:
        print('*** Missing required argument for driver_roms##.py: ' + a)
        sys.exit()
if (argsd['cpu_choice'] == 'cpu-g2') and (argsd['group_choice'] not in ['macc','coenv']):
    print('*** Need a --group_choice for this --cpu_choice.')
    sys.exit()
if (argsd['cpu_choice'] == 'compute') and (argsd['group_choice'] != 'macc'):
    print('*** Need --group_choice macc for --cpu_choice compute.')
    sys.exit()
if (argsd['run_type'] == 'backfill') and (argsd['ds0'] == None):
    print('*** Need at least -0 YYYY.MM.DD for -r backfill')
    sys.exit()
if (argsd['kopah_user'] is None) and (argsd['to_kopah'] == True):
    print('-kuser, kopah_user, is blank. You need to enter your kopah user name for macc storage. \n' \
    'This is likely your username on klone. \n' \
    'If you do not have macc storage credentials check with Kate or Parker.')
    sys.exit()
if (argsd['kopah_user'] == 'pmacc') and (os.environ.get('USER') not in ('parker', 'pmacc', 'kmhewett','katehewett')):
    print('Check with Kate or Parker on kopah storage credentials. Do not send to pmacc without checking first. Thanks!')
    sys.exit()

# Override some flags when testing to_kopah
if args.test_to_kopah:
    args.run_dot_in = False
    args.run_roms = False

# get Ldir
Ldir = Lfun.Lstart(gridname=args.gridname, tag=args.tag, ex_name=args.ex_name)




