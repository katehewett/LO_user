#!/usr/bin/env python3
import os
import sys
import subprocess
from datetime import datetime

# ==============================================================================
# Resolves $HOME directly using native OS paths
HOME_DIR = os.path.expanduser("~")

if str(HOME_DIR) == '/Users/katehewett':
    WORKING_DIR = os.path.join(HOME_DIR, "Documents", "LO_user", "test")
elif (str(HOME_DIR) == '/mmfs1/home/kmhewett'):
    WORKING_DIR = os.path.join(HOME_DIR, "LO_user", "test")

LOG_FILE = os.path.join(WORKING_DIR, "s5_transfer.log")
SOURCE_FILE = os.path.join(WORKING_DIR, "test_file.txt")
DEST_BUCKET = "s3://liveocean-kmhewett/test/"

# ==============================================================================
# Ensure the log directory exists
os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)

# Open log file in append mode ('a')
with open(LOG_FILE, "a") as log:
    # Record Start Time
    log.write("==========================================\n")
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log.write(f"Job Started: {start_time}\n")
    log.flush()  # Force write to disk

    # ==============================================================================
    # get the s5-macc details in python 
    env = os.environ.copy()       # copy current system environment variables

    # Map your existing MACC keys to the AWS variables s5cmd looks for
    macc_key = env.get("MACC_KEY")
    macc_secret = env.get("MACC_SECRET")

    if not macc_key or not macc_secret:
        log.write("CRITICAL ERROR: MACC_KEY or MACC_SECRET is missing from environment!\n")
        log.write("Make sure they are 'export'ed in your ~/.bash_profile and you ran 'source ~/.bash_profile'.\n")
        log.write("==========================================\n")
        sys.exit("Error: Environment variables missing. Check s5_transfer.log for details.")

    # Map variables for s5cmd
    env["AWS_ACCESS_KEY_ID"] = macc_key
    env["AWS_SECRET_ACCESS_KEY"] = macc_secret
    
    # Build the exact command execution list
    cmd = ["s5cmd", "cp", SOURCE_FILE, DEST_BUCKET]

    try:
        # Run s5cmd, merging stderr into stdout, and pipe all output directly to the log file
        result = subprocess.run(
            cmd, 
            env=env, #pass custom env with keys
            stderr=subprocess.STDOUT, 
            text=True
        )
        exit_code = result.returncode
    except FileNotFoundError:
        log.write("CRITICAL ERROR: 's5cmd' executable was not found in your system PATH.\n")
        exit_code = 127
    except Exception as e:
        log.write(f"CRITICAL ERROR: Execution failed: {str(e)}\n")
        exit_code = 1

    # Record Exit/Error Code
    log.write(f"s5cmd Exit Code: {exit_code}\n")

    # Record Stop Time
    stop_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log.write(f"Job Stopped: {stop_time}\n")
    log.write("==========================================\n")