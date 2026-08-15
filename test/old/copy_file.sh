#!/bin/bash
#testing s5cmd for trouble shooting forecast issues August 11-15th 2026

# Tells Bash to allow aliases inside this script.
shopt -s expand_aliases

# Imports your active aliases into the script environment.
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi

# set paths
WORKING_DIR="$HOME/LO_user/test"
LOG_FILE="${WORKING_DIR}/s5_transfer.log"
SOURCE_FILE="${WORKING_DIR}/test_file.txt"
DEST_BUCKET="s3://liveocean-kmhewett/test/"

# Ensure the log directory exists
mkdir -p "$(dirname "$LOG_FILE")"

# Record Start Time
echo "==========================================" >> "$LOG_FILE"
echo "Job Started: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"

# Execute s5cmd command and capture stderr
s5-hewett cp "$SOURCE_FILE" "$DEST_BUCKET" >> "$LOG_FILE" 2>&1
EXIT_CODE=$?

# Record Exit/Error Code
echo "s5cmd Exit Code: $EXIT_CODE" >> "$LOG_FILE"

# Record Stop Time
echo "Job Stopped: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
echo "==========================================" >> "$LOG_FILE"