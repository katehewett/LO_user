#!/bin/bash

# code to create lowpass files for new LO with hourly data 2017 (Aurora's roms out), 2024 and 2025 (parker's roms out) (hourly data)

LOe=/dat1/parker/LO/extract/lowpass

python3 $LOe/extract_lowpass.py -gtx cas7_t2_x11b -ro 1 -0 2024.01.01 -1 2024.12.31 -Nproc 4 > lowpass_log2024.log
python3 $LOe/extract_lowpass.py -gtx cas7_t2_x11b -ro 1 -0 2025.01.01 -1 2025.12.31 -Nproc 4 > lowpass_log2025.log

python3 $LOe/extract_lowpass.py -gtx cas7_t1_x11b -ro 5 -0 2017.01.01 -1 2017.12.31 -Nproc 4 > lowpass_log2017.log 