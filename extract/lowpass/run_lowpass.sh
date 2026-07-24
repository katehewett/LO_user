#!/bin/bash

# code to create lowpass files for 2017, 2024 and 2025 (hourly data)

LOe=/dat1/parker/LO/extract/lowpass

python3 $LOe/extract_lowpass.py -gtx cas7_t1_x11ab -ro 1 -0 2017.01.01 -1 2017.12.31 -Nproc 4 > lowpass_log2017.log
python3 $LOe/extract_lowpass.py -gtx cas7_t1_x11ab -ro 1 -0 2024.01.01 -1 2024.12.31 -Nproc 4 > lowpass_log2024.log
python3 $LOe/extract_lowpass.py -gtx cas7_t1_x11ab -ro 1 -0 2025.01.01 -1 2025.12.31 -Nproc 4 > lowpass_log2025.log