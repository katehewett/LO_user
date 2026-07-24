#!/bin/bash

# These exports were completed for UBC (Susan A and Doug L et al) for
# 2012 - 2023 
# will repeat 2017 and first 8 months of 2024, with a special lowpass extraction 
# had already started 2013 - 2017

LOe=/dat1/kmhewett/LO/extract/box

python3 $LOe/extract_box_nocat.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2018.01.01 -1 2018.12.31 -job ubc0 > ubc0_2018.log 
python3 $LOe/extract_box_nocat.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2019.01.01 -1 2019.12.31 -job ubc0 > ubc0_2019.log 
python3 $LOe/extract_box_nocat.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2020.01.01 -1 2020.12.31 -job ubc0 > ubc0_2020.log 
python3 $LOe/extract_box_nocat.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2021.01.01 -1 2021.12.31 -job ubc0 > ubc0_2021.log 
python3 $LOe/extract_box_nocat.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2022.01.01 -1 2022.12.31 -job ubc0 > ubc0_2022.log 
python3 $LOe/extract_box_nocat.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2023.01.01 -1 2023.12.31 -job ubc0 > ubc0_2023.log 
python3 $LOe/extract_box_nocat.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2024.01.01 -1 2024.12.31 -job ubc0 > ubc0_2024.log 
python3 $LOe/extract_box_nocat.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2012.10.07 -1 2012.12.31 -job ubc0 > ubc0_2012.log 