# code to extract UBC 
# These exports were completed for UBC (Susan A and Doug L et al) for
# 2012 - 2023 
# will repeat 2017 and first 8 months of 2024, with a special lowpass extraction 

#LOe=/dat1/kmhewett/LO/extract/box

import os
os.environ["LOe"] = '/dat1/kmhewett/LO/extract/box'

python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2012.10.07 -1 2012.12.31 -lt average -job ubc0 > ubc_2012.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2013.01.01 -1 2013.12.31 -lt average -job ubc0 > ubc_2013.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2014.01.01 -1 2014.12.31 -lt average -job ubc0 > ubc_2014.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2015.01.01 -1 2015.12.31 -lt average -job ubc0 > ubc_2015.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2016.01.01 -1 2016.12.31 -lt average -job ubc0 > ubc_2016.log

python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2017.01.01 -1 2017.12.31 -lt average -job ubc0 > ubc_2017.log

python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2018.01.01 -1 2018.12.31 -lt average -job ubc0 > ubc_2018.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2019.01.01 -1 2019.12.31 -lt average -job ubc0 > ubc_2019.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2020.01.01 -1 2020.12.31 -lt average -job ubc0 > ubc_2020.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2021.01.01 -1 2021.12.31 -lt average -job ubc0 > ubc_2021.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2022.01.01 -1 2022.12.31 -lt average -job ubc0 > ubc_2022.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2023.01.01 -1 2023.12.31 -lt average -job ubc0 > ubc_2023.log

python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2024.01.01 -1 2024.08.24 -lt average -job ubc0 > ubc_2024.log