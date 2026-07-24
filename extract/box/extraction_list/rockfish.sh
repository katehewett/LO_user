#!/bin/bash

# code to extract YOY_Rockfish_Recruit 
# These exports were completed for Dayv Lowry and Adam Obaza 
# asked for 2015 - 2025 rockfish recruit work 

LOe=/dat1/kmhewett/LO/extract/box

python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2015.01.01 -1 2015.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2015.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2016.01.01 -1 2016.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2016.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2017.01.01 -1 2017.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2017.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2018.01.01 -1 2018.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2018.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2019.01.01 -1 2019.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2019.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2020.01.01 -1 2020.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2020.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2021.01.01 -1 2021.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2021.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2022.01.01 -1 2022.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2022.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2023.01.01 -1 2023.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2023.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2024.01.01 -1 2024.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2024.log
python3 $LOe/extract_box_chunks.py -gtx cas7_t1_x11ab -ro 1 -0 2025.01.01 -1 2025.12.31 -lt average -job YOY_Rockfish_Recruit > rockfish_2025.log