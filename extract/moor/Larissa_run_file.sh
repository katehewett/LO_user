#!/bin/bash

# These exports were completed for work for Larissa. so much. data.
# 2015 - 2025

LOe=/dat1/kmhewett/LO/extract/moor

python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2014.01.01 -1 2015.12.31 -job OCNMS_jobs -get_all True > Diaz_ocnms1.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2016.01.01 -1 2017.12.31 -job OCNMS_jobs -get_all True > Diaz_ocnms2.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2018.01.01 -1 2019.12.31 -job OCNMS_jobs -get_all True > Diaz_ocnms3.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2020.01.01 -1 2021.12.31 -job OCNMS_jobs -get_all True > Diaz_ocnms4.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2022.01.01 -1 2023.12.31 -job OCNMS_jobs -get_all True > Diaz_ocnms5.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2024.01.01 -1 2025.12.31 -job OCNMS_jobs -get_all True > Diaz_ocnms6.log

python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2014.09.30 -1 2015.12.31 -job ORCA_jobs -get_all True > Diaz_orca1.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2016.01.01 -1 2017.12.31 -job ORCA_jobs -get_all True > Diaz_orca2.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2018.01.01 -1 2019.12.31 -job ORCA_jobs -get_all True > Diaz_orca3.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2020.01.01 -1 2021.12.31 -job ORCA_jobs -get_all True > Diaz_orca4.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2022.01.01 -1 2023.12.31 -job ORCA_jobs -get_all True > Diaz_orca5.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2024.01.01 -1 2025.12.31 -job ORCA_jobs -get_all True > Diaz_orca6.log

python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2014.09.30 -1 2015.12.31 -job RCA_jobs -get_all True > Diaz_rca1.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2016.01.01 -1 2017.12.31 -job RCA_jobs -get_all True > Diaz_rca2.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2018.01.01 -1 2019.12.31 -job RCA_jobs -get_all True > Diaz_rca3.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2020.01.01 -1 2021.12.31 -job RCA_jobs -get_all True > Diaz_rca4.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2022.01.01 -1 2023.12.31 -job RCA_jobs -get_all True > Diaz_rca5.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2024.01.01 -1 2025.12.31 -job RCA_jobs -get_all True > Diaz_rca6.log

python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2014.09.30 -1 2015.12.31 -job CEA_jobs -get_all True > Diaz_cea1.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2016.01.01 -1 2017.12.31 -job CEA_jobs -get_all True > Diaz_cea2.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2018.01.01 -1 2019.12.31 -job CEA_jobs -get_all True > Diaz_cea3.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2020.01.01 -1 2021.12.31 -job CEA_jobs -get_all True > Diaz_cea4.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2022.01.01 -1 2023.12.31 -job CEA_jobs -get_all True > Diaz_cea5.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2024.01.01 -1 2025.12.31 -job CEA_jobs -get_all True > Diaz_cea6.log

python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2018.03.01 -1 2019.09.22 -job anemone_jobs -get_all True > Diaz_anemone.log

python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2014.09.30 -1 2015.12.31 -job Hatchery_jobs -get_all True > Diaz_Hatchery_jobs1.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2016.01.01 -1 2017.12.31 -job Hatchery_jobs -get_all True > Diaz_Hatchery_jobs2.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2018.01.01 -1 2019.12.31 -job Hatchery_jobs -get_all True > Diaz_Hatchery_jobs3.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2020.01.01 -1 2021.12.31 -job Hatchery_jobs -get_all True > Diaz_Hatchery_jobs4.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2022.01.01 -1 2023.12.31 -job Hatchery_jobs -get_all True > Diaz_Hatchery_jobs5.log
python3 $LOe/multi_mooring_driver.py -gtx cas7_t1_x11ab -ro 1 -lt average -0 2024.01.01 -1 2025.12.31 -job Hatchery_jobs -get_all True > Diaz_Hatchery_jobs6.log
