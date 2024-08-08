# README for extract_corrosive_volume_rev1.py 
# and extract_corrosive_box.py

#### These scripts use LO history (or lowpass) files and run PCO2SYS to find aragonite saturation state (ARAG). 
#### Using (3) threshold values for aragonite saturation state (ARAG<0.5; <1; <1.7) we flag zones with "corrosive" water. 
#### The thickness of corrosive water for each threshold are saved along with cell area (DA), which are used to get a 
#### corrosive volume at each lat,lon.  

- The script extract_corrosive_volume_rev1.py calls get_one_corrosive_volume_rev1.py. Execution of this script will save corrosive volumes at three thresholds for the entire LO domain using the grid, dates, and list-type entered as command line arguments. The driver points to history (or lowpass) files, which are are fed to get_one_corrosive_volume_rev1.py. At that step, ARAG is calculated, and corrosive dz's are found for each of the 3 thresholds. Because PCO2SYS needs to get one layer at a time (for large arrays) there is a slight slow down in processing time when compared to get_hypoxic_volume.py. For jobs where the whole domain is not required, consider using extract_corrosive_box.py.

- The script extract_corrosive_box.py also requires a job name to be entered in the command line. The job name can be entered by a user "jobs_list.py", or an existing job can be called. This is similar to extract_box or extract_box_chunks.py. First, the driver will use ncks to clip varaibles to the lat/lon listed in jobs_list. Temporary box files are created, which are fed to get_one_corrosive_volume_rev2.py. ARAG is then calculated, and corrosive dz's are found for each of the 3 thresholds. This script was developed because it allowed the runs to process faster with a smaller domain, and because the code was developed for an OA indicators project focused on the WA/OR shelf, we clipped the remainder of the LO domain. 

- The final step in both drivers (_volume.rev1.py and _box.py) is to concatonate all the temporary volume files created in get_one_corrosive_volume.py; save a final compressed .nc file; and then cleanup/delete the temporary (box and) volume files. 







