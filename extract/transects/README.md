# README for the transect code

### This is code to do transect extractions. It is based off the LO/extract/tef code 

#### WORKFLOW OVERVIEW
**NOTE: when copied LO/extract/tef folder, the scripts hadn't been updated to use xarray --> need to update from what it uses here (netCDF4 to xarray**

Testing using NHL transect data (NH-1 thru NH-45) 
'extract_sections.py' was used to form 'extract_transect.py'

---
#### PREAMBLE

in tef, Parker has: 
`LO/pgrid/tef_section_maker.py` is a tool to define sections, if needed.  This is a GUI where you define section endpoints (forced to be either perfectly E-W or N-S), name the section, and define which direction is "landward" (i.e. which way is positive transport).  You can define one or more sections in a session.  At the end, when you push DONE it produces screen output suitable for pasting into new or existing lists in `tef_fun.get_sect_df()`.

`tef_fun.py` is an important module for this code. It includes `tef_fun.get_sect_df()` which returns a DataFrame of all the section names and their lon,lat endpoints.

here... 'LO_user/extract/transects/' has NHL stations entered in transect_fun

---
#### EXTRACTION CODE

`extract_transects.py` should(still testing?) create a NetCDF file for a section with arrays tracer values on the section, arranged as (t, z, x-or-y). Using command line arguments you can change the run, the day range, the sections to extract, and the variables extracted. Typically this will be done on a remote machine, like perigee, although the defaults are designed to work with model output I have saved on my mac.

**NOTE**: this code runs multiple subprocess instances of `extract_section_one_time.py`, (set by the Nproc command line argument which has a default value of 10). This significantly speeds things up, but it tends to occupy the machine, e.g. if you use -Nproc 20 on perigee you are using all the cores and may slow down other jobs.

**NOTE**: this code also automatically runs the two subsequent steps, `process_sections.py` and `bulk_calc.py`.  These can also be run as stand-alone (use -test True when running `extract_sections.py`) to facilitate debugging.

**PERFORMANCE**: For a full year of cas6, with -get_bio True and all NPZDOC variables this takes 5 hours on perigee (6.5 hours when including all steps).

The **command line arguments** are defined in `LO/lo_tools/lo_tools/extract_argfun.py`, with the usual requirements for gtagex, and beginning and end dates.  You also specify:
- -ro, --roms_out_num, is a integer specifying where to look for the ROMS history profiles
- -get_bio is a Boolean.  If True then you get the list in tef_fun.vn_list.  If False (the default) then vn_list = ['salt'].
- -sect_name is a string to specify the sections to get, either a single one such as ai1, or all of them using "-sect_name all" in the command line.  The full list is in tef_fun.get_sect_df().

Input: ROMS history files over some date range, e.g. [*] = 2017.01.01_2017.12.31

Output: LO_output/extract/[gtagex]/tef/extractions_[*]/[sect name].nc where:

Variables in the NetCDF files:
- salt is hourly salinity in each cell (t, z, x-or-y) [same for all other variables]
- q is hourly transport in each cell (t, z, x-or-y)
- vel is velocity in each cell (t, z, x-or-y) positive to East or North
- DA is the area of each cell (t, z, x-or-y) hence: q = vel * DA
- z0 is the average z-position of cell centers (assumes SSH=0), useful for plotting
- DA0 is the average cross-sectional area of each cell (assumes SSH=0)
- h is depth on the section (x-or-y) positive down
- zeta is SSH on the section (t, x-or-y) positive up
- ocean_time is a vector of time in seconds since (typically) 1/1/1970.

---
