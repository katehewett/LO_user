"""
Module of functions to create job definitions for a box extraction.
"""

def get_vol(job, Lon, Lat):
    vn_list = 'h,pm,pn,mask_rho,salt,temp,rho,oxygen,zeta' # default list
    # specific jobs
    if job == 'LO_oxygen_WA':
        aa = [-126, -122.5, 46, 49]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,rho,oxygen,zeta'
        
    return aa, vn_list
