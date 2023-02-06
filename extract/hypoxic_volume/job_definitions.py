"""
Module of functions to create job definitions for a box extraction.
"""

def get_box(job, Lon, Lat):
    vn_list = 'h,f,pm,pn,mask_rho,salt,temp,rho,zeta,u,v,ubar,vbar' # default list
    # specific jobs
    if job == 'LO_oxygen_WA':
        aa = [-126, -122.5, 46, 49]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,oxygen'
        
    return aa, vn_list
