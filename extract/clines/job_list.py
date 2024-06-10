"""
Module of functions to create job definitions for a cline extraction.

kh: LO_user 
"""

def get_cbox(job, Lon, Lat):
    vn_list = 'h,f,pm,pn,mask_rho,salt,temp,zeta' # default list
    # specific jobs
    if job == 'shelf_box':
        aa = [-127.5, -123.5, 43, 50]
        vn_list = 'h,f,pm,pn,mask_rho,salt,temp,zeta'

    return aa, vn_list
