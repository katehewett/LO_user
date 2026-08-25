"""
Module of functions to create job definitions for a box test extraction.

kh: LO_user 
"""

def get_cbox(job, Lon, Lat):
    vn_list = 'h,f,pm,pn,mask_rho,salt,temp,zeta' # default list
    # specific jobs
    if job == 'LO_domain':
        aa = [Lon[0], Lon[-1], Lat[0], Lat[-1]]
        vn_list = 'h,f,pm,pn,lat_rho,lon_rho,mask_rho,salt,temp,NH4,NO3,oxygen,TIC,alkalinity,zeta'
    if job == 'PNW_shelf':
        aa = [-128, -123.5, 42, 50]
        vn_list = 'h,f,pm,pn,lat_rho,lon_rho,mask_rho,salt,temp,oxygen,zeta'
    return aa, vn_list
