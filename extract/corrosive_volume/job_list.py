"""
Module of functions to create job definitions for a corrosive box extraction.

kh: LO_user 
"""

def get_cbox(job, Lon, Lat):
    vn_list = 'h,f,pm,pn,mask_rho,salt,temp,zeta' # default list
    # specific jobs
    if job == 'OA_indicators':
        aa = [-125.5, -123.5, 42.75, 48.75]
        vn_list = 'h,f,pm,pn,lat_rho,lon_rho,mask_rho,salt,temp,alkalinity,TIC,zeta'

    return aa, vn_list
