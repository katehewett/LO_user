"""
Module of functions to create job definitions for a corrosive box extraction.

kh: LO_user 
"""

def get_cbox(job, Lon, Lat):
    vn_list = 'h,pm,pn,zeta,mask_rho,salt,temp,zeta' # default list
    # specific jobs
    if job == 'OA_indicators':
        aa = [-125.5, -123.5, 42.75, 48.75]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,alkalinity,TIC'

    return aa, vn_list
