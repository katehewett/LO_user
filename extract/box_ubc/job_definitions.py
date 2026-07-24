"""
Module of functions to create job definitions for a box extraction.
This is a special case for UBC in response to a request for data from Susan/Doug at UBC

kh: LO_user 
"""

def get_box(job, Lon, Lat):
    vn_list = 'h,f,pm,pn,mask_rho,salt,temp,zeta,u,v,ubar,vbar' # default list
    # specific jobs
    if job == 'ubc0':
        aa = [-125.016452048434, -124.494612925929, 48.312, 48.7515055163539]
        # old version
        # vn_list = ('h,f,pm,pn,mask_rho,salt,temp,zeta,NO3,phytoplankton,'
        #         + 'zooplankton,detritus,Ldetritus,oxygen,TIC,alkalinity')
        # new version
        vn_list = ('h,f,pm,pn,mask_rho,salt,temp,zeta,NO3,NH4,phytoplankton,'
                + 'zooplankton,SdetritusN,LdetritusN,oxygen,TIC,alkalinity')
    elif job == 'tester':
        aa = [-125, -124, 47, 49]  
        vn_list = 'h,pm,pn,mask_rho,salt,zeta'

    return aa, vn_list
