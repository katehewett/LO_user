"""
Module of functions to create job definitions for a box extraction.

kh: LO_user 
"""

def get_box(job, Lon, Lat):
    vn_list = 'h,f,pm,pn,mask_rho,salt,temp,zeta,u,v,ubar,vbar' # default list
    # specific jobs
    if job == 'sequim0':
        aa = [-123.15120787, -122.89090010, 48.07302111, 48.19978336]
    elif job == 'taiping_hc':
        aa = [-122.66394, -122.61417, 47.93171, 47.94398]
        vn_list = 'h,f,pm,pn,mask_rho,salt,temp,zeta,oxygen,phytoplankton,NO3'
    elif job == 'PS':
        # 3 MB per save (26 GB/year for hourly)
        aa = [-123.5, -122.05, 47, 49]
    elif job == 'garrison':
        aa = [-129.9, -122.05, 42.1, 51.9]
        vn_list = 'h,f,pm,pn,mask_rho,salt,temp,oxygen'
    elif job == 'full':
        aa = [Lon[0], Lon[-1], Lat[0], Lat[-1]]
        vn_list = 'h,f,pm,pn,mask_rho,salt,temp,oxygen'
    elif job == 'liu_wind':
        aa = [-123.5, -122.05, 47, 48.5]
        vn_list = 'h,mask_rho,Uwind,Vwind'
    elif job == 'liu_ps':
        aa = [-123.3, -122.2, 47, 49]
        vn_list = 'h,f,pm,pn,mask_rho,salt,temp,Uwind,Vwind,shflux'
    elif job == 'surface0':
        aa = [Lon[0], Lon[-1], Lat[0], Lat[-1]]
        # For reasons I do not understand, this gets zeta even when it is not on the
        # list.  I will put it here to be explicit.
        vn_list = 'h,f,pm,pn,mask_rho,salt,temp,zeta'
    elif job == 'surface1':
        # For Samantha 2021.12.06
        aa = [Lon[0], Lon[-1], Lat[0], Lat[-1]]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,sustr,svstr,zeta' 
    elif job == 'ubc0':
        aa = [-125.016452048434, -124.494612925929, 48.312, 48.7515055163539]
        vn_list = ('h,f,pm,pn,mask_rho,salt,temp,zeta,NO3,phytoplankton,'
                + 'zooplankton,detritus,Ldetritus,oxygen,TIC,alkalinity,CaCO3,rho')
    elif job == 'cox':
        aa = [-123.204529, -122.728532, 48.393771, 48.726895]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,phytoplankton'
    elif job == 'jerry0':
        aa = [-122.52, -122.40, 47.40, 47.85]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,w,zeta,u,v,ubar,vbar,Uwind,Vwind'
    elif job == 'pisces0':
        aa = [-127, -124, 46, 48]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,zeta'
    elif job == 'desanto':
        aa = [-125.028, -124.8993, 45.2581, 45.3481]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,zeta,u,v'
    elif job == 'gheibi':
        aa = [-123.15, -122.84, 48.68, 48.95]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,zeta,u,v'
    elif job == 'byrd':
        aa = [Lon[0], Lon[-1], Lat[0], Lat[-1]]
        vn_list = 'h,mask_rho,Uwind,Vwind,u,v'
    elif job == 'barbanell':
        aa = [-125.5, -122.1, 47, 50.3]
        vn_list = 'h,mask_rho,temp'
    elif job == 'bass':
        aa = [-126.5, -124, 48.4, 49]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,oxygen'
    elif job == 'bass2':
        aa = [-125.5, -122.1, 47, 49]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,oxygen'
    elif job == 'harcourt':
        aa = [-125.6, -124.2, 46.6, 47.2]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,oxygen,zeta,u,v,w,Uwind,Vwind'
    elif job == 'NHL1_65':
        aa = [-124.0999985, -125.5999985,   45.65169907,   43.65169907]
        vn_list = 'h,f,pm,pn,mask_rho,salt,temp,rho,oxygen,zeta,u,v' 
    elif job == 'NHL1_200':
        aa = [-124.0999985, -128.7720032,   45.65169907,   43.65169907]
        vn_list = 'h,f,pm,pn,mask_rho,salt,temp,rho,oxygen,zeta,u,v'
    elif job == 'LO_oxygen_WA':
        aa = [-126, -122.5, 46, 49]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,oxygen'
    elif job == 'Full_plusO2Chem':
        aa = [-129.9798, -122.00184, 42.00675, 52] # using the Lon[0], Lon[-1], Lat[0], Lat[-1] doesn't work... - fix
        vn_list = 'h,pm,pn,mask_rho,zeta,salt,temp,oxygen,TIC,alkalinity'
    elif job == 'Freeman_OCNMS':
        aa = [-126.2, -123.8, 46.6, 48.8]  
        vn_list = 'h,pm,pn,mask_rho,zeta,salt,temp,phytoplankton,zooplankton,NO3,oxygen,TIC,alkalinity,Uwind,Vwind' 
    elif job == 'NHL_transect':
        aa = [-124.7, -124.1, 44.64, 44.66]  
        vn_list = 'h,pm,pn,mask_rho,zeta,u,v,w,salt,temp,AKs,AKv,sustr,svstr,bustr,bvstr,Uwind,Vwind'
    elif job == 'NHL_bio':
        aa = [-124.7, -124.1, 44.64, 44.66]  
        vn_list = 'h,pm,pn,mask_rho,zeta,NO3,NH4,phytoplankton,zooplankton,SdetritusN,LdetritusN,SdetritusC,LdetritusC,oxygen,alkalinity,TIC,rho'
    elif job == 'shelf_box':
        aa = [-127.5, -123.5, 43, 50]
        vn_list = 'h,f,pm,pn,zeta,mask_rho,salt,temp'
    elif job == 'OA_indicators':
        aa = [-125.5, -123.5, 42.75, 48.75]
        vn_list = 'h,pm,pn,mask_rho,salt,temp,alkalinity,TIC'
    elif job == 'OA_indicators_phys':
        aa = [-125.5, -123.5, 42.75, 48.75]
        vn_list = 'h,pm,pn,mask_rho,zeta,u,v,sustr,svstr,bustr,bvstr,Uwind,Vwind'
    elif job == 'aakritiv': # student extraction for surface nitrate
        aa = [-125.3, -123.8, 46.5, 48]
        vn_list = 'h,pm,pn,mask_rho,NO3'
    elif job == 'CE_transect':
        aa = [-124.731, -124.348, 47.353, 47.353]  
        vn_list = 'h,pm,pn,mask_rho,zeta,u,v,w,salt,temp,AKv,oxygen,alkalinity,TIC'
    elif job == 'CE_transect_surf':
        aa = [-124.731, -124.348, 47.353, 47.353]  
        vn_list = 'h,pm,pn,mask_rho,zeta,sustr,svstr,Uwind,Vwind'
    elif job == 'CE_transect_bio':
        aa = [-124.731, -124.348, 47.353, 47.353]  
        vn_list = 'h,pm,pn,mask_rho,zeta,NO3,NH4,phytoplankton,zooplankton,SdetritusN,LdetritusN,SdetritusC,LdetritusC'
    elif job == 'tester':
        aa = [-125, -124, 47, 49]  
        vn_list = 'h,pm,pn,mask_rho,salt,zeta'

    return aa, vn_list
