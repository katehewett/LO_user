# get indices for plotting // ilon and ilat will be the index of the mlon mlat inputs::
#mlon = -124.502
#mlat = 45.135

mlon = -124.06
mlat = 45.135

alat = G['lat_rho']
alon = G['lon_rho'] 

Lon = alon[0,:]
Lat = alat[:,0]

# error checking
if (mlon < Lon[0]) or (mlon > Lon[-1]):
    print('ERROR: lon out of bounds ')
    sys.exit()
if (mlat < Lat[0]) or (mlat > Lat[-1]):
    print('ERROR: lat out of bounds ')
    sys.exit()
# get indices
ilon = zfun.find_nearest_ind(Lon, mlon)
ilat = zfun.find_nearest_ind(Lat, mlat)

'''NT, NZ, NR, NC = np.shape(dsb.salt)
zeta = dsb.zeta.values
h = dsb.h.values

z_rho, z_w = zrfun.get_z(h, zeta, S) # (30, 1111, 356) (31, 1111, 356)
NW = z_w.shape[0]

#### isolate the one profile to test method #### 
z_w = z_w[:,ilat,ilon]
z_rho = z_rho[:,ilat,ilon]

tempC = dsb.temp.values[:,:,ilat,ilon]
SP = dsb.salt.values[:,:,ilat,ilon]'''

plat = alat[ilat,ilon]
plon = alon[ilat,ilon]