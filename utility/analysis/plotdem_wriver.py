import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('classic')
import matplotlib as mpl
import xarray as xr
import basinpaint
from scipy.stats import gaussian_kde

# Before to run you need to wrap fortran libraries (basinpaint.f90): python -m numpy.f2py -c basinpaint.f90 -m basinpaint 
FILENAME = '../../output/Po_test.static_fields.nc'

# Longitude and Latitude indexes of the watershed you want to plot (you can take these values using 
# ncview moving with the mouse cursor over the point of the river you want)
RLON = 297
RLAT = 115

# Min and Maximum values for the colorbar
CMIN = -800
CMAX = 3000

# This threshold (units: drained points) is used to 
# discriminates the river points you want to plot over 
# the DEM  
THRE = 200

data = xr.open_dataset(FILENAME)
nlat=data['dem'].shape[0]
nlon=data['dem'].shape[1]
demval=data['dem'][:,:]
fdmval=data['fdm'][:,:]
lusval=data['lus'][:,:]
draval=data['dra'][:,:]
demftype=demval.values.astype('float32',order='F')
draftype=draval.values.astype('float32',order='F')
fdmftype=fdmval.values.astype('float64',order='F')
lusftype=lusval.values.astype('float64',order='F')
dembas=basinpaint.basinpaint(0,RLON,RLAT,demftype,fdmftype,lusftype,nlat,nlon)
drabas=basinpaint.basinpaint(THRE,RLON,RLAT,draftype,fdmftype,lusftype,nlat,nlon)
demval.values=dembas
dravalm = np.ma.masked_where(drabas < THRE , drabas)

print dravalm[1,1]
pltdem = plt.imshow(demval,origin='lower',cmap='terrain',vmin=CMIN,vmax=CMAX)
plt.colorbar(extend='both')
pltdra = plt.imshow(dravalm, interpolation='none',origin='lower', cmap='Blues',vmin=0,vmax=10)
plt.show()
