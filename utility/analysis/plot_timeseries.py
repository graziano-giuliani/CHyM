import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('classic')
import matplotlib as mpl
import xarray as xr
import basinpaint
from scipy.stats import gaussian_kde
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)

FILENAME1 = '/home/netapp-clima/users/fdi_sant/CHYM/last5Work/chym-esp_v609/output/Po_test0009_2008.nc_por_cat'
FILENAME2 = '/home/netapp-clima/users/fdi_sant/CHYM/last5Work/chym-esp_v609/output/Po_test0035_1.2e-7_2008.nc_por_cat'
FILENAME3 = '/home/netapp-clima/users/fdi_sant/CHYM/last5Work/chym-esp_v609/output/Po_test015_2e-8_2008.nc_por_cat'

# Longitude and Latitude indexes of the watershed you want to plot
RLON1 = 328
RLAT1 = 112
RLON2 = 77
RLAT2 = 30
RLON3 = 19
RLAT3 = 8

data1 = xr.open_dataset(FILENAME1)
data2 = xr.open_dataset(FILENAME2)
data3 = xr.open_dataset(FILENAME3)
porval1=data1['por'][:,RLAT1,RLON1]
porval2=data2['por'][:,RLAT2,RLON2]
porval3=data3['por'][:,RLAT3,RLON3]
time=data1['time']
x = np.linspace(0, time.shape[0], time.shape[0])
fig, ax = plt.subplots(1)
pl1 = ax.plot(time, porval1)
pl2 = ax.plot(time, porval2)
#pl3 = ax.plot(time, porval3)
formatter = DateFormatter('%m/%d/%y-%H')
ax.xaxis.set_major_formatter(formatter)
ax.xaxis.set_tick_params(rotation=30, labelsize=10)
ax.xaxis.set_major_locator(plt.MaxNLocator(10))
#plt.legend((pl1[0],pl2[0],pl3[0]),('0.009','0.035','0.15'),loc='upper left')
plt.legend((pl1[0],pl2[0]),('0.009','0.035'),loc='upper left')
plt.title(porval1.long_name)
plt.show()
