"""
Package imports
"""
import xarray as xr
import numpy as np
import pandas as pd
import csv
from datetime import datetime
import glob

"""
Variable declarations and user settings
"""
#Number of years in the calculation: 1980 to 2010
num_years = 31

#Number of days in MJJAS period
num_days_MJJAS = 153

#Number of days in the JJA period
num_days_JJA = 92

p_level = 500

#Output pathname
era5_path = "/glade/campaign/collections/rda/data/ds633.0/e5.oper.an.pl/"

#Westernmost longitude for subsetting. Choose -180 to 180
wlon = 360-180

#Easternmost latitude for subsetting. Choose -180 to 180
elon = 360-40

#Southernmost latitude for subsetting. Choose -90 to 90
slat = 10

#Northernmost latitude for subsetting.Choose -90 to 90
nlat = 75

#Percentiles to be calculated
#percentiles = [90, 95, 98, 99]
percentiles = [98]

files = sorted(glob.glob(era5_path + "198?05/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "198?06/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "198?07/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "198?08/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "198?09/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "199?05/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "199?06/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "199?07/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "199?08/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "199?09/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "200?05/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "200?06/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "200?07/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "200?08/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "200?09/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "201005/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "201006/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "201007/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "201008/e5.oper.an.pl.128_129_z*.nc") \
               + glob.glob(era5_path + "201009/e5.oper.an.pl.128_129_z*.nc") \
              )
print(len(files))
filename = "era5_500mb_height_mean.nc"
"""
open one file to grab the grid shape
"""
ds_t = xr.open_dataset(files[0])
ds_t = ds_t.where((ds_t['latitude']<nlat) & (ds_t['latitude']>slat) &  (ds_t['longitude']>wlon) & (ds_t['longitude']<elon), drop = True)
lons = ds_t['longitude'].values
lats = ds_t['latitude'].values

def process_file(file):
    ds = xr.open_dataset(file)
    times = ds['time'].values
    ds_t = ds.sel(time = times[0], level = p_level)
    print(ds_t['time'].values)
    ds_t = ds_t.where((ds_t['latitude'] < nlat) & (ds_t['latitude'] > slat) & (ds_t['longitude'] > wlon) & (ds_t['longitude'] < elon), drop=True)
    return ds_t['Z'].values / 9.8, ds_t['time'].values
"""
Percentile calculation
"""
all_data = []
all_times = []
all_data = []
all_times = []
year_it = 0 #iterator
for i, file in enumerate(files): #Go through every file, open it, and get the times
    print(file)
    val, times = process_file(file)
    all_data.append(val)
    all_times.append(times)

all_data_np = np.array(all_data)
all_data_np = np.reshape(all_data_np, [31,153, lats.shape[0], lons.shape[0]]) #SHAPE: YEARS, DAYS, LAT, LON
print(all_data_np.shape)
all_times = np.array(all_times)

rolling_percentile = np.zeros([92, lats.shape[0], lons.shape[0]]) #Array to store the percentiles in
it = 0
for i in range(32,124): # Iterate over the JJA period, excluding May and September
    start_index = i - 15  # Index of the starting date
    end_index = i + 15  # Index of the ending date
    # Extract the temperature data for the given period
    period_data = all_data_np[:,start_index:end_index + 1, :, :]
    period_data = period_data.reshape([period_data.shape[0] * 31, lats.shape[0], lons.shape[0]])
    rolling_percentile[it,:,:] = np.mean(period_data,axis = 0)
    it += 1

filename = "500mb_height_mean.nc"
dates = pd.date_range(start=pd.Timestamp('2000-06-01'), end=pd.Timestamp('2000-08-31'), periods = 92) #Using "dummy" dates
rolling_per_da = xr.DataArray(rolling_percentile, coords=[dates, lats, lons], dims=['time', 'lat', 'lon'], name = "500mb_mean")
rolling_per_da.to_netcdf(path=filename, mode = 'w')
