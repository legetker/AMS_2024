"""
Code used to calculate 30-day rolling percentile values for the JJA period based 
on a daily ERA5 field. For instance, I might use the code to calculate 98th percentile
max daily temperatures over the period 1980-2010. I am using a rolling 30-day distribution.
To calculate the percentile for June 15th, I would use dates from June 1stto June 30th 
over the entire 1980-2010 period, then use all these values to calculate a normal distribution.
Author: Lauren Getker, NCSU, designed for use on Cheyenne
"""

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

#Output pathname
pathname = "/glade/work/lgetker/daily_percentiles_NA/"  

#Westernmost longitude for subsetting. Choose -180 to 180
wlon = 360-180

#Easternmost latitude for subsetting. Choose -180 to 180
elon = 360-40

#Southernmost latitude for subsetting. Choose -90 to 90
slat = 10

#Northernmost latitude for subsetting.Choose -90 to 90
nlat = 75

#Percentiles to be calculated
percentiles = [90, 95, 98, 99]
#percentiles = [98]

#Variables to use
variables = sorted(["max_t", "min_t", "max_hi", "min_hi", "max_hi_or_temp", "min_hi_or_temp"])
#variables = ['max_hi']
for v in range(len(variables)):
    files = sorted(glob.glob("/glade/work/lgetker/daily_" + variables[v] + "/198*MJJAS*") \
                        + glob.glob("/glade/work/lgetker/daily_" + variables[v] + "/199*MJJAS*") \
                        + glob.glob("/glade/work/lgetker/daily_" + variables[v] + "/200*MJJAS*") \
                        + glob.glob("/glade/work/lgetker/daily_" + variables[v] + "/2010*MJJAS*"))
    print(files)
    for p in range(len(percentiles)):
        #Output filename
        filename = "daily_" + str(percentiles[p]) + "th_per_" + variables[v] + "_NA.nc"
        """
        open one file to grab the grid shape
        """
        ds_t = xr.open_dataset(files[0])
        ds_t = ds_t.where((ds_t['lat']<nlat) & (ds_t['lat']>slat) &  (ds_t['lon']>wlon) & (ds_t['lon']<elon), drop = True)
        sample_lat = len(ds_t['lat'].values)
        sample_lon = len(ds_t['lon'].values)

        """
        Percentile calculation
        """
        all_data = np.zeros([num_days_MJJAS, num_years, sample_lat, sample_lon]) #Array to store all data over the 1980-2010 period
        all_times = []
        year_it = 0 #iterator
        for i in range(len(files)): #Go through every file, open it, and get the times
            ds = xr.open_dataset(files[i])
            times = ds['time'].values
            day_it = 0 #iterator
            for j in range(len(times)): #Go through every time and read the data into the "all_data" array
                ds_t = ds.sel(time = times[j])
                lons = ds_t['lon'].values
                lats = ds_t['lat'].values
                if variables[v] == "max_t":
                    val = ds_t['2t_max']
                elif variables[v] == "min_t":
                    val = ds_t['2t_min']
                elif variables[v] == "max_hi":
                    val = ds_t['2hi_max']                
                elif variables[v] == "min_hi":
                    val = ds_t['2hi_min']
                elif variables[v] == "max_hi_or_temp" or variables[v] == "min_hi_or_temp":
                    val = ds_t['temp_or_hi']
                all_data[day_it, year_it, :, :] = val
                all_times.append(times[j])
                day_it += 1
            year_it += 1

        rolling_percentile = np.zeros([num_days_JJA, sample_lat, sample_lon]) #Array to store the percentiles in
        it = 0
        for i in range(32,124): # Iterate over the JJA period, excluding May and September
            start_index = i - 15  # Index of the starting date
            end_index = i + 15  # Index of the ending date
            # Extract the temperature data for the given period
            period_data = all_data[start_index:end_index + 1, :, :, :]
            period_data = period_data.reshape([period_data.shape[0] * 31, sample_lat, sample_lon])
            rolling_percentile[it,:,:] = np.percentile(period_data, percentiles[p], axis = 0)
            it += 1

        dates = pd.date_range(start=pd.Timestamp('2000-06-01'), end=pd.Timestamp('2000-08-31'), periods = 92) #Usincd ..g "dummy" dates
        rolling_per_da = xr.DataArray(rolling_percentile, coords=[dates, lats, lons], dims=['time', 'lat', 'lon'], name= variables[v] + "_" + str(percentiles[p]))
        rolling_per_da.to_netcdf(path=pathname + filename, mode = 'w')

