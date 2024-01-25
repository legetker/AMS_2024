"""
The purpose of this code is to calculate gridded daily max heat indices based on 
hourly ERA5 files, and then save the fields out as netCDF files.In this case I 
am subsetting the data over the North American continent.
Author: Lauren Getker, NCSU, designed for use on Cheyenne
"""

"""
Package imports
"""
import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime
import metpy
import metpy.calc as mpcalc
from metpy.units import units
import glob

"""
USER SETTINGS
"""
#Path to ERA5 data
era5_path = "/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/"

#Westernmost longitude for subsetting. Choose -180 to 180
wlon = 360-180

#Easternmost latitude for subsetting. Choose -180 to 180
elon = 360-40

#Southernmost latitude for subsetting. Choose -90 to 90
slat = 10

#Northernmost latitude for subsetting.Choose -90 to 90
nlat = 75

#Number of years
years = 2022 - 1940

#Days in MJJAS season
days_in_season = 153

#Output pathname
pathname = "/glade/work/lgetker/daily_min_t/"  

#Output filename
filename = "_daily_mins_MJJAS_NA.nc"

"""
Pathnames
Grabbing MJJAS ERA5 files from the Research Data Archive
"""
may_t = sorted(glob.glob(era5_path + "????05/e5.oper.an.sfc.128_167_2t.ll025sc.*0100_*.nc"))
june_t = sorted(glob.glob(era5_path + "????06/e5.oper.an.sfc.128_167_2t.ll025sc.*0100_*.nc"))
july_t = sorted(glob.glob(era5_path + "????07/e5.oper.an.sfc.128_167_2t.ll025sc.*0100_*.nc"))
august_t = sorted(glob.glob(era5_path + "????08/e5.oper.an.sfc.128_167_2t.ll025sc.*0100_*.nc"))
sep_t = sorted(glob.glob(era5_path + "????09/e5.oper.an.sfc.128_167_2t.ll025sc.*0100_*.nc"))

may_td = sorted(glob.glob(era5_path + "????05/e5.oper.an.sfc.128_168_2d.ll025sc.*0100_*.nc"))
june_td = sorted(glob.glob(era5_path + "????06/e5.oper.an.sfc.128_168_2d.ll025sc.*0100_*.nc"))
july_td = sorted(glob.glob(era5_path + "????07/e5.oper.an.sfc.128_168_2d.ll025sc.*0100_*.nc"))
august_td = sorted(glob.glob(era5_path + "????08/e5.oper.an.sfc.128_168_2d.ll025sc.*0100_*.nc"))
sep_td = sorted(glob.glob(era5_path + "????09/e5.oper.an.sfc.128_168_2d.ll025sc.*0100_*.nc"))

all_t = sorted(may_t + june_t + july_t + august_t + sep_t)
all_td = sorted(may_td + june_td + july_td + august_td + sep_td)

"""
open one file to grab the grid shape
"""
ds_t = xr.open_dataset(all_t[0])
ds_t = ds_t.where((ds_t['latitude']<nlat) & (ds_t['latitude']>slat) &  (ds_t['longitude']>wlon) & (ds_t['longitude']<elon), drop = True)
sample_lat = len(ds_t['latitude'].values)
sample_lon = len(ds_t['longitude'].values)

"""
Max value calculations
"""
day_it = 0 #iterator
month_it = 0 #iterator
daily_times = [] #list of times
maxes = np.zeros([days_in_season, sample_lat, sample_lon]) #array of max values for each day
maxes.fill(999)
for i in range(len(all_t)): #go through each monthly data file, subset the data and get the grid
    ds_t = xr.open_dataset(all_t[i])
    ds_t = ds_t.where((ds_t['latitude']<nlat) & (ds_t['latitude']>slat) &  (ds_t['longitude']>wlon) & (ds_t['longitude']<elon), drop = True)
    
    ds_td = xr.open_dataset(all_td[i])
    ds_td = ds_td.where((ds_td['latitude']<nlat) & (ds_td['latitude']>slat) &  (ds_td['longitude']>wlon) & (ds_td['longitude']<elon), drop = True)

    times = ds_t['time'].values
    lats = ds_t['latitude'].values
    lons = ds_t['longitude'].values
    
    iterator = 0
    for j in range(len(times)): #go through each timestep, checking if a given grid cell hit its max heat index
        if iterator <= 22: #check if we have gone through a whole day
            pd_time =  pd.to_datetime(times[j])
            ds_t_time = ds_t.sel(time = times[j], method = 'nearest')
            ds_td_time = ds_td.sel(time = times[j], method = 'nearest')
            
            t = ds_t_time['VAR_2T'].values * units('K')
            t = t.to(units.degC)
            td = ds_td_time['VAR_2D'].values * units('K')
            td = td.to(units.degC)
            # Using metpy to calculate heat index
            rh = metpy.calc.relative_humidity_from_dewpoint(t, td) 
            hi = mpcalc.heat_index(t, rh, mask_undefined=False)
            hi = hi.to(units.degC)
            
            hi = t.magnitude  #CHANGE

            maxes[day_it,:,:] = np.where(np.less(hi, maxes[day_it,:,:]), hi, maxes[day_it, :,:]) #check if the value at each cell exceeds the max/min
            iterator += 1
        else:
            daily_times.append(times[j])
            iterator = 0
            day_it += 1
    month_it += 1
    if month_it == 5: #save out yearly files at netcdf when we get through a single year worth of data
        daily_times_np = np.array(daily_times)
        year = pd_time.year
        xr_maxes = xr.DataArray(data = maxes, 
                                coords=[daily_times_np, lats, lons], 
                                dims=['time', 'lat', 'lon'],
                                attrs=dict(
                                    description="Daily min 2m temperature",
                                    units="degC")
                                )
        xr_maxes = xr_maxes.rename("2t_min")
        xr_maxes.to_netcdf(pathname + str(year) + filename, mode = 'w', format = 'NETCDF4')

        maxes = np.zeros([days_in_season, sample_lat, sample_lon])
        maxes.fill(999)
        month_it = 0
        day_it = 0
        daily_times = []
