"""
Authors: Jacob Radford, CIRA; and Lauren Getker, NCSU
"""
#Package imports
import cartopy
import pygrib as pg
import numpy as np
import cv2
from geopy.distance import geodesic
import xarray as xr
import pandas as pd
import glob
import matplotlib.pylab as plt
"""
USER SETTINGS
"""
#Percentile value used as a threshold
percentile = 98

#Variable to be used as a threshold: temp, hi, or hit
value = "hit"

#Choose a threshold method: "max_only", "min_only", "max_and_floor", "max_min_and_floor", or "floor_only"
method = "max_min_and_floor"

#Only identify heatwaves over the continents?
use_lsm = True

#Minimum floor value, F
floor_value = 32

#Path to ERA5 surface variables
era5_path = "/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/"

#Westernmost longitude for subsetting. Choose -180 to 180
wlon = 360-180

#Easternmost latitude for subsetting. Choose -180 to 180
elon = 360-40

#Southernmost latitude for subsetting. Choose -90 to 90
slat = 10

#Northernmost latitude for subsetting.Choose -90 to 90
nlat = 75

#File location of max and min temp/hi/HIT values
max_min_path = '/glade/work/lgetker/'

#File location of gridded percentile values
per_path = "/glade/work/lgetker/daily_percentiles_NA/"

#Minimum heat wave area (square kilometers)
areaThreshold = 200000

#Output pathname
output_path = "/glade/work/lgetker/detected_hws/"

#Output filename
#filename = "hw_83yr_" + str(percentile) + "_" + method + "_" + str(areaThreshold) +  "_data.nc"
#filename = "floor_only_" + str(areaThreshold) + "_32hi_data.nc"
filename = "98_per_" + str(areaThreshold) + "_hit_min_and_floor_data.nc"

#Convert longitudes to 0-360 
if elon < 0:
    elon = elon + 360
if wlon < 0:
    wlon = wlon + 360

if value == "hit":
    all_t = sorted(glob.glob(max_min_path +"daily_max_hi_or_temp/*")) #daily max hit values
    t_min = sorted(glob.glob(max_min_path +"daily_min_hi_or_temp/*")) #daily min hit values
    max_per_ds = xr.open_dataset(per_path + "daily_" + str(percentile) + "th_per_max_hi_or_temp_NA.nc") #max percentile
    min_per_ds = xr.open_dataset(per_path + "daily_" + str(percentile) + "th_per_min_hi_or_temp_NA.nc") #min percentile
elif value == "temp":
    all_t = sorted(glob.glob(max_min_path + "daily_max_t/*")) #daily max hit values
    t_min = sorted(glob.glob(max_min_path +"daily_min_t/*")) #daily min hit values
    max_per_ds = xr.open_dataset(per_path + "daily_" + str(percentile) + "th_per_max_t_NA.nc") #max percentile
    min_per_ds = xr.open_dataset(per_path + "daily_" + str(percentile) + "th_per_min_t_NA.nc") #min percentile
elif value == "hi":
    all_t = sorted(glob.glob(max_min_path + "daily_max_hi/*")) #daily max hit values
    t_min = sorted(glob.glob(max_min_path +"daily_min_hi/*")) #daily min hit values
    max_per_ds = xr.open_dataset(per_path + "daily_" + str(percentile) + "th_per_max_hi_NA.nc") #max percentile
    min_per_ds = xr.open_dataset(per_path + "daily_" + str(percentile) + "th_per_min_hi_NA.nc") #min percentile    
else:
    sys.exit("Incorrect value name.")

"""
Get landmask for subsetting data
"""
landmask = xr.open_dataset("/glade/work/lgetker/e5.oper.invariant.128_172_lsm.ll025sc.1979010100_1979010100.nc")
landmask = landmask.where((landmask['latitude']<nlat) & (landmask['latitude']>slat) &  (landmask['longitude']>wlon) & (landmask['longitude']<elon), drop = True)
lsm_times = landmask['time'].values
landmask = landmask.sel(time = lsm_times[0])
landmask = landmask['LSM'].values

all_dates = []
centroids_lat = []
centroids_lon = []
areas = []
masks_to_save = []
times_it = 0
it = 0
for j in range(len(all_t)):  #Iterate through each year and open files
    ds_t_full = xr.open_dataset(all_t[j])
    ds_t_full = ds_t_full.where((ds_t_full['lat']<nlat) & (ds_t_full['lat']>slat) &  (ds_t_full['lon']>wlon) & (ds_t_full['lon']<elon), drop = True)
    
    ds_min_full = xr.open_dataset(t_min[j])
    ds_min_full = ds_min_full.where((ds_min_full['lat']<nlat) & (ds_min_full['lat']>slat) &  (ds_min_full['lon']>wlon) & (ds_min_full['lon']<elon), drop = True)
    
    times = ds_t_full['time'].values
    for i in range(len(times)): #Iterate through each day
        time =  pd.to_datetime(times[i])
        generic_time = pd.Timestamp(year = 2000, month = time.month, day = time.day) #Percentiles are saved with a "dummy" time of 2000
        if time.month != 6 and time.month != 7 and time.month != 8: #check that it's a JJA date
            continue
        times_it += 1
        if value == "temp":
            temperatureThreshold = max_per_ds.sel(time = generic_time)['max_t_' + str(percentile)].values
        elif value == "hi":
            temperatureThreshold = max_per_ds.sel(time = generic_time)['max_hi_' + str(percentile)].values
        elif value == "hit":
            temperatureThreshold = max_per_ds.sel(time = generic_time)['max_hi_or_temp_' + str(percentile)].values 
            
        if value == "temp":
            minThreshold = min_per_ds.sel(time = generic_time)['min_t_' + str(percentile)].values
        elif value == "hi":
            minThreshold = min_per_ds.sel(time = generic_time)['min_hi_' + str(percentile)].values
        elif value == "hit":
            minThreshold = min_per_ds.sel(time = generic_time)['min_hi_or_temp_' + str(percentile)].values    
        
        lats = ds_t_full['lat'].values
        lons = ds_t_full['lon'].values
        
        ds_t = ds_t_full.sel(time = times[i], method = 'nearest')
        ds_min = ds_min_full.sel(time = times[i], method = 'nearest')

        if value == "hit":
            t_vals = ds_t['temp_or_hi']
            min_vals = ds_min['temp_or_hi'] 
        elif value == "temp":
            t_vals = ds_t['2t_max']
            min_vals = ds_min['2t_min']
        elif value == "hi":
            t_vals = ds_t['2hi_max']
            min_vals = ds_min['2hi_min'] 
        else:
            sys.exit("Incorrect value name.")
        
        # Convert longitudes from 0-360 range to -180-180 range
        lons[lons > 180] -= 360
        lons_sorted = lons
        lats_sorted = lats

        max_max_t = np.max(t_vals)
        min_min_t = np.min(t_vals)
        
        #Some arrays for checking masks
        all_zeros = np.zeros_like(t_vals)
        nan_arr = np.zeros_like(t_vals)
        nan_arr[nan_arr == 0] = np.nan
        ones_arr = np.zeros_like(t_vals)
        ones_arr[ones_arr == 0] = 1

        #Create a copy of the temeprature values array, then apply masks as specified by the user
        temperatureCopy = t_vals.values
        if method == "min_only" or method == "max_min_and_floor":
            print("applying min")
            temperatureCopy[min_vals<minThreshold] = np.nan
        if use_lsm == True:
            print("using lsm")
            temperatureCopy[landmask == 0] = np.nan
        if method == "max_only" or method == "max_min_and_floor":
            print("applying max")
            temperatureCopy[temperatureCopy<temperatureThreshold] = np.nan
        if method == "max_and_floor" or method == "max_min_and_floor" or method == "floor_only":
            print("applying floor")
            temperatureCopy[temperatureCopy <= floor_value] = np.nan

        #Threshold at temperature threshold, but this time set values above threshold to 1 and everything else to nan
        binaryThreshold = np.where((np.isnan(temperatureCopy)), all_zeros, ones_arr) 
        binaryThreshold[binaryThreshold==0] = np.nan
        #This section is the main stuff and can take a little while
        #This changes array data type to 8-bit integer
        #It also converts nan values back to 0
        binaryThreshold = np.array(binaryThreshold.astype('uint8'))

        #OpenCV to identify the contours
        contours, _ = cv2.findContours(binaryThreshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        #Array to store the masks in
        masks = np.zeros([len(contours) + 1, binaryThreshold.shape[0], binaryThreshold.shape[1]])
        mask_it = 0
        for contour in contours:
            #This identifies all of the points contained in a contour
            mask = np.zeros_like(binaryThreshold)
            mask_copy = np.zeros_like(mask)
            cv2.drawContours(mask,[contour],0,255,-1)
            pts = cv2.findNonZero(mask)

            #Go through all of the points within the contour and calculate the points area
            contourarea = 0
            all_lats = []
            all_lons =[]
            for pt in pts:
                x = pt[0][0]
                y = pt[0][1]

                #Latitude and longitude of center of grid cell
                ptlat = lats_sorted[y]
                all_lats.append(ptlat)
                ptlon = lons_sorted[x]
                all_lons.append(ptlon)

                #Latitudes of top and bottom of grid cell
                latplus = ptlat + 0.125
                latminus = ptlat - 0.125

                #Longitudes of left and right of grid cell
                lonplus = ptlon + 0.125
                lonminus = ptlon - 0.125

                #Geodesic distance between top and bottom and left and right side of grid cells
                gridptwidth = geodesic((ptlat,lonplus), (ptlat,lonminus)).kilometers
                gridptheight = geodesic((latplus,ptlon), (latminus,ptlon)).kilometers

                #Calculate area of grid cell
                gridptarea = gridptwidth*gridptheight

                #Add the grid cell's area to the contour area
                contourarea+=gridptarea

            #Check if the contour area is greater than some threshold
            #If so, keep the mask
            if (contourarea>areaThreshold):
                filteredData = np.copy(t_vals)
                filteredData[mask == 0] = np.nan

                masks_to_save.append(filteredData)
                areas.append(contourarea)
                centroids_lat.append(np.mean(all_lats))
                centroids_lon.append(np.mean(all_lons))
                all_dates.append(time)
                mask_it += 1

all_dates_np = np.array(all_dates)
masks_np = np.array(masks_to_save)
areas_np = np.array(areas)
centroids_lat_np = np.array(centroids_lat)
centroids_lon_np = np.array(centroids_lon)

print(masks_np.shape)
print(all_dates_np.shape)
print(centroids_lat_np.shape)
print(areas_np.shape)

ds = xr.Dataset(
    {
        'climate_data': (['time', 'lat', 'lon'], masks_np),
        'area': (['time'], areas_np),
        'centroid_lats' : (['time'], centroids_lat_np),
        'centroid_lons' : (['time'], centroids_lon_np)
    },
    coords={'time': all_dates_np, 'lat': lats, 'lon': lons},
    attrs = dict(
            description="Heat wave masks",
            units="degC")
    )
ds.to_netcdf(output_path + filename, mode = 'w', format = 'NETCDF4')
