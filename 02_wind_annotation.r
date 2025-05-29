#Script for downloading and preparing the wind data for Safi et al 2025.
#Each GPS point will be associated with the wind speed at that location in the closest hour
#Elham Nourani, PhD. elham.nourani@unil.ch

#input: GPS data is downloaded from Movebnak(https://www.movebank.org); Wind data is downloaded from the Copernicus Data Store(https://cds.climate.copernicus.eu/)
#output: "your_path/your_downloaded_gps_data.rds", "your_path/your_annotated_gps.rds"

library(tidyverse)
library(lubridate)
library(terra)
library(sf)
library(ncdf4)
library(parallel)

#---------------------------------------------------------------------------------
## Step 0: Download all GPS data                                             #####
#---------------------------------------------------------------------------------

#Make sure you have set up your credentials following Move2 tutorials.
#Data used in this study was downloaded on 15.04.2024. Use this timestamp as the cutoff point when downloading the data to reproduce the results of this study.

#creds <- movebank_store_credentials(username = "your_user_name", rstudioapi::askForPassword())

movebank_retrieve(study_id = 2201086728, entity_type= "tag_type")

gps <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "gps",   #download data for all individuals 
                         entity_type = "event",  attributes = "all",
                         timestamp_end = as.POSIXct("2024-04-15 00:00:00"))

#write to file. Will be used in 03b_data_prep_days.r
saveRDS(gps, file = "your_path/your_downloaded_gps_data.rds")

#---------------------------------------------------------------------------------
## Step 1: prepare tracking data                                             #####
#---------------------------------------------------------------------------------

all_gps_apr <- gps %>%  #2022-09-25 07:56:41.0000 to 2024-04-15 10:02:33.0000
  drop_na(individual_local_identifier, location_lat) %>% #remove NA individuals and NA locations.
  mutate(yr = year(timestamp),
         mn = month(timestamp),
         dy = day(timestamp),
         hr = hour(timestamp),
         unique_hr = paste(yr,mn,dy,hr, sep = "_"),
         closest_hr = round(timestamp, units = "hours") %>% as.character()) %>% 
  as.data.frame()

#split up the data based on year
all_gps_ls <- split(all_gps_apr, all_gps_apr$yr)

#---------------------------------------------------------------------------------
## Step 2: extract wind u and v for each point                               #####
#---------------------------------------------------------------------------------

#This step was done when CDS API was not functioning, so no code is provided for data download. 
#The data were downloaded manually from the website (see download request at the end of script) and one netcdf file was saved per year (2022,2023,2024)

#list the downloaded files
nc_files <- list.files("path_to_your_files", pattern = ".nc", full.names = T)

#decide on an output directory to save the annotated files
output_path <- "path_to_your_output_dir/"

#go over each nc file
lapply(all_gps_ls, function(x){
  
  #extract the year
  yr <- x$yr[1]
  
  ########### read in wind data for this year #####
  #open nc data for this year
  u <- nc_files[grep(yr, nc_files)] %>% 
    rast("u")
  
  v <- nc_files[grep(yr, nc_files)] %>% 
    rast("v")
  
  #I couldn't extract time using the terra::time() function, so I'm extracting them from the layer names
  times <- names(u) %>%
    str_split("time=") %>%
    map_chr(2)
  
  #reference time according to the nc file
  ref_date <- ymd_hms("1970-01-01 00:00:00")
  
  #convert the hours into date + hour
  timestamp <- ref_date + seconds(times)
  
  #rename the layers based on the converted timestamps. these are the same for both u and v. but keep the u and v in the names for clarity when extracting values
  #from the two layers that would otherwise have the same name (ie. the timestamp)
  names(u) <- paste0("u_900_", timestamp)
  names(v) <- paste0("v_900_", timestamp)
  
  ########### split the tracking data into unique hours #####
  x_ls <- split(x, x$closest_hr)
  
  # Define the number of cores to use
  num_cores <- detectCores() - 10 #run on 2 cores (my machine has 12. adjust based on your machine). keep in mind the available RAM when deciding on the number of cores.
  
  wind_this_yr <- mclapply(x_ls, function(y){ #for each hour
    
    #Extract the unique hour
    unique_hr <- y$closest_hr[1]
    
    #Check whether there is any wind data for this hour
    if (any(str_detect(names(u), unique_hr))) {
      # Extract the corresponding rasters
      wind <- c(
        u[[str_which(names(u), unique_hr)]],
        v[[str_which(names(v), unique_hr)]]
      )
      
      # Convert tracking data to a SpatVector
      y_vec <- vect(y, geom = c("location_long", "location_lat"), crs = "EPSG:4326")
      
      # Extract values for each point and append directly to y
      extracted_wind <- extract(x = wind, y = y_vec, method = "bilinear", bind = FALSE, ID = FALSE)
      colnames(extracted_wind) <- c("u_900", "v_900")
      
      # Append extracted values to y
      y_df <- y %>%
        bind_cols(as.data.frame(extracted_wind))
      
    } else {
      # If there are no matching wind data for this hour, create y_df with NA values
      y_df <- y %>%
        mutate(u_900 = as.numeric(NA),
               v_900 = as.numeric(NA))
    }
    
    rm(wind, y)
    
    y_df
    
  }, mc.cores = num_cores) %>% 
    bind_rows()
  
  saveRDS(wind_this_yr, file = paste0(output_path,"gps_annotated_", yr, ".rds"))
  
})

#---------------------------------------------------------------------------------
## Step 3: put all files together and calculate wind speed                   #####
#---------------------------------------------------------------------------------

ann_ls <- list.files("path_to_your_output_dir/", full.names = T) %>% 
  map(readRDS) %>% 
  map(bind_rows) %>% 
  bind_rows() %>% 
  select(1:48) %>% 
  mutate(wind_speed =  sqrt(u_900^2 + v_900^2))

#save to file. This will be used in 03a_data_prep_bursts.r
saveRDS(ann_ls, file = "your_path/your_annotated_gps.rds")

#----------------------------------------------------- 
# ## This step was done by downloading the data from the gui with the following request. The same request can be made through the API

# dataset = "reanalysis-era5-pressure-levels"
# request = {
#   "product_type": ["reanalysis"],
#   "variable": [
#     "u_component_of_wind",
#     "v_component_of_wind"
#   ],
#   "year": ["2022"],
#   "month": [
#     "07", "08", "09", "10",
#     "11", "12"
#   ],
#   "day": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12",
#     "13", "14", "15",
#     "16", "17", "18",
#     "19", "20", "21",
#     "22", "23", "24",
#     "25", "26", "27",
#     "28", "29", "30",
#     "31"
#   ],
#   "time": [
#     "03:00", "04:00", "05:00",
#     "06:00", "07:00", "08:00",
#     "09:00", "10:00", "11:00",
#     "12:00", "13:00", "14:00",
#     "15:00", "16:00", "17:00",
#     "18:00"
#   ],
#   "pressure_level": ["900"],
#   "data_format": "netcdf",
#   "download_format": "zip",
#   "area": [64, -11, -27, 36]
# }
# 
# client = cdsapi.Client()
# client.retrieve(dataset, request).download()
# 
# 
# import cdsapi
# 
# dataset = "reanalysis-era5-pressure-levels"
# request = {
#   "product_type": ["reanalysis"],
#   "variable": [
#     "u_component_of_wind",
#     "v_component_of_wind"
#   ],
#   "year": ["2023"],
#   "month": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12"
#   ],
#   "day": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12",
#     "13", "14", "15",
#     "16", "17", "18",
#     "19", "20", "21",
#     "22", "23", "24",
#     "25", "26", "27",
#     "28", "29", "30",
#     "31"
#   ],
#   "time": [
#     "04:00", "05:00", "06:00",
#     "07:00", "08:00", "09:00",
#     "10:00", "11:00", "12:00",
#     "13:00", "14:00", "15:00",
#     "16:00", "17:00"
#   ],
#   "pressure_level": ["900"],
#   "data_format": "netcdf",
#   "download_format": "zip",
#   "area": [64, -11, -27, 36]
# }
# 
# client = cdsapi.Client()
# client.retrieve(dataset, request).download()
# 
# import cdsapi
# 
# dataset = "reanalysis-era5-pressure-levels"
# request = {
#   "product_type": ["reanalysis"],
#   "variable": [
#     "u_component_of_wind",
#     "v_component_of_wind"
#   ],
#   "year": ["2024"],
#   "month": [
#     "01", "02", "03",
#     "04"
#   ],
#   "day": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12",
#     "13", "14", "15",
#     "16", "17", "18",
#     "19", "20", "21",
#     "22", "23", "24",
#     "25", "26", "27",
#     "28", "29", "30",
#     "31"
#   ],
#   "time": [
#     "06:00", "07:00", "08:00",
#     "09:00", "10:00", "11:00",
#     "12:00", "13:00", "14:00",
#     "15:00", "16:00", "17:00",
#     "18:00"
#   ],
#   "pressure_level": ["900"],
#   "data_format": "netcdf",
#   "download_format": "zip",
#   "area": [64, -11, -27, 36]
# }

