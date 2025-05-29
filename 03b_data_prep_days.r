#Script for preparing daily metrics for Safi et al 2025.
#This script calculating migration metrics for analysis in 04_data_analysis.r
#Elham Nouarni, PhD.  elham.nourani@unil.ch


#inputs: "your_path/your_downloaded_gps_data.rds" (from 02_wind_annotation.r), "updated_life_cycle_nov24.rds" (from Edmond repository), "EGM96_us_nga_egm96_15.tif" (from Edmond repository), "your_path/your_processed_acc_data.rds" (from 01_imu_processing.r)
#output: "your_path/your_hourly_migration_metrics.rds"

library(move2)
library(sf)
library(terra)

#---------------------------------------------------------------------------------
## Step 1: Add life stage                                                    #####
#---------------------------------------------------------------------------------

#GPS data was downloaded in 02_wind_annotation.r
gps <- readRDS("your_path/your_downloaded_gps_data.rds")

#life-cycle stages from Edmond repository
life_cycle <- readRDS("updated_life_cycle_nov24.rds")

#open gps data. make a unique ID for each ind-day combo. for speed, etc. calculations.
gps_no_winter <- gps %>%  #2022-09-25 07:56:41.0000 to 2024-04-15 10:02:33.0000
  drop_na(individual_local_identifier, location_lat) %>% #remove NA individuals and NA locations.
  mutate(unique_date = as.Date(timestamp),
         yr = year(timestamp),
         mn = month(timestamp),
         dy = day(timestamp),
         hr = hour(timestamp),
         unique_hr = paste(yr,mn,dy,hr, sep = "_"),
         closest_hr = round(timestamp, units = "hours") %>% as.character()) %>% 
  mutate(ind_day = paste0(individual_local_identifier, "_", unique_date)) %>% 
  full_join(life_cycle %>% select(individual_local_identifier, migration_start, migration_end, first_exploration), by = "individual_local_identifier") %>% 
  group_by(individual_local_identifier) %>% 
  rowwise() %>% 
  mutate(life_stage = case_when(
    between(unique_date, migration_start, migration_end) ~ "migration",
    between(unique_date, first_exploration, migration_start) ~ "post-fledging",
    unique_date < first_exploration ~ "pre-fledging",
    unique_date > migration_end ~ "wintering",
    TRUE ~ NA_character_
  )) %>% 
  ungroup() %>% 
  filter(!(life_stage %in% c("wintering", "pre-fledging"))) %>% 
  as.data.frame()

#---------------------------------------------------------------------------------
## Step 2: Calculate daily distance, and flight altitude                    #####
#---------------------------------------------------------------------------------

#open geoid layer for calculating flight altitude. This is included in the Edmond repository
geo <- rast("EGM96_us_nga_egm96_15.tif")

#filter the data for 1 hour. find days that have 12-14 hours of data. Calculate daily distance as the sum of hourly distances for the day.
gps_1hr <- gps_no_winter %>% 
  #remove the points at (0,0) ... there are 54 of them!!
  filter(!(location_lat == 0 & location_long == 0)) %>% 
  #only keep migration days
  filter(life_stage == "migration") %>%  
  st_as_sf(coords = c("location_long", "location_lat"), crs = "EPSG:4326") %>% 
  mutate(dt_1hr = round_date(timestamp, "1 hour")) %>%  #the hourly subset will be used just to pick days that have 12-14 hours of data. 
  group_by(individual_local_identifier, unique_date, dt_1hr) %>% 
  slice(1) %>% 
  ungroup() %>% 
  #calculate flight height
  extract(x = geo, y = ., method = "simple", bind = T) %>% 
  st_as_sf() %>% 
  mutate(height_msl = height_above_ellipsoid %>% as.numeric() - geoid_undulation) %>% 
  group_by(individual_local_identifier, unique_date) %>% 
  #summarize(n_hrs_in_day = n()) #most days have 12-14 hours of data
  filter(n() >= 12, .preserve = T) %>%  #keep days with >= 12 hours of data
  arrange(timestamp, .by_group = T) %>% 
  mutate(hrly_step_length = if_else(row_number() == 1, NA, st_distance(geometry, lag(geometry), by_element = TRUE) %>% as.numeric() / 1000), #kilometers
         daily_distance = sum(hrly_step_length, na.rm = T),
         daily_avg_speed = mean(hrly_step_length, na.rm = T),
         daily_avg_altitude = mean(height_msl, na.rm = T),  #make sure to remove negative values when modeling
         daily_max_altitude = max(height_msl, na.rm = T)) %>% 
  ungroup() 

#---------------------------------------------------------------------------------
## Step 3: Calculate daily VeDBA                                             #####
#---------------------------------------------------------------------------------

#open acc data prepared in 01b_imu_processing.r which includes the mean odba and vedba for each row (1.2 seconds) of data

acc_g <- readRDS("your_path/your_processed_acc_data.rds") %>% 
  mutate(unique_date = as.Date(timestamp),
         ind_day = paste0(individual_local_identifier, "_", as.character(unique_date)))

#extract days of migration based on gps data processed in the previous step
migr_days <-gps_1hr %>% 
  distinct(ind_day)

#filter for migration days, calculate vedba, summarize for each day
acc_migr <- acc_g %>% 
  filter(ind_day %in% migr_days$ind_day) %>%  #only keep acc data for migration days 
  mutate(dt_1hr = round_date(timestamp, "1 hour")) %>%  #the hourly subset will be used just to pick days that have 12-14 hours of data. 
  group_by(individual_local_identifier, unique_date, dt_1hr) %>% 
  mutate(hrly_mean_vedba = mean(VeDBA, na.rm = T),
         hrly_max_vedba = max(VeDBA, na.rm = T),
         hrly_min_vedba = min(VeDBA, na.rm = T),
         hrly_IQR_vedba = quantile(VeDBA, prob = 0.75) - quantile(VeDBA, prob = 0.25)) %>%  #calculate interquartile range... this includes non-flight data tooo.... so maybe not so useful
  ungroup() %>% 
  group_by(individual_local_identifier, unique_date) %>% 
  mutate(daily_mean_vedba = mean(VeDBA, na.rm = T),
         daily_max_vedba = max(VeDBA, na.rm = T),
         daily_min_vedba = min(VeDBA, na.rm = T),
         daily_IQR_vedba = quantile(VeDBA, prob = 0.75) - quantile(VeDBA, prob = 0.25)) %>%  #calculate interquartile range... this includes non-flight data tooo.... so maybe not so useful
  ungroup() %>% 
  group_by(individual_local_identifier, unique_date, dt_1hr) %>% 
  slice(1) %>%#just keep one row per hour 
  ungroup()

#append to the gps data summarized for migration period
migr_hourly <- gps_1hr %>% #this is an sf object. keep the lat and long
  mutate(location_long = st_coordinates(.)[,1],
         location_lat = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  select(individual_id, deployment_id, tag_id, study_id, individual_local_identifier, tag_local_identifier, individual_taxon_canonical_name, location_long, location_lat,
         ind_day, migration_start, migration_end, first_exploration, life_stage, dt_1hr,  unique_date, height_msl, hrly_step_length, daily_distance, daily_avg_speed, daily_avg_altitude,
         daily_max_altitude) %>% 
  full_join(acc_migr %>% dplyr::select(individual_local_identifier, ind_day, dt_1hr, unique_date, hrly_mean_vedba, hrly_max_vedba, hrly_min_vedba, hrly_IQR_vedba,
                                       daily_mean_vedba, daily_max_vedba, daily_min_vedba, daily_IQR_vedba)) %>% 
  as.data.frame()

#write to file. This file will be used in 04_data_analysis.r Step 2 to analyse how lateralization impacts daily migration performance
saveRDS(migr_hourly, file = "your_path/your_hourly_migration_metrics.rds")

