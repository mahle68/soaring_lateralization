# soaring_lateralization
Scripts for reproducing the results of Safi et al 2025


The repository includes R scripts for reproducing the results of Safi et al 2025 **"Ontogenetic Loss of Lateralisation Improves Migratory Flight Performance"**

The following scripts are included:

**00_imu_functions.r** includes the functions for coverting quaternion data to pitch, yaw, and roll angles. Plus additional functions specifically for processing data from biologging devices manufactured by e-obs GmbH

**01_imu_processing.r** includes the workflow for processing accelerometry and quaternion data

**02_wind_annotation.r** includes the workflow for downloading gps data, downloading the corresponding wind data, and annotating the gps data with wind speed

**03a_data_prep_bursts.r** includes the workflow for calculating the laterality index

**03b_data_prep_days.r** includes the workflow for calculating daily migration metrics using the GPS (speed, flight altitude) and accelerometry (VeDBA)

**04_data_analysis.r** includes the workflow for analyzing the data at two scales: the 8-second burst (is lateralization more likely in difficult flight conditions) and the daily scale (does lateralization impact migration performance)

**05_plot_Fig2.r** includes the code to generate the two panels of Figure 2

**000_bonus_gps_processing.r** includes the workflow for downloading GPS data and using the high-resolution bursts to classify the flight types. This is not a necessary step for reproducing the results of the study, but can be useful for data exploration. 

The following datasets are made available in the Edmond repository xyz

**"meta_data.rds"** includes the meta-data for each individual, including the timing of start of different life stages

**"imu_wind_laterality_bursts.rds"** includes the dataset ready for modeling following Step 1 of 04_data_analysis.r 

**"EGM96_us_nga_egm96_15.tif"** geoid layer used for calculating flight altitude in Step 2 of 03b_data_prep_days.r 

**"migration_perfromance_days.rds"** contains daily migration metrics. ready for modeling following Step 2 of 04_data_analysis.r 

The raw tracking data is stored on Movebank.or under DOI: xyz
