# script for analyzing the data for Safi et al 2025
# Data analysis is done at two scales: 8-second bursts (year-round) and daily (only during migration)
# Elham Nourani, PhD. elham.nourani@unil.ch

#inputs: "thinned_laterality_w_gps_wind_all_filters2_public_prep.rds" (from 03a_data_prep_bursts.r & available from the Edmond repository), "your_path/your_hourly_migration_metrics.rds" (from 03b_data_prep_days.r)
#outputs: Fig 3 all panels. Extended Data Figure 1. Extended Data Table 1

library(tidyverse)
library(corrr)
library(INLA)
library(gridExtra)
library(patchwork)
library(terra)
library(xtable)


#---------------------------------------------------------------------------------
## Step 1: Is laterality more likely when the task is difficult?             #####
#---------------------------------------------------------------------------------

#read in filtered data. this is not filtered for days since tagging
filtered_w_LI <- readRDS("thinned_laterality_w_gps_wind_all_filters2_public_prep.rds")


#### ----------------------- filter for day since tagging and z-transform
quantile(filtered_w_LI$days_since_tagging, probs = 0.9) #284 ... 

circling_data <- filtered_w_LI %>% 
  #filter for days since tagging. only keep the 
  filter(days_since_tagging < 284) %>% 
  #make sure to do all the filtering before scaling the variables!!!!
  mutate_at(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
              "abs_cum_yaw", "wind_speed", "days_since_tagging", "location_lat_closest_gps_raw"),
            list(z = ~scale(.))) %>% 
  as.data.frame()

#### ----------------------- look at multi-collinearity

circling_data %>% 
  dplyr::select(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
                  "abs_cum_yaw", "days_since_tagging", "wind_speed", "location_lat_closest_gps_raw")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #correlated variables: sd of roll and pitch, sd of yaw and cumulative yaw, age and latitude

#### ----------------------- model: binomial logistic regression with inla

#create new data: to make predictions for values that fall on a regular grid for visualization of interaction terms

#to make sure the predictions cover the parameter space, create a dataset with all possible combinations. The max of the variables might be outliers, so use the 99% quantile instead
grd_pitch_yaw <- expand.grid(x = seq(from = quantile(circling_data$mean_pitch_mean, .01, na.rm = T), to = max(circling_data$mean_pitch_mean, na.rm = T),  length.out = 50),
                             y = seq(from = min(circling_data$abs_cum_yaw, na.rm = T), to = quantile(circling_data$abs_cum_yaw, .99, na.rm = T),  length.out = 50)) %>% 
  rename(mean_pitch_mean = x,
         abs_cum_yaw = y)  %>% 
  mutate(wind_speed = attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "pitch_yaw")

grd_wind_yaw <- expand.grid(x = seq(from = min(circling_data$wind_speed,na.rm = T), to = max(circling_data$wind_speed, na.rm = T),  length.out = 50),
                            y = seq(from = min(circling_data$abs_cum_yaw, na.rm = T), to = quantile(circling_data$abs_cum_yaw, .99, na.rm = T),  length.out = 50)) %>% 
  rename(wind_speed = x,
         abs_cum_yaw = y)  %>% 
  mutate(mean_pitch_mean = attr(circling_data[,colnames(circling_data) == "mean_pitch_mean_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "wind_yaw")

grd_wind_pitch <- expand.grid(x = seq(from = min(circling_data$wind_speed,na.rm = T), to = max(circling_data$wind_speed, na.rm = T),  length.out = 50),
                              y = seq(from = quantile(circling_data$mean_pitch_mean, .01, na.rm = T), to = max(circling_data$mean_pitch_mean, na.rm = T),  length.out = 50)) %>% 
  rename(wind_speed = x,
         mean_pitch_mean = y)  %>% 
  mutate(abs_cum_yaw = attr(circling_data[,colnames(circling_data) == "abs_cum_yaw_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "wind_pitch")


#merge all together
grd_all <- bind_rows(grd_pitch_yaw, grd_wind_yaw, grd_wind_pitch) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation for consistency
  mutate(mean_pitch_mean_z = (mean_pitch_mean - attr(circling_data[,colnames(circling_data) == "mean_pitch_mean_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "mean_pitch_mean_z"],'scaled:scale'),
         abs_cum_yaw_z = (abs_cum_yaw - attr(circling_data[,colnames(circling_data) == "abs_cum_yaw_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "abs_cum_yaw_z"],'scaled:scale'),
         wind_speed_z = (wind_speed - attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:scale'),
         days_since_tagging_z = (days_since_tagging - attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:scale'),
         individual_local_identifier = sample(circling_data$individual_local_identifier, nrow(.), replace = T) %>% as.factor(),
         laterality_bi = NA) 


#bin the age variable and append the new data
data <- circling_data %>% 
  dplyr::select(c(intersect(colnames(grd_all), colnames(.)), "life_stage", "location_lat_closest_gps_raw_z")) %>% 
  full_join(grd_all) %>% 
  mutate(individual_local_identifier2 = individual_local_identifier, #repeat individual ID column to be used in the model formula for random effects
         individual_local_identifier3 = individual_local_identifier,
         age_group = inla.group(days_since_tagging_z, n = 100, method = "quantile"), #age will be included as a smooth term
         lat_group = inla.group(location_lat_closest_gps_raw_z, n = 10, method = "quantile"))  


#re-order life stage, so that post-fledging is the reference level
data$life_stage <- factor(data$life_stage, levels = c("post-fledging", "migration", "wintering"))

#### build the model -----------------------
m_inla <- inla(laterality_bi ~ 1 + mean_pitch_mean_z * abs_cum_yaw_z * wind_speed_z + 
                 f(individual_local_identifier, mean_pitch_mean_z, model = "iid") +  #random effect on the slope 
                 f(individual_local_identifier2, abs_cum_yaw_z, model = "iid") + #random effect on the slope
                 f(individual_local_identifier3, wind_speed_z, model = "iid") + #random effect on the slope
                 f(age_group, model = "rw1"), #smooth term for age
               data = data, family = "binomial", #run the model as a binomial logistic regression
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

#### model evaluation -----------------------

#model validation metrics
eval <- data.frame(CPO = mean(m_inla$cpo$cpo, na.rm = T), # 0.51
                   Mlik = as.numeric(m_inla$mlik[,1][2])) # -12951.97

#### coefficients plot (Fig 3 a) -----------------------------------------------------------------------------

# posterior means of coefficients
graph <- as.data.frame(summary(m_inla)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

#graph$Model<-i
graph$Factor <- rownames(graph)

#remove weeks since dispersal
VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

#plot the coefficients
X11(width = 3.2, height = 2)
(coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#0d0887", size = 1.5)  +
    labs(x = "Estimate", y = "") +
    scale_y_discrete(labels = rev(c("Intercept", "Average pitch", "Absolute total yaw",
                                    "Wind speed", "Average pitch: Absolute total yaw", "Average pitch: Wind speed",
                                    "Absolute total yaw: Wind speed", "Average pitch : Absolute total yaw: \nWind speed "))) +
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#0d0887", linewidth = 0.5) +
    ggtitle("a") +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"), #top, right, bottom, left
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          panel.grid.minor = element_line(color = "white"),
          plot.title = element_text(face = "bold"), # make title bold
          axis.title.x = element_text(margin = margin(t = 2))) #increase distance between x-axis values and title
)


#### plot the smooth term (Fig 3 b) -----------------------------------------------------------------------------

# average start and end of migration... to add to the plot
population_migr_period <- filtered_w_LI %>% 
  filter(life_stage == "migration") %>% 
  group_by(individual_local_identifier) %>% 
  arrange(start_timestamp, .by_group = T) %>% 
  summarize(migration_start = min(days_since_tagging),
            migration_end = max(days_since_tagging)) %>% 
  ungroup() %>% 
  summarize(avg_migration_start = mean(migration_start) %>% round(), 
            avg_migration_end = mean(migration_end) %>% round()) %>% 
  rename(xmin = avg_migration_start,
         xmax = avg_migration_end) %>% 
  mutate(ymin = -0.35,
         ymax = 0.48) %>% #format in case I want a shaded rectangle for the migration period.
  as.data.frame()


# Extract the summary of the smooth term 
smooth_effects <- m_inla$summary.random$age_group %>% 
  #back transform the age values
  mutate(age = ID * attr(data[,colnames(data) == "days_since_tagging_z"],'scaled:scale') +
           attr(data[,colnames(data) == "days_since_tagging_z"],'scaled:center') )

# Plot the smooth term....
X11(width = 3.42, height = 2)
(s <- ggplot(smooth_effects, aes(x = age, y = mean)) +
    geom_line(color = "#0d0887", linewidth = 0.4) +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "#0d0887", alpha = 0.12) +
    geom_hline(yintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_vline(xintercept = c(population_migr_period$xmin, population_migr_period$xmax),
               linetype = "dotted",  color = "black", linewidth = 0.5, show.legend = TRUE) +
    labs(y = "Effect Size", x = "Days since tagging") +
    ggtitle("b") +
    ylim(-0.35, 0.48) +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"), #top, right, bottom, left
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          plot.title = element_text(face = "bold"), # make title bold
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 2)))
)

#combine the two plots to make the multi-panel plot for the coefficients and the smooth term
X11(width = 6.7, height = 2)
model_output_p <- grid.arrange(coefs, s, nrow = 1, widths = c(0.6, 0.4))

#### individual-specific coefficients plot (Extended Data Fig 1) -----------------------------------------------------------------------------

#extract laterality for each individual. to use for coloring. 
handedness <- circling_data %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(individual_local_identifier, laterality_dir_ind)

#extract random effects, ID is for the unique individuals
#mean_pitch_mean_z
random_effects_pitch <- m_inla$summary.random$individual_local_identifier

#abs_cum_yaw_z
random_effects_yaw <- m_inla$summary.random$individual_local_identifier2

#wind_speed_z
random_effects_wind <- m_inla$summary.random$individual_local_identifier3

pitch <- random_effects_pitch %>% 
  mutate(coef = mean + graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Upper),
         variable = "mean_pitch_mean_z")

yaw <- random_effects_yaw %>% 
  mutate(coef = mean + graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Upper),
         variable = "abs_cum_yaw_z")

wind <- random_effects_wind %>% 
  mutate(coef = mean + graph %>% filter(Factor == "wind_speed_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "wind_speed_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "wind_speed_z") %>% pull(Upper),
         variable = "wind_speed_z")

three_vars <- bind_rows(pitch, yaw, wind) %>%
  full_join(handedness, by = c("ID" = "individual_local_identifier")) %>% 
  mutate(v_line = case_when(
    variable == "mean_pitch_mean_z" ~ graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate),
    variable == "abs_cum_yaw_z" ~ graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate),
    variable == "wind_speed_z" ~ graph %>% filter(Factor == "wind_speed_z") %>% pull(Estimate)
  ))


#re-order the individuals based on the laterality index
three_vars$laterality_dir_ind <- factor(three_vars$laterality_dir_ind, levels = c("left_handed", "ambidextrous", "right_handed"))
three_vars$ID <- reorder(three_vars$ID, desc(three_vars$laterality_dir_ind))


X11(width = 7, height = 8)
(coefs_inds <- ggplot(three_vars, aes(x = coef, y = ID, color = laterality_dir_ind)) +
    geom_vline(data = filter(three_vars, variable == "mean_pitch_mean_z"), 
               aes(xintercept = graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", linewidth = 0.5) + 
    geom_vline(data = filter(three_vars, variable == "abs_cum_yaw_z"), 
               aes(xintercept = graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", linewidth = 0.5) +  
    geom_vline(data = filter(three_vars, variable == "wind_speed_z"), 
               aes(xintercept = graph %>% filter(Factor == "wind_speed_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", linewidth = 0.5) +  
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = lower, xmax = upper), size = 0.8, position = position_dodge(width = .7)) +
    scale_color_manual(values = c("left_handed" = "#9c179e" , "ambidextrous" = "#fb9f3a", "right_handed" =  "#0d0887"),
                       labels = c("Left", "No bias", "Right"),
                       name = "Lateralisation") +
    scale_y_discrete(labels = levels(three_vars$ID)) +
    labs(x = "Estimate", y = "Individual ID") +
    theme_classic() +
    theme(text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 2)),
          strip.text = element_text(size = 10)) + # Increase panel name size
    facet_wrap(~ variable, scales = "free_x", labeller = as_labeller(c(
      "mean_pitch_mean_z" = "Average pitch",
      "abs_cum_yaw_z" = "Absolute total yaw",
      "wind_speed_z" = "Wind speed"
    )) # Separate panels for each variable
    ))

#### plot the interaction terms (Fig 3 c-e) -----------------------------------------------------------------------------

## yaw:pitch----------------------------------------------------------------
#extract information for rows that had NAs as response variables
na_rows <- which(data$interaction == "pitch_yaw")

preds <- data.frame(yaw = data[na_rows,"abs_cum_yaw"],
                    pitch = data[na_rows,"mean_pitch_mean"],
                    preds = m_inla$summary.fitted.values[na_rows,"mean"]) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot

pred_py <- preds %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = preds)) +
  scale_fill_gradientn(#colors = c("#440154FF", "#39568CFF" ,"#20A387FF", "#83e22b", "#FDE725FF"),
    #colors = c("#0d0887", "#7e03a8" ,"#cc4778", "#f89540", "#f0f921"),
    colors = c("#0d0887", "#7e03a8", "white" ,"#d5546e", "#fdae32", "#f0f921"),
    values = c(0, 0.3, 0.5, 0.6, 0.7, 1),
    limits = c(0, 1),
    na.value = "white",
    name = "Probability\n of laterality") +
  guides(fill = guide_colourbar(title.vjust = 2.5)) + # the legend title needs to move up a bit
  labs(x = "Absolute total yaw", y = "Average pitch") +
  ggtitle("c") +
  theme_classic() +
  theme(plot.margin = margin(0, 6, 0, 0, "pt"),
        text = element_text(size = 8),
        legend.direction="vertical",
        legend.position = "right",
        legend.key.width=unit(.2,"cm"),
        legend.key.height=unit(.4,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(face = "bold"), # make title bold
        panel.grid.minor = element_line(color = "white"),
        axis.title.x = element_text(margin = margin(t = 2))) + # increase distance between x-axis values and title
  scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
  scale_y_continuous(expand = c(0, 0))


## wind:yaw----------------------------------------------------------------
#extract information for rows that had NAs as response variables
na_rows <- which(data$interaction == "wind_yaw")

preds <- data.frame(yaw = data[na_rows,"abs_cum_yaw"],
                    wind = data[na_rows,"wind_speed"],
                    preds = m_inla$summary.fitted.values[na_rows,"mean"]) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot

pred_wy <- preds %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = preds)) +
  scale_fill_gradientn(colors = c("#0d0887", "#7e03a8", "white" ,"#d5546e", "#fdae32", "#f0f921"),
                       values = c(0, 0.3, 0.5, 0.6, 0.7, 1),
                       limits = c(0, 1),
                       na.value = "white",
                       name = "Probability\n of laterality") +
  guides(fill = guide_colourbar(title.vjust = 2.5)) + # the legend title needs to move up a bit
  labs(x = "Absolute total yaw", y = "Wind speed") +
  ggtitle("d") +
  theme_classic() +
  theme(plot.margin = margin(0, 6, 0, 0, "pt"),
        text = element_text(size = 8),
        legend.direction="vertical",
        legend.position = "right",
        legend.key.width=unit(.2,"cm"),
        legend.key.height=unit(.4,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(face = "bold"), # make title bold
        panel.grid.minor = element_line(color = "white"),
        axis.title.x = element_text(margin = margin(t = 2))) + # increase distance between x-axis values and title
  scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
  scale_y_continuous(expand = c(0, 0))



##wind:pitch----------------------------------------------------------------

#extract information for rows that had NAs as response variables
na_rows <- which(data$interaction == "wind_pitch")

preds <- data.frame(pitch = data[na_rows,"mean_pitch_mean"],
                    wind = data[na_rows,"wind_speed"],
                    preds = m_inla$summary.fitted.values[na_rows,"mean"]) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot

pred_wp <- preds %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = preds)) +
  scale_fill_gradientn(colors = c("#0d0887", "#7e03a8", "white" ,"#d5546e", "#fdae32", "#f0f921"),
                       values = c(0, 0.3, 0.5, 0.6, 0.7, 1),
                       limits = c(0, 1),
                       na.value = "white",
                       name = "Probability\n of laterality") +
  guides(fill = guide_colourbar(title.vjust = 2.5)) + # the legend title needs to move up a bit
  labs(x = "Average pitch", y = "Wind speed") +
  ggtitle("e") +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        text = element_text(size = 8),
        legend.direction="vertical",
        legend.position = "right",
        legend.key.width=unit(.2,"cm"),
        legend.key.height=unit(.4,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(face = "bold"), # make title bold
        panel.grid.minor = element_line(color = "white"),
        axis.title.x = element_text(margin = margin(t = 2))) + # increase distance between x-axis values and title
  scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
  scale_y_continuous(expand = c(0, 0))


#combine the three plots for the coefficients and the interaction term-------------- horizontal
X11(width = 6.7, height = 2)
combined <- pred_py + pred_wy + pred_wp & theme(legend.position = "right")
(p <- combined + plot_layout(guides = "collect", nrow = 1))

#combine panels into a multi-panel figure in gimp:
#https://theplosblog.plos.org/2019/12/multi-panel-figures-using-gimp-to-combine-individual-images-for-use-in-plos-articles/

#---------------------------------------------------------------------------------
## Step 2: Does laterality impact migration performance?                     #####
#---------------------------------------------------------------------------------

#### data prep -----------------------------------------------------------------------------

#Calculate the daily summaries for pitch, yaw, and wind. Then append GPS- and ACC- derived variables (from 03b_data_prep_days.r)

##open migration data from 03b_data_prep_days.r
migr_hrly <- readRDS("your_path/your_hourly_migration_metrics.rds")

#extract days of migration
migr_days <- migr_hrly %>% 
  distinct(ind_day)

#calculate hourly and daily summaries for wind using the imu, then append migration dataframe (cricling_data was created in Step 1)
migr_hrly_w <- circling_data %>% 
  mutate(ind_day =  paste0(individual_local_identifier, "_", as.character(unique_date))) %>% 
  filter(ind_day %in% migr_days$ind_day) %>%  #subset for days in the migration data 
  mutate(dt_1hr = round_date(start_timestamp, "1 hour")) %>%  #assign unique hour
  group_by(individual_local_identifier, unique_date, dt_1hr) %>%  
  mutate(hrly_mean_wind_speed = mean(wind_speed, na.rm = T), #summarize wind for every hour. I have a unique wind value for each hour anyway
         hrly_mean_cum_yaw = mean(abs_cum_yaw, na.rm = T),
         hrly_max_cum_yaw = max(abs_cum_yaw, na.rm = T),
         hrly_max_mean_pitch = max(mean_pitch_mean, na.rm = T),
         hrly_mean_mean_pitch = mean(mean_pitch_mean, na.rm = T),) %>% 
  ungroup() %>% 
  group_by(individual_local_identifier, unique_date) %>% #summarize wind for every hour
  mutate(daily_max_wind = max(wind_speed, na.rm = T),
         daily_mean_wind = mean(wind_speed, na.rm = T),
         daily_max_cum_yaw = max(abs_cum_yaw, na.rm = T),
         daily_mean_cum_yaw = mean(abs_cum_yaw, na.rm = T),
         daily_max_mean_pitch = max(mean_pitch_mean, na.rm = T),
         daily_mean_mean_pitch = mean(mean_pitch_mean, na.rm = T)) %>% 
  ungroup() %>% 
  select("individual_local_identifier", "days_since_tagging", "unique_date", "life_stage", "laterality_bank_day", "laterality_dir_day", "laterality_bank_stage", "laterality_dir_stage", 
         "laterality_bank_ind","laterality_dir_ind", "location_lat_closest_gps_raw_z", "ind_day", "dt_1hr", "hrly_mean_wind_speed", "hrly_mean_cum_yaw", "hrly_max_cum_yaw", "hrly_max_mean_pitch", "hrly_mean_mean_pitch" ,
         "daily_max_wind", "daily_mean_wind", "daily_max_cum_yaw", "daily_mean_cum_yaw", "daily_max_mean_pitch","daily_mean_mean_pitch") %>% 
  group_by(individual_local_identifier, unique_date, dt_1hr) %>% 
  slice(1) %>% #just keep one row per hour. 
  #bind to the daily migration data. get rid of days in the migr_hrly data that don't have laterality data
  left_join(migr_hrly) %>% 
  as.data.frame()


migr_daily_w <- migr_hrly_w %>% #ww: wind
  select(-contains("hrly")) %>%  #remove hourly data
  group_by(individual_local_identifier, unique_date) %>%  #keep one row per day. the daily data is repeated for every row of that day in the hourly dataset
  slice(1) %>% 
  ungroup() %>% 
  mutate(laterality_dir_day = as.character(laterality_dir_day), #convert to character, so that"ambidextrous" is the reference level
         laterality_bi_day = ifelse(laterality_dir_day == "ambidextrous", 0, 1)) %>% 
  as.data.frame()

### Model migration performance as a function of laterality -----------------------------------------------------

#z-transform the variables
data_m <- migr_daily_w %>%
  mutate(daily_abs_LI = abs(laterality_bank_day),
         laterality_bi_ind = ifelse(laterality_dir_ind == "ambidextrous", "ambidextrous", "handed"),
         laterality_bi_stage = ifelse(laterality_dir_stage == "ambidextrous", "ambidextrous", "handed")) %>% #because this is migration data only, there should be one assignment per individual
  #make sure to do all the filtering before scaling the variables!!!!
  mutate(across(contains("daily") , ~scale(.), .names = "{.col}_z"),
         days_since_tagging_z = scale(days_since_tagging)) %>%
  as.data.frame()

#check for autocorrelation
data_m %>% 
  dplyr::select(ends_with("_z")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #correlated: mean and max wind, max cum yaw and mean cum yaw, daily distance an average speed, mean vedba and IQR vedba; max and avg altitude;
# sd of vertical speed with avg altitude (0.57) and max altitude (0.64)

#model separately for daily_mean_vedba, daily_max_altitude, daily_distance, daily_mean_cum_yaw, daily_mean_mean_pitch, daily_max_altitude.
#create a character vector to use in the model formula
response_vars <- c( "daily_distance", "daily_max_altitude",
                    "daily_mean_vedba", "daily_mean_cum_yaw", 
                    "daily_mean_mean_pitch")

#prepare character vector to be used for naming the plot panels
response_names <- c(
  expression(atop(bold("f") * phantom("                             "), "Daily distance (km)")),
  expression(atop(bold("g") * phantom("                                  "), "Max flight altitude (m)")),
  expression(atop(bold("h") * phantom("                        "), "Avg VeDBA (g)")), 
  expression(atop(bold("i") * phantom("                               "), "Avg total yaw (deg)")), 
  expression(atop(bold("j") * phantom("                        "), "Avg pitch (deg)"))
)

#for each response variable, build a linear model, extract coefficients in latex format, plot the coefficient plot
plots_ls <- lapply(1:length(response_vars), function(response){
  
  #model formula
  formula <- paste0(response_vars[response], " ~ 1 + daily_max_wind_z + daily_abs_LI_z + f(individual_local_identifier, model = 'iid')") %>% formula()
  
  #run the model
  m_inla <- inla(formula,
                 data = data_m,
                 control.compute = list(cpo = TRUE),
                 control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted
  
  # extract coefficients 
  graph <- as.data.frame(summary(m_inla)$fixed)
  colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
  
  graph$Factor <- rownames(graph)
  
  VarOrder <- rev(unique(graph$Factor))
  VarNames <- VarOrder
  
  graph$Factor <- factor(graph$Factor, levels = VarOrder)
  levels(graph$Factor) <- VarNames
  
  #export the output as a latex table ----------- uncomment this section to write the tables to file
  # Convert to LaTeX table
  # latex_table <- xtable(graph)
  # 
  # # Specify the file path for each response variable
  #  file_path <- paste0("your_path/latex_table_", response_vars[response], "2.txt")
  # 
  # # Open a connection to the file
  # sink(file_path)
  # 
  # # Print the LaTeX code to the file
  # print(latex_table, type = "latex", include.rownames = FALSE)
  # 
  # # Close the connection to the file
  # sink()
  
  #plot the coefficients -----------
  
  #remove intercept for better visualization. The intercepts will be reported in the table
  graph <- graph[-1,]
  
  ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#0d0887", size = 1.5)  +
    labs(x = if (response == 5) "Estimate" else "", 
         y = "") +
    scale_y_discrete(labels = if (response %in% c(1,4)) rev(c("Intercept", "Max. wind speed", "Abs. Laterality Index")) else c("", "")) + 
    #scale_y_discrete(labels = if (response %in% c(1,4)) rev(c("Intercept", "Max. wind speed", "Abs. Laterality Index", "Days since tagging")) else c("", "", "")) + 
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#0d0887", linewidth = 0.5) +
    ggtitle(response_names[response]) +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"),
          text = element_text(size = 8),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 2))) #increase distance between x-axis values and title
})

#put all plots together
X11(width = 6.7, height = 2.5)
combined <- reduce(plots_ls[1:5], `+`)
(p <- combined + plot_layout(ncol = 3))


