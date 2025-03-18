# ----------------
# Run calibrarion
# ----------------

rm(list = ls())

library(broom)
library(ggplot2)
library(gstat)
library(interp)
library(lubridate)
library(raster)
library(usmap)
library(tidyverse)
library(usmap)

output_prefix <- "residential"   # "residential" or "commercial"
load_type <- "heating" # "cooling" or "heating"


if ( load_type == "cooling" ) {
    dem_column <- "cooling_demand"
    i_months_ssn <- c(6, 7, 8)
} else { 
    dem_column <- "heat_demand"
    i_months_ssn <- c(1, 2, 12)
}   

# paths must be adjusted
data_dir <- "._heating implementation_/_data processing_/_res_/t_th_22C_adv_hum"
# data_ts_fl <- "ts_residential_nrel_heating_t_th_22_2018_ref_adv_humid_eff.csv"
# data_ts_fl <- "ts_residential_nrel_cooling_t_th_22_2018_ref_adv_humid_eff.csv"
data_ts_fl <- paste0(
    "ts_",
    output_prefix,
    "_nrel_",
    load_type,
    # "_t_th_22_2018_ref_adv_humid_eff.csv"
    "_t_th_22_2018_ref_adv_humid_eff.csv"
    )

# load data -------------------------------------------------------------------

load_df <- read.csv(file.path(data_dir, data_ts_fl))
load_df["time_date"] = ymd(unlist(load_df["time"]))

# build regression of the normal dd-demand on the measured one ----------------

# there are states with zero cooling load (if CDDs are used)
if ( load_type == "cooling" ) {
    aggre_mean_load1 <- load_df[, c("state", "cooling_demand", "load_daily")] %>%
        group_by(state) %>%
        summarize_all(mean)
    non_zero_load_states <- unlist(
        aggre_mean_load1[aggre_mean_load1$cooling_demand > 0, "state"]
    )
    load_df2 <- load_df %>%
        dplyr::filter(state %in% non_zero_load_states)   
} else {
    load_df2 <- load_df
}

load_df_normal <- load_df2 %>%
    group_by(state) %>%
    mutate(
        dem_mes_normal = load_daily / max(load_daily),
        dem_dd_normal = heat_demand / max(heat_demand)
        # dem_dd_normal = cooling_demand / max(cooling_demand)
    )

test_normal_clean <- load_df_normal[
    , c("state", "dem_mes_normal", "dem_dd_normal")
    ] %>%
    dplyr::filter(

        # # to fix the treshold issue of the cooling demand
        # dem_dd_normal > 0.1 &
        # dem_dd_normal < 0.9

        dem_mes_normal > 0.1 &
        dem_mes_normal < 0.9
        # dem_mes_normal > 0.05 &
        # dem_mes_normal < 0.95        
    ) %>%
    nest_by() %>% 
    mutate(model1 = list(lm(dem_mes_normal ~ dem_dd_normal, data = data)),
           across(starts_with("model"),  ~ list(Predict = predict(., data)),
            .names = "{.col}_Predict")) %>% 
    select(-model1)  %>%
    ungroup %>% 
    unnest(c(data, model1_Predict))     

test_pl2 <- test_normal_clean %>%
    ggplot(
        aes(
            x = dem_mes_normal,
            y = model1_Predict,
            color = state,
            group = state
        )
    )+
    geom_line()+
    coord_cartesian(xlim = c(0, 1.2), ylim = c(0, 1.2))+
    facet_wrap(vars(state))
  ggsave(
    paste0(
        "regr_explor_",
        output_prefix,
        "_nrel_",
        load_type,
        "_map_t_th22C_cut25perc.pdf"
        )
    )    

res_normal_clean <- test_normal_clean %>%
    mutate(
        mae = abs(model1_Predict - dem_mes_normal),
        mape = abs(model1_Predict - dem_mes_normal) / dem_mes_normal
    )

res_error_normal_clean <- res_normal_clean %>%
    group_by(state) %>%
    summarise_all(median) %>%
    mutate(
        aggr_dem_err = (dem_mes_normal - model1_Predict) / dem_mes_normal
    )

# plot the errors -------------------------------------------------------------    
state <- map_data("state")
names_map <- state.abb
names(names_map) <- tolower(state.name)
state["state"] <- names_map[state$region]

state_data <- state %>% 
    left_join(res_error_normal_clean) 

title_txt <- paste0(
    "Demand error: ", load_type, ", ", output_prefix
)
ggplot(data=state_data, aes(x=long, y=lat, fill=mape, group=group)) + 
  geom_polygon(color = "gray10") + 

  # heating_v5, resid (filtering by measurments)
  scale_fill_distiller(palette="Oranges", direction=0, limits=c(0.07, 0.31)) +  

  # # heating_v5, comm (filtering by measurments)
  # scale_fill_distiller(palette="Oranges", direction=0, limits=c(0.1, 0.27)) +

  # # heating_v5, comm (filtering by measurments)
  # scale_fill_distiller(palette="Oranges", direction=0, limits=c(0.1, 0.27)) +  

  # # cooling (filtering by cdd-deriven)
  # scale_fill_distiller(palette="Oranges", direction=0, limits=c(0, 0.25)) +  

  # # heating (filtering by measurments)
  # scale_fill_distiller(palette="Oranges", direction=0, limits=c(0, 0.25)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  ggtitle(title_txt) + 
  coord_fixed(1.3)
  ggsave(
    paste0(
        # "sensit_",
        "clean_lm_mape_",
        output_prefix,
        "_nrel_",
        load_type,
        "_map_t_th22C_cut15perc.pdf"
        )
    )  

ggplot(data=state_data, aes(x=long, y=lat, fill=aggr_dem_err, group=group)) + 
  geom_polygon(color = "gray10") + 

  # heating_v5, resid (filtering by measurments)
  scale_fill_distiller(palette="PiYG", direction=0, limits=c(-0.27, 0.04)) +
  # # heating_v5, comm (filtering by measurments)
  # scale_fill_distiller(palette="PiYG", direction=0, limits=c(-0.26, 0.02)) +
  # # cooling (filtering by cdd-deriven demand)
  # scale_fill_distiller(palette="PiYG", direction=0, limits=c(-0.6, 0.06)) +
  # # heating (filtering by measurments)
  # scale_fill_distiller(palette="PiYG", direction=0, limits=c(-0.3, 0.10)) +

  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  ggtitle(title_txt) + 
  coord_fixed(1.3)
  ggsave(
    paste0(
        # "sensit_",
        "clean_lm_aggr_err_",
        output_prefix,
        "_nrel_",
        load_type,
        "_map_t_th22C.pdf"
        )
    )  
