# -------------------------------------------------
# Preprocess ReStock/ComStock data by states
# (attach coordinates of a centroid for each state)
# -------------------------------------------------


rm(list=ls())


load_type <- "cooling" # "cooling" or "heating"
buildinds_type <- "residential" # "commercial" or "residential"

# paths must be adjusted
data_dir <- "./_heating implementation_/_data processing_/_res_/timeseries"
data_fl <- paste0(buildinds_type, "_", load_type ,"_by_states.csv")

inputs_dir <- "./_heating implementation_/_workflow_/data_inputs"
us_centroids_fl <- "us_states_geo.csv"
us_pop_fl <- "./_heating implementation_/suppl_data/gpw_v4_population_density_rev11_2000_2.5m.tif"

res_dir <- "./_heating implementation_/_data processing_/_res_"

library(ggplot2)
library(gstat)
library(interp)
library(lubridate)
library(raster)
library(usmap)
library(tidyverse)
library(usmap)


# load data -------------------------------------------------------------------

load_df <- read.csv(file.path(data_dir, data_fl))
load_df["time"] = ymd_hms(unlist(load_df["timestamp"]))

# https://inkplant.com/code/state-latitudes-longitudes
states_centers <- read.csv(
    file.path(inputs_dir, us_centroids_fl),
    sep=",",
    header=TRUE,    
    stringsAsFactors=FALSE
)

# Hawai and Alaska are out of scope for the analysis
# states_to_consider <- us_map()$abbr[!(us_map()$abbr %in% c("AK", "HI"))]
us_map_df <- data.frame(
    State=us_map()$full[!(us_map()$abbr %in% c("AK", "HI"))],
    Code=us_map()$abbr[!(us_map()$abbr %in% c("AK", "HI"))]
)

geo_df <- merge(states_centers, us_map_df)

# prepare heat demand data ----------------------------------------------------
heat_wide_df <- pivot_wider(
    load_df[c("in.state", "overall_resload", "time")],
    names_from="in.state",
    values_from="overall_resload"
    )

heat_loc_df <- geo_df
heat_loc_df$y = heat_loc_df[["Latitude"]]
heat_loc_df$x = heat_loc_df[["Longitude"]]
coordinates(heat_loc_df) = ~x+y

# x            <- st_sf(a = 1:2, b = 2:3, geom=st_sfc(st_point(1:2),st_point(3:4)))
# ts_data      <- list()

ts_data <- list()
ts_data[[1]] <- heat_wide_df
heat_loc_df$ts_data <- ts_data

# # TODO fix CRS --------------------------------------------------------------
# proj4string(grid_east_germany) <- proj4string(dwd_east_sp)
# head(grid_east_germany)

# -----------------------------------------------------------------------------
res_wide <- load_df %>% 
    pivot_wider(
        names_from = `in.state`,
        values_from = overall_resload
    )
res_wide <- res_wide %>%
    mutate(
        day=day(time),
        month=month(time),
        year=year(time),
    ) %>%
    mutate(
        time_date=ymd(paste(year, month, day, sep="-"))
    )    

res_daily <- res_wide[, -1] %>%
    group_by(
        time_date
    ) %>%
    summarise(across(AL:WY, sum))

# geo_df <- geo_df
# geo_df$y = geo_df[["Latitude"]]
# geo_df$x = geo_df[["Longitude"]]
# coordinates(geo_df) = ~x+y    
# res_daily.loc <- geo_df

res_daily_long <- res_daily %>% 
    pivot_longer(
        cols=AL:WY,
        names_to="state",
        values_to="load_daily"
    )

geo_df <- merge(states_centers, us_map_df) %>%
    rename(
        state_name=State,
        state=Code,
        x=Longitude,
        y=Latitude
    )

res_daily_long_sp <- merge(res_daily_long, geo_df, all=TRUE)

res_daily_long_sp$time_date <- ymd(res_daily_long_sp$time_date)
res_daily_long_sp_clean <- res_daily_long_sp %>%
    group_by(state) %>%
    arrange(time_date)

write.csv(
    res_daily_long_sp_clean,
    file.path(
        res_dir,
        paste0(buildinds_type, "_nrel_", load_type, "_2018_ref.csv")
    ),
    row.names = FALSE)

# NAs must be dropped
res_daily_long_sp2 <- res_daily_long_sp[-((nrow(res_daily_long_sp) - 1):nrow(res_daily_long_sp)),]

# generate grid ---------------------------------------------------------------
xh <- 0.5
x_margin <- 5
grid_us_fine <- expand.grid(
  x = seq(
    from = round(min(res_daily_long_sp2$x)) - x_margin,
    to = round(max(res_daily_long_sp2$x)) + x_margin,
    by = xh
  ),
  y = seq(
    from = round(min(res_daily_long_sp2$y)) - x_margin,
    to = round(max(res_daily_long_sp2$y)) + x_margin,
    by = xh
  )
)
coordinates(grid_us_fine) <- ~ x + y
class(grid_us_fine)

gridded(grid_us_fine) <- TRUE
class(grid_us_fine)

# graphical check
plot(grid_us_fine,
  main = paste("Centroids by States\n and Interpolation Grid"),
  col = "grey",
  cex.main = 0.9
)
# plot(east_germany_states_sp, add = TRUE, border = "red")
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "blue")

# read population data --------------------------------------------------------
pop_raster <- raster(us_pop_fl)

mainland_raster <- raster::crop(pop_raster, extent(grid_us_fine))
pop_raster <- resample(
    mainland_raster, raster(grid_us_fine), method="bilinear"
    )
pop_grid_df <- raster::as.data.frame(pop_raster, xy=TRUE) %>%
    rename(pop=gpw_v4_population_density_rev11_2000_30_sec)

# create a popullation dataframe by states ------------------------------------
data(statepop)
statepop_df <- statepop %>%
    rename(
        state=abbr,
        pop=pop_2022
    )
geo_df2 <- merge(geo_df, statepop_df)

# evaluate an interpolated field ----------------------------------------------
# average by the whole period
coordinates(res_daily_long_sp2) <- ~ x + y
idw_gstat <- gstat(
  formula = load_daily ~ 1, # intercept-only model
  data = res_daily_long_sp2
)

heating_load_field <- predict(
  object = idw_gstat,
  newdata = grid_us_fine
)

dev.new()
plot(heating_load_field)
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")

# pdf("heat_load_2018_v3.pdf")
pdf("cooling_load_2018_v3.pdf")
# pdf("cooling_load_2018_v3.pdf")
plot(heating_load_field)
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")
dev.off()

# a single day ----------------------------------------------------------------
day_snapshot <- 1
month_snapshot <- 2

load_snapshot_df <- res_daily_long_sp %>%
  dplyr::filter(
    month(ymd(time_date)) %in% month_snapshot,
    day(ymd(time_date)) %in% day_snapshot
  )
coordinates(load_snapshot_df) <- ~ x + y

idw_gstat_test <- gstat(
  formula = load_daily ~ 1, # intercept-only model
  data = load_snapshot_df
)

heating_load_field_test <- predict(
  object = idw_gstat_test,
  newdata = grid_us_fine
)

dev.new()
plot(heating_load_field_test)
title(
  main = paste(day_snapshot, month_snapshot, "2018", sep="/")
)
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")

pdf(
    paste0(
        # "heat_load_2018_",
        "cooling_load_2018_",        
        paste(day_snapshot, month_snapshot, "2018", sep="-"),
        ".pdf")
)
# pdf("cooling_load_2018_v3.pdf")
plot(heating_load_field_test)
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")
dev.off()

# bi-linear interpolation -----------------------------------------------------
day_snapshot <- 1
month_snapshot <- 2

test_df <- res_daily_long_sp %>%
  dplyr::filter(
    year(ymd(time_date)) %in% 2018,
    month(ymd(time_date)) %in% month_snapshot,
    day(ymd(time_date)) %in% day_snapshot
  )

xh <- 0.5
x_margin <- 5
grid_us_test <- expand.grid(
  x = seq(
    from = round(min(res_daily_long_sp2$x)) - x_margin,
    to = round(max(res_daily_long_sp2$x)) + x_margin,
    by = xh
  ),
  y = seq(
    from = round(min(res_daily_long_sp2$y)) - x_margin,
    to = round(max(res_daily_long_sp2$y)) + x_margin,
    by = xh
  )
)  

test_interp <- interpp(
    x=test_df$x,
    y=test_df$y,
    z=test_df$load_daily,
    xo=grid_us_test$x,
    yo=grid_us_test$y, 
    linear = TRUE,
    extrap = TRUE, duplicate = "mean", dupfun = NULL,
    deltri = "shull"
)

test_interp_df <- data.frame(
    x=test_interp$x,
    y=test_interp$y,
    z=test_interp$z
    )

# p2 <- ggplot(data = test_interp_df, aes(x, y)) +
#   geom_raster(aes(fill = z)) +
#   scale_fill_distiller(palette = "Spectral", na.value = NA) + 
#   theme_classic() +
#   ggtitle("Linear interpolation")


usa <- map_data('usa')
p3 <- ggplot(data = usa, aes(x=long, y=lat, group=group)) +
  # geom_sf(color = "black", fill = NA) +
  geom_polygon(color='gray90')+
  geom_raster(
    data = test_interp_df, aes(x = x, y = y, fill = z),
    inherit.aes = FALSE
    ) +
  geom_raster(
    data = test_interp_df, aes(x = x, y = y, fill = z),
    inherit.aes = FALSE
    ) +  
  geom_point(
    data = geo_df,
    aes(x = x, y = y),
    color="white",
    inherit.aes = FALSE
    ) +
  # scale_fill_distiller(palette = "Spectral", na.value = NA)
  # "bpy.colors"
  # scale_fill_distiller(palette = "magma", na.value = NA)
  scale_fill_viridis_c(option = "C")
plot(p3)

# ggsave(
#     paste0(
#         "heating_load_2018_linpp_",
#         paste(day_snapshot, month_snapshot, "2018", sep="-"),
#         ".pdf")
# )
ggsave(
    paste0(
        "cooling_load_2018_linpp_",
        paste(day_snapshot, month_snapshot, "2018", sep="-"),
        ".pdf")
)



# grid_us_fine_load2 <- predict(
#   object = idw_temp,
#   newdata = grid_us_fine2
# )

# dev.new()
# plot(grid_us_fine_load2)
# plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")

# # vary parameters of the model ------------------------------------------------

# # neighbors <- nrow(res_daily_long_sp2)
# neighbors <- 50
# beta <- 2

# # idw_temp23 <- gstat(
# #   formula = log(load_daily) ~ 1, # intercept-only model
# #   data = res_daily_long_sp2,
# #   locations = ~x+y,
# #   nmax = neighbors,
# #   set = list(idp = beta)
# # )

# coordinates(res_daily_long_sp2) <- ~ x + y
# class(res_daily_long_sp2)

# idw_temp23 <- gstat(
#   formula = log(load_daily) ~ 1, # intercept-only model
#   data = res_daily_long_sp2,
#   # locations = ~x+y,
#   nmax = neighbors,
#   set = list(idp = beta)
# )

# # > range(res_daily_long_sp2$x)
# # [1] -122.07094  -69.38193
# # > range(res_daily_long_sp2$y)
# # [1] 27.76628 47.52891
# xy <- expand.grid(-115:-60, 25:50)
# names(xy) <- c("x","y")
# gridded(xy) = ~x+y

# g.dummy <- gstat(formula = z~1, dummy = TRUE, beta = 0,
#     model = vgm(1,"Exp",15), nmax = 10) 
# yy <- predict(g.dummy, xy, nsim = 4)


# idw_temp23 <- gstat(
#   formula = load_daily ~ 1, # intercept-only model
#   data = res_daily_long_sp2
# )
# yy <- predict(idw_temp23, xy)
# dev.new()
# plot(yy)
# plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")



# grid_us_fine_load23 <- predict(
#   object = idw_temp23,
#   newdata = xy
# )
# dev.new()
# plot(grid_us_fine_load23 )
# plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")

# pdf("HDD_2018.pdf")
# plot(grid_us_fine_load23 )
# plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")
# dev.off()

# # data_sp <- stConstruct(
# #     res_daily_long_sp2,
# #     space=c("x", "y"),
# #     time="time_date"
# # )

# # stConstruct(x, space, time, SpatialObj = NULL, TimeObj = NULL,
# # crs = CRS(as.character(NA)), interval, endTime)




# # neighbors <- length(heat_loc_df)
# # beta <- 2

# # idw_temp <- gstat(
# #   formula = load_daily ~ 1, # intercept-only model
# #   # locations=res_daily.loc,
# #   data = data_sp,
# #   nmax = neighbors,
# #   set = list(idp = beta)
# # )


# # i <- idw(log(overall_resload) ~ 1, res_daily.loc, grid_us, idp = 0.5, data=test)


# # idw_temp <- gstat(
# #   formula = overall_resload ~ 1, # intercept-only model
# #   locations=test.loc,
# #   data = test,
# #   nmax = neighbors,
# #   set = list(idp = beta)
# # )





# # # data list; each element is a list with the formula, locations, data, nvars, beta,
# # # etc., for a variable

# # idw_temp <- gstat(
# #   formula = ts_data ~ 1, # intercept-only model
# #   data = heat_loc_df,
# #   nmax = neighbors,
# #   set = list(idp = beta)
# # )

# # grid_us_load <- predict(
# #   object = idw_temp,
# #   newdata = grid_us
# # )


