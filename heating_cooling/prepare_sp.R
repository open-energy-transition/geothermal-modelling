# Following the tutorial by Freie Uni Berlin
# https://www.geo.fu-berlin.de/en/v/soga-r/Advances-statistics/Geostatistics/Inverse-Distance-Weighting-IDW/IDW-Interpolation-of-Weather-Data/index.html

rm(list=ls())

data_dir <- "./"
data_fl <- "residential_load_by_states.csv"

library(ggplot2)
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
    "us_states_geo.csv",
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
#data_df <- merge(test, load_df)


# # generate grid -------------------------------------------------------------
# # set the grid cell size in meter
# xh <- 100000

xh <- 5
grid_us <- expand.grid(
  x = seq(
    from = round(min(geo_df$Longitude)),
    to = round(max(geo_df$Longitude)),
    by = xh
  ),
  y = seq(
    from = round(min(geo_df$Latitude)),
    to = round(max(geo_df$Latitude)),
    by = xh
  )
)
coordinates(grid_us) <- ~ x + y
class(grid_us)

gridded(grid_us) <- TRUE
class(grid_us)

# grid_us <- expand.grid(
#   x = seq(
#     from = round(extent_us@xmin),
#     to = round(extent_us@xmax),
#     by = xh
#   ),
#   y = seq(
#     from = round(extent_us@ymin),
#     to = round(extent_us@ymax),
#     by = xh
#   )
# )

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

# graphical check -------------------------------------------------------------
plot(grid_us,
  main = paste("Centroids by States\n and Interpolation Grid"),
  col = "grey",
  cex.main = 0.9
)
# plot(east_germany_states_sp, add = TRUE, border = "red")
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "blue")

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
    summarise(across(AL:WI, sum))

geo_df <- geo_df
geo_df$y = geo_df[["Latitude"]]
geo_df$x = geo_df[["Longitude"]]
coordinates(geo_df) = ~x+y    
res_daily.loc <- geo_df

res_daily_long <- res_daily %>% 
    pivot_longer(
        cols=AL:WI,
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
res_daily_long_sp2 <- res_daily_long_sp[-((nrow(res_daily_long_sp) - 1):nrow(res_daily_long_sp)),]


idw_temp <- gstat(
  formula = load_daily ~ 1, # intercept-only model
  data = res_daily_long_sp2,
  locations = ~x+y,
  nmax = neighbors,
  set = list(idp = beta)
)

grid_us_load <- predict(
  object = idw_temp,
  newdata = grid_us
)

# try a finer mesh ------------------------------------------------------------
xh <- 1
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

grid_us_fine_load <- predict(
  object = idw_temp,
  newdata = grid_us_fine
)

plot(grid_us_fine_load)
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")

# try a very fine mesh ------------------------------------------------------------
xh <- 0.1
x_margin <- 5
grid_us_fine2 <- expand.grid(
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
coordinates(grid_us_fine2) <- ~x+y
class(grid_us_fine2)

gridded(grid_us_fine2) <- TRUE
class(grid_us_fine2)

grid_us_fine_load2 <- predict(
  object = idw_temp,
  newdata = grid_us_fine2
)

dev.new()
plot(grid_us_fine_load2)
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")

# vary parameters of the model ------------------------------------------------

# neighbors <- nrow(res_daily_long_sp2)
neighbors <- 50
beta <- 2

# idw_temp23 <- gstat(
#   formula = log(load_daily) ~ 1, # intercept-only model
#   data = res_daily_long_sp2,
#   locations = ~x+y,
#   nmax = neighbors,
#   set = list(idp = beta)
# )

coordinates(res_daily_long_sp2) <- ~ x + y
class(res_daily_long_sp2)

idw_temp23 <- gstat(
  formula = log(load_daily) ~ 1, # intercept-only model
  data = res_daily_long_sp2,
  # locations = ~x+y,
  nmax = neighbors,
  set = list(idp = beta)
)

# > range(res_daily_long_sp2$x)
# [1] -122.07094  -69.38193
# > range(res_daily_long_sp2$y)
# [1] 27.76628 47.52891
xy <- expand.grid(-115:-60, 25:50)
names(xy) <- c("x","y")
gridded(xy) = ~x+y

g.dummy <- gstat(formula = z~1, dummy = TRUE, beta = 0,
    model = vgm(1,"Exp",15), nmax = 10) 
yy <- predict(g.dummy, xy, nsim = 4)


idw_temp23 <- gstat(
  formula = load_daily ~ 1, # intercept-only model
  data = res_daily_long_sp2
)
yy <- predict(idw_temp23, xy)
dev.new()
plot(yy)
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")



grid_us_fine_load23 <- predict(
  object = idw_temp23,
  newdata = xy
)
dev.new()
plot(grid_us_fine_load23 )
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")

pdf("HDD_2018.pdf")
plot(grid_us_fine_load23 )
plot(heat_loc_df, add = TRUE, pch = 19, cex = 0.5, col = "white")
dev.off()

# data_sp <- stConstruct(
#     res_daily_long_sp2,
#     space=c("x", "y"),
#     time="time_date"
# )

# stConstruct(x, space, time, SpatialObj = NULL, TimeObj = NULL,
# crs = CRS(as.character(NA)), interval, endTime)




# neighbors <- length(heat_loc_df)
# beta <- 2

# idw_temp <- gstat(
#   formula = load_daily ~ 1, # intercept-only model
#   # locations=res_daily.loc,
#   data = data_sp,
#   nmax = neighbors,
#   set = list(idp = beta)
# )


# i <- idw(log(overall_resload) ~ 1, res_daily.loc, grid_us, idp = 0.5, data=test)


# idw_temp <- gstat(
#   formula = overall_resload ~ 1, # intercept-only model
#   locations=test.loc,
#   data = test,
#   nmax = neighbors,
#   set = list(idp = beta)
# )





# # data list; each element is a list with the formula, locations, data, nvars, beta,
# # etc., for a variable

# idw_temp <- gstat(
#   formula = ts_data ~ 1, # intercept-only model
#   data = heat_loc_df,
#   nmax = neighbors,
#   set = list(idp = beta)
# )

# grid_us_load <- predict(
#   object = idw_temp,
#   newdata = grid_us
# )


