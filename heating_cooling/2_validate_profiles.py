# ------------------------------------------
# Generate validation dataframes by states:
# 1) aggregation by days and load users
# 2) evaluation of HDDs/CDDs for each state
# (taken for it's centroid)
# 3) adding meteo-corrections following
# ------------------------------------------

# TODO:
# 1) add a correction to the wind speed (an updated cutout is needed with v_10m)
# 2) add a blending effect for cooling


import atlite
import os
import pypsa

import geopandas as gpd
import matplotlib.pyplot as plt
import xarray as xr

import datetime as dt
import numpy as np
import pandas as pd

STATES = ['AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DC', 'DE', 'FL', 'GA', 'IA',
       'ID', 'IL', 'IN', 'KS', 'KY', 'LA', 'MA', 'MD', 'ME', 'MI', 'MN',
       'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NJ', 'NM', 'NV', 'NY',
       'OH', 'OK', 'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VA',
       'VT', 'WA', 'WI', 'WV', 'WY']

output_prefix = "commercial"   # "residential" or "commercial"
load_type = "cooling" # "cooling" or "heating"

t_threshold = 22

# folders management ----------------------------------------------------------

# paths must be adjusted
cutout_path = "./pypsa-earth/cutouts"
cutout_fl = "us-2018-era5.nc"

ref_data_dir = "./_heating implementation_/_data processing_/_res_/timeseries"
ref_data_fl = output_prefix + "_nrel_" + load_type + "_2018_ref.csv" 

res_dir = "./_heating implementation_/_data processing_/_res_" 

# calculation of climate parameters -------------------------------------------
def RH(t, td):
    """
    t is the ambient air temperature [C]
    td is the dew point [C]
    RH is the relative humidity [%]
    """
    b = 17.625
    c = 243.04
    a = ((c * b) * (td - t)) / ((c + t) * (c + td))
    f = 100 * np.exp(a)
    return f

# https://www.conservationphysics.org/atmcalc/atmoclc1.html
def p_s(t):
    """
    p_s is the water vapor saturation pressure [Pa]
    t is ambient air temperature [C]
    """
    p_s = 610.78 * np.exp( t / ( t + 238.3 ) * 17.2694 )
    return p_s

# kg water vapour / kg dry air = 0.62 *10-5 *p
def AH(t, RH):
    """
    t is the ambient air temperature [C]
    RH is the relative humidity [-]
    """
    p = RH * p_s(t)
    ah = p * 0.62e-5
    return ah        

# pure PyPSA-Eur approach -----------------------------------------------------
def convert_heat_demand1(ds, threshold, a, constant, hour_shift):
    # Temperature is in Kelvin; take daily average
    T = ds["temperature"]
    T = T.assign_coords(
        time=(T.coords["time"] + np.timedelta64(dt.timedelta(hours=hour_shift)))
    )
    T = T.resample(time="1D").mean(dim="time") 

    threshold += 273.15
    heat_demand = a * (threshold - T)

    heat_demand = heat_demand.clip(min=0.0)

    return (constant + heat_demand).rename("heat_demand")


def heat_demand1(cutout, threshold=15.0, a=1.0, constant=0.0, hour_shift=0.0, **params):
    return cutout.convert_and_aggregate(
        convert_func=convert_heat_demand2,
        threshold=threshold,
        a=a,
        constant=constant,
        hour_shift=hour_shift,
        **params,
    )     

def convert_cooling_demand1(ds, threshold, a, constant, hour_shift):
    # Temperature is in Kelvin; take daily average
    T = ds["temperature"]
    T = T.assign_coords(
        time=(T.coords["time"] + np.timedelta64(dt.timedelta(hours=hour_shift)))
    )
    T = T.resample(time="1D").mean(dim="time") 

    threshold += 273.15
    cooling_demand = a * (T - threshold)

    cooling_demand = cooling_demand.clip(min=0.0)

    return (constant + cooling_demand).rename("cooling_demand")


def cooling_demand1(cutout, threshold=23.0, a=1.0, constant=0.0, hour_shift=0.0, **params):
    return cutout.convert_and_aggregate(
        convert_func=convert_cooling_demand1,
        threshold=threshold,
        a=a,
        constant=constant,
        hour_shift=hour_shift,
        **params,
    )

# adding meteo-corrections styled according to demand-ninja -------------------
def convert_heat_demand2(ds, threshold, a, constant, hour_shift):
    # Temperature is in Kelvin; take daily average
    T = ds["temperature"]
    T = T.assign_coords(
        time=(T.coords["time"] + np.timedelta64(dt.timedelta(hours=hour_shift)))
    )
    T = T.resample(time="1D").mean(dim="time")

    #(0.012 * cutout_atl.data["influx_toa"][1]).plot()
    dT_radiation = 0.012 * ds["influx_toa"]
    dT_radiation = dT_radiation.assign_coords(
        time=(dT_radiation.coords["time"] + np.timedelta64(dt.timedelta(hours=hour_shift)))
    )
    dT_radiation = dT_radiation.resample(time="1D").mean(dim="time")    

    threshold += 273.15
    heat_demand = a * (threshold - T - dT_radiation)

    heat_demand = heat_demand.clip(min=0.0)

    return (constant + heat_demand).rename("heat_demand")


def heat_demand2(cutout, threshold=15.0, a=1.0, constant=0.0, hour_shift=0.0, **params):
    return cutout.convert_and_aggregate(
        convert_func=convert_heat_demand2,
        threshold=threshold,
        a=a,
        constant=constant,
        hour_shift=hour_shift,
        **params,
    )     

def convert_heat_demand3(ds, threshold, a, constant, hour_shift):
    # Temperature is in Kelvin; take daily average
    T = ds["temperature"]
    T = T.assign_coords(
        time=(T.coords["time"] + np.timedelta64(dt.timedelta(hours=hour_shift)))
    )
    T = T.resample(time="1D").mean(dim="time")

    #(0.012 * cutout_atl.data["influx_toa"][1]).plot()
    dT_radiation = 0.012 * ds["influx_toa"]
    dT_radiation = dT_radiation.assign_coords(
        time=(dT_radiation.coords["time"] + np.timedelta64(dt.timedelta(hours=hour_shift)))
    )
    dT_radiation = dT_radiation.resample(time="1D").mean(dim="time")

    T_apparent = T + dT_radiation

    # T_smooth = T_apparent.rolling(time=3, min_periods=2, center=True).mean()
    T_smooth = T_apparent.rolling(time=5, min_periods=2, center=True).mean()

    ## numbagg >= 0.2.1 is required for rolling_exp but currently version {_NUMBAGG_VERSION} is installed
    #T_smooth = T.rolling_exp(time=3).mean()

    threshold += 273.15
    heat_demand = a * (threshold - T_smooth)

    heat_demand = heat_demand.clip(min=0.0)

    return (constant + heat_demand).rename("heat_demand")


def heat_demand3(cutout, threshold=15.0, a=1.0, constant=0.0, hour_shift=0.0, **params):
    return cutout.convert_and_aggregate(
        convert_func=convert_heat_demand3,
        threshold=threshold,
        a=a,
        constant=constant,
        hour_shift=hour_shift,
        **params,
    )

def convert_heat_demand5(ds, threshold, a, constant, hour_shift):
    # Temperature is in Kelvin; take daily average
    T = ds["temperature"]
    T = T.assign_coords(
        time=(T.coords["time"] + np.timedelta64(dt.timedelta(hours=hour_shift)))
    )
    T = T.resample(time="1D").mean(dim="time")

    # 0.012 ± 0.008 °C W−1 m−2
    dT_radiation = 0.012 * ds["influx_toa"]
    dT_radiation = dT_radiation.assign_coords(
        time=(dT_radiation.coords["time"] + np.timedelta64(dt.timedelta(hours=hour_shift)))
    )
    dT_radiation = dT_radiation.resample(time="1D").mean(dim="time")

    Td = ds["dewpoint temperature"]
    phi = RH((T - 273.15), (Td - 273.15))

    # 0.050 ± 0.065 °C (g kg−1)−1
    dT_apparent_phi = 50 * AH((T - 273.15), phi/100)

    T_apparent = T + dT_radiation + dT_apparent_phi

    T_smooth = T_apparent.rolling(time=3, min_periods=2, center=True).mean()

    ## numbagg >= 0.2.1 is required for rolling_exp but currently version {_NUMBAGG_VERSION} is installed
    #T_smooth = T.rolling_exp(time=3).mean()

    threshold += 273.15
    heat_demand = a * (threshold - T_smooth)

    heat_demand = heat_demand.clip(min=0.0)

    return (constant + heat_demand).rename("heat_demand")


def heat_demand5(cutout, threshold=15.0, a=1.0, constant=0.0, hour_shift=0.0, **params):
    return cutout.convert_and_aggregate(
        convert_func=convert_heat_demand5,
        threshold=threshold,
        a=a,
        constant=constant,
        hour_shift=hour_shift,
        **params,
    )

def convert_cooling_demand5(ds, threshold, a, constant, hour_shift):
    # Temperature is in Kelvin; take daily average
    T = ds["temperature"]
    T = T.assign_coords(
        time=(T.coords["time"] + np.timedelta64(dt.timedelta(hours=hour_shift)))
    )
    T = T.resample(time="1D").mean(dim="time") 

    # 0.012 ± 0.008 °C W−1 m−2
    dT_radiation = 0.012 * ds["influx_toa"]
    dT_radiation = dT_radiation.assign_coords(
        time=(dT_radiation.coords["time"] + np.timedelta64(dt.timedelta(hours=hour_shift)))
    )
    dT_radiation = dT_radiation.resample(time="1D").mean(dim="time")

    Td = ds["dewpoint temperature"]
    phi = RH((T - 273.15), (Td - 273.15))

    # 0.050 ± 0.065 °C (g kg−1)−1
    dT_apparent_phi = 50 * AH((T - 273.15), phi/100)

    T_apparent = T + dT_radiation + dT_apparent_phi

    T_smooth = T_apparent.rolling(time=3, min_periods=2, center=True).mean()

    ## numbagg >= 0.2.1 is required for rolling_exp but currently version {_NUMBAGG_VERSION} is installed
    #T_smooth = T.rolling_exp(time=3).mean()


    threshold += 273.15
    cooling_demand = a * (T_smooth - threshold)

    cooling_demand = cooling_demand.clip(min=0.0)

    return (constant + cooling_demand).rename("cooling_demand")


def cooling_demand5(cutout, threshold=23.0, a=1.0, constant=0.0, hour_shift=0.0, **params):
    return cutout.convert_and_aggregate(
        convert_func=convert_cooling_demand1,
        threshold=threshold,
        a=a,
        constant=constant,
        hour_shift=hour_shift,
        **params,
    )



# read inputs -----------------------------------------------------------------

print("ts_" + output_prefix + "_nrel_" + load_type + "_t_th_" + str(t_threshold) + "_2018_ref.csv")

# reference data NREL database
ref_df = pd.read_csv(os.path.join(ref_data_dir, ref_data_fl))

# a cutout for the reference year
cutout_atl = atlite.Cutout(os.path.join(cutout_path, cutout_fl))

# calculate heat demand using HDD/CDD representation --------------------------

if load_type == "cooling":
    test1_ts = convert_cooling_demand1(
        ds=cutout_atl.data, threshold=t_threshold, a=1, constant=0, hour_shift=0
    )
else:
    test1_ts = convert_heat_demand1(
        ds=cutout_atl.data, threshold=t_threshold, a=1, constant=0, hour_shift=0
    )

# validation run --------------------------------------------------------------

valid_ts_list = [None] * len(STATES)
valid_corr_df = pd.DataFrame(
    {
        "state": STATES,
        "corr_with_model": None
    }
)

for i, st in enumerate(STATES):
    x_state = ref_df[ref_df.state == st]["x"].iloc[0]
    y_state = ref_df[ref_df.state == st]["y"].iloc[0]

    mod_demand_ts_df = (
        test1_ts
        .sel(x=x_state, y=y_state, method='nearest')
        .to_dataframe()
    )

    ref_demand_ts = (
        #ref_df[(pd.to_datetime(ref_df.time_date).dt.month == 1) & (ref_df.state == st) & (pd.to_datetime(ref_df.time_date).dt.year == 2018)]
        ref_df[(ref_df.state == st) & (pd.to_datetime(ref_df.time_date).dt.year == 2018)]
        .sort_values("time_date")[["time_date", "load_daily", "state"]]
    )

    ref_demand_ts.time_date = ref_demand_ts.time_date + " 00:00:00"
    ref_demand = ref_demand_ts.set_index(
        pd.DatetimeIndex(ref_demand_ts["time_date"], freq="1D")
    )
    valid_df = mod_demand_ts_df.join(ref_demand)
    valid_ts_list[i] = valid_df

    
    if load_type == "cooling":
        corr_val = valid_df["cooling_demand"].corr(valid_df["load_daily"])
        mean_mod_dem = mod_demand_ts_df["cooling_demand"].mean()
    else:
        corr_val = valid_df["heat_demand"].corr(valid_df["load_daily"])
        mean_mod_dem = mod_demand_ts_df["heat_demand"].mean()

    mean_nrel_dem = ref_demand.load_daily.mean()
    mean_min_ratio_nrel = ref_demand.load_daily.mean() / ref_demand.load_daily.min()
    min_nrel_dem = ref_demand.load_daily.min()

    valid_corr_df.loc[i, "state"] = st
    valid_corr_df.loc[i, "corr_with_model"] = corr_val

    valid_corr_df.loc[i, "mean_mod_dem"] = mean_mod_dem
    valid_corr_df.loc[i, "mean_nrel_dem"] = mean_nrel_dem
    valid_corr_df.loc[i, "mean_min_ratio_nrel"] = mean_min_ratio_nrel
    valid_corr_df.loc[i, "min_nrel_dem"] = min_nrel_dem

    print(st, round(corr_val, 2))

valid_corr_df.to_csv(
    os.path.join(
        res_dir,
        output_prefix + "_nrel_" + load_type + "_t_th_" + str(t_threshold) + "_2018_baseline.csv" 
    )
)

valid_ts_df = pd.concat(valid_ts_list)
valid_ts_df.to_csv(
    os.path.join(
        res_dir,
        "ts_" + output_prefix + "_nrel_" + load_type + "_t_th_" + str(t_threshold) + "_2018_baseline.csv" 
    )
) 


