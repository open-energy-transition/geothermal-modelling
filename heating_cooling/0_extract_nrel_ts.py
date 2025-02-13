# ---------------------------------------
# Extract ReStock/ComStock data by states
# ---------------------------------------

# paths must be adjusted
resstock_dir = "./_geothermal_/data/_NREL_profiles_/ResStock/"
comstock_dir = "./_geothermal_/data/_NREL_profiles_/ComStock/"

output_dir = "./_heating implementation_/_data processing_/_res_/"

import os
import pypsa

import geopandas as gpd
import pandas as pd

# define constants ------------------------------------------------------------
STATES=["AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", 
        "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", 
        "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", 
        "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", 
        "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", 
        "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY"]

HEATING_COLS_RESSTOCK = [
    "out.electricity.heating.energy_consumption.kwh",
    "out.fuel_oil.heating.energy_consumption.kwh",
    "out.natural_gas.heating.energy_consumption.kwh",
    "out.propane.heating.energy_consumption.kwh",
]

HEATING_COLS_COMSTOCK = [
    "out.electricity.heating.energy_consumption.kwh",
    "out.natural_gas.heating.energy_consumption.kwh",
    "out.other_fuel.heating.energy_consumption.kwh",
]

COOLING_COLS_RESSTOCK = [
    "out.electricity.cooling.energy_consumption.kwh",
    # "out.electricity.cooling_fans_pumps.energy_consumption.kwh"
]

COOLING_COLS_COMSTOCK = [
    "out.district_cooling.cooling.energy_consumption.kwh",
    "out.electricity.cooling.energy_consumption.kwh",
    "out.other_fuel.cooling.energy_consumption.kwh"
]

# parameters of the run -------------------------------------------------------

load_type = "heating" # "cooling" or "heating"

output_prefix = "residential"   # "residential" or "commercial"

if (output_prefix == "residential"):
    data_dir = resstock_dir
else:
    data_dir = comstock_dir 

if ((output_prefix == "residential") & (load_type == "heating")):
        load_type_columns = HEATING_COLS_RESSTOCK
if ((output_prefix == "commercial") & (load_type == "heating")):
        load_type_columns = HEATING_COLS_COMSTOCK
if ((output_prefix == "residential") & (load_type == "cooling")):
        load_type_columns = COOLING_COLS_RESSTOCK 
if ((output_prefix == "commercial") & (load_type == "cooling")):
        load_type_columns = COOLING_COLS_COMSTOCK                      

load_output_suff = "_" + load_type + "_by_states.csv"

# define functions ------------------------------------------------------------
def extract_state_load(state_path, data_cols=COOLING_COLS_RESSTOCK):
    state_fls = os.listdir(state_path)
    load_dfs_list = [None] * len(state_fls)

    for i, fls in enumerate(state_fls):
        load_df = pd.read_csv(
            os.path.join(state_fls_path, fls)
        )
        load_df["overall_resload"] = load_df[data_cols].sum(axis=1)
        load_df["heat_consumer"] = fls
        print(load_df[["upgrade", "in.state", "timestamp", "overall_resload"]].head(2))

        load_dfs_list[i] = load_df[["upgrade", "in.state", "timestamp", "overall_resload"]]

    state_resid_df = pd.concat(load_dfs_list)
    consolid_profile = (
        state_resid_df[["timestamp", "in.state", "overall_resload"]]
        .groupby(["timestamp"])
        .sum()
    )

    consolid_profile["in.state"] = consolid_profile["in.state"].str[0:2]

    return consolid_profile 

# run extraction --------------------------------------------------------------

# assuming that heat consumers are same for all the states
state_fls_path = os.path.join(data_dir, STATES[0])
state_fls = os.listdir(state_fls_path)
load_dfs_list = [None] * len(state_fls)

# it's convenient to put all the data into a list
states = STATES
dfs_list = [None] * len(states)

for i, state in enumerate(states):
    print(output_prefix)
    print(i, state)
    state_fls_path = os.path.join(data_dir, STATES[i])
    print(state_fls_path)
    state_df = extract_state_load(state_fls_path, data_cols=load_type_columns)
    dfs_list[i] = state_df

state_resid_df = pd.concat(dfs_list)

state_resid_df.to_csv(output_dir + output_prefix + load_output_suff)