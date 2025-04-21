import logging
import pandas as pd
import geopandas as gpd
import numpy as np
import os
import pathlib
import pypsa
import sys

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("energyplus_aggregate.log"),
        logging.StreamHandler()
    ]
)

def get_state_id(
        state_fl_name,
        file_suff,
        abbrev_df,
        pumas_df,
):
    state_abbr = state_fl_name.replace(file_suff, "").upper()
    state_full_name = (
        abbrev_df[abbrev_df.state == state_abbr]
        .full_name
        .to_list()[0]
    )
    state_id = pumas_df[pumas_df.State == state_full_name].STATEFIP.unique()[0]
    state_geoid = str(state_id)

    return state_geoid

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("preprocess_demand_data")

    default_path = pathlib.Path(__file__).parent.parent.parent

    state_heat_dir = pathlib.Path(default_path, snakemake.input.state_heat_dir)
    state_cool_dir = pathlib.Path(default_path, snakemake.input.state_cool_dir)
    shapes_path = pathlib.Path(default_path, snakemake.input.shapes_path)
    puma_path = pathlib.Path(default_path, snakemake.input.puma_path)
    states_path = pathlib.Path(default_path, snakemake.input.states_path)

    heat_demand_path = pathlib.Path(default_path, snakemake.output.heat_demand)
    cool_demand_path = pathlib.Path(default_path, snakemake.output.cool_demand)

    # gis information is processed for the model and PUMA regions
    # to match demand time series for the model bus regions
    model_gdf = gpd.read_file(shapes_path)
    puma_gdf = gpd.read_file(puma_path).to_crs(model_gdf.crs)

    # abbreviation vs full names transformaiton is needed
    # to deal with transformaitons of PUMA codes
    states_abbr_df = pd.read_csv(states_path)

    puma_centroid = puma_gdf.copy()
    puma_centroid.geometry = puma_centroid.geometry.centroid
    puma_centroid_merged = gpd.sjoin_nearest(puma_centroid, model_gdf, how="left")

    # consolidating load profiles
    heating_state_fls = list(state_heat_dir.iterdir())
    cooling_state_fls = list(state_cool_dir.iterdir())

    # Alaska and Hawai are not included into the power model
    if "ak.csv" in heating_state_fls:
        heating_state_fls.remove("ak.csv")
    if "hi.csv" in heating_state_fls:
        heating_state_fls.remove("hi.csv")

    if "ak.csv" in cooling_state_fls:
        cooling_state_fls.remove("ak.csv")
    if "hi.csv" in cooling_state_fls:
        cooling_state_fls.remove("hi.csv")

    # a single time-series dataframe is needed to look-up for each PUMA -----------
    heating_ts_national_list = [None] * len(heating_state_fls)
    cooling_ts_national_list = [None] * len(cooling_state_fls)

    logger.info("Build a consolidated national-wide load dataframe")
    for i, st_fl in enumerate(heating_state_fls):
        logger.info("Heating consolidation for " + st_fl)
        state_heat_df = pd.read_csv(
            os.path.join(state_heat_dir, st_fl),
        ).set_index("time")

        # the column names should correspond to GEOID to make further lookup work
        state_geoid = get_state_id(
            st_fl,
            file_suff=".csv",
            abbrev_df=states_abbr_df,
            pumas_df=puma_centroid_merged,
        )
        state_heat_df.columns = [state_geoid + col for col in state_heat_df.columns]
        heating_ts_national_list[i] = state_heat_df

    heating_ts_national_df = pd.concat(heating_ts_national_list, axis=1)

    # cooling requires a special treatment
    for i, st_fl in enumerate(cooling_state_fls):
        logger.info("Cooling consolidation for " + st_fl)
        state_cool_df = pd.read_csv(
            os.path.join(state_cool_dir, st_fl),
        ).set_index("time")

        # the column names should correspond to GEOID to make further lookup work
        state_geoid = get_state_id(
            st_fl,
            file_suff=".csv",
            abbrev_df=states_abbr_df,
            pumas_df=puma_centroid_merged,
        )
        state_cool_df.columns = [state_geoid + col for col in state_cool_df.columns]

        cooling_ts_national_list[i] = state_cool_df

    cooling_ts_national_df = pd.concat(cooling_ts_national_list, axis=1)

    # time-series for each PUMA should be aggregated ------------------------------
    load_buses = puma_centroid_merged.name.unique()
    pumas_heating_list = [None] * len(load_buses)
    pumas_cooling_list = [None] * len(load_buses)

    logger.info("Aggregate by PUMAs")
    for i, bus in enumerate(load_buses):
        logger.info("Load aggreagtion for a bus: " + bus)
        bus_pumas = puma_centroid_merged[puma_centroid_merged.name == bus]["GEOID"]

        # some PUMAs can be missed from the time-series data columns
        pumas_in_heating_data = heating_ts_national_df.columns.intersection(
            bus_pumas
        ).to_list()
        pumas_in_cooling_data = cooling_ts_national_df.columns.intersection(
            bus_pumas
        ).to_list()

        pumas_heating_df = pd.DataFrame(index=heating_ts_national_df.index)
        pumas_heating_df[bus] = heating_ts_national_df[pumas_in_heating_data].sum(
            axis=1
        )
        pumas_heating_list[i] = pumas_heating_df

        pumas_cooling_df = pd.DataFrame(index=cooling_ts_national_df.index)
        pumas_cooling_df[bus] = cooling_ts_national_df[pumas_in_cooling_data].sum(
            axis=1
        )
        pumas_cooling_list[i] = pumas_cooling_df

    heating_load_aggreg_df = pd.concat(pumas_heating_list, axis=1)
    cooling_load_aggreg_df = pd.concat(pumas_cooling_list, axis=1)

    heating_load_aggreg_df.to_csv(heat_demand_path)
    cooling_load_aggreg_df.to_csv(cool_demand_path)
