import logging
import pandas as pd
import geopandas as gpd
import pathlib
import numpy as np

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("energyplus_aggregate.log"), logging.StreamHandler()],
)

# TODO Put the parameters into the config
SHARE_WATER_SH_DEMAND = 0.20
RATIO_SERV_TO_RESID = 1
SIMPLIFIED_WRMWATER = True
FILTER_ISLAND_PUMAS = False
DATA_IS_SCALED = (
    True  # post processing was done externally to workflow, presented as MWh
)

# TODO Revise the scaling part accounting for the projections
# Assuming that heating & cooling loads are ~25% of the electricity consumption
RESID_HEATING_TOTAL = 0.25 * 3.3e9  # [MWh]
RESID_COOLING_TOTAL = 0.25 * 3.3e9  # [MWh]

SERVICE_HEATING_TOTAL = 0.25 * 3.3e9  # [MWh]
SERVICE_COOLING_TOTAL = 0.25 * 3.3e9  # [MWh]

THERM_LOAD_SECTORS = [
    "residential heating",
    "services heating",
    "residential water",
    "services water",
    "residential cooling",
    "services cooling"
]


def get_state_id(
    state_fl_name,
    file_suff,
    abbrev_df,
    pumas_df,
):
    state_abbr = str(state_fl_name).replace(file_suff, "").upper()
    state_full_name = abbrev_df[abbrev_df.state == state_abbr].full_name.to_list()[0]
    state_id = pumas_df[pumas_df.State == state_full_name].STATEFIP.unique()[0]
    state_geoid = str(state_id)

    return state_geoid


def consolidate_pumas(
    data_path,
    states_abbr_df,
    puma_centroid_merged,
):
    # consolidating load profiles
    data_state_fls = list(data_path.iterdir())

    data_state_fls_clean = [
        fl_path
        for fl_path in data_state_fls
        if fl_path.name not in ["ak.csv", "hi.csv", "AK.csv", "HI.csv"]
    ]

    # To remove hidden files listed in the directory
    data_state_fls_clean = [
        fl_path for fl_path in data_state_fls if not fl_path.name.startswith(".")
    ]

    # a single time-series dataframe is needed to look-up for each PUMA -----------
    data_ts_national_list = [None] * len(data_state_fls_clean)

    logger.info("Build a consolidated national-wide load dataframe")
    for i, st_fl_path in enumerate(data_state_fls_clean):
        logger.info("Consolidation for " + str(st_fl_path.name))
        state_heat_df = pd.read_csv(
            pathlib.Path(st_fl_path),
        ).set_index("time")

        # the column names should correspond to GEOID to make further lookup work
        state_geoid = get_state_id(
            st_fl_path.name,
            file_suff=".csv",
            abbrev_df=states_abbr_df,
            pumas_df=puma_centroid_merged,
        )

        state_heat_df.columns = [state_geoid + col for col in state_heat_df.columns]
        data_ts_national_list[i] = state_heat_df

    data_ts_national_df = pd.concat(data_ts_national_list, axis=1)

    return data_ts_national_df

def find_onshore_pumas(puma_centroid_merged):
    if "distances" in puma_centroid_merged.columns:
        i_pumas_inside_onshore = puma_centroid_merged.distances == 0
        codes_pumas_inside_onshore = (
            puma_centroid_merged[i_pumas_inside_onshore].STATEFIP 
            + puma_centroid_merged[i_pumas_inside_onshore].PUMA
        )
    return codes_pumas_inside_onshore

def filter_by_island_pumas(puma_centroid_merged, data_df):
    onshore_codes = find_onshore_pumas(puma_centroid_merged)
    onshore_data_cols = resstock_heating_ts_national_df.columns.intersection(onshore_codes)
    
    return data_df[onshore_data_cols]

def lookup_bus_pumas(data_ts_national_df, bus, bus_pumas):
    # some PUMAs can be missed from the time-series data columns
    pumas_in_data = data_ts_national_df.columns.intersection(bus_pumas).to_list()
    bus_pumas_df = pd.DataFrame(index=data_ts_national_df.index)
    bus_pumas_df[bus] = data_ts_national_df[pumas_in_data].sum(axis=1)
    return bus_pumas_df


def add_level_column(df, level_name="residential space"):
    df.columns = pd.MultiIndex.from_product([[level_name], list(df.columns)])
    return df

# Implies that years are not duplicated each feature 
# & `scenario` column is present
# TODO add sanity-checks
def extract_growth_rates(
        scenario,
        growth_df,
        year,
        feature
    ):
    growth_rate = (
        growth_rates_df.query("(year == @year) & feature == @feature")[scenario]
    )
    return growth_rate

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("aggregate_energyplus")
        
    scenario = snakemake.params.thermal_scenario
    heat_year = snakemake.params.thermal_proj_year
    
    default_path = pathlib.Path(__file__).parent.parent.parent

    state_resstock_heat_dir = pathlib.Path(
        default_path, snakemake.input.state_resstock_heat_dir
    )
    state_resstock_wrmwater_dir = pathlib.Path(
        default_path, snakemake.input.state_resstock_wrmwater_dir
    )
    state_resstock_cool_dir = pathlib.Path(
        default_path, snakemake.input.state_resstock_cool_dir
    )

    state_comstock_heat_dir = pathlib.Path(
        default_path, snakemake.input.state_comstock_heat_dir
    )
    state_comstock_wrmwater_dir = pathlib.Path(
        default_path, snakemake.input.state_comstock_wrmwater_dir
    )
    state_comstock_cool_dir = pathlib.Path(
        default_path, snakemake.input.state_comstock_cool_dir
    )

    shapes_path = pathlib.Path(default_path, snakemake.input.shapes_path[0])
    puma_path = pathlib.Path(default_path, snakemake.input.puma_path)
    states_path = pathlib.Path(default_path, snakemake.input.states_path)

    heat_demand_path = pathlib.Path(default_path, snakemake.output.heat_demand_path[0])
    cool_demand_path = pathlib.Path(default_path, snakemake.output.cool_demand_path[0])
    growth_rates_path = pathlib.Path(default_path, snakemake.input.growth_path)

    # gis information is processed for the model and PUMA regions
    # to match demand time series for the model bus regions
    model_gdf = gpd.read_file(shapes_path)
    puma_gdf = gpd.read_file(puma_path).to_crs(model_gdf.crs)

    # abbreviation vs full names transformaiton is needed
    # to deal with transformaitons of PUMA codes
    states_abbr_df = pd.read_csv(states_path)

    puma_centroid = puma_gdf.copy()
    puma_centroid.geometry = puma_centroid.geometry.centroid
    puma_centroid_merged = gpd.sjoin_nearest(
        puma_centroid, model_gdf, how="left", distance_col="distances"
    )

    growth_rates_df = pd.read_csv(growth_rates_path)

    growth_rates_dict = {
        x: extract_growth_rates(
            scenario="reference case",
            growth_df=growth_rates_df,
            year=heat_year,
            feature=x
        ) for x in THERM_LOAD_SECTORS
    }

    # consolidating load profiles
    resstock_heating_ts_national_df = consolidate_pumas(
        data_path=state_resstock_heat_dir,
        states_abbr_df=states_abbr_df,
        puma_centroid_merged=puma_centroid_merged,
    )

    resstock_cooling_ts_national_df = consolidate_pumas(
        data_path=state_resstock_cool_dir,
        states_abbr_df=states_abbr_df,
        puma_centroid_merged=puma_centroid_merged,
    )

    comstock_heating_ts_national_df = consolidate_pumas(
        data_path=state_comstock_heat_dir,
        states_abbr_df=states_abbr_df,
        puma_centroid_merged=puma_centroid_merged,
    )

    comstock_cooling_ts_national_df = consolidate_pumas(
        data_path=state_comstock_cool_dir,
        states_abbr_df=states_abbr_df,
        puma_centroid_merged=puma_centroid_merged,
    )

    if FILTER_ISLAND_PUMAS:
        # assuming that islands are not necessarily connected 
        # with the mainland power system
        resstock_heating_ts_national_df = filter_by_island_pumas(
            puma_centroid_merged, data_df=resstock_heating_ts_national_df
        ) 
        resstock_cooling_ts_national_df = filter_by_island_pumas(
            puma_centroid_merged, data_df=resstock_cooling_ts_national_df
        )

        comstock_heating_ts_national_df = filter_by_island_pumas(
            puma_centroid_merged, data_df=comstock_heating_ts_national_df
        )  
        comstock_cooling_ts_national_df = filter_by_island_pumas(
            puma_centroid_merged, data_df=comstock_cooling_ts_national_df
        )

    # Scaling is needed to account for the effect of scenarios
    # on the thermal load
    resstock_heating_ts_national_df = (
        (growth_rates_dict["residential heating"].iloc[0]) * resstock_heating_ts_national_df
    )
    resstock_cooling_ts_national_df = (
        (growth_rates_dict["residential cooling"].iloc[0]) * resstock_cooling_ts_national_df
    )
    comstock_heating_ts_national_df = (
        (growth_rates_dict["services heating"].iloc[0]) * comstock_heating_ts_national_df
    )
    comstock_cooling_ts_national_df = (
        (growth_rates_dict["services cooling"].iloc[0]) * comstock_cooling_ts_national_df
    )       

    # time-series for each PUMA should be aggregated ------------------------------
    load_buses = puma_centroid_merged.name.unique()

    resstock_pumas_heating_list = [None] * len(load_buses)
    comstock_pumas_heating_list = [None] * len(load_buses)

    pumas_cooling_list = [None] * len(load_buses)

    logger.info("Aggregate by PUMAs")
    for i, bus in enumerate(load_buses):
        logger.info("Load aggreagtion for a bus: " + bus)
        bus_pumas = puma_centroid_merged[puma_centroid_merged.name == bus]["GEOID"]

        # resstock ------------------------------------------------------------
        resstock_pumas_heating_df = lookup_bus_pumas(
            data_ts_national_df=resstock_heating_ts_national_df,
            bus=bus,
            bus_pumas=bus_pumas,
        )
        resstock_pumas_heating_list[i] = resstock_pumas_heating_df

        resstock_pumas_cooling_df = lookup_bus_pumas(
            data_ts_national_df=resstock_cooling_ts_national_df,
            bus=bus,
            bus_pumas=bus_pumas,
        )

        # comstock ------------------------------------------------------------
        comstock_pumas_heating_df = lookup_bus_pumas(
            data_ts_national_df=comstock_heating_ts_national_df,
            bus=bus,
            bus_pumas=bus_pumas,
        )

        # Comstock has 15 minutes granularity and restock 1 hour,
        # Grouping all comstock data to hourly granularity
        comstock_pumas_heating_df = comstock_pumas_heating_df.groupby(
            np.arange(len(comstock_pumas_heating_df.index)) // 4
        ).sum()
        comstock_pumas_heating_df.index = resstock_pumas_cooling_df.index

        comstock_pumas_heating_list[i] = comstock_pumas_heating_df

        comstock_pumas_cooling_df = lookup_bus_pumas(
            data_ts_national_df=comstock_cooling_ts_national_df,
            bus=bus,
            bus_pumas=bus_pumas,
        )

        # Comstock has 15 minutes granularity and restock 1 hour,
        # Grouping all comstock data to hourly granularity
        comstock_pumas_cooling_df = comstock_pumas_cooling_df.groupby(
            np.arange(len(comstock_pumas_cooling_df.index)) // 4
        ).sum()
        comstock_pumas_cooling_df.index = resstock_pumas_cooling_df.index

        # for cooling we don't distinguish between the residential and services sector
        pumas_cooling_list[i] = resstock_pumas_cooling_df + comstock_pumas_cooling_df

    resstock_heating_load_aggreg_df = pd.concat(resstock_pumas_heating_list, axis=1)
    comstock_heating_load_aggreg_df = pd.concat(comstock_pumas_heating_list, axis=1)

    cooling_load_aggreg_df = pd.concat(pumas_cooling_list, axis=1)

    if SIMPLIFIED_WRMWATER:
        # A temporally solution for warm water
        resstock_wrmwater_load_aggreg_df = pd.DataFrame(
            index=resstock_heating_load_aggreg_df.index,
            data=[SHARE_WATER_SH_DEMAND * resstock_heating_load_aggreg_df.sum(axis=0)]
            * len(resstock_heating_load_aggreg_df.index),
        )
        comstock_wrmwater_load_aggreg_df = pd.DataFrame(
            index=comstock_heating_load_aggreg_df.index,
            data=[SHARE_WATER_SH_DEMAND * comstock_heating_load_aggreg_df.sum(axis=0)]
            * len(comstock_heating_load_aggreg_df.index),
        )
    else:
        resstock_wrmwater_ts_national_df = consolidate_pumas(
            data_path=state_comstock_wrmwater_dir,
            states_abbr_df=states_abbr_df,
            puma_centroid_merged=puma_centroid_merged,
        )
        comstock_wrmwater_ts_national_df = consolidate_pumas(
            data_path=state_comstock_wrmwater_dir,
            states_abbr_df=states_abbr_df,
            puma_centroid_merged=puma_centroid_merged,
        )

        if FILTER_ISLAND_PUMAS:
            # Assuming that islands are not nessecarily connected 
            # with the mainland power system
            resstock_wrmwater_ts_national_df = filter_by_island_pumas(
                puma_centroid_merged, data_df=resstock_wrmwater_ts_national_df
            )
            comstock_wrmwater_ts_national_df = filter_by_island_pumas(
                puma_centroid_merged, data_df=comstock_wrmwater_ts_national_df
            )        

        # Scaling is needed to account for the effect of scenarios
        # on the thermal load
        resstock_wrmwater_ts_national_df = (
            (growth_rates_dict["residential water"].iloc[0]) * resstock_wrmwater_ts_national_df
        )
        
        comstock_wrmwater_ts_national_df = (
            (growth_rates_dict["services water"].iloc[0]) * comstock_wrmwater_ts_national_df
        )
        
        resstock_pumas_wrmwater_list = [None] * len(load_buses)
        comstock_pumas_wrmwater_list = [None] * len(load_buses)

        for i, bus in enumerate(load_buses):
            logger.info("Load aggreagtion for a bus: " + bus)
            bus_pumas = puma_centroid_merged[puma_centroid_merged.name == bus]["GEOID"]
            # resstock ------------------------------------------------------------
            resstock_pumas_wrmwater_df = lookup_bus_pumas(
                data_ts_national_df=resstock_wrmwater_ts_national_df,
                bus=bus,
                bus_pumas=bus_pumas,
            )
            resstock_pumas_wrmwater_list[i] = resstock_pumas_wrmwater_df
            # comstock ------------------------------------------------------------
            comstock_pumas_wrmwater_df = lookup_bus_pumas(
                data_ts_national_df=comstock_wrmwater_ts_national_df,
                bus=bus,
                bus_pumas=bus_pumas,
            )
            comstock_pumas_wrmwater_list[i] = comstock_pumas_wrmwater_df

        resstock_wrmwater_load_aggreg_df = pd.concat(
            resstock_pumas_wrmwater_list, axis=1
        )
        comstock_wrmwater_load_aggreg_df = pd.concat(
            comstock_pumas_wrmwater_list, axis=1
        )

    if not DATA_IS_SCALED:
        # Scaling is needed to be consistent with the overall energy balanse
        resstock_heating_scale = (
            RESID_HEATING_TOTAL / resstock_heating_load_aggreg_df.sum().sum()
        )
        comstock_heating_scale = (
            SERVICE_HEATING_TOTAL / comstock_heating_load_aggreg_df.sum().sum()
        )

        resstock_heating_load_aggreg_df = (
            resstock_heating_scale * resstock_heating_load_aggreg_df
        )
        comstock_heating_load_aggreg_df = (
            comstock_heating_scale * comstock_heating_load_aggreg_df
        )

    # A multi-index dataframe is needed
    # 1) heating load has multiple components
    add_level_column(df=resstock_heating_load_aggreg_df, level_name="residential space")
    add_level_column(df=comstock_heating_load_aggreg_df, level_name="services space")
    add_level_column(
        df=resstock_wrmwater_load_aggreg_df, level_name="residential water"
    )
    add_level_column(df=comstock_wrmwater_load_aggreg_df, level_name="services water")

    heating_overall_load = pd.concat(
        [
            resstock_wrmwater_load_aggreg_df,
            resstock_heating_load_aggreg_df,
            comstock_wrmwater_load_aggreg_df,
            comstock_heating_load_aggreg_df,
        ],
        axis=1,
    )

    if not DATA_IS_SCALED:
        # 2) cooling
        cooling_scale = (
            RESID_COOLING_TOTAL + SERVICE_COOLING_TOTAL
        ) / cooling_load_aggreg_df.sum().sum()
        cooling_load_aggreg_df = cooling_scale * cooling_load_aggreg_df

    add_level_column(df=cooling_load_aggreg_df, level_name="space")

    # Year adjustments in the data to match the snapshot year
    snapshot_year = int(snakemake.params.snapshot_start[:4])
    data_year = pd.to_datetime(heating_overall_load.index).year.unique()[0]
    year_offset = data_year - snapshot_year
    heating_overall_load.index = pd.to_datetime(
        heating_overall_load.index
    ) - pd.DateOffset(years=year_offset)
    cooling_load_aggreg_df.index = pd.to_datetime(
        cooling_load_aggreg_df.index
    ) - pd.DateOffset(years=year_offset)

    heating_overall_load.to_csv(heat_demand_path)
    cooling_load_aggreg_df.to_csv(cool_demand_path)
