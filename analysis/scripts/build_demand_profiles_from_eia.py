"""
Converts annual spatial utility demand data into time series profiles

Relevant Settings
-----------------

.. code:: yaml

    US:
        demand_modelling:
            demand_year:


.. seealso::
    Documentation of the configuration file ``config.usa_baseline.yaml``


Inputs
------

- ``shapes/country_shapes.geojson``

- ``analysis/gdrive_data/data/electricity_demand_data/demand_data/table_10_EIA_utility_sales.xlsx``

- ``analysis/gdrive_data/data/electricity_demand_data/demand_data/Electric_Retail_Service_Territories.geojson``

- ``analysis/gdrive_data/data/shape_files/gadm41_USA_1.json``

- ``analysis/gdrive_data/data/electricity_demand_data/use_es_capita.xlsx``

- ``analysis/gdrive_data/data/electricity_demand_data/HS861 2010-.xlsx``


Outputs
-------

- ``resources/demand_profiles.csv``


Description
-----------
"""

import pathlib
import datetime as dt
import pandas as pd
import geopandas as gpd
from _helpers_usa import get_colors
import numpy as np
import pypsa
import os


def parse_inputs(default_path, distance_crs):
    """
    Load all input data

    Parameters
    ----------
    default_path: str
        current total directory path
    distance_crs: CRS code
        Co-ordinate reference system to retrieve coordinate in 'm'

    Returns
    -------
    df_ba_demand: pandas dataframe
        Balancing Authority demand profiles
    gdf_ba_shape: geopandas dataframe
        Balancing Authority shapes
    df_utility_demand: geopandas dataframe
        Output of preprocess_demand_data - demand mapped to utility level shapes and holes
    pypsa_network: pypsa
        network to obtain pypsa bus information
    """
    BA_demand_path1 = pathlib.Path(default_path, snakemake.input.BA_demand_path1[0])
    BA_demand_path2 = pathlib.Path(default_path, snakemake.input.BA_demand_path2[0])
    df_ba_demand1 = pd.read_csv(BA_demand_path1, index_col="period")
    df_ba_demand2 = pd.read_csv(BA_demand_path2, index_col="period")
    df_ba_demand = df_ba_demand1._append(df_ba_demand2[3:])
    # df_ba_demand = df_ba_demand.replace(0, np.nan)
    # df_ba_demand = df_ba_demand.dropna(axis=1)

    balancing_authority_shapefile = pathlib.Path(
        default_path, snakemake.input.BA_shape_path
    )
    gdf_ba_shape = gpd.read_file(balancing_authority_shapefile)
    gdf_ba_shape = gdf_ba_shape.to_crs(distance_crs)

    utility_demand_path = pathlib.Path(
        default_path, snakemake.input.utility_demand_path
    )
    df_utility_demand = gpd.read_file(utility_demand_path)
    # Todo: index_right no longer an available in the upstream packages - needs to be modified
    df_utility_demand.rename(columns={"index_right": "index_right_1"}, inplace=True)

    df_utility_demand = df_utility_demand.to_crs(distance_crs)

    pypsa_network_path = pathlib.Path(
        default_path, snakemake.input.pypsa_network_path[0]
    )
    pypsa_network = pypsa.Network(pypsa_network_path)

    return df_ba_demand, gdf_ba_shape, df_utility_demand, pypsa_network


def build_demand_profiles(
    df_utility_demand, df_ba_demand, gdf_ba_shape, pypsa_network, geo_crs, distance_crs
):
    """
    Build spatiotemporal demand profiles

    Parameters
    ----------
    df_utility_demand: geopandas dataframe
        Output of preprocess_demand_data - demand mapped to utility level shapes and holes
    df_ba_demand: pandas dataframe
        Balancing Authority demand profiles
    gdf_ba_shape: geopandas dataframe
        Balancing Authority shapes
    pypsa_network: pypsa
        network to obtain pypsa bus information
    geo_crs: CRS code
        Co-ordinate reference system to retrieve coordinates in 'geometrical degrees'
    distance_crs: CRS code
        Co-ordinate reference system to retrieve coordinate in 'm'

    Returns
    -------
    df_demand_bus_timeshifted: pandas dataframe
        bus-wise demand profiles
    """
    # Obtaining the centroids of the Utility demands
    df_utility_centroid = df_utility_demand.copy()
    df_utility_centroid.geometry = df_utility_centroid.geometry.centroid

    # Removing those shapes which do not have a representation in the demand timeseries data
    demand_columns = df_ba_demand.columns
    demand_columns = [x.split("_")[1] for x in demand_columns if x.endswith("_D")]
    drop_shape_rows = list(set(gdf_ba_shape["EIAcode"].tolist()) - set(demand_columns))
    gdf_ba_shape_filtered = gdf_ba_shape[~gdf_ba_shape["EIAcode"].isin(drop_shape_rows)]
    gdf_ba_shape_filtered.geometry = gdf_ba_shape_filtered.geometry.centroid
    gdf_ba_shape_filtered["color_ba"] = get_colors(len(gdf_ba_shape_filtered))

    df_utility_centroid = gpd.sjoin_nearest(
        df_utility_centroid, gdf_ba_shape_filtered, how="left"
    )
    # Todo: index_right no longer an available in the upstream packages - needs to be modified
    df_utility_centroid.rename(columns={"index_right": "index_right_2"}, inplace=True)

    # temporal scaling factor
    df_utility_centroid["temp_scale"] = df_utility_centroid.apply(
        lambda x: df_ba_demand[f"E_{x['EIAcode']}_D"].sum()
        / 1e3
        / x["Sales (Megawatthours)"],
        axis=1,
    )

    # Mapping demand utilities to nearest PyPSA bus
    df_reqd = pypsa_network.buses.query('carrier == "AC"')
    pypsa_gpd = gpd.GeoDataFrame(
        df_reqd, geometry=gpd.points_from_xy(df_reqd.x, df_reqd.y), crs=geo_crs
    )
    pypsa_gpd = pypsa_gpd.to_crs(distance_crs)
    pypsa_gpd["color"] = get_colors(len(pypsa_gpd))

    df_utility_centroid = gpd.sjoin_nearest(df_utility_centroid, pypsa_gpd, how="left")
    # Todo: index_right no longer an available in the upstream packages - needs to be modified
    df_utility_centroid.rename(columns={"index_right": "PyPSA_bus"}, inplace=True)

    df_demand_bus = pd.DataFrame(
        index=pd.to_datetime(df_ba_demand.index),
        columns=df_reqd.index.tolist(),
        data=0.0,
    )

    for col in df_demand_bus.columns:
        utility_rows = df_utility_centroid.query("PyPSA_bus == @col")
        for i in np.arange(0, len(utility_rows)):
            row = utility_rows.iloc[i]
            if not np.isnan(row["temp_scale"]):
                demand_data = df_ba_demand[f"E_{row['EIAcode']}_D"].copy()
                demand_data /= row["temp_scale"]
                demand_data /= 1000  # Converting to MWh
                df_demand_bus[col] += demand_data.tolist()

    # The EIA profiles start at 6:00:00 hours on 1/1 instead of 00:00:00 hours - rolling over the time series to start at 00:00 hours
    df_demand_bus_timeshifted = df_demand_bus[-9:-3]._append(df_demand_bus[:-9])
    df_demand_bus_timeshifted = df_demand_bus_timeshifted[:8760]

    df_demand_bus_timeshifted = df_demand_bus_timeshifted.reset_index().rename(
        columns={"period": "time"}
    )
    df_demand_bus_timeshifted.set_index("time", inplace=True)
    df_demand_bus_timeshifted.index = pypsa_network.snapshots

    return df_demand_bus_timeshifted


def modify_pypsa_network_demand(df_demand_profiles, pypsa_network, pypsa_network_path):
    # To-do: Change to using pypsa_network.snapshot_weightings
    time_resolution = 8760 / len(pypsa_network.snapshots)
    # Groupby time resolution and then convert from kWh -> kW
    df_demand_profiles.index = np.arange(0, len(df_demand_profiles.index))
    df_demand_profiles_agg = (
        df_demand_profiles.groupby(
            df_demand_profiles.index // int(time_resolution)
        ).sum()
        / time_resolution
    )
    df_demand_profiles_agg = df_demand_profiles_agg[0 : len(pypsa_network.snapshots)]
    df_demand_profiles_agg.index = pypsa_network.snapshots
    pypsa_network.loads_t.p_set = df_demand_profiles_agg
    pypsa_network.export_to_netcdf(pypsa_network_path)


def read_scaling_factor(demand_scenario, horizon, default_path):
    """
    Reads scaling factor for future projections
    Parameters
    ----------
    demand_scenario: str
        Future demand projection scenario
    horizon: int
        Horizon for demand projection
    Returns
    -------
    scaling_factor: pandas dataframe
        Scaling factor of demand projection scenario
    """
    horizon = (
        2024 if horizon == 2025 else horizon
    )  # select 2024 demand projection for 2025 horizon
    foldername = os.path.join(default_path, snakemake.input.demand_projections_path)
    filename = f"Scaling_Factor_{demand_scenario}_Moderate_{horizon}_by_state.csv"
    scaling_factor = pd.read_csv(os.path.join(foldername, filename), sep=";")
    scaling_factor["time"] = pd.to_datetime(scaling_factor["time"])
    # logger.info(f"Read {filename} for scaling the demand for {horizon}.")
    return scaling_factor


def scale_demand_profiles(df_demand_profiles, pypsa_network, scaling_factor):
    """
    Scales demand profiles for each state based on the NREL EFS demand projections
    Parameters
    ----------
    df_demand_profiles: pandas dataframe
        Hourly demand profiles for buses of base network
    pypsa_network: netcdf file
        base.nc network
    scaling_factor: pandas dataframe
        Hourly scaling factor per each state
    Returns
    -------
    scaled_demand_profiles: pandas dataframe
         Scaled demand profiles based on the demand projections
    """
    # read gadm file
    gadm_shape = gpd.read_file(snakemake.input.gadm_shape)

    # create geodataframe out of x and y coordinates of buses
    buses_gdf = gpd.GeoDataFrame(
        pypsa_network.buses,
        geometry=gpd.points_from_xy(pypsa_network.buses.x, pypsa_network.buses.y),
        crs=snakemake.params.geo_crs,
    ).reset_index()

    gadm_centroid = gadm_shape.to_crs(snakemake.params.geo_crs)
    gadm_centroid.geometry = gadm_centroid.geometry.centroid

    # map gadm shapes to each bus
    spatial_gadm_bus_mapping = (
        buses_gdf.sjoin_nearest(gadm_centroid, how="left")
        .set_index("Bus")["ISO_1"]
        .str.replace("US-", "")
    )

    # breakpoint()

    # convert demand_profiles from wide to long format
    df_demand_long = df_demand_profiles.melt(
        ignore_index=False, var_name="Bus", value_name="demand"
    )

    # map Bus IDs to State Codes
    df_demand_long["region_code"] = df_demand_long["Bus"].map(spatial_gadm_bus_mapping)

    # merge with scaling_factor DataFrame based on region_code and time
    scaling_factor["time"] = scaling_factor["time"].apply(
        lambda t: t.replace(year=df_demand_long.index[0].year)
    )
    # breakpoint()
    df_demand_long = df_demand_long.reset_index().rename(columns={"snapshot": "time"})
    df_scaled = df_demand_long.merge(
        scaling_factor, on=["region_code", "time"], how="left"
    )
    del scaling_factor

    # multiply demand by scaling factor
    df_scaled["scaling_factor"] = df_scaled["scaling_factor"].fillna(1)
    df_scaled["scaled_demand"] = df_scaled["demand"] * df_scaled["scaling_factor"]

    # pivot back to original wide format
    scaled_demand_profiles = df_scaled.pivot(
        index="time", columns="Bus", values="scaled_demand"
    )
    scaled_demand_profiles = scaled_demand_profiles[
        sorted(scaled_demand_profiles.columns)
    ]
    # logger.info(f"Scaled demand based on scaling factor for each state.")

    return scaled_demand_profiles


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("build_demand_profiles_from_eia")

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "demand_modelling")
    plot_path = pathlib.Path(default_path, "analysis", "plots", "demand_modelling")
    output_demand_profile_path = pathlib.Path(
        default_path, snakemake.output.demand_profile_path
    )
    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(plot_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(output_demand_profile_path).parent.mkdir(parents=True, exist_ok=True)

    geo_crs = snakemake.params.geo_crs
    distance_crs = snakemake.params.distance_crs
    demand_scenario = snakemake.params.demand_scenario
    demand_horizon = snakemake.params.demand_horizon

    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(
        log_path, f"output_build_demand_profiles_{today_date[:10]}.txt"
    )
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    (df_ba_demand, gdf_ba_shape, df_utility_demand, pypsa_network) = parse_inputs(
        default_path, distance_crs
    )

    df_demand_profiles = build_demand_profiles(
        df_utility_demand,
        df_ba_demand,
        gdf_ba_shape,
        pypsa_network,
        geo_crs,
        distance_crs,
    )

    if demand_horizon > 2020:
        scaling_factor = read_scaling_factor(
            demand_scenario, demand_horizon, default_path
        )
        df_demand_profiles = scale_demand_profiles(
            df_demand_profiles, pypsa_network, scaling_factor
        )

    df_demand_profiles.to_csv(output_demand_profile_path)

    log_output_file.close()
