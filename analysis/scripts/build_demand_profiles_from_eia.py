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


def parse_inputs(default_path, distance_crs):
    """
    Load all input data

    Parameters
    ----------
    default_path: str
        current total directory path
    distance_crs: CRS code
        Coordinate reference system to retrieve coordinate in 'm'

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
        Coordinate reference system to retrieve coordinates in 'geometrical degrees'
    distance_crs: CRS code
        Coordinate reference system to retrieve coordinate in 'm'

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

    df_demand_profiles.to_csv(output_demand_profile_path)

    log_output_file.close()
