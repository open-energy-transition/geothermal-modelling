import pathlib
import datetime as dt
import pandas as pd
import geopandas as gpd
from _helpers_usa import get_colors
import numpy as np

def parse_inputs(default_path):
    BA_demand_path1 = pathlib.Path(default_path, snakemake.input.BA_demand_path1)
    BA_demand_path2 = pathlib.Path(default_path, snakemake.input.BA_demand_path2)
    df_ba_demand1 = pd.read_csv(BA_demand_path1, index_col="period")
    df_ba_demand2 = pd.read_csv(BA_demand_path2, index_col="period")
    df_ba_demand = df_ba_demand1._append(df_ba_demand2)
    df_ba_demand = df_ba_demand.replace(0, np.nan)
    df_ba_demand = df_ba_demand.dropna(axis=1)

    balancing_authority_shapefile = pathlib.Path(default_path, snakemake.input.BA_shape_path)
    gdf_ba_shape = gpd.read_file(balancing_authority_shapefile)
    gdf_ba_shape = gdf_ba_shape.to_crs(3857)

    utility_demand_path = pathlib.Path(default_path, snakemake.input.Utiltiy_demand_path)
    df_utility_demand = gpd.read_file(utility_demand_path)
    df_utility_demand.rename(columns={"index_right": "index_right_1"}, inplace=True)

    df_utility_demand = df_utility_demand.to_crs(3857)

    return df_ba_demand, gdf_ba_shape, df_utility_demand

def build_demand_profiles(df_utility_demand, df_ba_demand, gdf_ba_shape):
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
    df_utility_centroid.rename(columns={"index_right": "index_right_2"}, inplace=True)

    # temporal scaling factor
    df_utility_centroid["temp_scale"] = df_utility_centroid.apply(
        lambda x: df_ba_demand[f"E_{x['EIAcode']}_D"].sum()
        / 1e3
        / x["Sales (Megawatthours)"],
        axis=1,
    )

    # Mapping demand utilities to nearest PyPSA bus
    df_reqd = df_n.buses.query('carrier == "AC"')
    pypsa_gpd = gpd.GeoDataFrame(
        df_reqd, geometry=gpd.points_from_xy(df_reqd.x, df_reqd.y), crs=4326
    )
    pypsa_gpd = pypsa_gpd.to_crs(3857)
    pypsa_gpd["color"] = get_colors(len(pypsa_gpd))

    df_utility_centroid = gpd.sjoin_nearest(df_utility_centroid, pypsa_gpd, how="left")
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

    return df_demand_bus


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("build_demand_profiles_from_eia")

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "demand_modelling")
    plot_path = pathlib.Path(default_path, "analysis", "plots", "demand_modelling")
    output_path = pathlib.Path(default_path, snakemake.output.demand_profile_path)
    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(plot_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)

    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(
        log_path, f"output_demand_modelling_{today_date[:10]}.txt"
    )
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    (
        df_ba_demand, 
        gdf_ba_shape, 
        df_utility_demand 

    ) = parse_inputs(default_path)

    df_demand_profiles = build_demand_profiles(df_utility_demand, df_ba_demand, gdf_ba_shape)

    df_demand_profiles.to_csv(output_path)

    log_output_file.close()
