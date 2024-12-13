import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import pathlib
import numpy as np
import pypsa
import random

if __name__ ==  '__main__':

    # Path 
    base_path = pathlib.Path(__file__).parent.parent.parent
    BA_demand_path1 = pathlib.Path("/Volumes","Sermisha_harddisk","01_OET","14_geothermal_modelling","processed","EIA930_2021_Jan_Jun_opt.csv")
    BA_demand_path2 = pathlib.Path("/Volumes","Sermisha_harddisk","01_OET","14_geothermal_modelling","processed","EIA930_2021_Jul_Dec_opt.csv")
    BA_shape_path = pathlib.Path("/Volumes","Sermisha_harddisk","01_OET","14_geothermal_modelling","processed","Balancing_Authorities.geojson")
    Utility_demand_path = pathlib.Path("/Volumes/Sermisha_harddisk/01_OET/14_geothermal_modelling/Outputs/Demand_mapped_v2.geojson")

    # Load data
    df_ba_demand1 = pd.read_csv(BA_demand_path1, index_col='period')
    df_ba_demand2 = pd.read_csv(BA_demand_path2, index_col='period')
    df_ba_demand = df_ba_demand1._append(df_ba_demand2)
    df_ba_demand = df_ba_demand.replace(0,np.nan)
    df_ba_demand = df_ba_demand.dropna(axis=1)
    gdf_ba_shape = gpd.read_file(BA_shape_path)
    df_utility_demand = gpd.read_file(Utility_demand_path)
    df_utility_demand.rename(columns={'index_right':'index_right_1'},inplace=True)
    df_n = pypsa.Network("/Volumes/Sermisha_harddisk/01_OET/14_geothermal_modelling/geothermal-modelling/workflow/pypsa-earth/networks/US_2021/elec_s_100flex_ec_lcopt_Co2L-1H.nc")

    # Convert CRS of shape files
    df_utility_demand = df_utility_demand.to_crs(3857)
    gdf_ba_shape  = gdf_ba_shape.to_crs(3857)

    # Obtaining the centroids of the Utility demands
    df_utility_centroid = df_utility_demand.copy()
    df_utility_centroid.geometry = df_utility_centroid.geometry.centroid

    get_colors = lambda n: ["#%06x" % random.randint(0, 0xFFFFFF) for _ in range(n)]

    # Removing those shapes which do not have a representation in the demand timeseries data
    demand_columns = df_ba_demand.columns
    demand_columns = [x.split("_")[1] for x in demand_columns if x.endswith("_D")]
    drop_shape_rows = list(set(gdf_ba_shape['EIAcode'].tolist()) - set(demand_columns))
    gdf_ba_shape_filtered = gdf_ba_shape[~gdf_ba_shape['EIAcode'].isin(drop_shape_rows)]
    gdf_ba_shape_filtered.geometry = gdf_ba_shape_filtered.geometry.centroid
    gdf_ba_shape_filtered['color_ba'] = get_colors(len(gdf_ba_shape_filtered))

    df_utility_centroid = gpd.sjoin_nearest(df_utility_centroid, gdf_ba_shape_filtered, how='left')
    df_utility_centroid.rename(columns={'index_right':'index_right_2'},inplace=True)

    # temporal scaling factor
    df_utility_centroid['temp_scale'] = df_utility_centroid.apply(lambda x: df_ba_demand[f"E_{x['EIAcode']}_D"].sum() / 1e3 / x['Sales (Megawatthours)'], axis=1)
    
    # Mapping demand utilities to nearest PyPSA bus
    df_reqd = df_n.buses.query('carrier == "AC"')
    pypsa_gpd = gpd.GeoDataFrame(df_reqd, geometry=gpd.points_from_xy(df_reqd.x,df_reqd.y), crs=4326)
    pypsa_gpd = pypsa_gpd.to_crs(3857)
    pypsa_gpd['color'] = get_colors(len(pypsa_gpd))

    df_utility_centroid = gpd.sjoin_nearest(df_utility_centroid, pypsa_gpd, how='left')
    df_utility_centroid.rename(columns={'index_right':'PyPSA_bus'},inplace=True)

    df_demand_bus = pd.DataFrame(index=pd.to_datetime(df_ba_demand.index), columns=df_reqd.index.tolist(), data=0.0)

    for col in df_demand_bus.columns:
        utility_rows = df_utility_centroid.query("PyPSA_bus == @col")
        for i in np.arange(0,len(utility_rows)):
            row = utility_rows.iloc[i]
            if not np.isnan(row['temp_scale']):
                demand_data = df_ba_demand[f"E_{row['EIAcode']}_D"].copy()
                demand_data /= row['temp_scale'] 
                demand_data /= 1000 #Converting to MWh
                df_demand_bus[col] += demand_data.tolist()
