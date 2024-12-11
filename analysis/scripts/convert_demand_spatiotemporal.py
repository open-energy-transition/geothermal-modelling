import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import pathlib
import numpy as np
import pypsa

if __name__ ==  '__main__':

    # Path 
    base_path = pathlib.Path(__file__).parent.parent.parent
    BA_demand_path = pathlib.Path("/Volumes","Sermisha_harddisk","01_OET","14_geothermal_modelling","processed","EIA930_2024_Jan_Jun_opt.csv")
    BA_shape_path = pathlib.Path("/Volumes","Sermisha_harddisk","01_OET","14_geothermal_modelling","processed","Balancing_Authorities.geojson")
    Utility_demand_path = pathlib.Path("/Volumes/Sermisha_harddisk/01_OET/14_geothermal_modelling/Outputs/Demand_mapped_v2.geojson")

    # Load data
    df_ba_demand = pd.read_csv(BA_demand_path, index_col='period')
    gdf_ba_shape = gpd.read_file(BA_shape_path)
    df_utility_demand = gpd.read_file(Utility_demand_path)
    df_utility_demand.rename(columns={'index_right':'index_r'},inplace=True)
    df_n = pypsa.Network("/Volumes/Sermisha_harddisk/01_OET/14_geothermal_modelling/geothermal-modelling/workflow/pypsa-earth/networks/US_2021/elec_s_100flex_ec_lcopt_Co2L-1H.nc")

    # Obtaining the centroids of the Utility demands
    df_utility_centroid = df_utility_demand.copy()
    df_utility_centroid.geometry = df_utility_centroid.geometry.centroid

    # Removing those shapes which do not have a representation in the demand timeseries data
    demand_columns = df_ba_demand.columns
    demand_columns = [x.split("_")[1] for x in demand_columns]
    drop_shape_rows = list(set(gdf_ba_shape['EIAcode'].tolist()) - set(demand_columns))
    gdf_ba_shape_filtered = gdf_ba_shape[~gdf_ba_shape['EIAcode'].isin(drop_shape_rows)]

    # Convert CRS of shape files
    df_utility_demand = df_utility_demand.to_crs(3857)
    gdf_ba_shape_filtered = gdf_ba_shape_filtered.to_crs(3857)
    gdf_ba_shape  = gdf_ba_shape.to_crs(3857)
    df_utility_demand = gpd.sjoin_nearest(df_utility_demand, gdf_ba_shape_filtered, how='inner')
    
    # temporal scaling factor
    df_utility_demand['temp_scale'] = df_utility_demand.apply(lambda x: df_ba_demand[f"E_{x['EIAcode']}_D"].sum() / 1e3 / x['Sales (Megawatthours)'], axis=1)
    # # PyPSA AC bus data
    # pypsa_gpd = gpd.GeoDataFrame(df_reqd, geometry=gpd.points_from_xy(df_reqd.x,df_reqd.y))

    # for i in np.arange(0,len(gdf_ba_shape)):
    #     ba_row = gdf_ba_shape.iloc[i]
    #     ba = ba_row['EIAcode']
    #     print(ba)
    #     try:
    #         print(df_ba_demand[f"E_{ba}_D"])
    #     except:
    #         print("BA not found")
    #     input()