import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib

# Defining file paths
base_path = "../"
geojson_path = pathlib.Path(base_path,"data","demand_data","Electric_Retail_Service_Territories.geojson")
utility_demand_data = pathlib.Path(base_path,"data","demand_data","table_10_EIA_utility_sales.xlsx")
# plots_path = pathlib.Path("..","plots","demand_plots")

# Loading the data
gdf_shapes = gpd.read_file(geojson_path)
df_utility_demand = pd.read_excel(utility_demand_data, skiprows=2)

df_utility_demand['Entity'] = df_utility_demand['Entity'].apply(str.upper)
demand_dict = dict(df_utility_demand[['Entity','Sales (Megawatthours)']].values)

gdf_shapes['Sales'] = gdf_shapes['NAME'].map(demand_dict)