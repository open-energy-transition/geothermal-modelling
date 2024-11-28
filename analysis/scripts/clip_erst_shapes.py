import geopandas as gpd
from shapely.geometry import LineString, Point
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import pathlib
import pandas as pd
from shapely.validation import make_valid

if __name__ ==  '__main__':

    # Path 
    base_path = pathlib.Path(__file__).parent.parent.parent
    erst_path = pathlib.Path(base_path, "analysis", "gdrive_data", "data", "electricity_demand_data", "demand_data", "Electric_Retail_Service_Territories.geojson")
    gadm_usa_path = pathlib.Path(base_path, "analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json")
    demand_utility_path = pathlib.Path(base_path, "analysis", "gdrive_data","data", "electricity_demand_data", "demand_data","table_10_EIA_utility_sales.xlsx")
    erst_final_path = pathlib.Path(base_path, "analysis", "gdrive_data", "data", "electricity_demand_data", "demand_data", "ERST_overlay_demand.geojson")
    # Load data
    df_erst_gpd = gpd.read_file(erst_path)
    df_gadm_usa = gpd.read_file(gadm_usa_path)
    df_demand_utility = pd.read_excel(demand_utility_path,skiprows=2)

    # ERST overlay with GADM shapes for intersection
    df_erst_overlay = gpd.overlay(df_erst_gpd, df_gadm_usa, how='intersection')
    df_erst_overlay['STATE'] = df_erst_overlay.apply(lambda x: x['ISO_1'].split('-')[1], axis=1)
    df_erst_overlay.set_index(['NAME','STATE'], inplace=True)

    # Demand mapping to ERST shape file
    df_demand_utility['Entity'] = df_demand_utility['Entity'].apply(str.upper)
    df_demand_utility.rename(columns={'Entity':'NAME','State':'STATE'}, inplace=True)
    df_demand_utility.set_index(['NAME','STATE'], inplace=True)

    df_erst_overlay = df_erst_overlay.join(df_demand_utility)
    
    # Make geometry valid
    df_erst_overlay['geometry'] = df_erst_overlay['geometry'].apply(make_valid)

    df_erst_overlay.to_file(erst_final_path, driver='GeoJSON')

    
    # demand_dict = dict(df_demand_utility[['Entity','Sales (Megawatthours)']].values)

    # df_erst_gpd = df_erst_gpd.to_crs(3857)

    # split_threshold = 20 # Length of shapes in km beyond which the shape is split into two parts

    # df_filter = df_erst_gpd.loc[df_erst_gpd['Shape__Length'] > split_threshold]
    # df_eg = df_filter.iloc[0:1]

    # # xmin, ymin, xmax, ymax = df_eg.total_bounds
    # # length = (xmax - xmin) / 1e3
    # # breadth = (ymax - ymin) / 1e3
    # length = df_eg['Shape__Length'].iloc[0]

    # if length > split_threshold:
    #     splitx = -((xmin - xmax) / 2 - xmin)
    #     split_line = LineString([Point(splitx, ymin),Point(splitx, ymax)])
    #     lhs = split_line.buffer(splitx, single_sided = True)
    #     rhs = split_line.buffer(-splitx, single_sided = True)
    #     s1 = df_eg.difference(lhs)
    #     s2 = df_eg.difference(rhs)

    # if breadth > 20:
    #     splity = -((ymin - ymax) / 2 - ymin)
    #     split_line = LineString([Point(xmin, splity),Point(xmax, splity)])
    #     uhs = split_line.buffer(splity, single_sided = True)
    #     dhs = split_line.buffer(-splity, single_sided = True)
    #     s3 = df_eg.difference(uhs)
    #     s4 = df_eg.difference(dhs)

    # fig, ax = plt.subplots()
    # df_eg.plot(ax=ax)
    # s1.plot(ax=ax, color='red', alpha=0.5, hatch="|")
    # s2.plot(ax=ax, color='yellow', alpha=0.5, hatch="-")
    # s3.plot(ax=ax, color='green', alpha=0.5, hatch="\\")
    # s4.plot(ax=ax, color='orange', alpha=0.5, hatch="/")


