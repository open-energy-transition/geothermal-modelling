import shapely
import geopandas as gpd
import pandas as pd
import numpy as np
import pypsa
from shapely import Point

gdf = gpd.read_file("../data/gadm41_USA_shp/gadm41_USA_1.shp")

df_ppl = pd.read_csv("../../Data/custom_powerplants.csv")

def map_geo(gdf,point,row):
    for j in range(0,len(gdf)):
        reg = gdf.iloc[j]
        geo = reg["geometry"]
        if geo.contains(point):
            region_name = reg["NAME_1"]
            df_ppl.loc[i,"busmap"] = region_name
    print(i)
            
# df = pd.DataFrame(index=df_ppl.Name.tolist())
df_ppl["busmap"] = ""
for i in range(0,len(df_ppl)):
    print(i)
    row = df_ppl.iloc[i]
    point = Point(row["lon"],row["lat"])
    map_geo(gdf,point,row)

gdf['US ID'] = gdf.apply(lambda x: x['ISO_1'].split('-')[1],axis=1)
state_dict = dict(gdf[['NAME_1','US ID']].values)
df_ppl['state'] = df_ppl['busmap'].map(state_dict)
df_ppl.to_csv("../data/custom_powerplants_with_state.csv")