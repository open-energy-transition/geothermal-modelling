import geopandas as gpd
import pathlib
import matplotlib.pyplot as plt
import os


base_path = "../../"
geojson_path = pathlib.Path(base_path,"workflow","pypsa-earth","resources","US_2021","shapes","gadm_shapes_nosimplify.geojson")
geojson_simplify_path = pathlib.Path(base_path,"workflow","pypsa-earth","resources","US_2021","shapes","gadm_shapes_simplify.geojson")
plots_path = pathlib.Path("..","plots","GADM_plots")

os.makedirs(plots_path,exist_ok=True)

file_simplify = open(geojson_simplify_path)
df_simplify = gpd.read_file(file_simplify)

file = open(geojson_path)
df = gpd.read_file(file)

m = df.explore()
outfp = f"{plots_path}/GADM_not_simplified.html"
m.save(outfp)    

m = df_simplify.explore()
outfp = f"{plots_path}/GADM_simplified.html"
m.save(outfp)

m = df.difference(df_simplify).explore()
outfp = f"{plots_path}/GADM_difference.html"
m.save(outfp)    