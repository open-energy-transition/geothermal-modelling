import pathlib
import geopandas as gpd

def get_gadm_mapping():
    base_path = pathlib.Path(__file__).parent.parent.parent
    gadm_shapes_path = pathlib.Path(base_path, "analysis", "data", "gadm41_USA_1.json")
    gadm_gdp_usa = gpd.read_file(gadm_shapes_path)
    gadm_gdp_usa_state = gadm_gdp_usa[["GID_1", "ISO_1"]]
    gadm_gdp_usa_state["state"] = gadm_gdp_usa_state["ISO_1"].str[-2:]
    gadm_gdp_usa_state["GID_1_new"] = gadm_gdp_usa_state["GID_1"].str.replace("USA", "US")
    gadm_gdp_usa_state = gadm_gdp_usa_state[["GID_1_new", "state"]]
    usa_state_dict = dict(gadm_gdp_usa_state.values)

    return usa_state_dict