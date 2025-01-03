import geopandas as gpd
import numpy as np
import pandas as pd
import re


def get_gadm_mapping(gadm_shapes_path):
    gadm_gdp_usa = gpd.read_file(gadm_shapes_path)
    gadm_gdp_usa_state = pd.DataFrame()
    gadm_gdp_usa_state["GID_1_new"] = gadm_gdp_usa["GID_1"].str.replace("USA", "US")
    gadm_gdp_usa_state["state"] = gadm_gdp_usa["ISO_1"].str[-2:]
    usa_state_dict = dict(gadm_gdp_usa_state.values)
    return usa_state_dict


def extract_time_res(network_path):

    # Search for the pattern in the filename
    match = re.search(r'(\d+)(H|SEG)', str(network_path.name))

    if match:
        number = np.float64(match.group(1))
        unit = match.group(2)

        if unit == 'H':
            return number
        elif unit == 'SEG':
            return 8760 / number
    else:
        return None
