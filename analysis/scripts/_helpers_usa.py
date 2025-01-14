# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import geopandas as gpd
import numpy as np
import pandas as pd
import re
import yaml
import pathlib

def get_gadm_mapping(gadm_shapes_path):
    gadm_gdp_usa = gpd.read_file(gadm_shapes_path)
    gadm_gdp_usa_state = pd.DataFrame()
    gadm_gdp_usa_state["GID_1_new"] = gadm_gdp_usa["GID_1"].str.replace("USA", "US")
    gadm_gdp_usa_state["state"] = gadm_gdp_usa["ISO_1"].str[-2:]
    usa_state_dict = dict(gadm_gdp_usa_state.values)
    return usa_state_dict


def get_state_node(gadm_shapes_path, state):
    usa_state_dict = {v: k for k, v in get_gadm_mapping(gadm_shapes_path).items()}
    return usa_state_dict[state]


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

def config(config_path):
    path_config = pathlib.Path(config_path)
    with open(path_config) as file:
        config_dict = yaml.safe_load(file)
    return config_dict