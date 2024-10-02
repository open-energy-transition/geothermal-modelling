# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import argparse
import pathlib
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import pypsa
import cartopy.crs as ccrs

# Initial configurations
eia_name = "EIA"
pypsa_name = "PyPSA"

# Set the paths
base_path = pathlib.Path(__file__).parent.parent.parent
log_file_dir_path = pathlib.Path(base_path, "logs")
plot_dir_path = pathlib.Path(base_path, "analysis", "plots")
pypsa_earth_path = pathlib.Path(base_path, "workflow", "pypsa-earth")
#network_pypsa_earth_path = pathlib.Path(pypsa_earth_path, "results", "US_2021", "networks", "elec_s_10_ec_lcopt_Co2L-25H.nc")
network_pypsa_earth_path = pathlib.Path(pypsa_earth_path, "results", "US_2021", "networks", "elec_s_50_ec_lcopt_Co2L-24H.nc") # just to test it
eia_generation_reference_path = pathlib.Path(base_path.parent, "geothermal-modelling", "analysis", "data", "EIA_statewise_data", "use_all_phy.xlsx")
gadm_shapes_path = pathlib.Path(base_path, "analysis", "data", "gadm41_USA_1.json")

# Load data
network_pypsa_earth = pypsa.Network(network_pypsa_earth_path)
eia_generation_reference = pd.read_excel(eia_generation_reference_path, sheet_name="Data")
gadm_gdp_usa = gpd.read_file(gadm_shapes_path)

# Parse argument
year_list = [val for val in list(eia_generation_reference.columns) if val not in ["Data_Status", "State", "MSN"]]
parser = argparse.ArgumentParser()
parser.add_argument("--year", help="Year to consider for the comparison", default=2020, type=int, choices=year_list)
args = parser.parse_args()

# EIA select relevant year and relevant rows
keywords_list = ["WYTCP", "NUETP", "HYTCP", "GEEGP"]
eia_generation_reference_filtered = eia_generation_reference[eia_generation_reference["MSN"].isin(keywords_list)][["State", "MSN", args.year]]

# PyPSA-Earth mapping the generators bus to the state names
gadm_gdp_usa_state = gadm_gdp_usa[["GID_1", "ISO_1"]]
gadm_gdp_usa_state["state"] = gadm_gdp_usa_state["ISO_1"].str[-2:]
gadm_gdp_usa_state["GID_1_new"] = gadm_gdp_usa_state["GID_1"].str.replace("USA", "US")
gadm_gdp_usa_state = gadm_gdp_usa_state[["GID_1_new", "state"]]
usa_state_dict = dict(gadm_gdp_usa_state.values)
network_pypsa_earth.generators.index = network_pypsa_earth.generators.index.map(usa_state_dict)
network_pypsa_earth.generators.to_csv("gen.csv")





