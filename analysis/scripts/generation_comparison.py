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
network_pypsa_earth_path = pathlib.Path(pypsa_earth_path, "results", "US_2021", "networks", "elec_s_10_ec_lcopt_Co2L-25H.nc")
eia_generation_reference_path = pathlib.Path(base_path.parent, "geothermal-modelling", "analysis", "data", "EIA_statewise_data", "use_all_phy.xlsx")
gadm_shapes_path = pathlib.Path(base_path, "analysis", "data", "gadm41_USA_1.json")

# Load data
network_pypsa_earth = pypsa.Network(network_pypsa_earth_path)
eia_generation_reference = pd.read_excel(eia_generation_reference_path, sheet_name="Data")
gadm_shapes = gpd.read_file(gadm_shapes_path)

# Parse argument
#year_list = eia_generation_reference["Period"].unique()
parser = argparse.ArgumentParser()
#parser.add_argument("--year", help="Year to consider for the comparison", default=2020, type=int, hoices=year_list)
parser.add_argument("--year", help="Year to consider for the comparison", default=2020, type=int)
args = parser.parse_args()

keywords_list = ["WYTCP", "NUETP", "HYTCP", "GEEGP"]
eia_generation_reference_filtered = eia_generation_reference[eia_generation_reference["MSN"].isin(keywords_list)][["State", "MSN", args.year]]
network_pypsa_earth
#eia_generation_reference_filtered.to_csv('test.csv')




