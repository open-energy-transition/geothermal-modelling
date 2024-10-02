# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import argparse
import pathlib
import pandas as pd
import pypsa
import plotly.express as px
import _helpers_usa as h

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

# Load data
network_pypsa_earth = pypsa.Network(network_pypsa_earth_path)
eia_generation_reference = pd.read_excel(eia_generation_reference_path, sheet_name="Data")

# Parse argument
year_list = [val for val in list(eia_generation_reference.columns) if val not in ["Data_Status", "State", "MSN"]]
parser = argparse.ArgumentParser()
parser.add_argument("--year", help="Year to consider for the comparison", default=2020, type=int, choices=year_list)
args = parser.parse_args()

# EIA select relevant year and relevant rows
keywords_list = ["WYTCP", "NUETP", "HYTCP", "GEEGP"]
eia_generation_reference_filtered = eia_generation_reference[eia_generation_reference["MSN"].isin(keywords_list)][["State", "MSN", args.year]]
eia_generation_reference_filtered = eia_generation_reference_filtered.rename(columns={'MSN': 'carrier'})
eia_generation_reference_filtered.carrier = eia_generation_reference_filtered.carrier.map({"WYTCP": "onwind",
                                                                                    "NUETP": "nuclear",
                                                                                    "HYTCP": "hydro",
                                                                                    "GEEGP": "geothermal"})
eia_generation_reference_filtered.set_index(['carrier', 'State'] ,inplace=True)
eia_generation_reference_filtered[args.year] /= 1e3 #conversion from million kWh i.e., GWh to TWh

# Calculating the annual generation from each of the generators in MWh
time_res = 24
network_pypsa_earth.generators = network_pypsa_earth.generators.assign(gen_MWh=network_pypsa_earth.generators_t.p.sum() * time_res)

# PyPSA-Earth mapping the generators bus to the state names
usa_state_dict = h.get_gadm_mapping()
network_pypsa_earth.generators.bus = network_pypsa_earth.generators.bus.str.replace("_AC","")
network_pypsa_earth.generators.bus = network_pypsa_earth.generators.bus.str.replace("_DC","")
network_pypsa_earth.generators.bus = network_pypsa_earth.generators.bus.str.replace(" csp","")
network_pypsa_earth.generators.bus = network_pypsa_earth.generators.bus.map(usa_state_dict)

# Selecting relevant data from PyPSA results
pypsa_generation_filtered = network_pypsa_earth.generators.reset_index()[['carrier', 'bus', 'gen_MWh']]
pypsa_generation_filtered = pypsa_generation_filtered.rename(columns={'bus': 'State'})
pypsa_generation_filtered.set_index(['carrier', 'State'], inplace=True)
pypsa_generation_filtered['gen_MWh'] /= 1e6

# Joining EIA and PyPSA generations to compare them
comparison_generation = pypsa_generation_filtered.copy()
comparison_generation = comparison_generation.join(eia_generation_reference_filtered)
comparison_generation = comparison_generation.rename(columns={'gen_MWh': pypsa_name, args.year: eia_name})

# iterate through each carrier, filter and plot the comparison between PyPSA and EIA generation statewise
for car in ['nuclear', 'onwind', 'hydro', 'geothermal']:
    comparison_generation_filter = comparison_generation.reset_index().query('carrier == @car')
    comparison_generation_filter = comparison_generation_filter.set_index('State')[[pypsa_name, eia_name]]
    fig = px.bar(comparison_generation_filter,barmode='group')
    fig.update_layout(title=f"Carrier = {car}",yaxis_title='Generation (TWh)')
    fig.show()

# Compute error %s of the generation
comparison_generation['error %'] = comparison_generation.apply(lambda x: (x[eia_name]-x[pypsa_name])*100/x[eia_name], axis=1)