# coding=utf-8

##################################################
#                                                #
# Author: Fabrizio Finozzi	                     #
# Email: fabrizio.finozzi.business@gmail.com     #
# Version: 0.1                                   #
# Date: 23.07.2024                               #
#                                                #
# Created for Open Energy Transition GmbH        #
#                                                #
##################################################

#########################################################################
#                                                                       #
# IMPORTANT                                                             #
#                                                                       #
# This software is distributed without any warranty.                    #
#                                                                       #
# Neither the author nor Open Energy Transition GmbH                    #
# are liable for any damage caused directly or indirectly by the use    #
# or misuse of this software.                                           #
#                                                                       #
#########################################################################

#############################################################
#                                                           #
# CHANGELOG                                                 #
#                                                           #
# 0.1 - In progress                                         #
#                                                           #
#############################################################

import argparse
import datetime as dt
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import pypsa
import re
from pathlib import Path
import numpy as np

# Initial configurations
eia_name = "EIA"
pypsa_name = "PyPSA"

# Set the paths
base_path = pathlib.Path(__file__).parent.parent.parent
log_file_dir_path = pathlib.Path(base_path, "logs")
plot_dir_path = pathlib.Path(base_path, "analysis", "plots")
pypsa_earth_path = pathlib.Path(base_path, "workflow", "pypsa-earth")
network_path = pathlib.Path(pypsa_earth_path, "results", "US_2021", "networks", "elec_s_10_ec_lcopt_Co2L-24H.nc")
eia_generation_path = pathlib.Path(base_path, "analysis", "data", "generation_eia.csv")
generation_plot_path = pathlib.Path(plot_dir_path, "electricity_generation.png")
eia_capacity_path = pathlib.Path(base_path, "analysis", "data", "capacities_eia.xlsx")
installed_capacity_plot_path = pathlib.Path(plot_dir_path, "installed_capacity.png")


def extract_time_res(filename):
    # Convert the filename to a string (although you already did this, no need to use network_path.name)
    filename = str(filename)

    # Search for the pattern in the filename
    match = re.search(r'(\d+)(H|SEG)', filename)

    if match:
        number = np.float64(match.group(1))
        unit = match.group(2)

        if unit == 'H':
            return number
        elif unit == 'SEG':
            return 8760 / number
    else:
        return None


# Example usage
filepath = Path(network_path)
time_res = extract_time_res(filepath)
print(time_res)

# Read reference data
df_eia_generation = pd.read_csv(eia_generation_path, index_col="Unnamed: 0")
df_eia_capacity = pd.read_excel(eia_capacity_path, skiprows=1, index_col="Energy Source")

# Parse argument
year_list = df_eia_generation["Period"].unique()
parser = argparse.ArgumentParser()
parser.add_argument("--year", help="Year to consider for the comparison", default=2020, type=int, choices=year_list)
args = parser.parse_args()

# Ensure the logs directory exists
pathlib.Path(log_file_dir_path).mkdir(exist_ok=True)

# Ensure the plots directory exists
pathlib.Path(plot_dir_path).mkdir(exist_ok=True)

# print(base_path)
# print(pypsa_earth_path)
# print(network_path)
# print(eia_generation_path)
# print(eia_capacity_path)

# Open log_output_file
today_date = str(dt.datetime.now())
log_output_file = open(log_file_dir_path / f"output_pypsa_earth_analysis_{today_date[:10]}.txt", "w")

log_output_file.write("        \n")
log_output_file.write("        \n")
log_output_file.write("Import data \n")
print("Import data \n")

df_network = pypsa.Network(network_path)
df_eia_generation_year = df_eia_generation.loc[df_eia_generation["Period"] == args.year].squeeze()

s_max_pu = df_network.lines["s_max_pu"].unique()[0]

log_output_file.write("        \n")
log_output_file.write("        \n")
log_output_file.write(f"s_max_pu {s_max_pu} \n")
print(f"s_max_pu {s_max_pu} \n")

######################
# Installed capacity #
######################

log_output_file.write("        \n")
log_output_file.write("        \n")
log_output_file.write("Compare the installed capacity \n")
print("Comparing the installed capacity... \n")

# ---> Prepare the PyPSA results

print(df_network.storage_units)

df_pypsa_hydro_phs_capacity = df_network.storage_units.groupby("carrier").p_nom_opt.sum()
df_pypsa_capacity = pd.concat([df_network.generators.groupby("carrier").p_nom_opt.sum(), df_pypsa_hydro_phs_capacity])
df_pypsa_capacity.loc["wind"] = df_pypsa_capacity.loc[["offwind-ac", "offwind-dc", "onwind"]].sum()
df_pypsa_capacity = df_pypsa_capacity.drop(["offwind-ac", "offwind-dc", "onwind"])
df_pypsa_capacity /= 1000
df_pypsa_capacity = df_pypsa_capacity.round(2)
df_pypsa_capacity.name = pypsa_name

print("Installed capacity:\n", df_pypsa_capacity)

# ---> Prepare the EIA reference data
df_eia_capacity.index = df_eia_capacity.index.str.lower()
df_eia_capacity.loc["other biomass"] = df_eia_capacity.loc[["other biomass", "wood and wood-derived fuels"]].sum()
df_eia_capacity = df_eia_capacity.rename(index={"hydroelectric conventional": "hydro", "hydroelectric pumped storage": "PHS" , "solar photovoltaic": "solar", "natural gas": "CCGT", "petroleum": "oil", "other biomass": "biomass"})
df_eia_capacity = df_eia_capacity.drop(["estimated total solar", "solar thermal", "wood and wood-derived fuels", "other energy sources", "total", "small scale photovoltaic", "estimated total photovoltaic", "other gases",])
df_eia_capacity = df_eia_capacity.iloc[:-1]
df_eia_capacity.loc["solar", "Generator Nameplate Capacity"] = df_eia_capacity.loc["solar", "Net Summer Capacity"]
df_eia_capacity = df_eia_capacity["Generator Nameplate Capacity"]
df_eia_capacity.name = eia_name
df_eia_capacity /= 1000

# ---> Prepare comparison dataframe
df_compare_capacity = pd.concat([df_pypsa_capacity, df_eia_capacity], axis=1)
#df_compare_capacity["error"] = df_compare_capacity.apply(lambda x: abs(x[pypsa_name] - x[eia_name])*100 / x[eia_name], axis=1)

# ---> Plot
df_compare_capacity.plot(kind="bar")
plt.xlabel("Carriers")
plt.ylabel("Installed Capacity (GW)")
plt.title(f"Year = {args.year}")
plt.grid(linestyle="--")
plt.subplots_adjust(bottom=0.3)
plt.savefig(installed_capacity_plot_path, dpi=800)

##########################
# Electricity generation #
##########################

log_output_file.write("        \n")
log_output_file.write("        \n")
log_output_file.write("Compare the electricity generation \n")
print("\nComparing the electricity generation... \n")

# ---> Prepare the PyPSA results
df_network.storage_units = df_network.storage_units.assign(p=df_network.storage_units_t.p.sum() * time_res)
df_network.generators = df_network.generators.assign(p=df_network.generators_t.p.sum() * time_res)
df_pypsa_generation = pd.concat([df_network.generators.groupby("carrier").p.sum(), df_network.storage_units.groupby("carrier").p.sum()])
df_pypsa_generation.loc["wind"] = df_pypsa_generation.loc[["offwind-ac", "offwind-dc", "onwind"]].sum()
df_pypsa_generation = df_pypsa_generation.drop(["offwind-ac", "offwind-dc", "onwind"])
df_pypsa_generation /= 1e6
df_pypsa_generation = df_pypsa_generation.round(2)
df_pypsa_generation.name = pypsa_name

print("Electricity generation:\n", df_pypsa_generation)

# ---> Prepare the EIA reference data
pypsa_cols = ["Coal", "Natural Gas", "Other Gas", "Nuclear", "Hydro", "Estimated Total Solar", "PHS", "Petroleum", "Wind", "Other Waste Biomass", "Geothermal"]
rename_cols = {"estimated total solar": "solar", "other waste biomass": "biomass", "natural gas": "CCGT", "petroleum": "oil", "phs": "PHS"}

df_eia_generation_year = df_eia_generation_year[pypsa_cols]
df_eia_generation_year.index = df_eia_generation_year.index.str.lower()
df_eia_generation_year = df_eia_generation_year.rename(index=rename_cols)
df_eia_generation_year.name = eia_name
df_eia_generation_year = df_eia_generation_year.drop("other gas")

## ---> Prepare comparison dataframe
df_compare_generation = pd.concat([df_pypsa_generation, df_eia_generation_year], axis=1)
#df_compare_generation["error"] = df_compare_generation.apply(lambda x: abs(x[pypsa_name] - x[eia_name])*100 / x[eia_name], axis=1)

# ---> Plot
df_compare_generation.plot(kind="bar")
plt.xlabel("Carriers")
plt.ylabel("Electricity Generation (TWh)")
plt.title(f"Year = {args.year}")
plt.grid(linestyle="--")
plt.subplots_adjust(bottom=0.3)
plt.savefig(generation_plot_path, dpi=800)

print("\nMarginal costs of electricity:")
marginal_costs = pd.concat([df_network.generators.groupby("carrier").marginal_cost.first().sort_values(), df_network.storage_units.groupby("carrier").marginal_cost.first().sort_values()])
print(marginal_costs.round(3))

print("\nTotal electricity generation (2020):")
print("EIA: ", df_eia_generation_year.sum().round(2), "TWh")
print("PyPSA: ", df_pypsa_generation.sum().round(2), "TWh \n")

log_output_file.close()
