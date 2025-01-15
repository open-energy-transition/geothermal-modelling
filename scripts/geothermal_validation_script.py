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

import pathlib
import os
import subprocess
import datetime as dt

# Initial configurations

# Set the paths
base_path = pathlib.Path(__file__).parent
log_file_dir = base_path / "logs"

# Ensure the logs directory exists
pathlib.Path(log_file_dir).mkdir(exist_ok=True)

# Define paths
pypsa_earth_path = pathlib.Path(base_path, "workflow", "pypsa-earth")
network_path = pathlib.Path(
    pypsa_earth_path, "networks", "US_2021", "elec_s_10_ec_lcopt_Co2L-24H.nc"
)
config_path_geothermal = pathlib.Path(base_path, "Config", "config.usa_PE.yaml")
config_path_pypsa_earth = pathlib.Path(pypsa_earth_path, "config.yaml")
custom_power_plants_geothermal = pathlib.Path(
    base_path, "data", "custom_powerplants.csv"
)
custom_power_plants_pypsa_earth = pathlib.Path(
    pypsa_earth_path, "data", "custom_powerplants.csv"
)
north_america_geothermal = pathlib.Path(base_path, "data", "NorthAmerica.csv")
north_america_pypsa_earth = pathlib.Path(
    pypsa_earth_path, "data", "ssp2-2.6", "2030", "era5_2013", "NorthAmerica.csv"
)
costs_geothermal = pathlib.Path(base_path, "data", "costs.csv")
costs_pypsa_earth = pathlib.Path(pypsa_earth_path, "data", "costs.csv")

# Open log_output
today_date = str(dt.datetime.now())
log_output_file = open(log_file_dir / f"output_pypsa_earth_{today_date[:10]}.txt", "w")

log_output_file.write("        \n")
log_output_file.write("        \n")
log_output_file.write("Copy config file to PyPSA-Earth \n")
print("Copy config file to PyPSA-Earth \n")
subprocess.run(["cp", config_path_geothermal, config_path_pypsa_earth])

log_output_file.write("        \n")
log_output_file.write("        \n")
log_output_file.write("Copy custom_powerplants.csv to PyPSA-Earth \n")
print("Copy custom_powerplants.csv to PyPSA-Earth \n")
subprocess.run(["cp", custom_power_plants_geothermal, custom_power_plants_pypsa_earth])

os.chdir(pypsa_earth_path)

log_output_file.write(f"Switch to PyPSA-Earth submodule folder {pathlib.Path.cwd()} \n")
print(f"Switch to PyPSA-Earth submodule folder {pathlib.Path.cwd()} \n")

subprocess.run(["snakemake", "--unlock"])

# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Execute build_powerplants \n")
# print("Execute build_powerplants \n")
# subprocess.run(["snakemake", "-call", "build_powerplants", "--cores", "all", "--printshellcmds", "--configfile", "config.yaml"])
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Copy NorthAmerica.csv to PyPSA-Earth \n")
# print("Copy NorthAmerica.csv to PyPSA-Earth \n")
# subprocess.run(["cp", north_america_geothermal, north_america_pypsa_earth])
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Execute build_demand_profile \n")
# print("Execute build_demand_profile \n")
# subprocess.run(["snakemake", "-call", "resources/US_2021/demand_profiles.csv", "--cores", "all", "--printshellcmds", "--configfile", "config.yaml"])
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Build renewable profile for onwind \n")
# print("Execute add_electricity \n")
# subprocess.run(["snakemake", "-call", "resources/US_2021/renewable_profiles/profile_onwind.nc", "--cores", "all", "--printshellcmds", "--configfile", "config.yaml"])
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Build renewable profile for offwind-ac \n")
# print("Execute add_electricity \n")
# subprocess.run(["snakemake", "-call", "resources/US_2021/renewable_profiles/profile_offwind-ac.nc", "--cores", "all", "--printshellcmds", "--configfile", "config.yaml"])
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Build renewable profile for offwind-dc \n")
# print("Execute add_electricity \n")
# subprocess.run(["snakemake", "-call", "resources/US_2021/renewable_profiles/profile_offwind-dc.nc", "--cores", "all", "--printshellcmds", "--configfile", "config.yaml"])
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Build renewable profile for hydro \n")
# print("Execute add_electricity \n")
# subprocess.run(["snakemake", "-call", "resources/US_2021/renewable_profiles/profile_hydro.nc", "--cores", "all", "--printshellcmds", "--configfile", "config.yaml"])
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Build renewable profile for solar \n")
# print("Execute add_electricity \n")
# subprocess.run(["snakemake", "-call", "resources/US_2021/renewable_profiles/profile_solar.nc", "--cores", "all", "--printshellcmds", "--configfile", "config.yaml"])
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Build renewable profile for csp \n")
# print("Execute add_electricity \n")
# subprocess.run(["snakemake", "-call", "resources/US_2021/renewable_profiles/profile_csp.nc", "--cores", "all", "--printshellcmds", "--configfile", "config.yaml"])
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Copy costs.csv to PyPSA-Earth \n")
# print("Copy costs.csv to PyPSA-Earth \n")
# subprocess.run(["cp", costs_geothermal, costs_pypsa_earth])
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Execute add_electricity \n")
# print("Execute add_electricity \n")
# subprocess.run(["snakemake", "-call", "networks/US_2021/elec_s_10.nc" , "--cores", "all", "--printshellcmds", "--configfile", "config.yaml"])
#
# # log_output_file.write("        \n")
# # log_output_file.write("        \n")
# # log_output_file.write("Set extendable carriers to False \n")
# # print("Set extendable carriers to False \n")
#
# # network_to_modify = pypsa.Network(network_path)
# # network_to_modify.generators.loc[:, "p_nom_extendable"] = False
# # network_to_modify.export_to_netcdf(network_path)
#
# log_output_file.write("        \n")
# log_output_file.write("        \n")
# log_output_file.write("Execute solve_all_networks \n")
# print("Execute solve_all_networks \n")
# subprocess.run(["snakemake", "-call", "solve_all_networks", "--cores", "all", "--printshellcmds", "--configfile", "config.yaml"])

log_output_file.close()
