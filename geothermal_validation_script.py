# coding=utf-8

##################################################
#                                                #
# Author: Fabrizio Finozzi	                     #
# Email: fabrizio.finozzi.business@gmail.com     #
#                                                #
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
import pypsa
from pathlib import Path

# Initial configurations

# --> set the paths
base_path=pathlib.Path.cwd()
filestolog_path = base_path / "logs"

# Ensure the logs directory exists
os.makedirs(filestolog_path, exist_ok=True)

pypsa_earth_path=pathlib.Path(base_path, "workflow", "pypsa-earth")
network_path=pathlib.Path(pypsa_earth_path, "networks", "US_2021", "elec_s_10_ec_lcopt_Co2L-24H.nc")

# --> open log_output
todaydate = str(dt.datetime.now())
logoutputfile = open(filestolog_path / f'output_pypsa_earth_{todaydate[:10]}.txt', 'w')

os.chdir(pypsa_earth_path)

# Fetch new metadata from the repositories
logoutputfile.write(f"Switch to PyPSA-Earth submodule folder {pathlib.Path.cwd()} \n")
print(f"Switch to PyPSA-Earth submodule folder {pathlib.Path.cwd()}")

logoutputfile.write("        \n")
logoutputfile.write("        \n")
logoutputfile.write("Execute build_powerplants \n")
print("Execute build_powerplants")
subprocess.run(["snakemake", "-call", "build_powerplants", "--cores", "all", "--printshellcmds", "--configfile", "../../Config/config.usa_PE.yaml"])

logoutputfile.write("        \n")
logoutputfile.write("        \n")
logoutputfile.write("Execute solve_all_networks #1 \n")
print("Execute solve_all_networks #1")
subprocess.run(["snakemake", "-call", "solve_all_networks", "--cores", "all", "--printshellcmds", "--configfile", "../../Config/config.usa_PE.yaml"])


logoutputfile.write("        \n")
logoutputfile.write("        \n")
logoutputfile.write("Set extendable carriers to False \n")
print("Set extendable carriers to False")

network_to_modify = pypsa.Network(network_path)
network_to_modify.generators.loc[:, "p_nom_extendable"] = False
network_to_modify.export_to_netcdf(network_path)


logoutputfile.write("        \n")
logoutputfile.write("        \n")
logoutputfile.write("Execute solve_all_networks #2 \n")
print("Execute solve_all_networks #2")
subprocess.run(["snakemake", "-call", "solve_all_networks", "--cores", "all", "--printshellcmds", "--configfile", "../../Config/config.usa_PE.yaml"])


logoutputfile.write("        \n")
logoutputfile.write("        \n")
logoutputfile.write("Execute solve_all_networks \n")
print("Execute solve_all_networks")
subprocess.run(["snakemake", "-call", "solve_all_networks", "--cores", "all", "--printshellcmds", "--configfile", "../../Config/config.usa_PE.yaml"])

logoutputfile.close()