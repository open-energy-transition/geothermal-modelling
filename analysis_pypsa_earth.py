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
import pandas as pd
import datetime as dt
import pypsa

# Initial configurations

# Set the paths
base_path = pathlib.Path.cwd()
log_file_dir = base_path / "logs"

# Ensure the logs directory exists
pathlib.Path(log_file_dir).mkdir(exist_ok=True)

# Define paths
pypsa_earth_path = pathlib.Path(base_path, "workflow", "pypsa-earth")
network_path = pathlib.Path(pypsa_earth_path, "results", "US_2021", "networks", "elec_s_10_ec_lcopt_Co2L-24H.nc")
eia_data_path = pathlib.Path(base_path, "data", "generation_eia.csv")

# Open log_output
today_date = str(dt.datetime.now())
log_output_file = open(log_file_dir / f'output_pypsa_earth_analysis_{today_date[:10]}.txt', 'w')

log_output_file.write("        \n")
log_output_file.write("        \n")
log_output_file.write("Import data \n")
print("Import data \n")

df_network = pypsa.Network(network_path)
df_eia = pd.read_csv(eia_data_path, index_col="Unnamed: 0")

s_max_pu = df_network.lines["s_max_pu"].unique()[0]

log_output_file.write("        \n")
log_output_file.write("        \n")
log_output_file.write(f"s_max_pu {s_max_pu} \n")
print(f"s_max_pu {s_max_pu} \n")

log_output_file.close()
