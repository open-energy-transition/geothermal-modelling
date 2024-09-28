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
import pathlib
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import pypsa
import cartopy.crs as ccrs


def plot_network_comparison(pypsa_df, eia_df, voltage_class, pypsa_title, fig_name):
    pypsa_df.lines["line_width"] = 0.0
    pypsa_df.lines.loc[
        pypsa_df.lines["v_nom_class"] == voltage_class, "line_width"] = 1.0

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True,
                                   subplot_kw={"projection": ccrs.PlateCarree(), "projection": ccrs.PlateCarree()},
                                   figsize=(20, 4))
    pypsa_df.plot(line_widths=pypsa_df.lines["line_width"], bus_sizes=0,
                                  line_colors="blue", ax=ax1)
    pypsa_df.plot(line_widths=pypsa_df.lines["line_width"], bus_sizes=0,
                                  line_colors="white", ax=ax2)
    eia_df.loc[eia_df["VOLT_CLASS"] == voltage_class].plot(ax=ax2, color="orange")
    fig.suptitle("Comparison for voltage class: {}".format(voltage_class))
    ax1.set_aspect('equal')
    ax1.title.set_text(pypsa_title)
    ax2.set_aspect("equal")
    ax2.title.set_text("EIA network")
    fig.savefig(fig_name)


def plot_intersection(pypsa_df, eia_df, voltage_class, fig_name):
    base_network_pe_volt_class = pypsa_df.lines.loc[
        base_network_pypsa_earth.lines["v_nom_class"] == voltage_class]
    base_network_pe_volt_class = gpd.GeoDataFrame(base_network_pe_volt_class,
                                                  geometry=gpd.GeoSeries.from_wkt(base_network_pe_volt_class.geometry),
                                                  crs="EPSG:4326")
    base_network_pe_volt_class["union_geo"] = 0
    base_network_pe_volt_class = base_network_pe_volt_class.dissolve(by="union_geo")
    base_network_pe_volt_class = base_network_pe_volt_class.to_crs(3857)

    eia_base_network_volt_class = eia_df.loc[eia_df["VOLT_CLASS"] == voltage_class]
    eia_base_network_volt_class["union_geo"] = 0
    eia_base_network_volt_class = eia_base_network_volt_class.dissolve(by="union_geo")
    eia_base_network_volt_class = eia_base_network_volt_class.to_crs(3857)
    intersection_geometries = base_network_pe_volt_class.intersection(eia_base_network_volt_class)

    fig, ax = plt.subplots(figsize=(20, 4))
    ax.set_axis_off()

    base_network_pe_volt_class.plot(ax=ax, color="blue", kind="geo")
    eia_base_network_volt_class.plot(ax=ax, color="orange", kind="geo")
    intersection_geometries.plot(ax=ax, color="black", markersize=8)
    fig.savefig(fig_name)


# To perform the network topology comparison, execute python network_comparison --plot_network_topology
# To perform the network topology comparison, execute python network_comparison --plot_network_crossings
# To perform both actions, execute python network_comparison --plot_network_topology --plot_network_crossings
parser = argparse.ArgumentParser()
parser.add_argument("--plot_network_topology", help="Boolean: plot the network topology", action="store_true")
parser.add_argument("--plot_network_crossings", help="Boolean: plot the network crossings", action="store_true")
args = parser.parse_args()

# Initial configurations
eia_name = "EIA"
pypsa_name = "PyPSA"

# Set the paths
base_path = pathlib.Path(__file__).parent.parent.parent
log_file_dir_path = pathlib.Path(base_path, "logs")
plot_dir_path = pathlib.Path(base_path, "analysis", "plots")
pypsa_earth_path = pathlib.Path(base_path, "workflow", "pypsa-earth")
base_network_pypsa_earth_path = pathlib.Path(pypsa_earth_path, "networks", "US_2021", "base.nc")
base_network_pypsa_usa_path = pathlib.Path(base_path.parent, "pypsa-usa", "workflow", "resources", "Default", "usa", "elec_base_network.nc")
eia_base_network_path = pathlib.Path(base_path.parent, "US_Electric_Power_Transmission_Lines_5037807202786552385.geojson")
gadm_shapes_path = pathlib.Path(base_path, "analysis", "data", "gadm41_USA_1.json")

# Load data
base_network_pypsa_earth = pypsa.Network(base_network_pypsa_earth_path)
base_network_pypsa_usa = pypsa.Network(base_network_pypsa_usa_path)
eia_base_network = gpd.read_file(eia_base_network_path)
gadm_shapes = gpd.read_file(gadm_shapes_path)

# clean EIA data

# --> remove lines corresponding to voltage = -999999.0 kV
eia_base_network = eia_base_network.loc[eia_base_network["VOLTAGE"] != -999999.0]

# --> remove lines corresponding to voltage class 'DC'. All lines in the base.nc are AC
eia_base_network = eia_base_network.loc[eia_base_network["VOLT_CLASS"] != 'Dc']

# --> remove lines corresponding to voltage class 'Not Available'
eia_base_network = eia_base_network.loc[eia_base_network["VOLT_CLASS"] != 'Not Available']

# assign a voltage class to the pypsa-earth base.nc
base_network_pypsa_earth.lines["v_nom_class"] = base_network_pypsa_earth.lines["v_nom"]

v_nom_class_dict_pypsa_earth = {
    55.: 'Under 100',
    57.1: 'Under 100',
    60.: 'Under 100',
    66.: 'Under 100',
    69.: 'Under 100',
    70.: 'Under 100',
    88.: 'Under 100',
    92.: 'Under 100',
    100.: "100-161",
    115.: "100-161",
    120.: "100-161",
    125.: "100-161",
    138.: "100-161",
    160.: "100-161",
    161.: "100-161",
    220.: "220-287",
    230.: "220-287",
    287.: "220-287",
    345.: "345",
    500.: "500",
    765.: "735 And Above"
}

base_network_pypsa_earth.lines["v_nom_class"] = base_network_pypsa_earth.lines["v_nom_class"].replace(v_nom_class_dict_pypsa_earth)

# Comparison for the network topologies (PyPSA-Earth vs EIA)
# --> plot the EIA reference network and the PyPSA-Earth network
# --> plot the intersections between the networks

eia_voltage_classes = list(eia_base_network["VOLT_CLASS"].unique())

if args.plot_network_topology:
    for selected_voltage_class in eia_voltage_classes:
        fig_name_map = pathlib.Path(plot_dir_path, "network_comparison_pearth_for_voltage_class_{}.png".format(str(selected_voltage_class)))
        plot_network_comparison(base_network_pypsa_earth, eia_base_network, selected_voltage_class, "PyPSA-Earth base network", fig_name_map)
        fig_name_intersection = pathlib.Path(plot_dir_path, "network_comparison_intersection_{}.png".format(str(selected_voltage_class)))
        plot_intersection(base_network_pypsa_earth, eia_base_network, selected_voltage_class, fig_name_intersection)

# Comparison for the number of crossings (PyPSA-Earth vs EIA)


# --> plot the EIA reference network and the PyPSA-Earth network

eia_base_network_subset = eia_base_network[["ID", "TYPE", "VOLTAGE", "VOLT_CLASS", "SUB_1", "SUB_2", "SHAPE__Len", "geometry"]]
eia_base_network_subset["geometry"] = eia_base_network_subset["geometry"].boundary
spatial_join_gadm_eia = eia_base_network_subset.sjoin(gadm_shapes, how="inner")[["ID", "TYPE", "GID_1", "VOLTAGE", "VOLT_CLASS", "SUB_1", "SUB_2", "SHAPE__Len", "ISO_1", "NAME_1", "geometry"]]
count_df = (spatial_join_gadm_eia.groupby(["ID"])["ID"].count()-1).to_frame(name="crossings_count").reset_index()
spatial_join_gadm_eia = pd.merge(spatial_join_gadm_eia, count_df, how="inner", on="ID")
