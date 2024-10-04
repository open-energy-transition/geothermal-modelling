# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import argparse
import pathlib
import datetime as dt
import geopandas as gpd
import matplotlib.pyplot as plt
import pypsa
import cartopy.crs as ccrs
import sys
import pandas as pd


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


def plot_network_intersection(pypsa_df, eia_df, voltage_class, fig_name):
    base_network_pe_volt_class = pypsa_df.lines.loc[
        pypsa_df.lines["v_nom_class"] == voltage_class]
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


def plot_network_crossings(pypsa_df, eia_df, voltage_class, pypsa_title, fig_name):
    pass


# To perform the network topology comparison, execute python network_comparison --plot_network_topology
# To perform the network topology comparison, execute python network_comparison --plot_network_crossings
# To perform both actions, execute python network_comparison --plot_network_topology --plot_network_crossings
parser = argparse.ArgumentParser()
parser.add_argument("--plot_network_topology", help="Boolean: plot the network topology", action="store_true")
parser.add_argument("--plot_network_capacity", help="Boolean: plot the network capacity", action="store_true")
args = parser.parse_args()

# Initial configurations
eia_name = "EIA"
pypsa_name = "PyPSA"

# Set the paths
base_path = pathlib.Path(__file__).parent.parent.parent
log_file_dir_path = pathlib.Path(base_path, "analysis", "logs")
plot_dir_path = pathlib.Path(base_path, "analysis", "plots")
pypsa_earth_path = pathlib.Path(base_path, "workflow", "pypsa-earth")
base_network_pypsa_earth_path = pathlib.Path(pypsa_earth_path, "networks", "US_2021", "base.nc")
base_network_pypsa_usa_path = pathlib.Path(base_path.parent, "pypsa-usa", "workflow", "resources", "Default", "usa", "elec_base_network.nc")
eia_base_network_path = pathlib.Path(base_path.parent, "US_electric_transmission_lines_original.geojson")
gadm_shapes_path = pathlib.Path(base_path, "analysis", "data", "gadm41_USA_1.json")

###########
# Load data
###########
base_network_pypsa_earth = pypsa.Network(base_network_pypsa_earth_path)
base_network_pypsa_usa = pypsa.Network(base_network_pypsa_usa_path)
eia_base_network = gpd.read_file(eia_base_network_path)
gadm_shapes = gpd.read_file(gadm_shapes_path)

today_date = str(dt.datetime.now())
log_output_file = open(log_file_dir_path / f"output_network_comparison_{today_date[:10]}.txt", "w")




##########
# EIA data
##########

log_output_file.write("        \n")
log_output_file.write("        \n")
log_output_file.write(" Data preparation on the EIA input \n")
log_output_file.write(" --> shape of eia_base_network after reading it {} \n".format(eia_base_network.shape))

# add positions for the start- and end-points of the transmission lines

# --> Step 1 - Identify the lines with geometry type MultiLineString
eia_base_network["geometry_type"] = eia_base_network.geom_type
lines_with_multilinestring_geometry = eia_base_network[eia_base_network["geometry_type"]=="MultiLineString"]["OBJECTID_1"].values.tolist()
multilinestring_ratio = len(lines_with_multilinestring_geometry)/eia_base_network.shape[0]*100.0
if multilinestring_ratio > 10.0:
    sys.exit("The ratio of lines with geometry type MultiLineString is above 10%")
else:
    # exclude the lines with geometry MultiLineString. This is because the .boundary method yields (for this line) more than two boundary points
    eia_base_network = eia_base_network[~eia_base_network["OBJECTID_1"].isin(lines_with_multilinestring_geometry)]

log_output_file.write(" --> shape of eia_base_network after excluding multilinestrings {} \n".format(eia_base_network.shape))

# --> Step 2 - Compute the start- and end-points of a line
eia_base_network[["sub_0_coors", "sub_1_coors"]] = eia_base_network["geometry"].boundary.explode(index_parts=True).unstack()
log_output_file.write(" --> shape of eia_base_network after computing boundaries {} \n".format(eia_base_network.shape))

# --> Step 3 - Compute the start- and end-points of the line are located
eia_base_network_modified = eia_base_network.loc[:, ("OBJECTID_1", "sub_0_coors")]
eia_base_network_modified["geometry"] = eia_base_network_modified["sub_0_coors"]
log_output_file.write(" --> shape of eia_base_network before sub_0 spatial join {} \n".format(eia_base_network.shape))
spatial_join_gadm_eia_sub_0 = eia_base_network_modified.sjoin(gadm_shapes, how="left").loc[:, ("OBJECTID_1", "GID_1", "ISO_1")].rename(columns={"GID_1": "gid_sub_0", "ISO_1": "iso_sub_0"})
log_output_file.write(" --> shape of after sub_0 spatial join {} \n".format(spatial_join_gadm_eia_sub_0.shape))

eia_base_network_modified = eia_base_network.loc[:, ("OBJECTID_1", "sub_1_coors")]
eia_base_network_modified["geometry"] = eia_base_network_modified["sub_1_coors"]
log_output_file.write(" --> shape of eia_base_network before sub_1 spatial join {} \n".format(eia_base_network_modified.shape))
spatial_join_gadm_eia_sub_1 = eia_base_network_modified.sjoin(gadm_shapes, how="left").loc[:, ("OBJECTID_1", "GID_1", "ISO_1")].rename(columns={"GID_1": "gid_sub_1", "ISO_1": "iso_sub_1"})
log_output_file.write(" --> shape of after sub_1 spatial join {} \n".format(spatial_join_gadm_eia_sub_1.shape))

eia_base_network = pd.merge(eia_base_network, spatial_join_gadm_eia_sub_0, how="inner", on="OBJECTID_1")
eia_base_network = pd.merge(eia_base_network, spatial_join_gadm_eia_sub_1, how="inner", on="OBJECTID_1")
log_output_file.write(" --> shape of eia_base_network_after spatial joins {} \n".format(eia_base_network.shape))

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
        plot_network_intersection(base_network_pypsa_earth, eia_base_network, selected_voltage_class, fig_name_intersection)

# Comparison for the transmission capacities (PyPSA-Earth vs EIA)
if args.plot_network_capacity:
    pass
