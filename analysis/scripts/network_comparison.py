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
import plotly.express as px
import pypsa
import cartopy.crs as ccrs
import sys
import pandas as pd
import numpy as np


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
    eia_df.loc[eia_df["v_nom_class"] == voltage_class].plot(ax=ax2, color="orange")

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
    eia_base_network_volt_class = eia_df.loc[eia_df["v_nom_class"] == voltage_class]
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


def plot_network_crossings(pypsa_df, eia_df, color_dictionary, voltage_classes_list, output_base_path, plot_base_path):

    ###############
    # GADM shapes #
    ###############
    gadm_eia_crossings_df = eia_df.groupby(["state_0", "state_1", "v_nom_class"]).count().reset_index().loc[:,
                          ("state_0", "state_1", "v_nom_class", "ID")].rename(columns={"ID": "crossings"})
    gadm_eia_crossings_df["source"] = "EIA"

    gadm_pearth_crossings_df = pypsa_df.lines.groupby(["state_0", "state_1", "v_nom_class"]).count().reset_index().loc[:,
                          ("state_0", "state_1", "v_nom_class", "Line")].rename(columns={"Line": "crossings"})
    gadm_pearth_crossings_df["source"] = "PyPSA"

    gadm_pearth_crossings_parallel_df = pypsa_df.lines.groupby(
        ["state_0", "state_1", "v_nom_class"])["num_parallel"].sum().reset_index().loc[:,
                          ("state_0", "state_1", "v_nom_class", "num_parallel")].rename(columns={"num_parallel": "crossings"})
    gadm_pearth_crossings_parallel_df["source"] = "PyPSA_parallel"

    gadm_network_counts = pd.concat([gadm_eia_crossings_df, gadm_pearth_crossings_df, gadm_pearth_crossings_parallel_df])
    gadm_network_counts = gadm_network_counts.set_index(
        ["source", "state_0", "state_1", "v_nom_class"]
    ).unstack("source").droplevel(axis=1, level=0).reset_index()
    gadm_network_counts = gadm_network_counts[["state_0", "state_1", "v_nom_class", "PyPSA", "EIA", "PyPSA_parallel"]]

    # --> state crossings. The start and end states were assigned by means of a spatial join with the GADM shapes
    state_crossings_counts = gadm_network_counts.query("state_0 != state_1")

    state_crossings_counts_voltage = state_crossings_counts.groupby("v_nom_class")[["PyPSA", "EIA", "PyPSA_parallel"]].sum().reindex(
        ["Under 100", "100-161", "220-287", "345", "500", "735 And Above"]).reset_index()
    state_crossings_counts_voltage.to_csv(pathlib.Path(output_base_path, "gadm_state_crossings_counts_by_voltage.csv"))

    fig = px.bar(state_crossings_counts_voltage,
                 x="v_nom_class",
                 y=["PyPSA", "PyPSA_parallel", "EIA"],
                 barmode="group",
                 color_discrete_map=color_dictionary,
                 text_auto='.2s',
                 title="Number of transmission line crossings per voltage class"
                 ).update_layout(
        xaxis_title="Voltage class (kV)", yaxis_title="Number of transmission line crossings"
    )
    fig.update_traces(textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
    fig.write_image(pathlib.Path(plot_base_path, "gadm_state_crossings_counts_by_voltage.png"))

    state_crossings_counts["delta_PyPSA"] = (state_crossings_counts["PyPSA"] - state_crossings_counts["EIA"]) / state_crossings_counts["EIA"]*100.0
    state_crossings_counts["delta_PyPSA_parallel"] = (state_crossings_counts["PyPSA_parallel"] - state_crossings_counts["EIA"]) / state_crossings_counts["EIA"]*100.0
    state_crossings_counts["coalesce"] = state_crossings_counts[["state_0", "state_1"]].agg('-->'.join, axis=1)
    state_crossings_counts.to_csv(pathlib.Path(output_base_path, "gadm_state_crossings_counts.csv"))

    for voltage_class in voltage_classes_list:
        filtered_df = state_crossings_counts.loc[state_crossings_counts["v_nom_class"] == voltage_class]
        fig = px.scatter(filtered_df,
                         x="coalesce",
                         y=["delta_PyPSA", "delta_PyPSA_parallel"],
                         color_discrete_map=color_dictionary,
                         title="Voltage class: {}".format(voltage_class)
                         ).update_layout(
        xaxis_title="States", yaxis_title="Error (%)")
        fig.write_image(pathlib.Path(plot_base_path, "gadm_state_crossings_counts_for_voltage_{}.png".format(voltage_class)))

    # --> state lines. The start and end states were assigned by means of a spatial join with the GADM shapes
    state_lines_counts = gadm_network_counts.query("state_0 == state_1")
    state_lines_counts_by_voltage = state_lines_counts.groupby("v_nom_class")[["PyPSA", "EIA", "PyPSA_parallel"]].sum().reindex(
        ["Under 100", "100-161", "220-287", "345", "500", "735 And Above"]).reset_index()
    state_lines_counts_by_voltage.to_csv(pathlib.Path(output_base_path, "gadm_state_lines_counts_by_voltage.csv"))

    fig = px.bar(state_lines_counts_by_voltage,
                 x="v_nom_class",
                 y=["PyPSA", "PyPSA_parallel", "EIA"],
                 barmode="group",
                 color_discrete_map=color_dictionary,
                 text_auto='.2s',
                 title="Number of transmission line per state per voltage class"
                 ).update_layout(
        xaxis_title="Voltage class (kV)", yaxis_title="Number of state lines"
    )
    fig.update_traces(textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
    fig.write_image(pathlib.Path(plot_base_path, "gadm_state_lines_counts_by_voltage.png"))

    ##############
    # IPM shapes #
    ##############
    ipm_eia_crossings_df = eia_df.groupby(["ipm_region_0", "ipm_region_1", "v_nom_class"]).count().reset_index().loc[:,
                          ("ipm_region_0", "ipm_region_1", "v_nom_class", "ID")].rename(columns={"ID": "crossings"})
    ipm_eia_crossings_df["source"] = "EIA"

    ipm_pearth_crossings_df = pypsa_df.lines.groupby(["ipm_region_0", "ipm_region_1", "v_nom_class"]).count().reset_index().loc[:,
                          ("ipm_region_0", "ipm_region_1", "v_nom_class", "Line")].rename(columns={"Line": "crossings"})
    ipm_pearth_crossings_df["source"] = "PyPSA"

    ipm_pearth_crossings_parallel_df = pypsa_df.lines.groupby(
        ["ipm_region_0", "ipm_region_1", "v_nom_class"])["num_parallel"].sum().reset_index().loc[:,
                          ("ipm_region_0", "ipm_region_1", "v_nom_class", "num_parallel")].rename(columns={"num_parallel": "crossings"})
    ipm_pearth_crossings_parallel_df["source"] = "PyPSA_parallel"

    ipm_network_counts = pd.concat(
        [ipm_eia_crossings_df, ipm_pearth_crossings_df, ipm_pearth_crossings_parallel_df])
    ipm_network_counts = ipm_network_counts.set_index(
        ["source", "ipm_region_0", "ipm_region_1", "v_nom_class"]
    ).unstack("source").droplevel(axis=1, level=0).reset_index()
    ipm_network_counts = ipm_network_counts[["ipm_region_0", "ipm_region_1", "v_nom_class", "PyPSA", "EIA", "PyPSA_parallel"]]

    # --> ipm region crossings. The start and end ipm regions were assigned by means of a spatial join with the IPM shapes
    region_crossings_counts = ipm_network_counts.query("ipm_region_0 != ipm_region_1")

    region_crossings_counts_voltage = region_crossings_counts.groupby("v_nom_class")[
        ["PyPSA", "EIA", "PyPSA_parallel"]].sum().reindex(
        ["Under 100", "100-161", "220-287", "345", "500", "735 And Above"]).reset_index()
    region_crossings_counts_voltage.to_csv(pathlib.Path(output_base_path, "ipm_region_crossings_counts_by_voltage.csv"))

    fig = px.bar(region_crossings_counts_voltage,
                 x="v_nom_class",
                 y=["PyPSA", "PyPSA_parallel", "EIA"],
                 barmode="group",
                 color_discrete_map=color_dictionary,
                 text_auto='.2s',
                 title="Number of transmission line crossings per voltage class"
                 ).update_layout(
        xaxis_title="Voltage class (kV)", yaxis_title="Number of transmission line crossings"
    )
    fig.update_traces(textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
    fig.write_image(pathlib.Path(plot_base_path, "ipm_region_crossings_counts_by_voltage.png"))

    region_crossings_counts["delta_PyPSA"] = (region_crossings_counts["PyPSA"] - region_crossings_counts["EIA"]) / region_crossings_counts["EIA"] * 100.0
    region_crossings_counts["delta_PyPSA_parallel"] = (region_crossings_counts["PyPSA_parallel"] - region_crossings_counts[
        "EIA"]) / region_crossings_counts["EIA"] * 100.0
    region_crossings_counts["coalesce"] = region_crossings_counts[["ipm_region_0", "ipm_region_1"]].agg('-->'.join, axis=1)
    region_crossings_counts.to_csv(pathlib.Path(output_base_path, "ipm_region_crossings_counts.csv"))

    for voltage_class in voltage_classes_list:
        filtered_df = region_crossings_counts.loc[region_crossings_counts["v_nom_class"] == voltage_class]
        fig = px.scatter(filtered_df,
                         x="coalesce",
                         y=["delta_PyPSA", "delta_PyPSA_parallel"],
                         color_discrete_map=color_dictionary,
                         title="Voltage class: {}".format(voltage_class)
                         ).update_layout(
            xaxis_title="IPM Region", yaxis_title="Error (%)")
        fig.write_image(
            pathlib.Path(plot_base_path, "ipm_region_crossings_counts_for_voltage_{}.png".format(voltage_class)))


def plot_network_capacity(pypsa_df, color_dictionary, base_path, output_base_path, plot_base_path):
    transmission_capacities_path = pathlib.Path(base_path, "analysis", "data", "transmission_single_epaipm.csv")
    transmission_capacities_df = pd.read_csv(transmission_capacities_path).rename(columns={
        "region_from": "ipm_region_0",
        "region_to": "ipm_region_1",
        "firm_ttc_mw": "capacity (MW)"
    })
    transmission_capacities_df["source"] = "IPM"
    transmission_capacities_df = transmission_capacities_df.loc[:, ("ipm_region_0", "ipm_region_1", "capacity (MW)", "source")]

    transmission_capacities_df.to_csv("ipm_test_cap.csv")

    pypsa_df.lines["s_nom_num_parallel"] = pypsa_df.lines["s_nom"] * pypsa_df.lines["num_parallel"]
    #pypsa_df.lines["s_nom_num_parallel_pu"] = pypsa_df.lines["s_nom"] * pypsa_df.lines["num_parallel"] * pypsa_df.lines["s_max_pu"]
    ipm_region_pearth_transmission_capacities = pypsa_df.lines.query("ipm_region_0 != ipm_region_1").groupby(["ipm_region_0", "ipm_region_1"])["s_nom_num_parallel"].sum().reset_index().loc[:,
                          ("ipm_region_0", "ipm_region_1", "s_nom_num_parallel")].rename(columns={"s_nom_num_parallel": "capacity (MW)"})
    ipm_region_pearth_transmission_capacities["source"] = "PyPSA"

    capacity_df = pd.concat(
        [ipm_region_pearth_transmission_capacities, transmission_capacities_df])

    capacity_df = capacity_df.set_index(
        ["source", "ipm_region_0", "ipm_region_1"]
    ).unstack("source").droplevel(axis=1, level=0).reset_index()
    capacity_df["Error (%)"] = (capacity_df["PyPSA"]-capacity_df["IPM"])/capacity_df["IPM"]*100.0
    capacity_df["coalesce"] = capacity_df[["ipm_region_0", "ipm_region_1"]].agg('-->'.join, axis=1)
    capacity_df.to_csv(pathlib.Path(output_base_path, "ipm_pearth_capacities.csv"))

    fig = px.scatter(capacity_df,
                     x="coalesce",
                     y="Error (%)",
                     color_discrete_map=color_dictionary,
                     title="Transmission capacities"
                     ).update_layout(
        xaxis_title="IPM Region", yaxis_title="Error (%)")
    fig.write_image(
        pathlib.Path(plot_base_path, "transmission_capacities.png"))

    fig = px.box(capacity_df,
                     y="Error (%)",
                     title="Error Box Plot on Transmission capacities"
                     )
    fig.write_image(
        pathlib.Path(plot_base_path, "box_transmission_capacities.png"))


def parse_input_arguments():
    """
    Example:
    -) to perform the network topology comparison, execute python network_comparison.py --plot_network_topology
    -) to perform the network topology comparison, execute python network_comparison.py --plot_network_crossings
    -) to perform both actions, execute python network_comparison.py --plot_network_topology --plot_network_crossings

    Returns
    Args
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot_network_topology", help="Boolean: plot the network topology", action="store_true")
    parser.add_argument("--plot_network_crossings", help="Boolean: plot the network crossings", action="store_true")
    parser.add_argument("--plot_network_capacity", help="Boolean: plot the network capacity", action="store_true")
    return parser.parse_args()


def parse_inputs(base_path, log_file_dir_path):
    pypsa_earth_path = pathlib.Path(base_path, "workflow", "pypsa-earth")
    base_network_pypsa_earth_path = pathlib.Path(pypsa_earth_path, "networks", "US_2021", "base.nc")
    eia_base_network_path = pathlib.Path(base_path.parent, "US_electric_transmission_lines_original.geojson")
    gadm_shapes_path = pathlib.Path(base_path, "analysis", "data", "gadm41_USA_1.json")
    ipm_shapes_path = pathlib.Path(base_path, "analysis", "data", "ipm_v6_regions", "IPM_Regions_201770405.shp")

    ###########
    # Load data
    ###########
    base_network_pypsa_earth = pypsa.Network(base_network_pypsa_earth_path)
    eia_base_network = gpd.read_file(eia_base_network_path)
    gadm_shapes = gpd.read_file(gadm_shapes_path)
    ipm_shapes = gpd.read_file(ipm_shapes_path).to_crs("4326")
    today_date = str(dt.datetime.now())
    log_output_file = open(log_file_dir_path / f"output_network_comparison_{today_date[:10]}.txt", "w")

    ##########
    # EIA data
    ##########
    log_output_file.write("        \n")
    log_output_file.write("        \n")
    log_output_file.write(" Data preparation on the EIA base network \n")
    log_output_file.write(" --> shape of eia_base_network after reading it in {} \n".format(eia_base_network.shape))

    # add positions for the start- and end-points of the transmission lines

    # --> Identify the lines with geometry type MultiLineString
    eia_base_network["geometry_type"] = eia_base_network.geom_type
    lines_with_multilinestring_geometry = eia_base_network[eia_base_network["geometry_type"]=="MultiLineString"]["OBJECTID_1"].values.tolist()
    multilinestring_ratio = len(lines_with_multilinestring_geometry)/eia_base_network.shape[0]*100.0
    if multilinestring_ratio > 10.0:
        sys.exit("The ratio of lines with geometry type MultiLineString is above 10%")
    else:
        # exclude the lines with geometry MultiLineString. This is because the .boundary method yields (for this line) more than two boundary points
        eia_base_network = eia_base_network[~eia_base_network["OBJECTID_1"].isin(lines_with_multilinestring_geometry)]

    log_output_file.write(" --> shape of eia_base_network after excluding multilinestrings {} \n".format(eia_base_network.shape))

    # --> Compute the start- and end-points of a line
    eia_base_network[["sub_0_coors", "sub_1_coors"]] = eia_base_network["geometry"].boundary.explode(index_parts=True).unstack()
    log_output_file.write(" --> shape of eia_base_network after computing boundaries {} \n".format(eia_base_network.shape))

    # --> Determine where the start- and end-points of the line are located. In particular, we perform the spatial
    # joins with:
    #  -) the GADM shapes (level 1) to get the US state
    #  -) the IPM shapes to get the IPM region

    # Bus 0
    eia_base_network_modified = eia_base_network.loc[:, ("OBJECTID_1", "sub_0_coors")]
    eia_base_network_modified["geometry"] = eia_base_network_modified["sub_0_coors"]
    log_output_file.write(" --> shape of eia_base_network before sub_0 spatial join {} \n".format(eia_base_network.shape))
    spatial_join_gadm_eia_sub_0 = eia_base_network_modified.sjoin(gadm_shapes, how="left").loc[:, ("OBJECTID_1", "ISO_1")].rename(columns={"ISO_1": "state_0"})
    spatial_join_ipm_eia_sub_0 = eia_base_network_modified.sjoin(ipm_shapes, how="left").loc[:, ("OBJECTID_1", "IPM_Region")].rename(columns={"IPM_Region": "ipm_region_0"})
    spatial_join_eia_sub_0 = pd.merge(spatial_join_gadm_eia_sub_0, spatial_join_ipm_eia_sub_0, how="inner", on="OBJECTID_1")
    log_output_file.write(" --> shape of eia_base_network after sub_0 spatial join with gadm {} \n".format(spatial_join_gadm_eia_sub_0.shape))
    log_output_file.write(" --> shape of eia_base_network after sub_0 spatial join with ipm{} \n".format(spatial_join_ipm_eia_sub_0.shape))
    log_output_file.write(" --> shape of eia_base_network after sub_0 spatial join {} \n".format(spatial_join_eia_sub_0.shape))

    # Bus 1
    eia_base_network_modified = eia_base_network.loc[:, ("OBJECTID_1", "sub_1_coors")]
    eia_base_network_modified["geometry"] = eia_base_network_modified["sub_1_coors"]
    log_output_file.write(" --> shape of eia_base_network before sub_1 spatial join {} \n".format(eia_base_network_modified.shape))
    spatial_join_gadm_eia_sub_1 = eia_base_network_modified.sjoin(gadm_shapes, how="left").loc[:, ("OBJECTID_1", "ISO_1")].rename(columns={"ISO_1": "state_1"})
    spatial_join_ipm_eia_sub_1 = eia_base_network_modified.sjoin(ipm_shapes, how="left").loc[:, ("OBJECTID_1", "IPM_Region")].rename(columns={"IPM_Region": "ipm_region_1"})
    spatial_join_eia_sub_1 = pd.merge(spatial_join_gadm_eia_sub_1, spatial_join_ipm_eia_sub_1, how="inner", on="OBJECTID_1")
    log_output_file.write(" --> shape of eia_base_network after sub_1 spatial join with gadm {} \n".format(spatial_join_gadm_eia_sub_1.shape))
    log_output_file.write(" --> shape of eia_base_network after sub_1 spatial join with ipm{} \n".format(spatial_join_ipm_eia_sub_1.shape))
    log_output_file.write(" --> shape of eia_base_network after sub_1 spatial join {} \n".format(spatial_join_eia_sub_1.shape))

    # --> Inner join the results
    eia_base_network = pd.merge(eia_base_network, spatial_join_eia_sub_0, how="inner", on="OBJECTID_1")
    eia_base_network = pd.merge(eia_base_network, spatial_join_eia_sub_1, how="inner", on="OBJECTID_1")
    log_output_file.write(" --> shape of eia_base_network after the inner joins {} \n".format(eia_base_network.shape))

    # Clean the EIA data from lines with unnecessary voltages and voltage classes
    eia_base_network = eia_base_network.rename(columns={"VOLTAGE": "v_nom", "VOLT_CLASS": "v_nom_class"})

    # --> Remove lines corresponding to voltage = -999999.0 kV
    eia_base_network = eia_base_network.loc[eia_base_network["v_nom"] != -999999.0]
    log_output_file.write(" --> shape of eia_base_network after removing the lines with voltage -999999.0 kV {} \n".format(eia_base_network.shape))

    # --> Remove lines corresponding to voltage class 'DC'. All lines in the base.nc are AC
    eia_base_network = eia_base_network.loc[eia_base_network["v_nom_class"] != 'Dc']
    log_output_file.write(" --> shape of eia_base_network after removing the lines with voltage class 'Dc' {} \n".format(eia_base_network.shape))

    # --> Remove lines corresponding to voltage class 'Not Available'
    eia_base_network = eia_base_network.loc[eia_base_network["v_nom_class"] != 'Not Available']
    log_output_file.write(" --> shape of eia_base_network after removing the lines with voltage class 'Not Available' {} \n".format(eia_base_network.shape))

    #####################
    # PyPSA-Earth base.nc
    #####################

    log_output_file.write("        \n")
    log_output_file.write("        \n")
    log_output_file.write(" Data preparation on the PyPSA-Earth base network \n")

    # --> Determine where the start- and end-points of the line are located. In particular, we perform the spatial
    # joins with:
    #  -) the GADM shapes (level 1) to get the US state
    #  -) the IPM shapes to get the IPM region
    log_output_file.write(" --> shape of pypsa-earth base network after reading it in {} \n".format(base_network_pypsa_earth.lines.shape))

    # Bus 0
    base_network_pypsa_earth_geopandas_bus_0 = gpd.GeoDataFrame(base_network_pypsa_earth.lines, geometry=gpd.GeoSeries.from_wkt(base_network_pypsa_earth.lines.bus_0_coors), crs="EPSG:4326").reset_index()
    spatial_join_gadm_pearth_bus_0 = base_network_pypsa_earth_geopandas_bus_0.sjoin(gadm_shapes, how="left").loc[:, ("Line", "ISO_1")].rename(columns={"ISO_1": "state_0"})
    spatial_join_ipm_pearth_bus_0 = base_network_pypsa_earth_geopandas_bus_0.sjoin(ipm_shapes, how="left").loc[:, ("Line", "IPM_Region")].rename(columns={"IPM_Region": "ipm_region_0"})
    spatial_join_pearth_bus_0 = pd.merge(spatial_join_gadm_pearth_bus_0, spatial_join_ipm_pearth_bus_0, how="inner", on="Line")
    log_output_file.write(" --> shape of base.nc after sub_0 spatial join with gadm {} \n".format(spatial_join_gadm_pearth_bus_0.shape))
    log_output_file.write(" --> shape of base.nc after sub_0 spatial join with ipm{} \n".format(spatial_join_ipm_pearth_bus_0.shape))
    log_output_file.write(" --> shape of base.nc after sub_0 spatial join {} \n".format(spatial_join_pearth_bus_0.shape))

    # Bus 1
    base_network_pypsa_earth_geopandas_bus_1 = gpd.GeoDataFrame(base_network_pypsa_earth.lines, geometry=gpd.GeoSeries.from_wkt(base_network_pypsa_earth.lines.bus_1_coors), crs="EPSG:4326").reset_index()
    spatial_join_gadm_pearth_bus_1 = base_network_pypsa_earth_geopandas_bus_1.sjoin(gadm_shapes, how="left").loc[:, ("Line", "ISO_1")].rename(columns={"ISO_1": "state_1"})
    spatial_join_ipm_pearth_bus_1 = base_network_pypsa_earth_geopandas_bus_1.sjoin(ipm_shapes, how="left").loc[:, ("Line", "IPM_Region")].rename(columns={"IPM_Region": "ipm_region_1"})
    spatial_join_pearth_bus_1 = pd.merge(spatial_join_gadm_pearth_bus_1, spatial_join_ipm_pearth_bus_1, how="inner", on="Line")
    log_output_file.write(" --> shape of base.nc after sub_1 spatial join with gadm {} \n".format(spatial_join_gadm_pearth_bus_1.shape))
    log_output_file.write(" --> shape of base.nc after sub_1 spatial join with ipm{} \n".format(spatial_join_ipm_pearth_bus_1.shape))
    log_output_file.write(" --> shape of base.nc after sub_1 spatial join {} \n".format(spatial_join_pearth_bus_1.shape))

    base_network_pypsa_earth.lines = pd.merge(base_network_pypsa_earth.lines, spatial_join_pearth_bus_0, how="inner", on="Line")
    base_network_pypsa_earth.lines = pd.merge(base_network_pypsa_earth.lines, spatial_join_pearth_bus_1, how="inner", on="Line")
    log_output_file.write(" --> shape of pypsa-earth base network after the spatial joins {} \n".format(base_network_pypsa_earth.lines.shape))

    # --> Assign a voltage class to the pypsa-earth base.nc
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

    return eia_base_network, base_network_pypsa_earth


if __name__ == '__main__':

    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs")
    plot_path = pathlib.Path(default_path, "analysis", "plots")
    output_path = pathlib.Path(default_path, "analysis", "outputs")
    ccs_color_dict = {"EIA": "#FF8C00", "PyPSA": "#0000FF", "PyPSA_parallel": "#228B22", "delta_PyPSA": "#0000FF", "delta_PyPSA_parallel": "#228B22", "Error (%)": "#FF7F50"}

    network_eia_df, network_pypsa_df = parse_inputs(default_path, log_path)

    args = parse_input_arguments()

    eia_voltage_classes = list(network_eia_df["v_nom_class"].unique())

    if args.plot_network_topology:
        for selected_voltage_class in eia_voltage_classes:
            fig_name_map = pathlib.Path(plot_path, "network_comparison_pearth_for_voltage_class_{}.png".format(
                str(selected_voltage_class)))
            plot_network_comparison(network_pypsa_df, network_eia_df, selected_voltage_class,
                                    "PyPSA-Earth base network", fig_name_map)
            fig_name_intersection = pathlib.Path(plot_path, "network_comparison_intersection_{}.png".format(
                str(selected_voltage_class)))
            plot_network_intersection(network_pypsa_df, network_eia_df, selected_voltage_class,
                                      fig_name_intersection)

    # Comparison for the transmission crossings (PyPSA-Earth vs EIA) using:
    # -) the GADM shapes(level 1) for the US state comparison
    # -) the IPM shapes for the IPM region comparison
    if args.plot_network_crossings:
        plot_network_crossings(network_pypsa_df, network_eia_df, ccs_color_dict, eia_voltage_classes, output_path, plot_path)

    # Comparison for the transmission capacities (PyPSA-Earth vs IPM transmission capacities)
    if args.plot_network_capacity:
        plot_network_capacity(network_pypsa_df, ccs_color_dict, default_path, output_path, plot_path)
