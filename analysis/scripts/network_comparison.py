# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import cartopy.crs as ccrs
import datetime as dt
import geopandas as gpd
import matplotlib.pyplot as plt
import logging
import numpy as np
import pandas as pd
import pathlib
import plotly.express as px
import pypsa
import shapely as spl
import sys

logger = logging.getLogger(__name__)


def plot_network_topology_comparison(pypsa_df, eia_df, voltage_class, pypsa_title, fig_name):
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


def plot_network_topology_intersection(pypsa_df, eia_df, voltage_class, fig_name):
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
    gadm_network_counts = gadm_network_counts.loc[:, ("state_0", "state_1", "v_nom_class", "PyPSA", "EIA", "PyPSA_parallel")]

    # --> state crossings. The start and end states were assigned by means of a spatial join with the GADM shapes
    state_crossings_counts = gadm_network_counts.query("state_0 != state_1").copy()

    state_crossings_counts_voltage = state_crossings_counts.groupby("v_nom_class")[["PyPSA", "EIA", "PyPSA_parallel"]].sum().reindex(
        ["Under 100", "100-161", "220-287", "345", "500", "735 And Above"]).reset_index()
    state_crossings_counts_voltage.to_csv(pathlib.Path(output_base_path, "gadm_state_crossings_counts_by_voltage.csv"), index=False)

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
    state_crossings_counts.to_csv(pathlib.Path(output_base_path, "gadm_state_crossings_counts.csv"), index=False)

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
    state_lines_counts_by_voltage.to_csv(pathlib.Path(output_base_path, "gadm_state_lines_counts_by_voltage.csv"), index=False)

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
    region_crossings_counts = ipm_network_counts.query("ipm_region_0 != ipm_region_1").copy()

    region_crossings_counts_voltage = region_crossings_counts.groupby("v_nom_class")[
        ["PyPSA", "EIA", "PyPSA_parallel"]].sum().reindex(
        ["Under 100", "100-161", "220-287", "345", "500", "735 And Above"]).reset_index()
    region_crossings_counts_voltage.to_csv(pathlib.Path(output_base_path, "ipm_region_crossings_counts_by_voltage.csv"), index=False)

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
    region_crossings_counts.to_csv(pathlib.Path(output_base_path, "ipm_region_crossings_counts.csv"), index=False)

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


def plot_network_capacity_ipm(pypsa_df, ipm_shapes_gdf, color_dictionary, log_output_file, base_path, output_base_path, plot_base_path, model):

    transmission_capacities_path = pathlib.Path(base_path, "analysis", "gdrive_data", "data", "transmission_grid_data", "transmission_single_epaipm.csv")
    transmission_capacities_df = pd.read_csv(transmission_capacities_path).rename(columns={
        "region_from": "ipm_region_0",
        "region_to": "ipm_region_1",
        "nonfirm_ttc_mw": "capacity (MW)"
    })
    transmission_capacities_df["source"] = "IPM"
    transmission_capacities_df = transmission_capacities_df.loc[:, ("ipm_region_0", "ipm_region_1", "capacity (MW)", "source")]
    transmission_capacities_df[["ipm_region_0", "ipm_region_1"]] = np.sort(transmission_capacities_df[["ipm_region_0", "ipm_region_1"]])
    grouped_df = transmission_capacities_df.groupby(["ipm_region_0", "ipm_region_1"])["capacity (MW)"].max()
    transmission_capacities_df = pd.merge(transmission_capacities_df, grouped_df, on=["ipm_region_0", "ipm_region_1"])
    transmission_capacities_df = transmission_capacities_df.loc[:, ("ipm_region_0", "ipm_region_1", "capacity (MW)_y", "source")].rename(columns={"capacity (MW)_y": "capacity (MW)"})
    transmission_capacities_df = transmission_capacities_df.drop_duplicates(keep="first")

    pypsa_df["s_nom_num_parallel"] = pypsa_df["s_nom"] * pypsa_df["num_parallel"]
    pypsa_df["s_nom_num_parallel_pu"] = pypsa_df["s_nom"] * pypsa_df["num_parallel"] * pypsa_df["s_max_pu"]
    pypsa_df[["ipm_region_0", "ipm_region_1"]] = np.sort(pypsa_df[["ipm_region_0", "ipm_region_1"]])
    ipm_region_pypsa_transmission_capacities = pypsa_df.query("ipm_region_0 != ipm_region_1").groupby(["ipm_region_0", "ipm_region_1"])["s_nom_num_parallel"].sum().reset_index().loc[:,
                          ("ipm_region_0", "ipm_region_1", "s_nom_num_parallel")].rename(columns={"s_nom_num_parallel": "capacity (MW)"})
    ipm_region_pypsa_transmission_capacities["source"] = "PyPSA"

    capacity_df = pd.concat(
        [ipm_region_pypsa_transmission_capacities, transmission_capacities_df])

    capacity_df = capacity_df.set_index(
        ["source", "ipm_region_0", "ipm_region_1"]
    ).unstack("source").droplevel(axis=1, level=0)

    # remove values where IPM is zero or None
    capacity_df = capacity_df.query("IPM != 0.0 & ~IPM.isna() & ~PyPSA.isna()").reset_index()

    # compute percentage error
    capacity_df["Error wrt IPM (%)"] = (capacity_df["PyPSA"]-capacity_df["IPM"])/capacity_df["IPM"]*100.0
    capacity_df["Error wrt PyPSA (%)"] = (capacity_df["PyPSA"]-capacity_df["IPM"])/capacity_df["PyPSA"]*100.0
    capacity_df["factor_IPM_over_PyPSA"] = capacity_df["IPM"]/capacity_df["PyPSA"]
    capacity_df["coalesce"] = capacity_df[["ipm_region_0", "ipm_region_1"]].agg('-->'.join, axis=1)

    capacity_df.to_csv(pathlib.Path(output_base_path, f"{model}_ipm_capacities.csv"), index=False)

    log_output_file.write("====")
    log_output_file.write("{}: median error wrt IPM: {} \n".format(model, np.round(capacity_df["Error wrt IPM (%)"].median(), 2)))
    log_output_file.write("{}: mean error wrt IPM: {} \n".format(model, np.round(capacity_df["Error wrt IPM (%)"].mean(), 2)))
    log_output_file.write("{}: median error wrt PyPSA: {} \n".format(model, np.round(capacity_df["Error wrt PyPSA (%)"].median(), 2)))
    log_output_file.write("{}: mean error wrt PyPSA: {} \n".format(model, np.round(capacity_df["Error wrt PyPSA (%)"].mean(), 2)))
    log_output_file.write("====")
    logger.info("====")
    logger.info("{}: median error wrt IPM: {} \n".format(model, np.round(capacity_df["Error wrt IPM (%)"].median(), 2)))
    logger.info("{}: mean error wrt IPM: {} \n".format(model, np.round(capacity_df["Error wrt IPM (%)"].mean(), 2)))
    logger.info("{}: median error wrt PyPSA: {} \n".format(model, np.round(capacity_df["Error wrt PyPSA (%)"].median(), 2)))
    logger.info("{}: mean error wrt PyPSA: {} \n".format(model, np.round(capacity_df["Error wrt PyPSA (%)"].mean(), 2)))
    logger.info("====")

    fig = px.scatter(capacity_df,
                     x="coalesce",
                     y="Error wrt IPM (%)",
                     color_discrete_map=color_dictionary,
                     title="Transmission capacities"
                     ).update_layout(
        xaxis_title="IPM Region", yaxis_title="Error (%)")
    fig.write_image(
        pathlib.Path(plot_base_path, f"{model}_transmission_capacities_wrt_IPM.png"))

    fig = px.scatter(capacity_df,
                     x="coalesce",
                     y="Error wrt PyPSA (%)",
                     color_discrete_map=color_dictionary,
                     title="Transmission capacities"
                     ).update_layout(
        xaxis_title="IPM Region", yaxis_title="Error (%)")
    fig.write_image(
        pathlib.Path(plot_base_path, f"{model}_transmission_capacities_wrt_PyPSA.png"))

    fig = px.box(capacity_df,
                     y="Error wrt IPM (%)",
                     title="Error Box Plot on Transmission capacities"
                     )
    fig.write_image(
        pathlib.Path(plot_base_path, f"{model}_box_transmission_capacities_wrt_IPM.png"))

    fig = px.box(capacity_df,
                     y="Error wrt PyPSA (%)",
                     title="Error Box Plot on Transmission capacities"
                     )
    fig.write_image(
        pathlib.Path(plot_base_path, f"{model}_box_transmission_capacities_wrt_PyPSA.png"))

    ipm_shapes_gdf = ipm_shapes_gdf.to_crs("3857")
    ipm_shapes_gdf["ipm_region_centroid"] = ipm_shapes_gdf.centroid
    capacity_df["ipm_region_0_centroid"] = capacity_df["ipm_region_0"].map(
        dict(ipm_shapes_gdf[['IPM_Region', 'ipm_region_centroid']].values))
    capacity_df["ipm_region_1_centroid"] = capacity_df["ipm_region_1"].map(
        dict(ipm_shapes_gdf[['IPM_Region', 'ipm_region_centroid']].values))
    capacity_df["geometry"] = capacity_df.apply(
        lambda x: spl.LineString([x.ipm_region_0_centroid, x.ipm_region_1_centroid]), axis=1)
    capacity_df = capacity_df.drop(["ipm_region_0_centroid", "ipm_region_1_centroid"], axis=1)
    capacity_df_normal = capacity_df.query("factor_IPM_over_PyPSA <= 1.0")
    capacity_df_outlier = capacity_df.query("factor_IPM_over_PyPSA > 1.0")
    ipm_geo_data_normal = gpd.GeoDataFrame(capacity_df_normal, geometry=capacity_df_normal.geometry, crs="EPSG:3857")
    ipm_geo_data_outlier = gpd.GeoDataFrame(capacity_df_outlier, geometry=capacity_df_outlier.geometry, crs="EPSG:3857")
    ipm_geo_data_normal = ipm_geo_data_normal.explore(column='factor_IPM_over_PyPSA', cmap="jet", style_kwds={"weight": 5.0})
    ipm_geo_data_normal.save(pathlib.Path(plot_base_path, f"{model}_interactive_map_line_normal.html"))
    ipm_map_line_outlier = ipm_geo_data_outlier.explore(column='factor_IPM_over_PyPSA', cmap="Reds_r", style_kwds={"weight": 5.0})
    ipm_map_line_outlier.save(pathlib.Path(plot_base_path, f"{model}_interactive_map_line_outlier.html"))


def plot_network_capacity_reeds(pypsa_df, reeds_shapes_gdf, log_output_file, base_path, output_base_path, plot_base_path, model):

    transmission_capacities_path = pathlib.Path(base_path, "analysis", "gdrive_data", "data", "pypsa_usa", "transmission", "transmission_capacity_init_AC_ba_NARIS2024.csv")
    transmission_capacities_df = pd.read_csv(transmission_capacities_path).rename(columns={
        "r": "reeds_0",
        "rr": "reeds_1",
        "MW_f0": "export capacity (MW)",
        "MW_r0": "import capacity (MW)",
    })
    transmission_capacities_df["source"] = "reeds"
    transmission_capacities_df["Capacity (MW)"] = transmission_capacities_df[["export capacity (MW)", "import capacity (MW)"]].max(axis=1)
    transmission_capacities_df = transmission_capacities_df.loc[:, ("reeds_0", "reeds_1", "Capacity (MW)", "source")]
    transmission_capacities_df[["reeds_0", "reeds_1"]] = np.sort(transmission_capacities_df[["reeds_0", "reeds_1"]])
    transmission_capacities_df = transmission_capacities_df.drop_duplicates(keep="first")

    pypsa_df["s_nom_num_parallel"] = pypsa_df["s_nom"] * pypsa_df["num_parallel"]
    pypsa_df["s_nom_num_parallel_pu"] = pypsa_df["s_nom"] * pypsa_df["num_parallel"] * pypsa_df["s_max_pu"]
    pypsa_df[["reeds_0", "reeds_1"]] = np.sort(pypsa_df[["reeds_0", "reeds_1"]])
    reeds_pypsa_transmission_capacities = pypsa_df.query("reeds_0 != reeds_1").groupby(["reeds_0", "reeds_1"])["s_nom_num_parallel"].sum().reset_index().loc[:, ("reeds_0", "reeds_1", "s_nom_num_parallel")].rename(columns={"s_nom_num_parallel": "Capacity (MW)"})
    reeds_pypsa_transmission_capacities["source"] = "PyPSA"

    capacity_df = pd.concat(
        [reeds_pypsa_transmission_capacities, transmission_capacities_df]).reset_index()
    capacity_df = capacity_df.loc[:, ("reeds_0", "reeds_1", "Capacity (MW)", "source")]
    capacity_df = capacity_df.query("reeds_0 != 'nan' & reeds_1 != 'nan'")
    capacity_df = capacity_df.set_index(
        ["source", "reeds_0", "reeds_1"]
    ).unstack("source").droplevel(axis=1, level=0).reset_index()

    # remove values where IPM is zero or None
    capacity_df = capacity_df.query("reeds != 0.0 & PyPSA != 0.0 & ~reeds.isna() & ~PyPSA.isna()").reset_index()

    # compute percentage error
    capacity_df["Error wrt reeds (%)"] = (capacity_df["PyPSA"]-capacity_df["reeds"])/capacity_df["reeds"]*100.0
    capacity_df["Error wrt PyPSA (%)"] = (capacity_df["PyPSA"]-capacity_df["reeds"])/capacity_df["PyPSA"]*100.0
    capacity_df["factor_reeds_over_PyPSA"] = capacity_df["reeds"]/capacity_df["PyPSA"]

    capacity_df.to_csv(pathlib.Path(output_base_path, f"{model}_reeds_capacities.csv"), index=False)

    log_output_file.write("====")
    log_output_file.write("{}: median error wrt reeds: {} \n".format(model, np.round(capacity_df["Error wrt reeds (%)"].median(), 2)))
    log_output_file.write("{}: mean error wrt reeds: {} \n".format(model, np.round(capacity_df["Error wrt reeds (%)"].mean(), 2)))
    log_output_file.write("{}: median error wrt PyPSA: {} \n".format(model, np.round(capacity_df["Error wrt PyPSA (%)"].median(), 2)))
    log_output_file.write("{}: mean error wrt PyPSA: {} \n".format(model, np.round(capacity_df["Error wrt PyPSA (%)"].mean(), 2)))
    log_output_file.write("====")
    logger.info("====")
    logger.info("{}: median error wrt reeds: {} \n".format(model, np.round(capacity_df["Error wrt reeds (%)"].median(), 2)))
    logger.info("{}: mean error wrt reeds: {} \n".format(model, np.round(capacity_df["Error wrt reeds (%)"].mean(), 2)))
    logger.info("{}: median error wrt PyPSA: {} \n".format(model, np.round(capacity_df["Error wrt PyPSA (%)"].median(), 2)))
    logger.info("{}: mean error wrt PyPSA: {} \n".format(model, np.round(capacity_df["Error wrt PyPSA (%)"].mean(), 2)))
    logger.info("====")

    reeds_shapes_gdf = reeds_shapes_gdf.to_crs("3857")
    reeds_shapes_gdf["reeds_centroid"] = reeds_shapes_gdf.centroid
    capacity_df["reeds_0_centroid"] = capacity_df["reeds_0"].map(
        dict(reeds_shapes_gdf[['rb', 'reeds_centroid']].values))
    capacity_df["reeds_1_centroid"] = capacity_df["reeds_1"].map(
        dict(reeds_shapes_gdf[['rb', 'reeds_centroid']].values))
    capacity_df["geometry"] = capacity_df.apply(
        lambda x: spl.LineString([x.reeds_0_centroid, x.reeds_1_centroid]), axis=1)
    capacity_df = capacity_df.drop(["reeds_0_centroid", "reeds_1_centroid"], axis=1)
    capacity_df_normal = capacity_df.query("factor_reeds_over_PyPSA <= 1.0")
    capacity_df_outlier = capacity_df.query("factor_reeds_over_PyPSA > 1.0")
    reeds_geo_data_normal = gpd.GeoDataFrame(capacity_df_normal, geometry=capacity_df_normal.geometry, crs="EPSG:3857")
    reeds_geo_data_outlier = gpd.GeoDataFrame(capacity_df_outlier, geometry=capacity_df_outlier.geometry, crs="EPSG:3857")
    reeds_geo_data_normal = reeds_geo_data_normal.explore(column='factor_reeds_over_PyPSA', cmap="jet", style_kwds={"weight": 5.0})
    reeds_geo_data_normal.save(pathlib.Path(plot_base_path, f"{model}_reeds_interactive_map_line_normal.html"))
    reeds_map_line_outlier = reeds_geo_data_outlier.explore(column='factor_reeds_over_PyPSA', cmap="Reds_r", style_kwds={"weight": 5.0})
    reeds_map_line_outlier.save(pathlib.Path(plot_base_path, f"{model}_reeds_interactive_map_line_outlier.html"))


def place_line_boundaries(lines_dataframe, gadm_dataframe, ipm_dataframe, reeds_dataframe, log_output_file, id_column_name, lines_dataframe_name, network_used="other"):
    """
    The function spatially joins the boundaries of transmission line with the:
    -) GADM (level 1) shape files
    -) IPM regions shape files
    -) ReEDS zones shape files

    Returns
    The function returns a modified lines_dataframe with additional six extra columns:
    -) state_0: US state where the bus_0 of the transmission line is located
    -) ipm_region_0: IPM region where the bus_0 of the transmission line is located
    -) reeds_0: ReEDS region where the bus_0 of the transmission line is located
    -) state_1: US state where the bus_1 of the transmission line is located
    -) ipm_region_1: IPM region where the bus_1 of the transmission line is located
    -) reeds_1: ReEDS region where the bus_1 of the transmission line is located
    """

    # Spatially join Bus 0 with the GADM and IPM shapes.
    if network_used == "pypsa_earth":
        lines_dataframe_modified = gpd.GeoDataFrame(lines_dataframe, geometry=gpd.GeoSeries.from_wkt(lines_dataframe.bus_0_coors), crs="EPSG:4326").reset_index()
    else:
        lines_dataframe_modified = lines_dataframe.loc[:, (id_column_name, "sub_0_coors")]
        lines_dataframe_modified["geometry"] = lines_dataframe_modified["sub_0_coors"]
    log_output_file.write(" --> shape of {} before sub_0 spatial join {} \n".format(lines_dataframe_name, lines_dataframe.shape))
    logger.info(" --> shape of {} before sub_0 spatial join {} \n".format(lines_dataframe_name, lines_dataframe.shape))
    spatial_join_gadm_sub_0 = lines_dataframe_modified.sjoin(gadm_dataframe, how="left").loc[:, (id_column_name, "ISO_1")].rename(columns={"ISO_1": "state_0"})
    spatial_join_ipm_sub_0 = lines_dataframe_modified.sjoin(ipm_dataframe, how="left").loc[:, (id_column_name, "IPM_Region")].rename(columns={"IPM_Region": "ipm_region_0"})
    spatial_join_reeds_sub_0 = lines_dataframe_modified.sjoin(reeds_dataframe, how="left").loc[:, (id_column_name, "rb")].rename(columns={"rb": "reeds_0"})
    spatial_join_gadm_ipm_sub_0 = pd.merge(spatial_join_gadm_sub_0, spatial_join_ipm_sub_0, how="inner", on=id_column_name)
    spatial_join_sub_0 = pd.merge(spatial_join_gadm_ipm_sub_0, spatial_join_reeds_sub_0, how="inner", on=id_column_name)
    log_output_file.write(" --> shape of {} after sub_0 spatial join with gadm {} \n".format(lines_dataframe_name, spatial_join_gadm_sub_0.shape))
    log_output_file.write(" --> shape of {} after sub_0 spatial join with ipm{} \n".format(lines_dataframe_name, spatial_join_ipm_sub_0.shape))
    log_output_file.write(" --> shape of {} after sub_0 spatial join {} \n".format(lines_dataframe_name, spatial_join_sub_0.shape))
    logger.info(" --> shape of {} after sub_0 spatial join with gadm {} \n".format(lines_dataframe_name, spatial_join_gadm_sub_0.shape))
    logger.info(" --> shape of {} after sub_0 spatial join with ipm{} \n".format(lines_dataframe_name, spatial_join_ipm_sub_0.shape))
    logger.info(" --> shape of {} after sub_0 spatial join {} \n".format(lines_dataframe_name, spatial_join_sub_0.shape))

    # Spatially join Bus 1 with the GADM and IPM shapes.
    if network_used == "pypsa_earth":
        lines_dataframe_modified = gpd.GeoDataFrame(lines_dataframe, geometry=gpd.GeoSeries.from_wkt(lines_dataframe.bus_1_coors), crs="EPSG:4326").reset_index()
    else:
        lines_dataframe_modified = lines_dataframe.loc[:, (id_column_name, "sub_1_coors")]
        lines_dataframe_modified["geometry"] = lines_dataframe_modified["sub_1_coors"]
    log_output_file.write(" --> shape of {} before sub_1 spatial join {} \n".format(lines_dataframe_name, lines_dataframe_modified.shape))
    logger.info(" --> shape of {} before sub_1 spatial join {} \n".format(lines_dataframe_name, lines_dataframe_modified.shape))
    spatial_join_gadm_sub_1 = lines_dataframe_modified.sjoin(gadm_dataframe, how="left").loc[:, (id_column_name, "ISO_1")].rename(columns={"ISO_1": "state_1"})
    spatial_join_ipm_sub_1 = lines_dataframe_modified.sjoin(ipm_dataframe, how="left").loc[:, (id_column_name, "IPM_Region")].rename(columns={"IPM_Region": "ipm_region_1"})
    spatial_join_reeds_sub_1 = lines_dataframe_modified.sjoin(reeds_dataframe, how="left").loc[:, (id_column_name, "rb")].rename(columns={"rb": "reeds_1"})
    spatial_join_gadm_ipm_sub_1 = pd.merge(spatial_join_gadm_sub_1, spatial_join_ipm_sub_1, how="inner", on=id_column_name)
    spatial_join_sub_1 = pd.merge(spatial_join_gadm_ipm_sub_1, spatial_join_reeds_sub_1, how="inner", on=id_column_name)
    log_output_file.write(" --> shape of {} after sub_1 spatial join with gadm {} \n".format(lines_dataframe_name, spatial_join_gadm_sub_1.shape))
    log_output_file.write(" --> shape of {} after sub_1 spatial join with ipm{} \n".format(lines_dataframe_name, spatial_join_ipm_sub_1.shape))
    log_output_file.write(" --> shape of {} after sub_1 spatial join {} \n".format(lines_dataframe_name, spatial_join_sub_1.shape))
    logger.info(" --> shape of {} after sub_1 spatial join with gadm {} \n".format(lines_dataframe_name, spatial_join_gadm_sub_1.shape))
    logger.info(" --> shape of {} after sub_1 spatial join with ipm{} \n".format(lines_dataframe_name, spatial_join_ipm_sub_1.shape))
    logger.info(" --> shape of {} after sub_1 spatial join {} \n".format(lines_dataframe_name, spatial_join_sub_1.shape))

    # --> Inner join the results
    lines_dataframe = pd.merge(lines_dataframe, spatial_join_sub_0, how="inner", on=id_column_name)
    lines_dataframe = pd.merge(lines_dataframe, spatial_join_sub_1, how="inner", on=id_column_name)
    lines_dataframe["state_0"] = lines_dataframe["state_0"].astype(str)
    lines_dataframe["state_1"] = lines_dataframe["state_1"].astype(str)
    lines_dataframe["ipm_region_0"] = lines_dataframe["ipm_region_0"].astype(str)
    lines_dataframe["ipm_region_1"] = lines_dataframe["ipm_region_1"].astype(str)
    lines_dataframe["reeds_0"] = lines_dataframe["reeds_0"].astype(str)
    lines_dataframe["reeds_1"] = lines_dataframe["reeds_1"].astype(str)
    log_output_file.write(" --> shape of {} after the inner joins {} \n".format(lines_dataframe_name, lines_dataframe.shape))
    logger.info(" --> shape of {} after the inner joins {} \n".format(lines_dataframe_name, lines_dataframe.shape))

    return lines_dataframe


def parse_inputs(base_path, log_output_file):
    """
    The function parses the necessary inputs for the analysis
    """
    base_network_pypsa_earth_path = pathlib.Path(base_path, snakemake.input.base_network_pypsa_earth_path)
    base_network_pypsa_usa_path = pathlib.Path(base_path, snakemake.input.base_network_pypsa_usa_path)
    lines_osm_raw_path = pathlib.Path(base_path, snakemake.input.lines_osm_raw_path)
    lines_osm_clean_path = pathlib.Path(base_path, snakemake.input.lines_osm_clean_path)
    eia_base_network_path = pathlib.Path(base_path, snakemake.input.eia_base_network_path)
    gadm_shapes_path = pathlib.Path(base_path, snakemake.input.gadm_shapes_path)
    ipm_shapes_path = pathlib.Path(base_path, snakemake.input.ipm_shapes_path)
    reeds_shapes_path = pathlib.Path(base_path, snakemake.input.reeds_shapes_path)

    #############
    # Load data #
    #############
    base_network_pypsa_earth = pypsa.Network(base_network_pypsa_earth_path)
    base_network_pypsa_usa = pd.read_csv(base_network_pypsa_usa_path)
    base_network_pypsa_usa.loc[:, "Line"] = base_network_pypsa_usa.loc[:, "Line"].astype(str)
    base_network_pypsa_usa.loc[:, "bus1"] = base_network_pypsa_usa.loc[:, "bus1"].astype(str)
    lines_osm_raw = gpd.read_file(lines_osm_raw_path)
    lines_osm_clean = gpd.read_file(lines_osm_clean_path)
    eia_base_network = gpd.read_file(eia_base_network_path)
    gadm_shapes = gpd.read_file(gadm_shapes_path)
    ipm_shapes = gpd.read_file(ipm_shapes_path).to_crs("4326")
    reeds_shapes = gpd.read_file(reeds_shapes_path).to_crs("4326")

    ############
    # EIA data #
    ############
    log_output_file.write("        \n")
    log_output_file.write("        \n")
    log_output_file.write(" Data preparation on the EIA base network \n")
    log_output_file.write(" --> shape of eia_base_network after reading it in {} \n".format(eia_base_network.shape))
    logger.info("        \n")
    logger.info("        \n")
    logger.info(" Data preparation on the EIA base network \n")
    logger.info(" --> shape of eia_base_network after reading it in {} \n".format(eia_base_network.shape))

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
    logger.info(" --> shape of eia_base_network after excluding multilinestrings {} \n".format(eia_base_network.shape))

    # --> Compute the start- and end-points of a line
    eia_base_network[["sub_0_coors", "sub_1_coors"]] = eia_base_network["geometry"].boundary.explode(index_parts=True).unstack()
    log_output_file.write(" --> shape of eia_base_network after computing boundaries {} \n".format(eia_base_network.shape))
    logger.info(" --> shape of eia_base_network after computing boundaries {} \n".format(eia_base_network.shape))

    # --> Determine where the start- and end-points of the line are located. In particular, we perform the spatial
    # joins with:
    #  -) the GADM shapes (level 1) to get the US state
    #  -) the IPM shapes to get the IPM region

    eia_base_network = place_line_boundaries(eia_base_network, gadm_shapes, ipm_shapes, reeds_shapes, log_output_file, "OBJECTID_1", "eia_base_network")

    # Clean the EIA data from lines with unnecessary voltages and voltage classes
    eia_base_network = eia_base_network.rename(columns={"VOLTAGE": "v_nom", "VOLT_CLASS": "v_nom_class"})

    # --> Remove lines corresponding to voltage = -999999.0 kV
    eia_base_network = eia_base_network.loc[eia_base_network["v_nom"] != -999999.0]
    log_output_file.write(" --> shape of eia_base_network after removing the lines with voltage -999999.0 kV {} \n".format(eia_base_network.shape))
    logger.info(" --> shape of eia_base_network after removing the lines with voltage -999999.0 kV {} \n".format(eia_base_network.shape))

    # --> Remove lines corresponding to voltage class 'DC'. All lines in the base.nc are AC
    eia_base_network = eia_base_network.loc[eia_base_network["v_nom_class"] != 'Dc']
    log_output_file.write(" --> shape of eia_base_network after removing the lines with voltage class 'Dc' {} \n".format(eia_base_network.shape))
    logger.info(" --> shape of eia_base_network after removing the lines with voltage class 'Dc' {} \n".format(eia_base_network.shape))

    # --> Remove lines corresponding to voltage class 'Not Available'
    eia_base_network = eia_base_network.loc[eia_base_network["v_nom_class"] != 'Not Available']
    log_output_file.write(" --> shape of eia_base_network after removing the lines with voltage class 'Not Available' {} \n".format(eia_base_network.shape))
    logger.info(" --> shape of eia_base_network after removing the lines with voltage class 'Not Available' {} \n".format(eia_base_network.shape))

    #################
    # OSM lines raw #
    #################

    log_output_file.write("        \n")
    log_output_file.write("        \n")
    log_output_file.write(" Data preparation on the OSM lines raw \n")
    logger.info("        \n")
    logger.info("        \n")
    logger.info(" Data preparation on the OSM lines raw \n")

    lines_osm_raw[["sub_0_coors", "sub_1_coors"]] = lines_osm_raw["geometry"].boundary.explode(index_parts=True).unstack()
    log_output_file.write(" --> shape of lines_osm_raw after computing boundaries {} \n".format(lines_osm_raw.shape))
    logger.info(" --> shape of lines_osm_raw after computing boundaries {} \n".format(lines_osm_raw.shape))

    # --> Determine where the start- and end-points of the line are located. In particular, we perform the spatial
    # joins with:
    #  -) the GADM shapes (level 1) to get the US state
    #  -) the IPM shapes to get the IPM region

    lines_osm_raw = place_line_boundaries(lines_osm_raw, gadm_shapes, ipm_shapes, reeds_shapes, log_output_file, "id", "lines_osm_raw")

    ###################
    # OSM lines clean #
    ###################

    log_output_file.write("        \n")
    log_output_file.write("        \n")
    log_output_file.write(" Data preparation on the OSM lines clean \n")
    logger.info("        \n")
    logger.info("        \n")
    logger.info(" Data preparation on the OSM lines clean \n")

    lines_osm_clean[["sub_0_coors", "sub_1_coors"]] = lines_osm_clean["geometry"].boundary.explode(index_parts=True).unstack()
    log_output_file.write(" --> shape of lines_osm_clean after computing boundaries {} \n".format(lines_osm_clean.shape))
    logger.info(" --> shape of lines_osm_clean after computing boundaries {} \n".format(lines_osm_clean.shape))

    # --> Determine where the start- and end-points of the line are located. In particular, we perform the spatial
    # joins with:
    #  -) the GADM shapes (level 1) to get the US state
    #  -) the IPM shapes to get the IPM region

    lines_osm_clean = place_line_boundaries(lines_osm_clean, gadm_shapes, ipm_shapes, reeds_shapes, log_output_file, "line_id", "lines_osm_clean")

    #######################
    # PyPSA-Earth base.nc #
    #######################

    log_output_file.write("        \n")
    log_output_file.write("        \n")
    log_output_file.write(" Data preparation on the PyPSA-Earth base network \n")
    logger.info("        \n")
    logger.info("        \n")
    logger.info(" Data preparation on the PyPSA-Earth base network \n")

    # --> Determine where the start- and end-points of the lines are located. In particular, we perform the spatial
    # joins with:
    #  -) the GADM shapes (level 1) to get the US state
    #  -) the IPM shapes to get the IPM region
    log_output_file.write(" --> shape of pypsa-earth base network after reading it in {} \n".format(base_network_pypsa_earth.lines.shape))
    logger.info(" --> shape of pypsa-earth base network after reading it in {} \n".format(base_network_pypsa_earth.lines.shape))

    base_network_pypsa_earth.lines = place_line_boundaries(base_network_pypsa_earth.lines, gadm_shapes, ipm_shapes, reeds_shapes, log_output_file, "Line", "base.nc", "pypsa_earth")

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

    ###########################
    # PyPSA-USA lines_gis.csv #
    ###########################

    log_output_file.write("        \n")
    log_output_file.write("        \n")
    log_output_file.write(" Data preparation on the PyPSA-USA lines_gis \n")
    logger.info("        \n")
    logger.info("        \n")
    logger.info(" Data preparation on the PyPSA-USA lines_gis \n")

    base_network_pypsa_usa = gpd.GeoDataFrame(base_network_pypsa_usa, geometry=gpd.GeoSeries.from_wkt(base_network_pypsa_usa["WKT_geometry"]), crs="EPSG:4326").reset_index()
    base_network_pypsa_usa[["sub_0_coors", "sub_1_coors"]] = base_network_pypsa_usa["geometry"].boundary.explode(index_parts=True).unstack()
    log_output_file.write(" --> shape of pypsa-usa lines_gis after computing boundaries {} \n".format(base_network_pypsa_usa.shape))
    logger.info(" --> shape of pypsa-usa lines_gis after computing boundaries {} \n".format(base_network_pypsa_usa.shape))

    # --> Determine where the start- and end-points of the lines are located. In particular, we perform the spatial
    # joins with:
    #  -) the GADM shapes (level 1) to get the US state
    #  -) the IPM shapes to get the IPM region
    log_output_file.write(" --> shape of pypsa-usa lines_gis after reading it in {} \n".format(base_network_pypsa_usa.shape))
    logger.info(" --> shape of pypsa-usa lines_gis after reading it in {} \n".format(base_network_pypsa_usa.shape))

    base_network_pypsa_usa = place_line_boundaries(base_network_pypsa_usa, gadm_shapes, ipm_shapes, reeds_shapes, log_output_file, "Line", "base.nc")

    # --> Assign a voltage class to the pypsa-usa lines_gis.csv
    base_network_pypsa_usa["v_nom_class"] = base_network_pypsa_usa["v_nom"]

    v_nom_class_dict_pypsa_usa = {
        69.: 'Under 100',
        100.: "100-161",
        115.: "100-161",
        138.: "100-161",
        161.: "100-161",
        230.: "220-287",
        345.: "345",
        500.: "500",
        765.: "735 And Above"
    }

    base_network_pypsa_usa["v_nom_class"] = base_network_pypsa_usa["v_nom_class"].replace(v_nom_class_dict_pypsa_usa)

    return eia_base_network, base_network_pypsa_earth, base_network_pypsa_usa, lines_osm_raw, lines_osm_clean, ipm_shapes, reeds_shapes


if __name__ == '__main__':

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs")
    plot_path = pathlib.Path(default_path, "analysis", "plots")
    output_path = pathlib.Path(default_path, "analysis", "outputs")
    ccs_color_dict = {"EIA": "#FF8C00", "PyPSA": "#0000FF", "PyPSA_parallel": "#228B22", "delta_PyPSA": "#0000FF", "delta_PyPSA_parallel": "#228B22", "Error (%)": "#FF7F50"}
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(log_path, f"output_network_comparison_{today_date[:10]}.txt")
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    # parse the input files
    network_eia_df, network_pypsa_earth_df, network_pypsa_usa_df, osm_lines_raw, osm_lines_clean, ipm_region_shapes, reeds_network_shapes = parse_inputs(default_path, log_output_file)

    # output dataframes after pre-processing
    network_pypsa_earth_df.export_to_netcdf(pathlib.Path(output_path, "modified_pypsa_earth_base.nc"))
    network_pypsa_usa_df.to_csv(pathlib.Path(output_path, "modified_pypsa_usa_base.csv"))
    network_eia_df.to_csv(pathlib.Path(output_path, "modified_eia_base.csv"), index=False)
    osm_lines_raw.to_csv(pathlib.Path(output_path, "modified_osm_lines_raw.csv"), index=False)
    osm_lines_clean.to_csv(pathlib.Path(output_path, "modified_osm_lines_clean.csv"), index=False)

    eia_voltage_classes = list(network_eia_df["v_nom_class"].unique())

    if snakemake.params.plot_network_topology:
        for selected_voltage_class in eia_voltage_classes:
            fig_name_map = pathlib.Path(plot_path, "network_comparison_pearth_for_voltage_class_{}.png".format(
                str(selected_voltage_class)))
            plot_network_topology_comparison(network_pypsa_earth_df, network_eia_df, selected_voltage_class,
                                    "PyPSA-Earth base network", fig_name_map)
            fig_name_intersection = pathlib.Path(plot_path, "network_comparison_intersection_{}.png".format(
                str(selected_voltage_class)))
            plot_network_topology_intersection(network_pypsa_earth_df, network_eia_df, selected_voltage_class,
                                      fig_name_intersection)

    # Comparison for the transmission crossings (PyPSA-Earth vs EIA) using:
    # -) the GADM shapes(level 1) for the US state comparison
    # -) the IPM shapes for the IPM region comparison
    if snakemake.params.plot_network_crossings:
        plot_network_crossings(network_pypsa_earth_df, network_eia_df, ccs_color_dict, eia_voltage_classes, output_path, plot_path)

    # Comparison for the transmission capacities (PyPSA-Earth/PyPSA-USA vs IPM transmission capacities)
    if snakemake.params.plot_network_capacity_ipm:
        plot_network_capacity_ipm(network_pypsa_earth_df.lines, ipm_region_shapes, ccs_color_dict, log_output_file, default_path, output_path, plot_path, "pypsa_earth")
        plot_network_capacity_ipm(network_pypsa_usa_df, ipm_region_shapes, ccs_color_dict, log_output_file, default_path, output_path, plot_path, "pypsa_usa")

    # Comparison for the transmission capacities (PyPSA-Earth/PyPSA-USA vs reeds transmission capacities)
    if snakemake.params.plot_network_capacity_reeds:
        plot_network_capacity_reeds(network_pypsa_earth_df.lines, reeds_network_shapes, log_output_file, default_path, output_path, plot_path, "pypsa_earth")

    log_output_file.close()
