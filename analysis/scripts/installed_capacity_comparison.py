# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import pathlib
import datetime as dt
import pandas as pd
import pypsa
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib.patches import Patch
from _helpers_usa import get_state_node, config


def parse_inputs(base_path, alternative_clustering):
    """
    The function parses the necessary inputs for the analysis
    """
    if alternative_clustering:
        network_pypsa_earth_path = pathlib.Path(base_path, snakemake.input.pypsa_earth_network_path)
    else:
        network_pypsa_earth_path = pathlib.Path(base_path, snakemake.input.pypsa_earth_network_nonac_path)
    eia_installed_capacity_reference_path = pathlib.Path(base_path, snakemake.input.eia_installed_capacity_path)
    eia_state_temporal_installed_capacity_path = pathlib.Path(base_path, snakemake.input.eia_state_temporal_installed_capacity_path)
    gadm_shapes_path = pathlib.Path(base_path, snakemake.input.gadm_shapes_path)

    # load the data
    network_pypsa_earth = pypsa.Network(network_pypsa_earth_path)
    eia_installed_capacity_reference = pd.read_excel(eia_installed_capacity_reference_path, skiprows=1, index_col="Energy Source")
    eia_state_temporal_installed_capacity_reference = pd.read_excel(eia_state_temporal_installed_capacity_path, skiprows=1)
        
    return network_pypsa_earth, eia_installed_capacity_reference, eia_state_temporal_installed_capacity_reference, gadm_shapes_path


def add_legend_patches(ax, colors, labels, legend_kw):
    patches = [Patch(color=color, label=label) for color, label in zip(colors, labels)]
    ax.legend(handles=patches, **legend_kw)


def plot_capacity_spatial_representation(pypsa_network, plot_type, state_to_drop, plot_base_path, gadm_base_path):

    iso_code_to_omit = []
    for state in state_to_drop:
        iso_code_to_omit.append(get_state_node(gadm_base_path, state) + "_AC")

    carriers = pypsa_network.carriers.index.tolist()[:-2]
    buses = pypsa_network.buses.query('carrier == "AC" and ~(index in @iso_code_to_omit)')
    bus_index = buses.index.tolist()
    fig = plt.figure(figsize=(15, 12))
    ax = plt.axes(projection=ccrs.EqualEarth())
    pypsa_network.lines.loc[pypsa_network.lines.bus0.isin(iso_code_to_omit), "s_nom"] = 0
    pypsa_network.lines.loc[pypsa_network.lines.bus1.isin(iso_code_to_omit), "s_nom"] = 0

    if plot_type == "capacity":
        cap = pypsa_network.generators.query("carrier != 'load' and bus in @b", local_dict={"b": bus_index}).groupby(["bus", "carrier"]).p_nom.sum()
        cap = cap[~(cap.index.get_level_values('bus').isin(iso_code_to_omit))]

        # Filter carriers with positive capacity, and avoid KeyError by checking presence in index
        carriers_with_capacity = []
        for carrier in carriers:
            if carrier in cap.index.get_level_values('carrier'):
                total_capacity = cap.xs(carrier, level='carrier').sum()
                if total_capacity > 0:
                    carriers_with_capacity.append(carrier)

        pypsa_network.plot(
            ax=ax,
            bus_sizes=cap / 2e4,
            line_widths=pypsa_network.lines.s_nom / 1e4,
            line_colors='rosybrown',
            margin=0.25,
            bus_alpha=0.8,
            color_geomap=True,
            link_alpha=0,
            title="Spatial representation of the installed capacity (GW)"
        )

        # Add the legend with filtered carriers
        add_legend_patches(
            ax,
            colors=[pypsa_network.carriers.loc[carrier].color for carrier in carriers_with_capacity],
            labels=carriers_with_capacity,
            legend_kw=dict(frameon=False, bbox_to_anchor=(0, 1), title="Carriers")
        )

        plt.savefig(pathlib.Path(plot_base_path,  f"installed_capacity_spatial_representation.png"), dpi=800)

def rename_carrier(x):
    if x == 'ccgt':
        return 'CCGT'
    elif x == 'phs':
        return 'PHS'
    else:
        return x

def plot_capacity_state_by_state_comparison(pypsa_network, eia_reference, year_to_use, log_file, plot_base_path, gadm_shapes_path):
    """
    The function plots the state-by-state comparison between the EIA reference installed capacity data and the PyPSA-Earth network
    """

    log_file.write("        \n")
    log_file.write("        \n")
    log_file.write("Compare the state-by-state installed capacity \n")

    eia_installed_capacity_by_state_year = eia_reference.loc[eia_reference["Year"] == year_to_use]
    eia_installed_capacity_by_state_year = eia_installed_capacity_by_state_year[
        ["Year", "State Code", "Fuel Source", "Nameplate Capacity (Megawatts)", "Producer Type"]]
    rename_cols = {"Year": "year", "State Code": "state", "Fuel Source": "carrier",
                   "Nameplate Capacity (Megawatts)": "installed_capacity"}
    eia_installed_capacity_by_state_year = eia_installed_capacity_by_state_year.rename(columns=rename_cols)
    eia_installed_capacity_by_state_year = eia_installed_capacity_by_state_year.replace({"carrier": {
        "Hydroelectric": "hydro", "Solar Thermal and Photovoltaic": "solar",
        "Natural Gas": "gas", "Petroleum": "oil", "Wind": "wind", "Nuclear": "nuclear", "Geothermal": "geothermal",
        "Pumped Storage": "PHS", "Wood and Wood Derived Fuels": "biomass"}})
    eia_installed_capacity_by_state_year = eia_installed_capacity_by_state_year.loc[
        eia_installed_capacity_by_state_year["Producer Type"] == 'Total Electric Power Industry']
    eia_installed_capacity_by_state_year = eia_installed_capacity_by_state_year.loc[
        eia_installed_capacity_by_state_year['state'] != 'US']
 
    series_to_use = eia_installed_capacity_by_state_year.groupby(["carrier", "state"]).installed_capacity.sum()
    df = pd.DataFrame({"carrier_state": series_to_use.index, "installed_capacity": series_to_use.values})
    df[["carrier", "state"]] = df["carrier_state"].apply(lambda x: pd.Series(x))
    eia_installed_capacity_by_state_year = df[["carrier", "state", "installed_capacity"]].copy()
    eia_installed_capacity_by_state_year.carrier = eia_installed_capacity_by_state_year.carrier.str.lower()

    eia_installed_capacity_by_state_year['carrier'] = eia_installed_capacity_by_state_year['carrier'].apply(lambda x: rename_carrier(x))

    # Installed capacity groupby and summed by carrier and bus
    series_gen_to_use = pypsa_network.generators.groupby(["carrier", "bus"]).p_nom.sum()
    series_sto_to_use = pypsa_network.storage_units.groupby(["carrier","bus"]).p_nom.sum()
    series_to_use = series_gen_to_use._append(series_sto_to_use)
    df = pd.DataFrame({"carrier_gid": series_to_use.index, "installed_capacity": series_to_use.values})

def plot_capacity_country_comparison(pypsa_network, eia_reference, year_to_use, log_file, plot_base_path):
    """
    The function plots the countrywide comparison between the EIA reference generation data and the PyPSA-Earth network
    """

    log_file.write("        \n")
    log_file.write("        \n")
    log_file.write("Compare the countrywide installed capacity \n")

    # prepare the PyPSA results
    df_pypsa_hydro_phs_capacity = pypsa_network.storage_units.groupby("carrier").p_nom_opt.sum()
    df_pypsa_capacity = pd.concat(
        [pypsa_network.generators.groupby("carrier").p_nom_opt.sum(), df_pypsa_hydro_phs_capacity])
    df_pypsa_capacity.loc["wind"] = df_pypsa_capacity.loc[["offwind-ac", "offwind-dc", "onwind"]].sum()
    df_pypsa_capacity = df_pypsa_capacity.drop(["offwind-ac", "offwind-dc", "onwind"])
    df_pypsa_capacity /= 1000
    df_pypsa_capacity = df_pypsa_capacity.round(2)
    df_pypsa_capacity.name = pypsa_name

    log_output_file.write("Installed capacity: {} \n".format(df_pypsa_capacity))

    # prepare the EIA reference data
    eia_reference.index = eia_reference.index.str.lower()
    eia_reference.loc["other biomass"] = eia_reference.loc[["other biomass", "wood and wood-derived fuels"]].sum()
    df_eia_capacity = eia_reference.rename(
        index={"hydroelectric conventional": "hydro", "hydroelectric pumped storage": "PHS",
               "solar photovoltaic": "solar", "natural gas": "CCGT", "petroleum": "oil", "other biomass": "biomass"})
    df_eia_capacity = df_eia_capacity.drop(
        ["estimated total solar", "solar thermal", "wood and wood-derived fuels", "other energy sources", "total",
         "small scale photovoltaic", "estimated total photovoltaic", "other gases", ])
    df_eia_capacity = df_eia_capacity.iloc[:-1]
    df_eia_capacity.loc["solar", "Generator Nameplate Capacity"] = df_eia_capacity.loc["solar", "Net Summer Capacity"]
    df_eia_capacity = df_eia_capacity["Generator Nameplate Capacity"]
    df_eia_capacity.name = eia_name
    df_eia_capacity /= 1000

    # comparison dataframe
    df_compare_capacity = pd.concat([df_pypsa_capacity, df_eia_capacity], axis=1)

    df_compare_capacity.plot(kind="bar")
    plt.xlabel("Carriers")
    plt.ylabel("Installed Capacity (GW)")
    plt.title(f"Year = {year_to_use}")
    plt.grid(linestyle="--")
    plt.subplots_adjust(bottom=0.3)
    plt.savefig(pathlib.Path(plot_base_path, f"countrywide_installed_capacity.png"), dpi=800)


if __name__ == '__main__':

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "installed_capacity")
    plot_path = pathlib.Path(default_path, "analysis", "plots", "installed_capacity")
    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(plot_path).mkdir(parents=True, exist_ok=True)
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(log_path, f"output_generation_comparison_{today_date[:10]}.txt")
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")
    config_path = pathlib.Path(default_path, "configs","config.usa_PE.yaml")

    # initial configurations
    eia_name = "EIA"
    pypsa_name = "PyPSA"
    year_for_comparison = snakemake.params.year_for_comparison
    config_dict = config(config_path)
    alternative_clustering = config_dict["cluster_options"]["alternative_clustering"]

    network_pypsa_earth_df, eia_installed_capacity_df, eia_state_temporal_installed_capacity_df, gadm_path = parse_inputs(default_path, alternative_clustering)

    if snakemake.params.plot_country_comparison:
        plot_capacity_country_comparison(network_pypsa_earth_df, eia_installed_capacity_df, year_for_comparison, log_output_file, plot_path)

    if snakemake.params.plot_state_by_state_comparison:
        plot_capacity_state_by_state_comparison(network_pypsa_earth_df, eia_state_temporal_installed_capacity_df, year_for_comparison, log_output_file, plot_path, gadm_path)

    if snakemake.params.plot_spatial_representation:
        plot_capacity_spatial_representation(network_pypsa_earth_df, "capacity", snakemake.params.state_to_omit, plot_path, gadm_path)

    log_output_file.close()
