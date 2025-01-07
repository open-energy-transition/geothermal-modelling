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


def parse_inputs(base_path):
    """
    The function parses the necessary inputs for the analysis
    """
    network_pypsa_earth_path = pathlib.Path(base_path, snakemake.input.pypsa_earth_network_path)
    eia_installed_capacity_reference_path = pathlib.Path(base_path, snakemake.input.eia_installed_capacity_path)
    gadm_shapes_path = pathlib.Path(base_path, snakemake.input.gadm_shapes_path)

    # load the data
    network_pypsa_earth = pypsa.Network(network_pypsa_earth_path)
    eia_installed_capacity_reference = pd.read_excel(eia_installed_capacity_reference_path, skiprows=1, index_col="Energy Source")

    return network_pypsa_earth, eia_installed_capacity_reference, gadm_shapes_path


def plot_state_by_state_comparison(pypsa_network, eia_reference, year_to_use, log_file, plot_base_path, gadm_shapes_path):
    """
    The function plots the state-by-state comparison between the EIA reference installed capacity data and the PyPSA-Earth network
    """

    log_file.write("        \n")
    log_file.write("        \n")
    log_file.write("Compare the state-by-state installed capacity \n")


def plot_country_comparison(pypsa_network, eia_reference, year_to_use, log_file, plot_base_path):
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
    log_path = pathlib.Path(default_path, "analysis", "logs")
    plot_path = pathlib.Path(default_path, "analysis", "plots")
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(log_path, f"output_generation_comparison_{today_date[:10]}.txt")
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    # initial configurations
    eia_name = "EIA"
    pypsa_name = "PyPSA"
    year_for_comparison = snakemake.params.year_for_comparison

    network_pypsa_earth_df, eia_installed_capacity_df, gadm_path = parse_inputs(default_path)

    if snakemake.params.plot_country_comparison:
        plot_country_comparison(network_pypsa_earth_df, eia_installed_capacity_df, year_for_comparison, log_output_file, plot_path)

    if snakemake.params.plot_state_by_state_comparison:
        plot_state_by_state_comparison(network_pypsa_earth_df, eia_installed_capacity_df, year_for_comparison, log_output_file, plot_path, gadm_path)

    log_output_file.close()
