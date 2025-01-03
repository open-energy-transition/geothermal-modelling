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
import plotly.express as px
import matplotlib.pyplot as plt
from _helpers_usa import extract_time_res, get_gadm_mapping


def parse_inputs(base_path):
    """
    The function parses the necessary inputs for the analysis
    """
    network_pypsa_earth_path = pathlib.Path(base_path, snakemake.input.pypsa_earth_network_path)
    eia_country_generation_reference_path = pathlib.Path(base_path, snakemake.input.eia_country_generation_path)
    eia_state_generation_reference_path = pathlib.Path(base_path, snakemake.input.eia_state_generation_path)
    gadm_shapes_path = pathlib.Path(base_path, snakemake.input.gadm_shapes_path)

    # extract the time resolution from the PyPSA-Earth network file name
    extracted_time_resolution = extract_time_res(network_pypsa_earth_path)

    # load the data
    network_pypsa_earth = pypsa.Network(network_pypsa_earth_path)
    eia_country_generation_reference = pd.read_csv(eia_country_generation_reference_path, index_col="Unnamed: 0")
    eia_state_generation_reference = pd.read_excel(eia_state_generation_reference_path, sheet_name="Data")

    return network_pypsa_earth, eia_country_generation_reference, eia_state_generation_reference, gadm_shapes_path, extracted_time_resolution


def plot_state_by_state_comparison(pypsa_network, eia_reference, year_to_use, time_resolution, log_file, plot_base_path, gadm_shapes_path):
    """
    The function plots the state-by-state comparison between the EIA reference generation data and the PyPSA-Earth network
    """

    log_file.write("        \n")
    log_file.write("        \n")
    log_file.write("Compare the state-by-state electricity generation \n")

    year_list = [val for val in list(eia_reference.columns) if val not in ["Data_Status", "State", "MSN"]]
    if year_to_use not in year_list:
        raise Exception("{} is not allowed. Please use choose a year from the list: {}".format(year_to_use, year_list))

    # select the relevant year and relevant rows
    keywords_list = ["WYTCP", "NUETP", "HYTCP", "GEEGP"]
    eia_generation_reference_filtered = eia_reference[eia_reference["MSN"].isin(keywords_list)][["State", "MSN", year_to_use]]
    eia_generation_reference_filtered = eia_generation_reference_filtered.rename(columns={'MSN': 'carrier'})
    eia_generation_reference_filtered.carrier = eia_generation_reference_filtered.carrier.map({"WYTCP": "onwind",
                                                                                    "NUETP": "nuclear",
                                                                                    "HYTCP": "hydro",
                                                                                    "GEEGP": "geothermal"})

    eia_generation_reference_filtered.set_index(['carrier', 'State'] ,inplace=True)
    eia_generation_reference_filtered[year_to_use] /= 1e3 #conversion from million kWh i.e., GWh to TWh

    # calculate the annual generation from each of the generators in MWh
    log_file.write(" --> Compute the state-by-state annual electricity generation \n")
    pypsa_network.generators = pypsa_network.generators.assign(gen_MWh=pypsa_network.generators_t.p.sum() * time_resolution)

    # PyPSA-Earth mapping the generators bus to the state names
    log_file.write(" --> Determine the PyPSA-Earth mapping the generators bus to the state names \n")
    usa_state_dict = get_gadm_mapping(gadm_shapes_path)
    pypsa_network.generators.bus = pypsa_network.generators.bus.str.replace("_AC", "")
    pypsa_network.generators.bus = pypsa_network.generators.bus.str.replace("_DC", "")
    pypsa_network.generators.bus = pypsa_network.generators.bus.str.replace(" csp", "")
    pypsa_network.generators.bus = pypsa_network.generators.bus.map(usa_state_dict)

    # select the relevant data from PyPSA results
    log_file.write(" --> Selecting relevant data from PyPSA results \n")
    pypsa_generation_filtered = pypsa_network.generators.reset_index()[['carrier', 'bus', 'gen_MWh']]
    pypsa_generation_filtered = pypsa_generation_filtered.rename(columns={'bus': 'State'})
    pypsa_generation_filtered.set_index(['carrier', 'State'], inplace=True)
    pypsa_generation_filtered['gen_MWh'] /= 1e6

    # Join the EIA and PyPSA generations to compare them
    log_file.write(" --> Joining EIA and PyPSA generations to compare them \n")
    comparison_generation = pypsa_generation_filtered.copy()
    comparison_generation = comparison_generation.join(eia_generation_reference_filtered)
    comparison_generation = comparison_generation.rename(columns={'gen_MWh': pypsa_name, year_to_use: eia_name})

    # iterate through each carrier, filter and plot the comparison between PyPSA and EIA generation statewise
    for car in ['nuclear', 'onwind', 'hydro', 'geothermal']:
        comparison_generation_filter = comparison_generation.reset_index().query('carrier == @car')
        comparison_generation_filter = comparison_generation_filter.set_index('State')[[pypsa_name, eia_name]]
        fig = px.bar(comparison_generation_filter, barmode='group')
        fig.update_layout(title=f"Carrier = {car}",yaxis_title='Generation (TWh)')
        fig.write_html(pathlib.Path(plot_base_path, f"{car}_generation_comparison.html"))

    # compute error %s of the generation
    comparison_generation['error %'] = comparison_generation.apply(lambda x: (x[eia_name]-x[pypsa_name])*100/x[eia_name], axis=1)


def plot_country_comparison(pypsa_network, eia_reference, year_to_use, time_resolution, log_file, plot_base_path):
    """
    The function plots the countrywide comparison between the EIA reference generation data and the PyPSA-Earth network
    """

    log_file.write("        \n")
    log_file.write("        \n")
    log_file.write("Compare the countrywide electricity generation \n")

    # prepare the PyPSA results
    pypsa_network.storage_units = pypsa_network.storage_units.assign(p=pypsa_network.storage_units_t.p.sum() * time_resolution)
    pypsa_network.generators = pypsa_network.generators.assign(p=pypsa_network.generators_t.p.sum() * time_resolution)
    df_pypsa_generation = pd.concat(
        [pypsa_network.generators.groupby("carrier").p.sum(), pypsa_network.storage_units.groupby("carrier").p.sum()])
    df_pypsa_generation.loc["wind"] = df_pypsa_generation.loc[["offwind-ac", "offwind-dc", "onwind"]].sum()
    df_pypsa_generation = df_pypsa_generation.drop(["offwind-ac", "offwind-dc", "onwind"])
    df_pypsa_generation /= 1e6
    df_pypsa_generation = df_pypsa_generation.round(2)
    df_pypsa_generation.name = pypsa_name
    df_eia_generation_year = eia_reference.loc[eia_reference["Period"] == year_to_use].squeeze()

    log_output_file.write("Electricity generation: {} \n".format(df_pypsa_generation))

    # prepare the EIA reference data
    pypsa_cols = ["Coal", "Natural Gas", "Other Gas", "Nuclear", "Hydro", "Estimated Total Solar", "PHS", "Petroleum",
                  "Wind", "Other Waste Biomass", "Geothermal"]
    rename_cols = {"estimated total solar": "solar", "other waste biomass": "biomass", "natural gas": "CCGT",
                   "petroleum": "oil", "phs": "PHS"}

    df_eia_generation_year = df_eia_generation_year[pypsa_cols]
    df_eia_generation_year.index = df_eia_generation_year.index.str.lower()
    df_eia_generation_year = df_eia_generation_year.rename(index=rename_cols)
    df_eia_generation_year.name = eia_name
    df_eia_generation_year = df_eia_generation_year.drop("other gas")

    # prepare comparison dataframe
    df_compare_generation = pd.concat([df_pypsa_generation, df_eia_generation_year], axis=1)

    # plot
    df_compare_generation.plot(kind="bar")
    plt.xlabel("Carriers")
    plt.ylabel("Electricity Generation (TWh)")
    plt.title(f"Year = {year_to_use}")
    plt.grid(linestyle="--")
    plt.subplots_adjust(bottom=0.3)
    plt.savefig(pathlib.Path(plot_base_path, f"countrywide_electricity_generation.png"), dpi=800)

    marginal_costs = pd.concat([pypsa_network.generators.groupby("carrier").marginal_cost.first().sort_values(),
                                pypsa_network.storage_units.groupby("carrier").marginal_cost.first().sort_values()])

    log_file.write("====")
    log_file.write("\nMarginal costs of electricity: {}".format(marginal_costs.round(3)))

    log_file.write("\nTotal electricity generation ({}):".format(year_to_use))
    log_file.write("EIA: {} TWh \n".format(df_eia_generation_year.sum().round(2)))
    log_file.write("PyPSA: {} TWh \n".format(df_pypsa_generation.sum().round(2)))
    log_file.write("====")


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

    network_pypsa_earth_df, eia_country_generation_df, eia_state_generation_df, gadm_path, time_res = parse_inputs(default_path)

    if snakemake.params.plot_country_comparison:
        plot_country_comparison(network_pypsa_earth_df, eia_country_generation_df, year_for_comparison, 24, log_output_file, plot_path)

    if snakemake.params.plot_state_by_state_comparison:
        plot_state_by_state_comparison(network_pypsa_earth_df, eia_state_generation_df, year_for_comparison, time_res, log_output_file, plot_path, gadm_path)

