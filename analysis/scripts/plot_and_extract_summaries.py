"""
Plotting scripts to visualize model results

Relevant Settings
-----------------

.. code:: yaml

    US:
        summary:
            energy_carriers: ['electricity','heat']

.. seealso::
    Documentation of the configuration file ``config.usa_baseline.yaml``


Inputs
------

- ``workflow/pypsa-earth/results/postnetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc``


Outputs
-------

- ``analysis/plots/summary_plots``
- ``analysis/outputs/summary_outputs``


Description
-----------

To-do
------
1. Add plot of technology costs (sorted by merit-order for each energy carrier)
2. Ordering of technologies on the plot
3. Reflection of colors as in the config file
"""

import pypsa
import pathlib
import pandas as pd
import plotly.express as px
import datetime as dt
import plotly.graph_objects as go
from _helpers_usa import (
    get_component,
    get_component_list,
    get_energy_carriers_key,
    drop_carriers,
)
import numpy as np


def parse_inputs(default_path):
    """
    Load all input data required for preprocessing and computing utility level demand data

    Parameters
    ----------
    default_path: str
        current total directory path

    Returns
    -------
    pypsa_network: PyPSA model
    """
    pypsa_network_path = pathlib.Path(
        default_path, snakemake.input.pypsa_earth_results_path[0]
    )
    pypsa_network = pypsa.Network(pypsa_network_path)

    return pypsa_network


def get_capacities(pypsa_network, energy_carriers_array):
    """
    To fetch installed capacities across all technologies for a list of energy carriers

    Parameters
    ----------
    pypsa_network: PyPSA model
    energy_carriers_array: str list
        list of energy carriers

    Returns
    -------
    installed_capacities: pandas dataframe
        Installed capacities in GW
    """

    installed_capacities = pd.DataFrame()
    for energy_carriers in energy_carriers_array:
        energy_carriers_key = get_energy_carriers_key(energy_carriers)
        energy_carriers_buses = pypsa_network.buses.index[
            pypsa_network.buses.carrier.str.contains(energy_carriers_key)
        ]
        components = get_component_list(energy_carriers)
        # Removing stores as they contain only the nominal storage capacity in MWh and not the power capacity in MW
        if "stores" in components:
            components.remove("stores")
        energy_carriers_capacities = pd.DataFrame()
        for comp in components:
            df = get_component(pypsa_network, comp)
            if comp != "links":
                capacities = (
                    df.query("bus in @sec", local_dict={"sec": energy_carriers_buses})
                    .groupby("carrier")
                    .p_nom_opt.sum()
                    .div(1e3)
                )
                energy_carriers_capacities = pd.concat(
                    [energy_carriers_capacities, capacities]
                )
            else:
                capacities = (
                    df.query("bus1 in @sec", local_dict={"sec": energy_carriers_buses})
                    .groupby("carrier")
                    .p_nom_opt.sum()
                    * df.query(
                        "bus1 in @sec", local_dict={"sec": energy_carriers_buses}
                    )
                    .groupby("carrier")
                    .efficiency.mean()
                ).div(1e3)
                capacities.name = "p_nom_opt"
                energy_carriers_capacities = pd.concat(
                    [energy_carriers_capacities, capacities]
                )
        energy_carriers_capacities["carrier"] = energy_carriers
        installed_capacities = pd.concat(
            [installed_capacities, energy_carriers_capacities]
        )
    installed_capacities = drop_carriers(installed_capacities, "rows")
    return installed_capacities


def get_generations(pypsa_network, energy_carriers_array):
    """
    To fetch energy generations in TWh

    Parameters
    ----------
    pypsa_network: PyPSA model
    energy_carriers_array: str list
        List of energy carriers

    Returns
    -------
    energy_generations: pandas dataframe
        Energy generations across energy carriers in TWh
    """
    energy_generations = pd.DataFrame()
    for energy_carriers in energy_carriers_array:
        energy_carriers_key = get_energy_carriers_key(energy_carriers)
        energy_carriers_buses = pypsa_network.buses.index[
            pypsa_network.buses.carrier.str.contains(energy_carriers_key)
        ]
        components = get_component_list(energy_carriers)
        energy_carriers_generations = pd.DataFrame(index=pypsa_network.snapshots)
        for comp in components:
            df = get_component(pypsa_network, comp)
            if comp != "links":
                generators = df.query(
                    "bus in @sec", local_dict={"sec": energy_carriers_buses}
                )
                indices = generators.index
                reqd_carriers = generators.carrier.unique()
                generators_ts = get_component(pypsa_network, comp + "_t")
                generations = (
                    generators_ts[indices]
                    .groupby(
                        generators.query("carrier in @reqd_carriers").carrier, axis=1
                    )
                    .sum()
                    .div(1e3)
                )
                energy_carriers_generations = energy_carriers_generations.join(
                    generations, lsuffix="", rsuffix="_" + energy_carriers
                )
            else:
                generators = df.query(
                    "bus1 in @sec", local_dict={"sec": energy_carriers_buses}
                )
                indices = generators.index
                reqd_carriers = generators.carrier.unique()
                charger_carriers = df.query(
                    "carrier.str.contains(' charger') and bus0 in (@sec)",
                    local_dict={"sec": energy_carriers_buses},
                ).carrier.unique()

                reqd_carriers = np.append(reqd_carriers, charger_carriers)
                generators_ts = get_component(pypsa_network, comp + "_t")
                for car in reqd_carriers:
                    if "CHP" in car or "Fuel cell" in car:
                        intersecting_columns = indices.intersection(
                            generators_ts.p2.columns
                        )
                        generations = (
                            (generators_ts.p2[intersecting_columns])
                            .groupby(
                                generators.query("carrier == @car").carrier, axis=1
                            )
                            .sum()
                            .div(1e3)
                            * -1
                        )
                    elif "Fischer-Tropsch" in car:
                        intersecting_columns = indices.intersection(
                            generators_ts.p3.columns
                        )
                        generations = (
                            (generators_ts.p3[intersecting_columns])
                            .groupby(
                                generators.query("carrier == @car").carrier, axis=1
                            )
                            .sum()
                            .div(1e3)
                            * -1
                        )
                    elif " charger" in car:
                        generations = (
                            (generators_ts.p1.filter(like=car)).sum(axis=1).div(1e3)
                        )
                        generations.name = car
                        generations = pd.DataFrame(generations)
                    else:
                        generations = (
                            (generators_ts.p1[indices])
                            .groupby(
                                generators.query("carrier == @car").carrier, axis=1
                            )
                            .sum()
                            .div(1e3)
                            * -1
                        )

                    energy_carriers_generations = energy_carriers_generations.join(
                        generations, lsuffix="", rsuffix="_" + energy_carriers
                    )

        time_granularity = pypsa_network.snapshot_weightings.objective
        energy_carriers_generations = (
            energy_carriers_generations.mul(time_granularity, axis=0).div(1e3).sum()
        )  # in TWh
        energy_carriers_generations.name = "Energy_TWh"
        energy_carriers_generations = pd.DataFrame(energy_carriers_generations)
        energy_carriers_generations["carrier"] = energy_carriers

        energy_generations = energy_generations._append(energy_carriers_generations)

    energy_generations = drop_carriers(energy_generations, "rows")
    return energy_generations


def get_demands(pypsa_network, energy_carriers_array):
    """
    Fetch demands across energy carriers

    Parameters
    ----------
    pypsa_network: PyPSA model
    energy_carriers_array: str list
        List of energy carriers

    Returns
    -------
    demands_ts: pandas dataframe
        Time series of demands for each PyPSA carrier
        e.g., AC, rail transport electricity, agriculture electricity

    demands_grouped_ts: pandas dataframe
        Time series of demands aggregated over each energy carrier
        e.g., electricity, heat
    """

    demands_ts = pd.DataFrame(index=pypsa_network.snapshots)
    demands_grouped_ts = pd.DataFrame(index=pypsa_network.snapshots)
    for energy_carriers in energy_carriers_array:
        energy_carriers_key = get_energy_carriers_key(energy_carriers)
        energy_carriers_buses = pypsa_network.buses.index[
            pypsa_network.buses.carrier.str.contains(energy_carriers_key)
        ]
        load_energy_carriers_columns = pypsa_network.loads.query(
            "bus in @sec", local_dict={"sec": energy_carriers_buses}
        ).index
        load_energy_carriers_ts = (
            pypsa_network.loads_t.p[load_energy_carriers_columns]
            .groupby(pypsa_network.loads.carrier, axis=1)
            .sum()
            .div(1e3)
        )

        demands_grouped_ts[energy_carriers] = load_energy_carriers_ts.sum(axis=1)
        demands_ts = demands_ts.join(
            load_energy_carriers_ts, how="left", lsuffix="_x", rsuffix="_y"
        )
    return demands_ts, demands_grouped_ts


def get_generation_demands_by_energy_carriers(pypsa_network, energy_carrier):
    """
    Fetch demands across different technologies for each energy carrier

    Parameters
    ----------
    pypsa_network: PyPSA model
    energy_carrier: str
       energy carrier

    Returns
    -------

    demands_grouped_ts: pandas dataframe
        Time series of demands in TWh aggregated for a particular energy carrier
    energy_carrier_generations: pandas dataframe
        Time series of generations in TWh across each of the technologies associated to an energy carrier
    """

    energy_carrier_key = get_energy_carriers_key(energy_carrier)
    load_energy_carrier_buses = pypsa_network.buses.index[
        pypsa_network.buses.carrier.str.contains(energy_carrier_key)
    ]
    load_energy_carrier_columns = pypsa_network.loads.query(
        "bus in @sec", local_dict={"sec": load_energy_carrier_buses}
    ).index
    load_energy_carrier_ts = (
        pypsa_network.loads_t.p[load_energy_carrier_columns]
        .groupby(pypsa_network.loads.carrier, axis=1)
        .sum()
        .div(1e3)
    )

    demands_grouped_ts = load_energy_carrier_ts.sum(axis=1)

    components = get_component_list(energy_carrier)
    energy_carrier_generations = pd.DataFrame(index=pypsa_network.snapshots)
    gen_energy_carrier_buses = pypsa_network.buses.index[
        pypsa_network.buses.carrier.str.contains(energy_carrier_key)
    ]

    for comp in components:
        df = get_component(pypsa_network, comp)
        if comp != "links":
            generators = df.query(
                "bus in @sec", local_dict={"sec": gen_energy_carrier_buses}
            )
            indices = generators.index
            reqd_carriers = generators.carrier.unique()
            generators_ts = get_component(pypsa_network, comp + "_t")
            generations = (
                generators_ts[indices]
                .groupby(
                    generators.query(
                        "carrier in @car", local_dict={"car": reqd_carriers}
                    ).carrier,
                    axis=1,
                )
                .sum()
                .div(1e3)
            )
            energy_carrier_generations = energy_carrier_generations.join(
                generations, lsuffix="", rsuffix="_" + comp
            )
        else:
            generators = df.query(
                "bus1 in @sec", local_dict={"sec": gen_energy_carrier_buses}
            )
            indices = generators.index
            reqd_carriers = generators.carrier.unique()
            charger_carriers = df.query(
                "carrier.str.contains(' charger') and bus0 in (@sec)",
                local_dict={"sec": gen_energy_carrier_buses},
            ).carrier.unique()

            reqd_carriers = np.append(reqd_carriers, charger_carriers)

            generators_ts = get_component(pypsa_network, comp + "_t")

            for car in reqd_carriers:
                if "CHP" in car or "Fuel cell" in car:
                    intersecting_columns = indices.intersection(
                        generators_ts.p2.columns
                    )
                    generations = (
                        (generators_ts.p2[intersecting_columns])
                        .groupby(generators.query("carrier == @car").carrier, axis=1)
                        .sum()
                        .div(1e3)
                        * -1
                    )
                elif "Fischer-Tropsch" in car:
                    intersecting_columns = indices.intersection(
                        generators_ts.p3.columns
                    )
                    generations = (
                        (generators_ts.p3[intersecting_columns])
                        .groupby(generators.query("carrier == @car").carrier, axis=1)
                        .sum()
                        .div(1e3)
                        * -1
                    )
                elif " charger" in car:
                    generations = (
                        (generators_ts.p1.filter(like=car)).sum(axis=1).div(1e3)
                    )
                    generations.name = car
                    generations = pd.DataFrame(generations)

                else:
                    generations = (
                        (generators_ts.p1[indices])
                        .groupby(generators.query("carrier == @car").carrier, axis=1)
                        .sum()
                        .div(1e3)
                        * -1
                    )

                energy_carrier_generations = energy_carrier_generations.join(
                    generations, lsuffix="", rsuffix="_" + comp
                )
    # discharge_indices = [x for x in energy_carrier_generations.columns if 'discharge' in x]
    # energy_carrier_generations = energy_carrier_generations.drop(discharge_indices, axis=1)
    energy_carrier_generations = drop_carriers(energy_carrier_generations, "columns")
    return demands_grouped_ts, energy_carrier_generations


def installed_capacity_plots(pypsa_network, energy_carriers_array, plot_base_path):
    """
    Plot installed capacity in GW across energy carriers

    Parameters
    ----------
    pypsa_network: PyPSA model
    energy_carriers_array: str list
        List of energy carriers
    plot_base_path: str
        Plot directory path

    Returns
    --------
    installed_capacities: pandas dataframe
        Installed capacities in GW across energy carriers

    """

    installed_capacities = get_capacities(pypsa_network, energy_carriers_array)

    fig = px.bar(
        installed_capacities,
        y="p_nom_opt",
        color="carrier",
        barmode="stack",
        text_auto="0.2f",
        width=1200,
        height=800,
    )
    fig.update_layout(
        uniformtext_minsize=5,
        uniformtext_mode="show",
    )
    fig.update_traces(textposition="outside")
    fig.write_image(f"{plot_base_path}/installed_capacity_GW.png")

    fig = px.bar(
        installed_capacities,
        y="p_nom_opt",
        facet_row="carrier",
        color="carrier",
        barmode="stack",
        text_auto="0.2f",
        width=1200,
        height=800,
    )
    fig.update_yaxes(matches=None)
    fig.update_layout(
        uniformtext_minsize=5,
        uniformtext_mode="show",
    )
    fig.update_traces(textposition="outside")
    fig.for_each_yaxis(lambda y: y.update(title=""))
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))

    fig.add_annotation(
        x=-0.05,
        y=0.5,
        text="Installed capacity (GW)",
        textangle=-90,
        xref="paper",
        yref="paper",
    )
    fig.write_image(f"{plot_base_path}/installed_capacity_facet_energy_carriers_GW.png")

    return installed_capacities


def energy_generation_plots(pypsa_network, energy_carriers_array, plot_base_path):
    """
    Plot energy generation in TWh across energy carriers

    Parameters
    ----------
    pypsa_network: PyPSA model
    energy_carriers_array: str list
        List of energy carriers
    plot_base_path: str
        Plot directory path

    Returns
    --------
    energy_generations: pandas dataframe
        Energy generations in TWh across energy carriers

    """
    energy_generations = get_generations(pypsa_network, energy_carriers_array)

    fig = px.bar(
        energy_generations,
        y="Energy_TWh",
        color="carrier",
        barmode="stack",
        text_auto="0.2f",
        width=1500,
        height=800,
    )
    fig.update_layout(
        uniformtext_minsize=6,
        uniformtext_mode="show",
    )
    fig.update_traces(textposition="outside")
    fig.write_image(f"{plot_base_path}/energy_generation_TWh.png", scale=3)

    fig = px.bar(
        energy_generations,
        y="Energy_TWh",
        facet_row="carrier",
        color="carrier",
        barmode="stack",
        text_auto="0.2f",
        width=1500,
        height=900,
    )
    fig.update_yaxes(matches=None)
    fig.update_layout(
        uniformtext_minsize=6,
        uniformtext_mode="show",
    )
    fig.update_traces(textposition="outside")
    fig.for_each_yaxis(lambda y: y.update(title=""))
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))

    fig.add_annotation(
        x=-0.05, y=0.5, text="Energy (TWh)", textangle=-90, xref="paper", yref="paper"
    )
    fig.write_image(
        f"{plot_base_path}/energy_generation_facet_energy_carriers_TWh.png", scale=3
    )

    return energy_generations


def demand_plots(pypsa_network, energy_carriers_array, plot_base_path):
    """
    Plot demands in TWh across energy carriers

    Parameters
    ----------
    pypsa_network: PyPSA model
    energy_carriers_array: str list
        List of energy carriers
    plot_base_path: str
        Plot directory path

    Returns
    --------
    demands_ts: pandas dataframe
    demands_agg: pandas dataframe
    """
    demands_ts, demands_grouped_ts = get_demands(pypsa_network, energy_carriers_array)
    time_granularity = pypsa_network.snapshot_weightings.objective
    demands_agg = (
        demands_grouped_ts.mul(time_granularity, axis=0).div(1e3).sum()
    )  # in TWh
    demands_agg.name = "load_TWh"

    fig = px.bar(
        demands_agg,
        y="load_TWh",
        barmode="group",
        text_auto="0.2f",
        width=1200,
        height=800,
    )
    fig.update_layout(
        uniformtext_minsize=5,
        uniformtext_mode="show",
    )
    fig.update_traces(textposition="outside")
    fig.write_image(f"{plot_base_path}/demands_energy_carrierswise_TWh.png")

    return demands_ts, demands_agg


def energy_balance_plot(plot_base_path, output_path, energy_carriers_array):
    """
    Plot energy balance between generations and demands in TWh across energy carriers

    Parameters
    ----------
    energy_carriers_array: str list
        List of energy carriers
    plot_base_path: str
        Plot directory path

    """
    for energy_carriers in energy_carriers_array:
        demands, generations = get_generation_demands_by_energy_carriers(
            pypsa_network, energy_carriers
        )
        generations.round(2).to_csv(
            pathlib.Path(
                output_path, f"Sector_generations_ts_{energy_carriers}_in_TWh.csv"
            )
        )

        fig = px.area(generations.where(generations > 0))

        generations_neg = generations.where(generations < 0).dropna(axis=1)
        breakpoint()
        for col in generations_neg.columns:
            fig.add_trace(
                go.Scatter(
                    x=generations_neg.index,
                    y=generations_neg[col],
                    stackgroup="one",
                    name=col,
                )
            )

        fig.add_trace(
            go.Scatter(
                x=demands.index,
                y=demands,
                line_color="black",
                name=f"{energy_carriers} demand",
            )
        )
        fig.update_layout(yaxis_title="Energy in TWh", title=energy_carriers)
        fig.update_layout(
            uniformtext_minsize=5,
            # uniformtext_mode="show",
        )
        fig.write_image(f"{plot_base_path}/energy_balance_ts_{energy_carriers}.png")


def compare_generation_demand_agg_plot(
    energy_generations, demands_agg, energy_carriers_array, plot_base_path
):
    """
    Plot generations vs demands in TWh across energy carriers

    Parameters
    ----------
    energy_generations: pandas dataframe
        Aggregated generation in TWh for each energy carrier
    demands_agg: pandas dataframe
        Aggregated demand in TWh for each energy carrier
    energy_carriers_array: str list
        List of energy carriers
    plot_base_path: str
        Plot directory path

    """
    energy_agg = energy_generations.groupby("carrier").sum()
    df_comparison = pd.DataFrame(index=energy_carriers_array)
    df_comparison["generations"] = energy_agg
    df_comparison["demands"] = demands_agg

    fig = px.bar(df_comparison, barmode="group", text_auto="0.2f")
    fig.update_layout(
        uniformtext_minsize=5,
        uniformtext_mode="show",
        yaxis_title="Energy in TWh",
        title="Generation vs demand",
    )
    fig.update_traces(textposition="outside")
    fig.write_image(f"{plot_base_path}/generation_vs_demand.png")

    fig = px.bar(df_comparison.T, barmode="stack", text_auto="0.2f")
    fig.update_layout(
        uniformtext_minsize=5,
        uniformtext_mode="show",
        yaxis_title="Energy in TWh",
        title="Generation vs demand",
    )
    fig.update_traces(textposition="inside")
    fig.write_image(f"{plot_base_path}/generation_vs_demand_stacked.png")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("plot_summaries")

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "plot_summaries")
    plot_path = pathlib.Path(default_path, snakemake.output.plot_path)
    output_path = pathlib.Path(default_path, snakemake.output.output_path)

    energy_carriers_array = snakemake.params.energy_carriers

    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(
        log_path, f"output_plot_summaries_{today_date[:10]}.txt"
    )
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    pathlib.Path(plot_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)

    pypsa_network = parse_inputs(default_path)

    installed_capacities = installed_capacity_plots(
        pypsa_network, energy_carriers_array, plot_path
    )
    installed_capacities.round(2).to_csv(
        pathlib.Path(output_path, "Installed_Capacities_in_GW.csv")
    )

    energy_generations = energy_generation_plots(
        pypsa_network, energy_carriers_array, plot_path
    )
    energy_generations.round(2).to_csv(
        pathlib.Path(output_path, "Energy_generations_in_TWh.csv")
    )

    demands, demands_agg = demand_plots(pypsa_network, energy_carriers_array, plot_path)
    demands.round(2).to_csv(pathlib.Path(output_path, "Demands_in_TWh.csv"))
    demands_agg.round(2).to_csv(
        pathlib.Path(output_path, "Aggregated_demands_in_TWh.csv")
    )

    compare_generation_demand_agg_plot(
        energy_generations, demands_agg, energy_carriers_array, plot_path
    )

    energy_balance_plot(plot_path, output_path, energy_carriers_array)

    log_output_file.close()
