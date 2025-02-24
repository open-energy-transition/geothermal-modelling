import pypsa
import matplotlib.pyplot as plt
import pathlib
import pandas as pd
import plotly.express as px

def parse_inputs(default_path):
    pypsa_network_path = pathlib.Path(default_path, snakemake.input.pypsa_earth_results_path[0])
    pypsa_network = pypsa.Network(pypsa_network_path)

    return pypsa_network


def get_component(network, component):
    if component == 'generators':
        return network.generators
    elif component == 'links':
        return network.links
    elif component == 'stores':
        return network.stores
    elif component == 'storage_units':
        return network.storage_units
    elif component == 'generators_t':
        return network.generators_t.p
    elif component == 'links_t':
        return network.links_t.p1
    elif component == 'stores_t':
        return network.stores_t.p
    elif component == 'storage_units_t':
        return network.storage_units_t.p

def get_capacities(pypsa_network, sector_array):
    installed_capacities = pd.DataFrame()
    for sector in sector_array:
        sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(sector)]
        components = ['generators', 'links', 'storage_units']
        sector_capacities = pd.DataFrame()
        for comp in components:
            df = get_component(pypsa_network, comp)
            if comp != 'links':
                capacities = (df.query('bus in @sec', local_dict={'sec':sector_buses}).groupby('carrier').p_nom_opt.sum().div(1e3))
                sector_capacities = pd.concat([sector_capacities,capacities])
            else:
                capacities = (df.query('bus1 in @sec', local_dict={'sec':sector_buses}).groupby('carrier').p_nom_opt.sum() * df.query('bus1 in @sec', local_dict={'sec':sector_buses}).groupby('carrier').efficiency.mean()).div(1e3)
                capacities.name = 'p_nom_opt'
                sector_capacities = pd.concat([sector_capacities,capacities])
        sector_capacities['sector'] = sector
        installed_capacities = pd.concat([installed_capacities,sector_capacities])
    return installed_capacities

def get_generations(pypsa_network, sector_array):
    energy_generations = pd.DataFrame(index=pypsa_network.snapshots)
    for sector in sector_array:
        sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(sector)]
        components = ['generators', 'links', 'storage_units']
        sector_generations = pd.DataFrame(index=pypsa_network.snapshots)
        for comp in components:
            df = get_component(pypsa_network, comp)
            if comp != 'links':
                generators = (df.query('bus in @sec', local_dict={'sec':sector_buses}))
                indices = generators.index
                reqd_carriers = generators.carrier.unique()
                generators_ts = get_component(pypsa_network, comp+'_t')
                generations = generators_ts[indices].groupby(generators.query('carrier in @reqd_carriers').carrier, axis=1).sum().div(1e3)
                sector_generations = sector_generations.join(generations)
            else:
                generators = (df.query('bus1 in @sec', local_dict={'sec':sector_buses}))
                indices = generators.index
                reqd_carriers = generators.carrier.unique()
                generators_ts = get_component(pypsa_network, comp+'_t')
                generations = generators_ts[indices].groupby(generators.query('carrier in @reqd_carriers').carrier, axis=1).sum().div(1e3) * -1
                sector_generations = sector_generations.join(generations)

        discharge_indices = [x for x in sector_generations.columns if 'discharge' in x]
        sector_generations = sector_generations.drop(discharge_indices, axis=1)
        energy_generations = energy_generations.join(sector_generations, how='left',lsuffix='_x', rsuffix='_y')
    return energy_generations

def get_demands(pypsa_network, sector_array):
    demands_ts = pd.DataFrame(index=pypsa_network.snapshots)
    demands_grouped_ts = pd.DataFrame(index=pypsa_network.snapshots)
    for sector in sector_array:
        sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(sector)]
        load_sector_columns = pypsa_network.demands.query('bus in @sec', local_dict={'sec':sector_buses}).index
        load_sector_ts = pypsa_network.demands_t.p[load_sector_columns].groupby(pypsa_network.demands.carrier, axis=1).sum().div(1e3)
        
        demands_grouped_ts[sector] = load_sector_ts.sum(axis=1)
        demands_ts = demands_ts.join(load_sector_ts, how='left', lsuffix='_x', rsuffix='_y')

    return demands_ts, demands_grouped_ts

def get_generation_demands_by_sector(pypsa_network, load_sector, gen_sector):
    load_sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(load_sector)]
    load_sector_columns = pypsa_network.demands.query('bus in @sec', local_dict={'sec':load_sector_buses}).index
    load_sector_ts = pypsa_network.demands_t.p[load_sector_columns].groupby(pypsa_network.demands.carrier, axis=1).sum().div(1e3)
    
    demands_grouped_ts = load_sector_ts.sum(axis=1)

    components = ['generators', 'links', 'storage_units']
    sector_generations = pd.DataFrame(index=pypsa_network.snapshots)
    gen_sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(gen_sector)]

    for comp in components:
        df = get_component(pypsa_network, comp)
        if comp != 'links':
            generators = (df.query('bus in @sec', local_dict={'sec':gen_sector_buses}))
            indices = generators.index
            reqd_carriers = generators.carrier.unique()
            generators_ts = get_component(pypsa_network, comp+'_t')
            generations = generators_ts[indices].groupby(generators.query('carrier in @car', local_dict={'car':reqd_carriers}).carrier, axis=1).sum().div(1e3)
            sector_generations = sector_generations.join(generations)
        else:
            generators = (df.query('bus1 in @sec', local_dict={'sec':gen_sector_buses}))
            indices = generators.index
            reqd_carriers = generators.carrier.unique()
            generators_ts = get_component(pypsa_network, comp+'_t')
            generations = generators_ts[indices].groupby(generators.query('carrier in @car', local_dict={'car':reqd_carriers}).carrier, axis=1).sum().div(1e3) * -1
            sector_generations = sector_generations.join(generations)

    discharge_indices = [x for x in sector_generations.columns if 'discharge' in x]
    sector_generations = sector_generations.drop(discharge_indices, axis=1)

    return demands_grouped_ts, sector_generations




# def plot_electricity_balance(pypsa_network, renewable_carriers, conventional_carriers, storage_unit_carriers):
#     # ac_buses = pypsa_network.buses.query("carrier == 'AC'")
#     links = (
#         pypsa_network.links_t.p1.groupby(pypsa_network.links.carrier, axis=1).sum().div(1e3) * -1
#     )
#     conventional_energy = links[link_carriers]
#     p_by_carrier = conventional_energy

#     generators = pypsa_network.generators_t.p.groupby(pypsa_network.generators.carrier, axis=1).sum().div(1e3)
#     renewable_energy = generators[generator_carriers]
#     p_by_carrier = pd.concat([p_by_carrier,generators],axis=1)

#     if not pypsa_network.storage_units.empty:
#         sto = (
#             pypsa_network.storage_units_t.p.groupby(pypsa_network.storage_units.carrier, axis=1).sum().div(1e3)
#         )
#         p_by_carrier = pd.concat([p_by_carrier,sto],axis=1)

#         fig, ax = plt.subplots(figsize=(6, 3))

#     # color = p_by_carrier.columns.map(pypsa_network.carriers.color)
    
#     p_by_carrier.where(p_by_carrier > 0).plot.area(
#         ax=ax,
#         linewidth=0,
#         # color=color,
#     )

#     charge = p_by_carrier.where(p_by_carrier < 0).dropna(how="all", axis=1)

#     if not charge.empty:
#         charge.plot.area(
#             ax=ax,
#             linewidth=0,
#             # color=charge.columns.map(pypsa_network.carriers.color)
#         )

#     # (pypsa_network.loads_t.p_set.sum(axis=1).loc[time].div(1e3)).plot(ax=ax, c="k",linewidth = 0.1)
#     plt.legend(loc=(1.05, 0),fontsize=5,reverse=True)
#     ax.set_ylabel("GW")
#     # ax.set_ylim(-5, 30)
#     plt.subplots_adjust(right=0.8)
#     ax.grid(linestyle="--")
#     return p_by_carrier

def installed_capacity_plots(pypsa_network, plot_base_path):
    installed_capacities = get_capacities(pypsa_network,['AC','heat','water','H2','oil','gas','biomass'])

    fig = px.bar(installed_capacities, y='p_nom_opt', color='sector', barmode='group', height=800)
    fig.write_image(
        f"{plot_base_path}/installed_capacity_sectorwise_GW.png", scale=1.5
    )

def energy_generation_plots(pypsa_network, plot_base_path):
    energy_generations = get_generations(pypsa_network,['AC','heat','water','H2','oil','gas','biomass'])
    time_granularity = pypsa_network.snapshots.diff()[-1].total_seconds() / 3600
    energy_generations_agg = energy_generations.sum() * time_granularity / 1e3 #in TWh
    energy_generations_agg.name = 'energy_gen'

    fig = px.bar(energy_generations_agg, y='energy_gen', barmode='group', height=800)
    fig.write_image(
        f"{plot_base_path}/energy_generation_sectorwise_TWh.png", scale=1.5
    )
    return energy_generations, energy_generations_agg

def demand_plots(pypsa_network, plot_base_path):
    demands = get_demands(pypsa_network,['low voltage','heat','water','H2','oil','gas','biomass'])
    time_granularity = pypsa_network.snapshots.diff()[-1].total_seconds() / 3600
    demands_agg = demands.sum() * time_granularity / 1e3 #in TWh
    demands_agg.name = 'load in TWh'

    fig = px.bar(demands_agg, y='load in TWh', barmode='group', height=800)
    fig.write_image(
        f"{plot_base_path}/demands_sectorwise_TWh.png", scale=1.5
    )

    return demands, demands_agg


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("plot_summaries")

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "plot_summaries")
    plot_path = pathlib.Path(default_path, snakemake.output.plot_path)

    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(
        log_path, f"output_plot_summaries_{today_date[:10]}.txt"
    )
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    config_path = pathlib.Path(default_path, "configs", "config.usa_baseline.yaml")

    pathlib.Path(plot_path).mkdir(parents=True, exist_ok=True)

    pypsa_network = parse_inputs(default_path)

    installed_capacity_plots(pypsa_network, plot_path)

    energy_generations, energy_generations_agg = energy_generation_plots(pypsa_network, plot_path)

    demands, demands_agg = demand_plots(pypsa_network, plot_path)
    
    log_output_file.close()