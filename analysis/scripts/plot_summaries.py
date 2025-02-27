import pypsa
import matplotlib.pyplot as plt
import pathlib
import pandas as pd
import plotly.express as px
import datetime as dt
import plotly.graph_objects as go

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
        return network.links_t
    elif component == 'stores_t':
        return network.stores_t.p
    elif component == 'storage_units_t':
        return network.storage_units_t.p

def get_component_list(sector):
    sector_component_dict = {
        'electricity': ['generators','links','storage_units','stores'],
        'heat':['links','stores'],
        'H2':['links','stores'],
        'water':['links','stores'],
        'oil':['links','stores'],
        'biomass':['links','stores'],
        'gas':['links','stores'],
        'battery':['links','stores']
    }
    return sector_component_dict[sector]

def get_sector_key(sector):
    if sector == 'electricity':
        return 'AC|low voltage'
    else:
        return sector

def drop_carriers(df):
    drop_carriers_list = ['electricity distribution grid','H2 pipeline','H2 pipeline repurposed']
    for car in drop_carriers_list:
        if car in df.index.tolist():
            df = df.drop(index=drop_carriers_list)

    return df

def get_capacities(pypsa_network, sector_array):
    installed_capacities = pd.DataFrame()
    for sector in sector_array:
        sector_key = get_sector_key(sector)
        sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(sector_key)]
        components = get_component_list(sector)
        # Removing stores as they contain only the nominal storage capacity in MWh and not the power capacity in MW
        components.remove('stores')
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
    installed_capacities = drop_carriers(installed_capacities)
    return installed_capacities

def get_generations(pypsa_network, sector_array):
    energy_generations = pd.DataFrame()
    for sector in sector_array:
        sector_key = get_sector_key(sector)
        sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(sector_key)]
        components = get_component_list(sector)
        sector_generations = pd.DataFrame(index=pypsa_network.snapshots)
        for comp in components:
            df = get_component(pypsa_network, comp)
            if comp != 'links':
                generators = (df.query('bus in @sec', local_dict={'sec':sector_buses}))
                indices = generators.index
                reqd_carriers = generators.carrier.unique()
                generators_ts = get_component(pypsa_network, comp+'_t')
                generations = generators_ts[indices].groupby(generators.query('carrier in @reqd_carriers').carrier, axis=1).sum().div(1e3)
                sector_generations = sector_generations.join(generations, lsuffix="", rsuffix="_"+sector)
            else:
                generators = (df.query('bus1 in @sec', local_dict={'sec':sector_buses}))
                indices = generators.index
                reqd_carriers = generators.carrier.unique()
                generators_ts = get_component(pypsa_network, comp+'_t')
                for car in reqd_carriers:
                    if ('CHP' in car or 'Fuel cell' in car):
                        intersecting_columns = indices.intersection(generators_ts.p2.columns)
                        generations = (generators_ts.p2[intersecting_columns]).groupby(generators.query('carrier == @car').carrier, axis=1).sum().div(1e3) * -1
                    elif 'Fischer-Tropsch' in car:
                        intersecting_columns = indices.intersection(generators_ts.p3.columns)
                        generations = (generators_ts.p3[intersecting_columns]).groupby(generators.query('carrier == @car').carrier, axis=1).sum().div(1e3) * -1
                    else:
                        generations = (generators_ts.p1[indices]).groupby(generators.query('carrier == @car').carrier, axis=1).sum().div(1e3) * -1

                    sector_generations = sector_generations.join(generations, lsuffix="", rsuffix="_"+sector)

        # discharge_indices = [x for x in sector_generations.columns if 'discharge' in x]
        # sector_generations = sector_generations.drop(discharge_indices, axis=1)

        time_granularity = pypsa_network.snapshots.diff()[-1].total_seconds() / 3600
        sector_generations = sector_generations.sum() * time_granularity / 1e3 #in TWh
        sector_generations.name='Energy_TWh'
        sector_generations = pd.DataFrame(sector_generations)
        sector_generations['sector'] = sector
        
        energy_generations = energy_generations._append(sector_generations)
    
    energy_generations = drop_carriers(energy_generations)
    return energy_generations

def get_demands(pypsa_network, sector_array):
    demands_ts = pd.DataFrame(index=pypsa_network.snapshots)
    demands_grouped_ts = pd.DataFrame(index=pypsa_network.snapshots)
    for sector in sector_array:
        sector_key = get_sector_key(sector)
        sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(sector_key)]
        load_sector_columns = pypsa_network.loads.query('bus in @sec', local_dict={'sec':sector_buses}).index
        load_sector_ts = pypsa_network.loads_t.p[load_sector_columns].groupby(pypsa_network.loads.carrier, axis=1).sum().div(1e3)
        
        demands_grouped_ts[sector] = load_sector_ts.sum(axis=1)
        demands_ts = demands_ts.join(load_sector_ts, how='left', lsuffix='_x', rsuffix='_y')

    return demands_ts, demands_grouped_ts

def get_generation_demands_by_sector(pypsa_network, sector):
    sector_key = get_sector_key(sector)
    load_sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(sector_key)]
    load_sector_columns = pypsa_network.loads.query('bus in @sec', local_dict={'sec':load_sector_buses}).index
    load_sector_ts = pypsa_network.loads_t.p[load_sector_columns].groupby(pypsa_network.loads.carrier, axis=1).sum().div(1e3)
    
    demands_grouped_ts = load_sector_ts.sum(axis=1)

    components = get_component_list(sector)
    sector_generations = pd.DataFrame(index=pypsa_network.snapshots)
    gen_sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(sector_key)]

    for comp in components:
        df = get_component(pypsa_network, comp)
        if comp != 'links':
            generators = (df.query('bus in @sec', local_dict={'sec':gen_sector_buses}))
            indices = generators.index
            reqd_carriers = generators.carrier.unique()
            generators_ts = get_component(pypsa_network, comp+'_t')
            generations = generators_ts[indices].groupby(generators.query('carrier in @car', local_dict={'car':reqd_carriers}).carrier, axis=1).sum().div(1e3)
            sector_generations = sector_generations.join(generations, lsuffix="", rsuffix="_"+sector)
        else:
            generators = (df.query('bus1 in @sec', local_dict={'sec':gen_sector_buses}))
            indices = generators.index
            reqd_carriers = generators.carrier.unique()
            generators_ts = get_component(pypsa_network, comp+'_t')
            for car in reqd_carriers:
                if ('CHP' in car or 'Fuel cell' in car):
                    intersecting_columns = indices.intersection(generators_ts.p2.columns)
                    generations = (generators_ts.p2[intersecting_columns]).groupby(generators.query('carrier == @car').carrier, axis=1).sum().div(1e3) * -1
                elif 'Fischer-Tropsch' in car:
                    intersecting_columns = indices.intersection(generators_ts.p3.columns)
                    generations = (generators_ts.p3[intersecting_columns]).groupby(generators.query('carrier == @car').carrier, axis=1).sum().div(1e3) * -1
                else:
                    generations = (generators_ts.p1[indices]).groupby(generators.query('carrier == @car').carrier, axis=1).sum().div(1e3) * -1

                sector_generations = sector_generations.join(generations, lsuffix="", rsuffix="_"+sector)

    # discharge_indices = [x for x in sector_generations.columns if 'discharge' in x]
    # sector_generations = sector_generations.drop(discharge_indices, axis=1)
    sector_generations = drop_carriers(sector_generations)
    return demands_grouped_ts, sector_generations

def installed_capacity_plots(pypsa_network, sector_array, plot_base_path):
    installed_capacities = get_capacities(pypsa_network,sector_array)

    fig = px.bar(installed_capacities, y='p_nom_opt', color='sector', barmode='stack',text_auto='0.2f', width=1200, height=800)
    fig.update_layout(
        uniformtext_minsize=5,
        uniformtext_mode='show',
    )
    fig.update_traces(textposition='outside')
    fig.write_image(
        f"{plot_base_path}/installed_capacity_GW.png"
    )

    px.bar(installed_capacities,y='p_nom_opt',facet_row='sector', color='sector',barmode='stack',text_auto='0.2f', width=1200, height=1500)
    fig.update_yaxes(matches=None)
    fig.update_layout(
        uniformtext_minsize=5,
        uniformtext_mode='show',
    )
    fig.update_traces(textposition='outside')
    fig.write_image(
        f"{plot_base_path}/installed_capacity_facet_sector_GW.png"
    )

    return installed_capacities

def energy_generation_plots(pypsa_network, sector_array, plot_base_path):
    energy_generations = get_generations(pypsa_network,sector_array)
    
    fig = px.bar(energy_generations, y='Energy_TWh', color='sector',barmode='stack',text_auto='0.2f',width=1200, height=1500)
    fig.update_layout(
        uniformtext_minsize=5,
        uniformtext_mode='show',
    )
    fig.update_traces(textposition='outside')
    fig.write_image(
        f"{plot_base_path}/energy_generation_TWh.png"
    )

    fig = px.bar(energy_generations, y='Energy_TWh', facet_row='sector', color='sector',barmode='stack',text_auto='0.2f',width=1200, height=1500)
    fig.update_yaxes(matches=None)
    fig.update_layout(
        uniformtext_minsize=5,
        uniformtext_mode='show',
    )
    fig.update_traces(textposition='outside')
    fig.write_image(
        f"{plot_base_path}/energy_generation_facet_sector_TWh.png"
    )

    return energy_generations

def demand_plots(pypsa_network, sector_array, plot_base_path):
    demands_ts, demands_grouped_ts = get_demands(pypsa_network, sector_array)
    time_granularity = pypsa_network.snapshots.diff()[-1].total_seconds() / 3600
    demands_agg = demands_grouped_ts.sum() * time_granularity / 1e3 #in TWh
    demands_agg.name = 'load_TWh'

    fig = px.bar(demands_agg, y='load_TWh', barmode='group',text_auto='0.2f',width=1200, height=800)
    fig.update_layout(
        uniformtext_minsize=5,
        uniformtext_mode='show',
    )
    fig.update_traces(textposition='outside')
    fig.write_image(
        f"{plot_base_path}/demands_sectorwise_TWh.png"
    )

    return demands_ts, demands_agg

def energy_balance_plot(plot_base_path, sector_array):

    for sector in sector_array:
        demands, generations = get_generation_demands_by_sector(pypsa_network, sector)

        fig = px.area(generations.where(generations > 0))

        generations_neg = generations.where(generations < 0).dropna(axis=1)
        for col in generations_neg.columns:
            fig.add_trace(go.Scatter(x=generations_neg.index, y=generations_neg[col], stackgroup='one', name=col))

        fig.add_trace(go.Scatter(x=demands.index, y=demands, line_color='black', name=f"{sector} demand"))
        fig.update_layout(yaxis_title = 'Energy in TWh', title=sector)
        fig.write_image(
            f"{plot_base_path}/energy_balance_ts_{sector}.png"
        )

def compare_generation_demand_agg_plot(energy_generations, demands_agg, sector_array, plot_base_path):
    energy_agg = energy_generations.groupby('sector').sum()
    df_comparison = pd.DataFrame(index=sector_array)
    df_comparison['generations'] = energy_agg
    df_comparison['demands'] = demands_agg

    fig = px.bar(df_comparison,barmode='group')
    fig.update_layout(yaxis_title = 'Energy in TWh', title="Generation vs demand")
    fig.write_image(
        f"{plot_base_path}/generation_vs_demand.png"
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("plot_summaries")

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "plot_summaries")
    plot_path = pathlib.Path(default_path, snakemake.output.plot_path)

    sector_array = snakemake.params.sector_array
    
    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(
        log_path, f"output_plot_summaries_{today_date[:10]}.txt"
    )
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    pathlib.Path(plot_path).mkdir(parents=True, exist_ok=True)

    pypsa_network = parse_inputs(default_path)

    installed_capacities = installed_capacity_plots(pypsa_network, sector_array, plot_path)

    energy_generations = energy_generation_plots(pypsa_network, sector_array, plot_path)

    demands, demands_agg = demand_plots(pypsa_network, sector_array, plot_path)

    compare_generation_demand_agg_plot(energy_generations, demands_agg, sector_array, plot_path)

    energy_balance_plot(plot_path, sector_array)

    log_output_file.close()