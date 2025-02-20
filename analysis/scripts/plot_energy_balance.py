import pypsa
import matplotlib.pyplot as plt
import pathlib
import pandas as pd

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
    energy_generations = pd.DataFrame()
    for sector in sector_array:
        sector_buses = pypsa_network.buses.index[pypsa_network.buses.carrier.str.contains(sector)]
        components = ['generators_t', 'links_t', 'storage_units_t']
        sector_generations = pd.DataFrame()
        for comp in components:
            df = get_component(pypsa_network, comp)
            if comp != 'links':
                capacities = (df.query('bus in @sec', local_dict={'sec':sector_buses}).groupby('carrier').p_nom_opt.sum().div(1e3))
                sector_generations = pd.concat([sector_generations,capacities])
            else:
                capacities = (df.query('bus1 in @sec', local_dict={'sec':sector_buses}).groupby('carrier').p_nom_opt.sum() * df.query('bus1 in @sec', local_dict={'sec':sector_buses}).groupby('carrier').efficiency.mean()).div(1e3)
                capacities.name = 'p_nom_opt'
                sector_generations = pd.concat([sector_generations,capacities])
        sector_generations['sector'] = sector
        energy_generations = pd.concat([energy_generations,sector_generations])
    return energy_generations

def plot_electricity_balance(pypsa_network, renewable_carriers, conventional_carriers, storage_unit_carriers):
    # ac_buses = pypsa_network.buses.query("carrier == 'AC'")
    links = (
        pypsa_network.links_t.p1.groupby(pypsa_network.links.carrier, axis=1).sum().div(1e3) * -1
    )
    conventional_energy = links[link_carriers]
    p_by_carrier = conventional_energy

    generators = pypsa_network.generators_t.p.groupby(pypsa_network.generators.carrier, axis=1).sum().div(1e3)
    renewable_energy = generators[generator_carriers]
    p_by_carrier = pd.concat([p_by_carrier,generators],axis=1)

    if not pypsa_network.storage_units.empty:
        sto = (
            pypsa_network.storage_units_t.p.groupby(pypsa_network.storage_units.carrier, axis=1).sum().div(1e3)
        )
        p_by_carrier = pd.concat([p_by_carrier,sto],axis=1)

        fig, ax = plt.subplots(figsize=(6, 3))

    # color = p_by_carrier.columns.map(pypsa_network.carriers.color)
    
    p_by_carrier.where(p_by_carrier > 0).plot.area(
        ax=ax,
        linewidth=0,
        # color=color,
    )

    charge = p_by_carrier.where(p_by_carrier < 0).dropna(how="all", axis=1)

    if not charge.empty:
        charge.plot.area(
            ax=ax,
            linewidth=0,
            # color=charge.columns.map(pypsa_network.carriers.color)
        )

    # (pypsa_network.loads_t.p_set.sum(axis=1).loc[time].div(1e3)).plot(ax=ax, c="k",linewidth = 0.1)
    plt.legend(loc=(1.05, 0),fontsize=5,reverse=True)
    ax.set_ylabel("GW")
    # ax.set_ylim(-5, 30)
    plt.subplots_adjust(right=0.8)
    ax.grid(linestyle="--")
    return p_by_carrier

def installed_capacity_plots(pypsa_network, plot_base_path):
    installed_capacities = get_capacities(pypsa_network,['AC','heat','water','H2','oil','gas','biomass'])

    fig = px.bar(installed_capacities, y='p_nom_opt', color='sector', barmode='group', height=800)
    fig.write_image(
        f"{plot_base_path}/installed_capacity_eia_statewise.png", scale=1.5
    )

def energy_generation_by_sector():
    pass

pypsa_network_path = pathlib.Path("workflow","pypsa-earth","results","postnetworks","elec_s_100_ec_lcopt_Co2L-24H_144H_2030_0.071_AB_10export.nc")
pypsa_network = pypsa.Network(pypsa_network_path)

renewable_carriers = config['electricity']['renewable_carriers']
conventional_carriers = config['electricity']['conventional_carriers']
storage_unit_carriers = ['PHS', 'hydro']

p_by_carrier = plot_electricity_balance(pypsa_network, renewable_carriers, conventional_carriers, storage_unit_carriers)
