import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import pathlib
import geopandas as gpd
import sys

def disaggregation_v1(holes_mapped_intersect_filter, holes_centroid, df_utilities_grouped_state):
    holes_mapped_intersect_filter['Sales (Megawatthours)'] = 0
    df_final = pd.DataFrame()
    tot = 0
    for state in df_utilities_grouped_state.index:
        holes_state = holes_mapped_intersect_filter.query('State == @state')
        demand = df_utilities_grouped_state.loc[state]

        if len(holes_state) == 1:
            holes_mapped_intersect_filter.loc[holes_mapped_intersect_filter['GID_1'] == holes_state['GID_1'].values[0],'Sales (Megawatthours)'] = demand
            df_final = df_final._append(holes_mapped_intersect_filter.loc[holes_mapped_intersect_filter['GID_1'] == holes_state['GID_1'].values[0]], ignore_index=True)
            tot += demand
        elif len(holes_state) == 0:
            state_utilities = df_erst_gpd.query("STATE == @state")
            state_utilities['Sales (Megawatthours)'] = state_utilities['Sales (Megawatthours)']/state_utilities['Sales (Megawatthours)'].sum() * demand
            # df_final = df_final._append(state_utilities, ignore_index=True)
            df_erst_gpd.loc[state_utilities.index, 'Sales (Megawatthours)'] += state_utilities['Sales (Megawatthours)']
            tot += demand
        elif len(holes_state) > 1:
            holes_state['Sales (Megawatthours)'] = holes_state['pop'] / holes_state['pop'].sum() * demand
            df_final = df_final._append(holes_state, ignore_index=True)
            # tot += (holes_state['pop'] / holes_state['pop'].sum() * demand).sum()
            tot += demand

    return df_final

def disaggregation_v2(holes_mapped_intersect_filter, holes_centroid, df_utilities_grouped_state, df_demand_utility):
    holes_mapped_intersect_filter['Sales (Megawatthours)'] = 0
    df_final = pd.DataFrame()
    tot = 0
    for state in df_utilities_grouped_state.index:
        holes_state = holes_mapped_intersect_filter.query('State == @state')
        demand = df_utilities_grouped_state.loc[state]
        avg_state_demand = df_demand_utility.query('State == @state')['Sales (Megawatthours)'].mean()

        if len(holes_state) == 1:
            holes_mapped_intersect_filter.loc[holes_mapped_intersect_filter['GID_1'] == holes_state['GID_1'].values[0],'Sales (Megawatthours)'] = avg_state_demand
            df_final = df_final._append(holes_mapped_intersect_filter.loc[holes_mapped_intersect_filter['GID_1'] == holes_state['GID_1'].values[0]], ignore_index=True)
            tot += avg_state_demand
        elif len(holes_state) == 0:
            state_utilities = df_erst_gpd.query("STATE == @state")
            state_utilities['Sales (Megawatthours)'] = state_utilities['Sales (Megawatthours)']/state_utilities['Sales (Megawatthours)'].sum() * avg_state_demand
            # df_final = df_final._append(state_utilities, ignore_index=True)
            df_erst_gpd.loc[state_utilities.index, 'Sales (Megawatthours)'] += state_utilities['Sales (Megawatthours)']
            tot += avg_state_demand
        elif len(holes_state) > 1:
            holes_state['Sales (Megawatthours)'] = holes_state['pop'] / holes_state['pop'].sum() * avg_state_demand
            df_final = df_final._append(holes_state, ignore_index=True)
            # tot += (holes_state['pop'] / holes_state['pop'].sum() * avg_state_demand).sum()
            tot += avg_state_demand
    
    return df_final

def rescale_demands(df_final, df_demand_utility, df_utilities_grouped_state):
    df_demand_statewise = df_demand_utility.groupby('State')['Sales (Megawatthours)'].sum()
    for state in df_utilities_grouped_state.index:
        actual_state_demand = df_demand_statewise.loc[state]
        missing_demand = df_utilities_grouped_state.loc[state]
        df_filter_final = df_final.query('State == @state')
        assigned_utility_demand = df_filter_final['Sales (Megawatthours)'].sum()
        unmet_demand = missing_demand - assigned_utility_demand
        # print(state)
        # print(missing_demand)
        # print(assigned_utility_demand)
        # print(actual_state_demand)
        # input()
        # if assigned_utility_demand > 0:
        #     rescaling_factor = missing_demand / assigned_utility_demand
        if assigned_utility_demand != 0:
            rescaling_factor = actual_state_demand / assigned_utility_demand
        elif assigned_utility_demand == 0:
            rescaling_factor = actual_state_demand / missing_demand
        else:
            rescaling_factor = 1
        # print(rescaling_factor)
        # print(df_final.loc[df_final['State'] == state,'Sales (Megawatthours)'].sum() * rescaling_factor)
        df_final.loc[df_final['State'] == state,'Sales (Megawatthours)'] *= rescaling_factor

    return df_final

def calc_percentage_unmet_demand_by_state(df_calc, df_ref, df_error, text):
    df_calc_state = df_calc.groupby('State')['Sales (Megawatthours)'].sum()
    df_ref_state = df_ref.groupby('State')['Sales (Megawatthours)'].sum()
    df_error[text] = (df_ref_state - df_calc_state) * 100 / (df_ref_state)
    return df_error
 
if __name__ ==  '__main__':

    # Path 
    base_path = pathlib.Path(__file__).parent.parent.parent
    demand_utility_path = pathlib.Path(base_path, "analysis", "gdrive_data","data", "electricity_demand_data", "demand_data","table_10_EIA_utility_sales.xlsx")
    country_gadm_path = pathlib.Path(base_path, "workflow", "pypsa-earth", "resources", "US_2021", "shapes", "country_shapes.geojson")
    erst_path = pathlib.Path(base_path, "analysis", "data", "ERST_demand_GADM.geojson")
    gadm_usa_path = pathlib.Path(base_path, "analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json")
    pypsa_earth_scripts_path = pathlib.Path(base_path, "workflow", "pypsa-earth", "scripts")

    # Load data
    df_demand_utility = pd.read_excel(demand_utility_path,skiprows=2)
    df_erst_gpd = gpd.read_file(erst_path)
    df_country = gpd.read_file(country_gadm_path)
    df_gadm_usa = gpd.read_file(gadm_usa_path)

    # Add system paths
    sys.path.insert(0,str(pypsa_earth_scripts_path))
    import build_shapes

    # Printing % Demand grouped by ownership
    df_demand_own = df_demand_utility.groupby('Ownership')['Sales (Megawatthours)'].sum()
    print(df_demand_own * 100 / df_demand_own.sum())

    # Cumulative distribution function for demand binned by shape area
    df_cdf = df_erst_gpd[['Sales (Megawatthours)','Shape__Area']]
    bins = [0,1,5,10,15,20,35,50,100,1000]
    df_cdf['area_bin'] = pd.cut(df_cdf['Shape__Area'],bins)
    df_cdf = df_cdf.groupby('area_bin').sum()
    df_cdf['percent sales'] = df_cdf['Sales (Megawatthours)'] * 100 / df_cdf['Sales (Megawatthours)'].sum()
    px.bar(df_cdf.reset_index(),y='Sales (Megawatthours)')

    df_error = pd.DataFrame()
    # Initial error percentages of unmet demand
    df_error = calc_percentage_unmet_demand_by_state(df_erst_gpd.rename(columns={'STATE':'State'}), df_demand_utility, df_error, 'Initial')

    # Obtain holes in ERST shape files
    df_erst_gpd_dissolved = df_erst_gpd.dissolve()
    holes = df_country.difference(df_erst_gpd_dissolved)
    holes_exploded = holes.explode()
    holes_exploded = gpd.GeoDataFrame(geometry=holes.explode(),crs=df_erst_gpd.crs)
    holes_exploded['Area'] = holes_exploded.area
    # Filtering out holes with very small areas -
    area_threshold = 0.05
    holes_exploded_filter = holes_exploded.query('Area > @area_threshold')
    # holes_exploded_filter = holes_exploded.copy()
    holes_mapped = holes_exploded_filter.sjoin(df_gadm_usa)

    # Compute intersecting areas of holes and states
    # To filter out sjoin mapping to states where a tiny area of the hole is present in the GADM shape
    holes_mapped_intersect = holes_mapped.copy()
    holes_mapped_intersect['geometry'] = holes_mapped.apply(lambda x: x.geometry.intersection(df_gadm_usa.loc[df_gadm_usa['GID_1'] == x['GID_1']].iloc[0].geometry), axis=1)
    holes_mapped_intersect['area'] = holes_mapped_intersect.area
    holes_mapped_intersect_filter = holes_mapped_intersect.loc[holes_mapped_intersect['area'] > 1e-3]

    # Compute population and GDP
    holes_mapped_intersect_filter['GADM_ID'] = (np.arange(0,len(holes_mapped_intersect_filter),1))
    holes_mapped_intersect_filter['GADM_ID'] = holes_mapped_intersect_filter['GADM_ID'].astype('str')
    holes_mapped_intersect_filter['country'] = 'US'
    build_shapes.add_population_data(holes_mapped_intersect_filter,['US'],'standard',nprocesses=4)
    holes_mapped_intersect_filter['State'] = holes_mapped_intersect_filter.apply(lambda x: x['HASC_1'].split('.')[1],axis=1)
    
    # Missing utilities in ERST shape files
    missing_utilities = (list(set(df_demand_utility.Entity.str.upper()) - set(df_erst_gpd.NAME)))
    df_missing_utilities = df_demand_utility.query('Entity.str.upper() in @missing_utilities')
    df_utilities_grouped_state = df_missing_utilities.groupby('State')['Sales (Megawatthours)'].sum()

    # Compute centroid of the holes
    holes_centroid = holes_mapped_intersect_filter.copy()
    holes_centroid.geometry = holes_mapped_intersect_filter.geometry.centroid
    holes_centroid['State'] = holes_centroid.apply(lambda x: x['HASC_1'].split('.')[1],axis=1)

    version = 'v2'
    # v1 - disaggregation of unmet demand to holes based on population
    # v2 - demand at holes: average demand of the state where the hole is multiplied by the population
    # v2 - Rescaling demand of all utilities in a state to account for unmet demand in that particular state 
    if version == 'v1':
        df_final = disaggregation_v1(holes_mapped_intersect_filter, holes_centroid, df_utilities_grouped_state)
    elif version == 'v2':
        df_final = disaggregation_v2(holes_mapped_intersect_filter, holes_centroid, df_utilities_grouped_state, df_demand_utility)
    
    df_final = df_final._append(df_erst_gpd.rename(columns={'STATE':'State'}))

    # error percentages of unmet demand after assigning average demand to states
    df_error = calc_percentage_unmet_demand_by_state(df_final, df_demand_utility, df_error, 'Mid-way')

    if version == 'v2':
        df_final = rescale_demands(df_final, df_demand_utility, df_utilities_grouped_state)

        # Final error percentages of unmet demand after rescaling
        df_error = calc_percentage_unmet_demand_by_state(df_final, df_demand_utility, df_error, 'Final')

    fig = px.bar(df_error, barmode='group')
    fig.update_layout(yaxis_title='Unmet demand error %', xaxis_title='State')
    fig.write_image(f"../Plots/unmet_demand_error.png")

    # Adding population data for the final merged dataframe
    build_shapes.add_population_data(df_final,['US'],'standard',nprocesses=4)

    # Per-capita consumption
    df_per_capita['Calculated'] = df_final.groupby('State')['Sales (Megawatthours)'].sum() / df_final.groupby('State')['pop'].sum()

    geo_df_final = gpd.GeoDataFrame(df_final, geometry='geometry')
    geo_df_final['Sales (TWh)'] = geo_df_final['Sales (Megawatthours)'] / 1e6
    # geo_df_final['per capita'] = geo_df_final['Sales (Megawatthours)'] / geo_df_final['population']
    # Plot the GeoDataFrames
    m = geo_df_final.explore(column='Sales (TWh)',cmap='jet')
    m.save(f"../Plots/demand_filled_TWh_USA_{version}.html")

    df_erst_gpd['Sales (TWh)'] = df_erst_gpd['Sales (Megawatthours)'] / 1e6
    m = df_erst_gpd.explore(column='Sales (TWh)',cmap='jet')
    m.save(f"../Plots/demand_with_holes_TWh_USA_{version}.html")

    df_final.to_file(f"Demand_mapped_{version}.geojson",driver="GeoJSON")
