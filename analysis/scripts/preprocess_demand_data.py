import pathlib
import datetime as dt
import geopandas as gpd
import pandas as pd
import os
from shapely.validation import make_valid
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import random 

def compute_demand_disaggregation(holes_mapped_intersect_filter, holes_centroid, df_utilities_grouped_state, df_demand_utility, df_gadm_usa):
    holes_mapped_intersect_filter['Sales (Megawatthours)'] = 0
    df_final = pd.DataFrame()
    tot = 0
    for state in df_utilities_grouped_state.index:
        holes_state = holes_mapped_intersect_filter.query('State == @state')
        demand = df_utilities_grouped_state.loc[state]
        full_state_pop = df_gadm_usa.query('State == @state')['pop'].values[0]
        state_demand = df_demand_utility.query('State == @state')['Sales (Megawatthours)'].sum()
        avg_per_capita_demand = state_demand / full_state_pop

        if len(holes_state) > 0:
            holes_state['Sales (Megawatthours)'] = holes_state['pop'] * avg_per_capita_demand
            df_final = df_final._append(holes_state, ignore_index=True)
    
    return df_final

def calc_percentage_unmet_demand_by_state(df_calc, df_ref, df_error, text):
    df_calc_state = df_calc.groupby('State')['Sales (Megawatthours)'].sum()
    df_ref_state = df_ref.groupby('State')['Sales (Megawatthours)'].sum()
    df_error[text] = (df_ref_state - df_calc_state) * 100 / (df_ref_state)
    return df_error

def calc_per_capita_kWh_state(df_calc, df_gadm, df_per_capita_cons, text):
    df_calc_per_capita = df_calc.groupby('State')['Sales (Megawatthours)'].sum() * 1000 / df_gadm.groupby('State')['pop'].sum()
    df_per_capita_cons[text] = df_calc_per_capita
    return df_per_capita_cons

def rescale_demands(df_final, df_demand_utility, df_utilities_grouped_state):
    df_demand_statewise = df_demand_utility.groupby('State')['Sales (Megawatthours)'].sum()
    for state in df_utilities_grouped_state.index:
        actual_state_demand = df_demand_statewise.loc[state]
        missing_demand = df_utilities_grouped_state.loc[state]
        df_filter_final = df_final.query('State == @state')
        assigned_utility_demand = df_filter_final['Sales (Megawatthours)'].sum()
        unmet_demand = missing_demand - assigned_utility_demand
        if assigned_utility_demand != 0:
            rescaling_factor = actual_state_demand / assigned_utility_demand
        elif assigned_utility_demand == 0:
            rescaling_factor = actual_state_demand / missing_demand
        else:
            rescaling_factor = 1
        df_final.loc[df_final['State'] == state,'Sales (Megawatthours)'] *= rescaling_factor

    return df_final

# df_erst_gpd: geopandas dataframe -> ERST shape file
# df_gadm_usa: geopandas dataframe -> GADM shape file
# df_demand_utility: pandas dataframe -> Demand data from EIA
def overlay_demands(df_erst_gpd, df_gadm_usa, df_demand_utility):

    # ERST overlay with GADM shapes for intersection
    df_erst_overlay = gpd.overlay(df_erst_gpd, df_gadm_usa, how='intersection')
    df_erst_overlay['STATE'] = df_erst_overlay.apply(lambda x: x['ISO_1'].split('-')[1], axis=1)
    df_erst_overlay.set_index(['NAME','STATE'], inplace=True)

    # Demand mapping to ERST shape file
    df_demand_utility['Entity'] = df_demand_utility['Entity'].apply(str.upper)
    df_demand_utility.rename(columns={'Entity':'NAME','State':'STATE'}, inplace=True)
    df_demand_utility.set_index(['NAME','STATE'], inplace=True)

    df_erst_overlay = df_erst_overlay.join(df_demand_utility)
    
    # Make geometry valid
    df_erst_overlay['geometry'] = df_erst_overlay['geometry'].apply(make_valid)

    return df_erst_overlay

# total_demand: float -> Total Demand (Sales) as per the EIA data
# mapped_demand: float -> The demand mapped to ERST file at that particular stage 
def compute_missing_percentage(total_demand, mapped_demand):
    return np.round(((total_demand - mapped_demand) * 100 / total_demand), 2)

# df_map : geopandas -> geopandas file that is to be plotted on the map
# filename : str -> Name to be given to the saved plot (HTML file)
# color : boolean -> If True, each row of geometry is coloured based on the color assigned to it the gpd dataframe
# cmap : boolean -> If True, each row of geometry is coloured based on the magnitude of the value in the cmap_col specified
# cmap_col : str -> Column on which a color map is used to plot the geometry
def save_map(df_map, filename, color, cmap, cmap_col=''):
    if "VAL_DATE" in df_map.columns:
        df_map = df_map.drop("VAL_DATE",axis=1)
    if "SOURCEDATE" in df_map.columns:
        df_map = df_map.drop("SOURCEDATE",axis=1)
    if not color and not cmap:
        m = df_map.explore()
    elif color and not cmap:
        m = df_map.explore(color=df_map['color'].tolist())
    elif not color and cmap:
        m = df_map.explore(column=cmap_col, cmap='jet')
    m.save(os.path.join(plot_path,filename))


if __name__ == '__main__':

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs")
    plot_path = pathlib.Path(default_path, "analysis", "plots", "demand_modelling")
    os.makedirs(log_path, exist_ok=True)
    os.makedirs(plot_path, exist_ok=True)
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(log_path, f"output_preprocess_demand_{today_date[:10]}.txt")
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    # Load snakemake params
    demand_year = snakemake.params.demand_year
    holes_area_threshold = snakemake.params.holes_area_threshold
    nprocesses = snakemake.params.nprocesses
    log_output_file.write("Loading snakemake parameters \n")
    log_output_file.write(f"demand_year = {demand_year} \n")
    log_output_file.write(f"holes_area_threshold = {holes_area_threshold} \n")
    log_output_file.write("-------------------------------------------------\n")

    # Load data
    df_demand_utility = pd.read_excel(snakemake.input.demand_utility_path,skiprows=2)
    df_erst_gpd = gpd.read_file(snakemake.input.erst_path)
    df_country = gpd.read_file(snakemake.input.country_gadm_path)
    df_gadm_usa = gpd.read_file(snakemake.input.gadm_usa_path)
    df_eia_per_capita = pd.read_excel(snakemake.input.eia_per_capita_path, sheet_name='Total per capita', skiprows=2, index_col='State')
    df_eia_per_capita = df_eia_per_capita[demand_year]
    log_output_file.write("Reading input files completed \n")
    log_output_file.write("-------------------------------------------------\n")
    total_demand = df_demand_utility['Sales (Megawatthours)'].sum() / 1e6
    log_output_file.write(f"Total sales (TWh) as in EIA sales data: {total_demand} \n")

    # Add system paths
    pypsa_earth_scripts_path = pathlib.Path("workflow", "pypsa-earth", "scripts")
    sys.path.insert(0,str(pypsa_earth_scripts_path))
    import build_shapes

    get_colors = lambda n: ["#%06x" % random.randint(0, 0xFFFFFF) for _ in range(n)]

    # Make geometry valid
    df_erst_gpd['geometry'] = df_erst_gpd['geometry'].apply(make_valid)

    # Plot initial ERST shapefile 
    save_map(df_erst_gpd, filename="ERST_initial.html", color=False, cmap=False)
    log_output_file.write("Generated initial ERST shapes file \n")
    
    # Map demands to ERST shapefile
    df_erst_gpd = overlay_demands(df_erst_gpd, df_gadm_usa, df_demand_utility)
    demand_mapped = df_erst_gpd['Sales (Megawatthours)'].sum() / 1e6
    log_output_file.write(f"Total sales (TWh) as in ERST mapped data: {demand_mapped} \n")
    missing_demand_percentage = compute_missing_percentage(total_demand, demand_mapped)
    log_output_file.write(f"Missing sales (TWh) : {missing_demand_percentage} \n")

    # df_erst_gpd_centroid = df_erst_gpd.copy()
    # df_erst_gpd_centroid.geometry = df_erst_gpd_centroid.geometry
    save_map(df_erst_gpd, filename="ERST_overlay.html", color=False, cmap=True, cmap_col='Sales (Megawatthours)')
    log_output_file.write("Generated ERST shapes file with mapped demands \n")

    ## Obtain holes in ERST shape files

    # Merge all geometry into one geometry
    df_erst_gpd_dissolved = df_erst_gpd.dissolve()
    # Find the difference in geometry between country and the merged ERST geometry
    holes = df_country.difference(df_erst_gpd_dissolved)

    # Convert holes geomtry : multipolygons into polygons to create a separate row for each hole
    holes_exploded = holes.explode()
    holes_exploded = gpd.GeoDataFrame(geometry=holes.explode(),crs=df_erst_gpd.crs)
    save_map(holes_exploded, filename="Holes.html", color=False, cmap=False)
    log_output_file.write("Generated holes exploded map \n")
    holes_exploded = holes_exploded.to_crs(6372)
    holes_exploded['Area'] = holes_exploded.area / 1e6
    # Filtering out holes with very small areas (only hole areas larger than area_threshold considered)
    holes_exploded_filter = holes_exploded.query('Area > @holes_area_threshold')
    save_map(holes_exploded_filter, filename="Holes_considered.html", color=False, cmap=False)
    log_output_file.write(f"Generated holes greater than {holes_area_threshold} \n")
    holes_exploded_filter = holes_exploded_filter.to_crs(4326)

    df_gadm_usa['color'] = get_colors(len(df_gadm_usa))
    holes_mapped = holes_exploded_filter.sjoin(df_gadm_usa)
    save_map(holes_mapped, filename="Holes_mapped_GADM.html", color=True, cmap=False)
    log_output_file.write("Generated holes mapped to GADM \n")

    # # Compute intersecting areas of holes and states
    # # To filter out sjoin mapping to states where a tiny area of the hole is present in the GADM shape
    # holes_mapped_intersect = holes_mapped.copy()
    # holes_mapped_intersect['geometry'] = holes_mapped.apply(lambda x: x.geometry.intersection(df_gadm_usa.loc[df_gadm_usa['GID_1'] == x['GID_1']].iloc[0].geometry), axis=1)
    # holes_mapped_intersect['area'] = holes_mapped_intersect.to_crs(6372).area
    # holes_mapped_intersect_filter = holes_mapped_intersect.loc[holes_mapped_intersect['area'] > 200]
    # holes_mapped_intersect_filter['GADM_ID'] = (np.arange(0,len(holes_mapped_intersect_filter),1))
    # holes_mapped_intersect_filter['GADM_ID'] = holes_mapped_intersect_filter['GADM_ID'].astype('str')
    # holes_mapped_intersect_filter['country'] = 'US'
    # holes_mapped_intersect_filter['State'] = holes_mapped_intersect_filter.apply(lambda x: x['HASC_1'].split('.')[1],axis=1)    
    # save_map(holes_mapped, filename="Holes_intersect.html", color=False, cmap=False)

    build_shapes.add_population_data(holes_mapped_intersect_filter,['US'],'standard',nprocesses=nprocesses)
    df_gadm_usa['State'] = df_gadm_usa.apply(lambda x: x['ISO_1'].split('-')[1], axis=1)
    df_gadm_usa['country'] = 'US'
    df_gadm_usa['GADM_ID'] = df_gadm_usa['GID_1']
    build_shapes.add_population_data(df_gadm_usa,['US'],'standard',nprocesses=nprocesses)


    # Initial error percentages of unmet demand
    df_error = pd.DataFrame()
    df_per_capita_cons = pd.DataFrame()
    df_error = calc_percentage_unmet_demand_by_state(df_erst_gpd.rename(columns={'STATE':'State'}), df_demand_utility, df_error, 'Initial')
    df_per_capita_cons = calc_per_capita_kWh_state(df_erst_gpd.rename(columns={'STATE':'State'}), df_gadm_usa, df_per_capita_cons, 'Initial')

    # Missing utilities in ERST shape files
    print(df_demand_utility)
    missing_utilities = (list(set(df_demand_utility.reset_index().NAME) - set(df_erst_gpd.NAME)))
    df_missing_utilities = df_demand_utility.query('NAME in @missing_utilities')
    df_utilities_grouped_state = df_missing_utilities.groupby('State')['Sales (Megawatthours)'].sum()

    # Compute centroid of the holes
    holes_centroid = holes_mapped_intersect_filter.copy()
    holes_centroid.geometry = holes_mapped_intersect_filter.geometry.centroid
    holes_centroid['State'] = holes_centroid.apply(lambda x: x['HASC_1'].split('.')[1],axis=1)

    df_final = compute_demand_disaggregation(holes_mapped_intersect_filter, holes_centroid, df_utilities_grouped_state, df_demand_utility, df_gadm_usa)

    df_final = df_final._append(df_erst_gpd.rename(columns={'STATE':'State'}))

    # error percentages of unmet demand after assigning average demand to states
    df_error = calc_percentage_unmet_demand_by_state(df_final, df_demand_utility, df_error, 'Mid-way')
    df_per_capita_cons = calc_per_capita_kWh_state(df_final, df_gadm_usa, df_per_capita_cons, 'Mid-way')

    df_final = rescale_demands(df_final, df_demand_utility, df_utilities_grouped_state)

    # Final error percentages of unmet demand after rescaling
    df_error = calc_percentage_unmet_demand_by_state(df_final, df_demand_utility, df_error, 'Final')
    df_per_capita_cons = calc_per_capita_kWh_state(df_final, df_gadm_usa, df_per_capita_cons, 'Final')

    # Per-capita consumption
    df_per_capita = pd.DataFrame()
    df_per_capita['Calculated'] = df_final.groupby('State')['Sales (Megawatthours)'].sum() * 1000 / df_gadm_usa.groupby('State')['pop'].sum() #Per capita consumption in kWh
    df_per_capita = df_per_capita.join(df_eia_per_capita)
    df_per_capita.rename(columns={2021:'EIA'}, inplace=True)

    fig = px.bar(df_per_capita, barmode='group')
    fig.update_layout(yaxis_title='Per capita consumption (kWh)', xaxis_title='State')
    fig.write_image(f"{plot_path}/per_capita_consumption.png")

    df_per_capita['error'] = (df_per_capita['Calculated'] - df_per_capita['EIA']) * 100 / df_per_capita['EIA']
    df_per_capita['error'] = df_per_capita['error'].abs()

    fig = px.bar(df_per_capita, y = 'error')
    fig.update_layout(yaxis_title='Error %', xaxis_title='State')
    fig.write_image(f"{plot_path}/per_capita_error.png")

    geo_df_final = gpd.GeoDataFrame(df_final, geometry='geometry')
    geo_df_final['Sales (TWh)'] = geo_df_final['Sales (Megawatthours)'] / 1e6
    # geo_df_final['per capita'] = geo_df_final['Sales (Megawatthours)'] / geo_df_final['population']
    # Plot the GeoDataFrames
    geo_df_final = geo_df_final.drop(columns=['SOURCEDATE','VAL_DATE'],axis=1)
    m = geo_df_final.explore(column='Sales (TWh)',cmap='jet')
    m.save(f"{plot_path}/demand_filled_TWh_USA.html")

    df_erst_gpd['Sales (TWh)'] = df_erst_gpd['Sales (Megawatthours)'] / 1e6
    df_erst_gpd = df_erst_gpd.drop(columns=['SOURCEDATE','VAL_DATE'],axis=1)

    m = df_erst_gpd.explore(column='Sales (TWh)',cmap='jet')
    m.save(f"{plot_path}/demand_with_holes_TWh_USA.html")

    # df_final.to_file(f"Demand_mapped.geojson",driver="GeoJSON")