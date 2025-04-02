"""
Preprocesses utility level demand data from EIA

Relevant Settings
-----------------

.. code:: yaml

    US:
        demand_modelling:
            demand_year:
            holes_area_threshold:
            nprocesses:
            plotting:

.. seealso::
    Documentation of the configuration file ``config.usa_baseline.yaml``


Inputs
------

- ``shapes/country_shapes.geojson``

- ``analysis/gdrive_data/data/electricity_demand_data/EIA930_{demand_year}_Jan_Jun_opt.csv``

- ``analysis/gdrive_data/data/electricity_demand_data/EIA930_{demand_year}_Jul_Dec_opt.csv``

- ``analysis/gdrive_data/data/shape_files/Balancing_Authorities.geojson``

- ``analysis/outputs/demand_modelling/ERST_mapped_demand_centroids.geojson``

- ``networks/base.nc``

Outputs
-------

- ``analysis/outputs/demand_modelling/ERST_mapped_demand_centroids.geojson``


Description
-----------
"""

import pathlib
import datetime as dt
import geopandas as gpd
import pandas as pd
import os
from shapely.validation import make_valid
import numpy as np
import plotly.express as px
import sys
from _helpers_usa import get_colors


def parse_inputs(default_path, demand_year):
    """
    Load all input data required for preprocessing and computing utility level demand data

    Parameters
    ----------
    default_path: str
        current total directory path
    demand_year: int
        Year for which electrical demand has been modelled

    Returns
    -------
    df_demand_utility: pandas dataframe
        Utility level demand data from EIA
    df_erst_gpd: geopandas dataframe
        Utility level shapefile
    df_country: geopandas dataframe
        country level shapefile
    df_gadm_usa: geopandas dataframe
        GADM regions for USA
    df_eia_per_capita: pandas dataframe
        per capita electricity demand by state from EIA
    df_additional_demand_data: pandas dataframe
        excess demand at state level not covered in "df_demand_utility"
    """
    # Load data
    demand_utility_path = pathlib.Path(
        default_path, snakemake.input.demand_utility_path
    )
    erst_gpd_path = pathlib.Path(default_path, snakemake.input.erst_path)
    country_gadm_path = pathlib.Path(default_path, snakemake.input.country_gadm_path)
    gadm_shape_usa_path = pathlib.Path(default_path, snakemake.input.gadm_usa_path)
    eia_per_capita_filepath = pathlib.Path(
        default_path, snakemake.input.eia_per_capita_path
    )
    additional_demand_filepath = pathlib.Path(
        default_path, snakemake.input.additional_demand_data_path
    )

    df_demand_utility = pd.read_excel(demand_utility_path, skiprows=2)
    df_erst_gpd = gpd.read_file(erst_gpd_path)
    df_country = gpd.read_file(country_gadm_path)
    df_gadm_usa = gpd.read_file(gadm_shape_usa_path)
    df_eia_per_capita = pd.read_excel(
        eia_per_capita_filepath,
        sheet_name="Total per capita",
        skiprows=2,
        index_col="State",
    )
    # Per capita data available only until 2022
    try:
        df_eia_per_capita = df_eia_per_capita[demand_year]
    except ValueError:
        df_eia_per_capita = df_eia_per_capita[2022]

    # The utility level data adds upto 3294 TWh of sales in 2021 whilst at the state level data adds upto 3800 TWh
    # The function is used to add data from another EIA table containing sales at the state level
    df_additional_demand = pd.read_excel(additional_demand_filepath, skiprows=2)
    df_additional_demand = df_additional_demand.query("Year == @demand_year")
    df_additional_demand.rename(columns={"STATE": "State"}, inplace=True)
    df_additional_demand.set_index("State", inplace=True)
    df_additional_demand_data = df_additional_demand[
        "Megawatthours.4"
    ]  # Refers to total sales (residential+commercial+industrial+transport)
    df_additional_demand_data = df_additional_demand_data.drop("US")
    df_demand_grouped = df_demand_utility.groupby("State")[
        "Sales (Megawatthours)"
    ].sum()
    df_additional_demand_data = df_additional_demand_data - df_demand_grouped

    log_output_file.write("Reading input files completed \n")
    log_output_file.write("-------------------------------------------------\n")

    return (
        df_demand_utility,
        df_erst_gpd,
        df_country,
        df_gadm_usa,
        df_eia_per_capita,
        df_additional_demand_data,
    )


def compute_demand_disaggregation(
    holes_mapped_intersect_filter,
    holes_centroid,
    df_missing_data_state,
    df_demand_utility,
    df_gadm_usa,
):
    """
    Compute demand allocation for holes in each state based on per capita electrical energy consumption

    Parameters
    ----------
    holes_mapped_intersect_filter: geopandas dataframe
        shape geometries of holes mapped to respective states
    holes_centroid:
        centroids of holes
    df_missing_data_state: pandas dataframe
        missing demand data at state level with reference to EIA state level data
    df_demand_utility: pandas dataframe
        Utility level demand data from EIA
    df_gadm_usa: geopandas dataframe
        GADM regions for USA
    Returns
    -------
    df_final: geopandas dataframe
        Allocated portion of unmet demand to holes(geometries) based on per capita electrical demand in the state to which the hole is mapped
    """
    holes_mapped_intersect_filter["Sales (Megawatthours)"] = 0
    df_final = pd.DataFrame()
    for state in df_missing_data_state.index:
        holes_state = holes_mapped_intersect_filter.query(
            "State == @s", local_dict={"s": state}
        )
        full_state_pop = df_gadm_usa.query("State == @state")["pop"].values[0]
        state_demand = df_demand_utility.query("State == @state")[
            "Sales (Megawatthours)"
        ].sum()
        avg_per_capita_demand = state_demand / full_state_pop

        if len(holes_state) > 0:
            holes_state["Sales (Megawatthours)"] = (
                holes_state["pop"] * avg_per_capita_demand
            )
            df_final = df_final._append(holes_state, ignore_index=True)

    return df_final


def calc_percentage_unmet_demand_by_state(
    df_calc, df_ref, df_ref_additional, df_error, error_column, state_kwd
):
    """
    Calculate percentage unmet demand at various stages of the algorithm
    Parameters
    ----------
    df_calc: pandas dataframe
        Actual demand data
    df_ref: pandas dataframe
        Reference demand data
    df_ref_additional: pandas dataframe
        Additional demand data to be added to existing reference data
    df_error: pandas dataframe
        Error values in percentage at different stages of algorithm
    error_column: str
        column name
    state_kwd: str
        to use 'STATE' / 'State' as column name based on passed argument
    Returns
    -------
    df_error: pandas dataframe
        Error values in percentage at different stages of algorithm
    """
    df_calc_state = df_calc.groupby(state_kwd)["Sales (Megawatthours)"].sum()
    df_ref_state = (
        df_ref.groupby(state_kwd)["Sales (Megawatthours)"].sum() + df_ref_additional
    )
    df_error[error_column] = (df_ref_state - df_calc_state) * 100 / (df_ref_state)
    return df_error


def calc_per_capita_kWh_state(
    df_calc, df_gadm, df_per_capita_cons, column_name, state_kwd
):
    """
    Calculate per capita statewise electrical demand at various stages of the algorithm
    Parameters
    ----------
    df_calc: pandas dataframe
        Actual demand data
    df_gadm: geopandas dataframe
        GADM shapes
    df_per_capita_cons: pandas dataframe
        Per capital electricity demand at different stages of algorithm
    column_name: str
        column name
    state_kwd: str
        to use 'STATE' / 'State' as column name based on passed argument
    Returns
    -------
    df_per_capita_cons: pandas dataframe
        Per capita electricity demand at different stages of algorithm
    """
    df_calc_per_capita = (
        df_calc.groupby(state_kwd)["Sales (Megawatthours)"].sum()
        * 1000
        / df_gadm.groupby("State")["pop"].sum()
    )
    df_per_capita_cons[column_name] = df_calc_per_capita
    return df_per_capita_cons


def rescale_demands(df_final, df_demand_utility, df_additional_sales_data):
    """
    Rescale demand of all utilities in a given state to match reference value
    Parameters
    ----------
    df_final: pandas dataframe
        Actual demand data
    df_demand_utility: pandas dataframe
        Reference utility level demand data
    df_additional_sales_data: pandas dataframe
        Additional state level demand data not accounted in utility level

    Returns
    -------
    df_final: pandas dataframe
        Final rescaled demand values for all utilities and holes
    """
    df_demand_statewise = (
        df_demand_utility.groupby("State")["Sales (Megawatthours)"].sum()
        + df_additional_sales_data
    )
    df_final["rescaling_factor"] = 0
    for state in df_demand_statewise.index:
        actual_state_demand = df_demand_statewise.loc[state]
        df_filter_final = df_final.query("State == @state")
        assigned_utility_demand = df_filter_final["Sales (Megawatthours)"].sum()
        if assigned_utility_demand != 0:
            rescaling_factor = actual_state_demand / assigned_utility_demand
        else:
            rescaling_factor = 1
        df_final.loc[df_final["State"] == state, "Sales (Megawatthours)"] *= (
            rescaling_factor
        )
        df_final.loc[df_final["State"] == state, "rescaling_factor"] = rescaling_factor

    return df_final


def overlay_demands(df_erst_gpd, df_gadm_usa, df_demand_utility):
    """
    Map ERST utility shapes to utility level demand (sales) data
    Parameters
    ----------
    df_erst_gpd: geopandas dataframe
        Utilty level shapes
    df_gadm_usa: geopandas dataframe
        GADM shapes
    df_demand_utility: pandas dataframe
        Utility level demand data from EIA

    Returns
    -------
    df_erst_overlay: geopandas dataframe
        Utility shapes mapped to their respective demands
    """
    # ERST overlay with GADM shapes for intersection
    df_erst_overlay = gpd.overlay(df_erst_gpd, df_gadm_usa, how="intersection")
    df_erst_overlay["STATE"] = df_erst_overlay.apply(
        lambda x: x["ISO_1"].split("-")[1], axis=1
    )
    df_erst_overlay.set_index(["NAME", "STATE"], inplace=True)

    # Demand mapping to ERST shape file
    df_demand_utility["Entity"] = df_demand_utility["Entity"].apply(str.upper)
    df_demand_utility.rename(columns={"Entity": "NAME", "State": "STATE"}, inplace=True)
    df_demand_utility.set_index(["NAME", "STATE"], inplace=True)

    df_erst_overlay = df_erst_overlay.join(df_demand_utility)

    # Make geometry valid
    df_erst_overlay["geometry"] = df_erst_overlay["geometry"].apply(make_valid)

    return df_erst_overlay


# total_demand: float -> Total Demand (Sales) as per the EIA data
# mapped_demand: float -> The demand mapped to ERST file at that particular stage
def compute_missing_percentage(total_demand, mapped_demand):
    return np.round(((total_demand - mapped_demand) * 100 / total_demand), 2)


def save_map(df_map, filename, color, cmap, cmap_col=""):
    """
    Save US map as html file
    Parameters
    ----------
    df_map : geopandas
        geopandas file that is to be plotted on the map
    filename : str
        Name to be given to the saved plot (HTML file)
    color : boolean
        If True, each row of geometry is coloured based on the color assigned to it the gpd dataframe
    cmap : boolean
        If True, each row of geometry is coloured based on the magnitude of the value in the cmap_col specified
    cmap_col : str
        Column on which a color map is used to plot the geometry

    """
    if "VAL_DATE" in df_map.columns:
        df_map = df_map.drop("VAL_DATE", axis=1)
    if "SOURCEDATE" in df_map.columns:
        df_map = df_map.drop("SOURCEDATE", axis=1)
    if not color and not cmap:
        m = df_map.explore()
    elif color and not cmap:
        m = df_map.explore(color=df_map["color"].tolist())
    elif not color and cmap:
        m = df_map.explore(column=cmap_col, cmap="jet")
    m.save(os.path.join(plot_path, filename))


def map_demands_utilitywise(
    df_demand_utility,
    df_erst_gpd,
    df_country,
    df_gadm_usa,
    df_eia_per_capita,
    df_additional_demand_data,
    log_output_file,
    demand_year,
    plotting,
    geo_crs,
    area_crs,
):
    """
    Map ERST utility shapes to utility level demand (sales) data
    Identify holes in the country geometry
    Calculate population at utility level
    Map portion of missing demands to holes based on per capita electricity consumption
    Rescale demands

    Parameters
    ----------
    df_demand_utility: pandas dataframe
        Utility level demand data from EIA
    df_erst_gpd: geopandas dataframe
        Utilty level shapes
    df_country: geopandas dataframe
        Country level shapes
    df_gadm_usa: geopandas dataframe
        GADM shapes
    df_eia_per_capita: pandas dataframe
        per capita statewise electricity demands
    df_additional_demand_data: pandas dataframe:
        statewise demand data - accounts for some demand not included at utility level
    log_output_file: file object
        logger file
    demand_year: int
        year for which electrical demand is modelled
    plotting: boolean
        If set to true - html plots are generated else no
    distance_crs: CRS code
        Co-ordinate reference system to retrieve coordinate in 'm'
    area_crs: CRS code
        Co-ordinate reference system to retrieve coordinates in 'sq.m'
    Returns
    -------
    df_final: geopandas dataframe
        Utility shapes and holes mapped to their respective demands and rescaled to match reference values
    """

    total_demand = df_demand_utility["Sales (Megawatthours)"].sum() / 1e6
    log_output_file.write(f"Total sales (TWh) as in EIA sales data: {total_demand} \n")

    # Add system paths
    pypsa_earth_scripts_path = pathlib.Path("workflow", "pypsa-earth", "scripts")
    sys.path.insert(0, str(pypsa_earth_scripts_path))
    import build_shapes

    # Make geometry valid
    df_erst_gpd["geometry"] = df_erst_gpd["geometry"].apply(make_valid)

    if plotting:
        # Plot initial ERST shapefile
        save_map(df_erst_gpd, filename="ERST_initial.html", color=False, cmap=False)
        log_output_file.write("Generated initial ERST shapes file \n")

    # Map demands to ERST shapefile
    df_erst_gpd = overlay_demands(df_erst_gpd, df_gadm_usa, df_demand_utility)
    demand_mapped = df_erst_gpd["Sales (Megawatthours)"].sum() / 1e6
    log_output_file.write(
        f"Total sales (TWh) as in ERST mapped data: {demand_mapped} \n"
    )
    missing_demand_percentage = compute_missing_percentage(total_demand, demand_mapped)
    log_output_file.write(f"Missing sales (%) : {missing_demand_percentage} \n")

    if plotting:
        save_map(
            df_erst_gpd,
            filename="ERST_overlay.html",
            color=False,
            cmap=True,
            cmap_col="Sales (Megawatthours)",
        )
        log_output_file.write("Generated ERST shapes file with mapped demands \n")

    ## Obtain holes in ERST shape files

    # Merge all geometry into one geometry
    df_erst_gpd_dissolved = df_erst_gpd.dissolve()
    # Find the difference in geometry between country and the merged ERST geometry
    holes = df_country.difference(df_erst_gpd_dissolved)

    # Convert holes geometry : multipolygons into polygons to create a separate row for each hole
    holes_exploded = holes.explode()
    holes_exploded = gpd.GeoDataFrame(geometry=holes.explode(), crs=df_erst_gpd.crs)

    if plotting:
        save_map(holes_exploded, filename="Holes.html", color=False, cmap=False)
        log_output_file.write("Generated holes exploded map \n")
    holes_exploded = holes_exploded.to_crs(area_crs)
    holes_exploded["Area"] = holes_exploded.area / 1e6

    # Filtering out holes with very small areas (only hole areas larger than area_threshold considered)
    holes_exploded_filter = holes_exploded.query("Area > @holes_area_threshold")

    if plotting:
        save_map(
            holes_exploded_filter,
            filename="Holes_considered.html",
            color=False,
            cmap=False,
        )
        log_output_file.write(f"Generated holes greater than {holes_area_threshold} \n")
    holes_exploded_filter = holes_exploded_filter.to_crs(geo_crs)

    df_gadm_usa["color"] = get_colors(len(df_gadm_usa))
    holes_mapped = gpd.overlay(holes_exploded_filter, df_gadm_usa, how="intersection")

    if plotting:
        save_map(
            holes_mapped, filename="Holes_mapped_GADM.html", color=True, cmap=False
        )
        log_output_file.write("Generated holes mapped to GADM \n")

    # # Compute intersecting areas of holes and states
    holes_mapped_intersect = holes_mapped.copy()
    holes_mapped_intersect["geometry"] = holes_mapped.apply(
        lambda x: x.geometry.intersection(
            df_gadm_usa.loc[df_gadm_usa["GID_1"] == x["GID_1"]].iloc[0].geometry
        ),
        axis=1,
    )
    holes_mapped_intersect["area"] = holes_mapped_intersect.area
    holes_mapped_intersect_filter = holes_mapped_intersect.copy()

    holes_mapped_intersect_filter["GADM_ID"] = np.arange(
        0, len(holes_mapped_intersect_filter), 1
    )
    holes_mapped_intersect_filter["GADM_ID"] = holes_mapped_intersect_filter[
        "GADM_ID"
    ].astype("str")
    holes_mapped_intersect_filter["country"] = "US"
    holes_mapped_intersect_filter["State"] = holes_mapped_intersect_filter.apply(
        lambda x: x["HASC_1"].split(".")[1], axis=1
    )
    if plotting:
        save_map(holes_mapped, filename="Holes_intersect.html", color=False, cmap=False)

    build_shapes.add_population_data(
        holes_mapped_intersect_filter, ["US"], "standard", nprocesses=nprocesses
    )
    df_gadm_usa["State"] = df_gadm_usa.apply(lambda x: x["ISO_1"].split("-")[1], axis=1)
    df_gadm_usa["country"] = "US"
    df_gadm_usa["GADM_ID"] = df_gadm_usa["GID_1"]
    build_shapes.add_population_data(
        df_gadm_usa, ["US"], "standard", nprocesses=nprocesses
    )

    # Initial error percentages of unmet demand
    df_error = pd.DataFrame()
    df_per_capita_cons = pd.DataFrame()
    df_erst_gpd = df_erst_gpd.reset_index()
    df_demand_utility = df_demand_utility.reset_index()
    df_error = calc_percentage_unmet_demand_by_state(
        df_erst_gpd,
        df_demand_utility,
        df_additional_demand_data,
        df_error,
        "Initial",
        "STATE",
    )
    df_per_capita_cons = calc_per_capita_kWh_state(
        df_erst_gpd, df_gadm_usa, df_per_capita_cons, "Initial", "STATE"
    )

    # Missing utilities in ERST shape files
    df_demand_utility = df_demand_utility.reset_index()
    df_demand_utility.rename(columns={"STATE": "State"}, inplace=True)
    df_missing_utilities = df_demand_utility.query("NAME in @missing_utilities")
    df_utilities_grouped_state = df_missing_utilities.groupby("State")[
        "Sales (Megawatthours)"
    ].sum()

    # Compute centroid of the holes
    holes_centroid = holes_mapped_intersect_filter.copy()
    holes_centroid.geometry = holes_mapped_intersect_filter.geometry.centroid
    holes_centroid["STATE"] = holes_centroid.apply(
        lambda x: x["HASC_1"].split(".")[1], axis=1
    )

    df_final = compute_demand_disaggregation(
        holes_mapped_intersect_filter,
        holes_centroid,
        df_utilities_grouped_state,
        df_demand_utility,
        df_gadm_usa,
    )

    df_final = df_final._append(df_erst_gpd.rename(columns={"STATE": "State"}))

    # error percentages of unmet demand after assigning average demand to states
    df_error = calc_percentage_unmet_demand_by_state(
        df_final,
        df_demand_utility,
        df_additional_demand_data,
        df_error,
        "Mid-way",
        "State",
    )
    df_per_capita_cons = calc_per_capita_kWh_state(
        df_final, df_gadm_usa, df_per_capita_cons, "Mid-way", "State"
    )

    df_final = rescale_demands(df_final, df_demand_utility, df_additional_demand_data)

    # Final error percentages of unmet demand after rescaling
    df_error = calc_percentage_unmet_demand_by_state(
        df_final,
        df_demand_utility,
        df_additional_demand_data,
        df_error,
        "Final",
        "State",
    )
    df_per_capita_cons = calc_per_capita_kWh_state(
        df_final, df_gadm_usa, df_per_capita_cons, "Final", "State"
    )

    fig = px.bar(df_error, barmode="group")
    fig.update_layout(yaxis_title="Error %", xaxis_title="State")
    fig.write_image(f"{plot_path}/unmet_demand_error_stages.png")

    fig = px.bar(df_per_capita_cons, barmode="group")
    fig.update_layout(yaxis_title="Per capita consumption (kWh)", xaxis_title="State")
    fig.write_image(f"{plot_path}/per_capita_consumption_stages.png")

    # Per-capita consumption
    df_per_capita = pd.DataFrame()
    df_per_capita["Calculated"] = (
        df_final.groupby("State")["Sales (Megawatthours)"].sum()
        * 1000
        / df_gadm_usa.groupby("State")["pop"].sum()
    )  # Per capita consumption in kWh
    df_per_capita = df_per_capita.join(df_eia_per_capita)
    if demand_year in df_per_capita.columns:
        df_per_capita.rename(columns={demand_year: "EIA"}, inplace=True)
    else:
        df_per_capita.rename(columns={2022: "EIA"}, inplace=True)

    fig = px.bar(df_per_capita, barmode="group")
    fig.update_layout(yaxis_title="Per capita consumption (kWh)", xaxis_title="State")
    fig.write_image(f"{plot_path}/per_capita_consumption.png")

    df_per_capita["error"] = (
        (df_per_capita["Calculated"] - df_per_capita["EIA"])
        * 100
        / df_per_capita["EIA"]
    )
    df_per_capita["error"] = df_per_capita["error"].abs()

    fig = px.bar(df_per_capita, y="error")
    fig.update_layout(yaxis_title="Error %", xaxis_title="State")
    fig.write_image(f"{plot_path}/per_capita_error.png")

    geo_df_final = gpd.GeoDataFrame(df_final, geometry="geometry")
    geo_df_final["Sales (TWh)"] = geo_df_final["Sales (Megawatthours)"] / 1e6

    # Plot the GeoDataFrames
    if plotting:
        save_map(
            geo_df_final,
            filename="demand_filled_TWh_USA.html",
            color=False,
            cmap=True,
            cmap_col="Sales (TWh)",
        )
        log_output_file.write("Plotted demand_filled_TWh_USA \n ")

    # Final error in demand mapping
    demand_mapped = geo_df_final["Sales (TWh)"].sum()
    log_output_file.write(
        f"Total sales (TWh) as in ERST mapped data: {demand_mapped} \n"
    )
    missing_demand_percentage = compute_missing_percentage(total_demand, demand_mapped)
    log_output_file.write(f"Missing sales (%) : {missing_demand_percentage} \n")

    return df_final


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("preprocess_demand_data")

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "demand_modelling")
    plot_path = pathlib.Path(default_path, "analysis", "plots", "demand_modelling")
    os.makedirs(log_path, exist_ok=True)
    os.makedirs(plot_path, exist_ok=True)
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(
        log_path, f"output_preprocess_demand_{today_date[:10]}.txt"
    )
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    # Load snakemake params
    demand_year = snakemake.params.demand_year
    holes_area_threshold = snakemake.params.holes_area_threshold
    nprocesses = snakemake.params.nprocesses
    plotting = snakemake.params.plotting
    geo_crs = snakemake.params.geo_crs
    area_crs = snakemake.params.area_crs

    log_output_file.write("Loading snakemake parameters \n")
    log_output_file.write(f"demand_year = {demand_year} \n")
    log_output_file.write(f"holes_area_threshold = {holes_area_threshold} \n")
    log_output_file.write("-------------------------------------------------\n")

    (
        df_demand_utility,
        df_erst_gpd,
        df_country,
        df_gadm_usa,
        df_eia_per_capita,
        df_additional_demand_data,
    ) = parse_inputs(default_path, demand_year)

    df_final = map_demands_utilitywise(
        df_demand_utility,
        df_erst_gpd,
        df_country,
        df_gadm_usa,
        df_eia_per_capita,
        df_additional_demand_data,
        log_output_file,
        demand_year,
        plotting,
        geo_crs,
        area_crs,
    )

    df_final.to_file(snakemake.output.utility_demand_path, driver="GeoJSON")

    log_output_file.close()
