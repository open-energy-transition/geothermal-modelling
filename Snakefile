from snakemake.utils import min_version

min_version("6.0")

import sys
import pathlib

sys.path.append("workflow/pypsa-earth")


rule copy_custom_powerplants:
    input:
        "data/custom_powerplants_eia.csv",
    output:
        "workflow/pypsa-earth/data/custom_powerplants.csv",
    shell:
        "cp {input} {output}"


rule build_custom_powerplants:
    input:
        eia_generators_data_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "powerplant_data",
            "EIA_generators",
            "eia8602021",
            "3_1_Generator_Y2021.xlsx",
        ),
        eia_plants_data_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "powerplant_data",
            "EIA_generators",
            "eia8602021",
            "2___Plant_Y2021.xlsx",
        ),
        ror_custom_powerplants_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "powerplant_data",
            "custom_powerplants_ror.csv",
        ),
    output:
        output_file_name="custom_powerplants_eia.csv"
    script:
        "analysis/scripts/build_custom_powerplants.py"


rule retrieve_data:
    params:
        gdrive_url="https://drive.google.com/drive/folders/1sWDPC1EEzVtgixBb8C-OqZiEX3dmTOec",
        cookies_path=pathlib.Path(".cache", "gdown"),
        output_path=pathlib.Path("analysis", "gdrive_data", "data"),
        delta_months=5,
    script:
        "analysis/scripts/download_from_gdrive.py"


rule network_comparison:
    params:
        plot_network_topology=True,  # Boolean: plot the network topology
        plot_network_crossings=True,  # Boolean: plot the network crossings
        plot_network_capacity_ipm=True,  # Boolean: plot the network capacity for the PyPSA vs IPM case
        plot_network_capacity_reeds=True,  # Boolean: plot the network capacity for the PyPSA vs reeds case
    input:
        base_network_pypsa_earth_path=pathlib.Path(
            "workflow", "pypsa-earth", "networks", "US_2021", "base.nc"
        ),
        base_network_pypsa_usa_path=pathlib.Path(
            "analysis", "gdrive_data", "data", "pypsa_usa", "lines_gis.csv"
        ),
        eia_base_network_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "transmission_grid_data",
            "US_electric_transmission_lines_original.geojson",
        ),
        gadm_shapes_path=pathlib.Path(
            "analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json"
        ),
        ipm_shapes_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "shape_files",
            "ipm_v6_regions",
            "IPM_Regions_201770405.shp",
        ),
        lines_osm_raw_path=pathlib.Path(
            "workflow",
            "pypsa-earth",
            "resources",
            "US_2021",
            "osm",
            "raw",
            "all_raw_lines.geojson",
        ),
        lines_osm_clean_path=pathlib.Path(
            "workflow",
            "pypsa-earth",
            "resources",
            "US_2021",
            "osm",
            "clean",
            "all_clean_lines.geojson",
        ),
        reeds_shapes_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "pypsa_usa",
            "Reeds_Shapes",
            "rb_and_ba_areas.shp",
        ),
    script:
        "analysis/scripts/network_comparison.py"


rule installed_capacity_comparison:
    params:
        year_for_comparison=2021,  #Should this be 2021?
        plot_country_comparison=True,  # Boolean: plot the countrywide generation comparison
        plot_state_by_state_comparison=True,  # Boolean: plot the state-by-state generation comparison
        plot_spatial_representation=True,  # Boolean: plot the map with the installed capacity per node
        state_to_omit=["AK", "HI"],
    input:
        eia_installed_capacity_path=pathlib.Path(
            "analysis", "gdrive_data", "data", "powerplant_data", "capacities_eia.xlsx"
        ),
        eia_state_temporal_installed_capacity_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "powerplant_data",
            "existcapacity_annual.xlsx",
        ),
        eia_raw_reference_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "powerplant_data",
            "custom_powerplants_eia_with_state.csv",
        ),
        gadm_shapes_path=pathlib.Path(
            "analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json"
        ),
        pypsa_earth_network_path=pathlib.Path(
            "workflow",
            "pypsa-earth",
            "results",
            "US_2021",
            "networks",
            "elec_s_50_ec_lcopt_Co2L-200H.nc",
        ),
        #pypsa_earth_network_nonac_path = pathlib.Path() #will be added later
    script:
        "analysis/scripts/installed_capacity_comparison.py"


rule map_network_to_gadm:
    input:
        gadm_shapes_path=pathlib.Path(
            "analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json"
        ),
        pypsa_earth_network_path=pathlib.Path(
            "workflow",
            "pypsa-earth",
            "networks",
            "US_2021",
            "elec_s.nc",
        ),
    output:
        mapped_network_output_file_name="elec_s_gadm_mapped.nc",
    script:
        "analysis/scripts/map_network_to_gadm.py"


rule generation_comparison:
    params:
        year_for_comparison=2021,  #Should this be 2021?
        plot_country_comparison=True,  # Boolean: plot the countrywide generation comparison
        plot_state_by_state_comparison=True,  # Boolean: plot the state-by-state generation comparison
    input:
        eia_country_generation_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "electricity_generation_data",
            "generation_eia.csv",
        ),
        eia_state_generation_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "electricity_generation_data",
            "EIA_statewise_data",
            "use_all_phy.xlsx",
        ),
        gadm_shapes_path=pathlib.Path(
            "analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json"
        ),
        pypsa_earth_network_path=pathlib.Path(
            "workflow",
            "pypsa-earth",
            "results",
            "US_2021",
            "networks",
            "elec_s_10_ec_lcopt_Co2L-25H.nc",
        ),
    script:
        "analysis/scripts/generation_comparison.py"


rule preprocess_demand_data:
    params:
        demand_year=2021,
        holes_area_threshold=100,  # to ignore holes smaller than this area in sq.km
        nprocesses=1,
    input:
        demand_utility_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "electricity_demand_data",
            "demand_data",
            "table_10_EIA_utility_sales.xlsx",
        ),
        country_gadm_path=pathlib.Path(
            "workflow",
            "pypsa-earth",
            "resources",
            "US_2021",
            "shapes",
            "country_shapes.geojson",
        ),
        #erst_path = pathlib.Path("analysis", "gdrive_data", "data", "electricity_demand_data", "demand_data", "ERST_overlay_demand.geojson"),
        erst_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "electricity_demand_data",
            "demand_data",
            "Electric_Retail_Service_Territories.geojson",
        ),
        gadm_usa_path=pathlib.Path(
            "analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json"
        ),
        eia_per_capita_path=pathlib.Path(
            "analysis",
            "gdrive_data",
            "data",
            "electricity_demand_data",
            "use_es_capita.xlsx",
        ),
    script:
        "analysis/scripts/preprocess_demand_data.py"


# rule demand_modelling:
#   script:
