from snakemake.utils import min_version
min_version("6.0")

import sys
import pathlib

sys.path.append("workflow/pypsa-earth")

rule copy_custom_powerplants:
    input: 
        "data/custom_powerplants.csv"
    output:
        "workflow/pypsa-earth/data/custom_powerplants.csv"
    shell:
        "cp {input} {output}"

rule retrieve_data:
    params:
        gdrive_url = "https://drive.google.com/drive/folders/1LwSoZDtnyUx5ki9SmBvdlGW3QwWjs4rA?usp=drive_link",
        output_path = "analysis/gdrive_data/"
    script:
        "analysis/scripts/download_from_gdrive.py"


rule network_comparison:
    params:
        plot_network_topology=True, # Boolean: plot the network topology
        plot_network_crossings=True, # Boolean: plot the network crossings
        plot_network_capacity_ipm=True, # Boolean: plot the network capacity for the PyPSA vs IPM case
        plot_network_capacity_reeds=True # Boolean: plot the network capacity for the PyPSA vs reeds case
    input:
        base_network_pypsa_earth_path=pathlib.Path("workflow", "pypsa-earth", "networks", "US_2021", "base.nc"),
        base_network_pypsa_usa_path=pathlib.Path("analysis", "gdrive_data", "data", "pypsa_usa", "lines_gis.csv"),
        eia_base_network_path=pathlib.Path("analysis", "gdrive_data", "data", "transmission_grid_data", "US_electric_transmission_lines_original.geojson"),
        gadm_shapes_path=pathlib.Path("analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json"),
        ipm_shapes_path=pathlib.Path("analysis", "gdrive_data", "data", "shape_files", "ipm_v6_regions", "IPM_Regions_201770405.shp"),
        lines_osm_raw_path=pathlib.Path("workflow", "pypsa-earth", "resources", "US_2021", "osm", "raw", "all_raw_lines.geojson"),
        lines_osm_clean_path=pathlib.Path("workflow", "pypsa-earth", "resources", "US_2021", "osm", "clean", "all_clean_lines.geojson"),
        reeds_shapes_path=pathlib.Path("analysis", "gdrive_data", "data", "pypsa_usa", "Reeds_Shapes", "rb_and_ba_areas.shp")
    script:
        "analysis/scripts/network_comparison.py"


rule installed_capacity_comparison:
    params:
        year_for_comparison=2020, #Should this be 2021?
        plot_country_comparison=True, # Boolean: plot the countrywide generation comparison
        plot_state_by_state_comparison=False, # Boolean: plot the state-by-state generation comparison
    input:
        eia_installed_capacity_path=pathlib.Path("analysis", "gdrive_data", "data", "powerplant_data", "capacities_eia.xlsx"),
        gadm_shapes_path=pathlib.Path("analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json"),
        pypsa_earth_network_path=pathlib.Path("workflow", "pypsa-earth", "results", "US_2021", "networks", "elec_s_10_ec_lcopt_Co2L-25H.nc")
    script:
        "analysis/scripts/installed_capacity_comparison.py"


rule generation_comparison:
    params:
        year_for_comparison=2020, #Should this be 2021?
        plot_country_comparison=True, # Boolean: plot the countrywide generation comparison
        plot_state_by_state_comparison=True, # Boolean: plot the state-by-state generation comparison
    input:
        eia_country_generation_path=pathlib.Path("analysis", "gdrive_data", "data", "electricity_generation_data", "generation_eia.csv"),
        eia_state_generation_path=pathlib.Path("analysis", "gdrive_data", "data", "electricity_generation_data", "EIA_statewise_data", "use_all_phy.xlsx"),
        gadm_shapes_path=pathlib.Path("analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json"),
        pypsa_earth_network_path=pathlib.Path("workflow", "pypsa-earth", "results", "US_2021", "networks", "elec_s_10_ec_lcopt_Co2L-25H.nc")
    script:
        "analysis/scripts/generation_comparison.py"

rule preprocess_demand_data:
    params:
        demand_year = 2021,
        holes_area_threshold = 0.05, # to ignore holes smaller than this area
        nprocesses = 1
    input:
        demand_utility_path = pathlib.Path("analysis", "gdrive_data", "data", "electricity_demand_data", "demand_data","table_10_EIA_utility_sales.xlsx"),
        country_gadm_path = pathlib.Path("workflow", "pypsa-earth", "resources", "US_2021", "shapes", "country_shapes.geojson"),
        erst_path = pathlib.Path("analysis", "gdrive_data", "data", "electricity_demand_data", "demand_data", "ERST_overlay_demand.geojson"),
        gadm_usa_path = pathlib.Path("analysis", "gdrive_data", "data", "shape_files", "gadm41_USA_1.json"),
        eia_per_capita_path = pathlib.Path("analysis", "gdrive_data", "data", "electricity_demand_data", "use_es_capita.xlsx")
    script:
        "analysis/scripts/preprocess_demand_data.py"

#rule demand_modelling:
 #   script:



