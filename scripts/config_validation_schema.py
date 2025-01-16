# coding=utf-8

##################################################
#                                                #
# Author: Aisling Pigott                         #
# Email: aisling.pigott@openenergytransition.org #
# Version: 0.1                                   #
# Date: 23.07.2024                               #
#                                                #
# Created for Open Energy Transition GmbH        #
#                                                #
##################################################

#########################################################################
#                                                                       #
# IMPORTANT                                                             #
#                                                                       #
# This software is distributed without any warranty.                    #
#                                                                       #
# Neither the author nor Open Energy Transition GmbH                    #
# are liable for any damage caused directly or indirectly by the use    #
# or misuse of this software.                                           #
#                                                                       #
#########################################################################

#############################################################
#                                                           #
# CHANGELOG                                                 #
#                                                           #
# 0.1 - In progress                                         #
#                                                           #
#############################################################

from schema import Schema, And, Or, Use, Optional, SchemaError
import yaml
import os

SUPPORTED_PROJECTIONS = ["EPSG:4326", "EPSG:3857", "ESRI:54009"]
SUPPORTED_AC_LINE_TYPES = 
    [
        "243-AL1/39-ST1A 20.0",
        "Al/St 240/40 2-bundle 220.0",
        "Al/St 240/40 3-bundle 300.0",
        "Al/St 240/40 4-bundle 380.0",
        "Al/St 240/40 4-bundle 380.0",
        "Al/St 560/50 4-bundle 750.0"
    ]
SUPPORED_DC_LINE_TYPES = [
        "HVDC XLPE 1000"
    ]


datetime_schema = Schema(lamda d: datetime.datetime.strptime(d, "%Y-%m-%d"))

schema = Schema(
    {
        Optional('version'):And(str), 
        Optional('tutorial'):And(bool), 
        Optional('logging'):And(dict), 
        Optional('countries'):And(list), 
        Optional('enable'):And(dict), 
        Optional('custom_rules'):And([Use(open)]), 
        Optional('run'):run_schema, 
        Optional('scenario'):scenario_schema, 
        Optional('snapshots'):snapshot_schema, 
        Optional('crs'):crs_schema, 
        Optional('augmented_line_connection'):And(dict), 
        Optional('cluster_options'):cluster_schema, 
        Optional('build_shape_options'):build_shape_schema, 
        Optional('clean_osm_data_options'):clean_osm_data_schema, 
        Optional('build_osm_network'):And(dict), 
        Optional('base_network'):And(dict), 
        Optional('load_options'):And(dict), 
        Optional('electricity'):And(dict), 
        Optional('lines'):lines_schema, 
        Optional('links'):links_schema, 
        Optional('transformers'):transformer_schema, 
        Optional('atlite'):atlite_schema, 
        Optional('renewable'):renewable_schema, 
        Optional('costs'):costs_schema, 
        Optional('monte_carlo'):And(dict), 
        Optional('solving'):And(dict), 
        Optional('plotting'):And(dict)
    }
)

run_schema = {
        "name": And(str),
        "shared_cutouts": And(bool)
    }

cluster_schema = {
        'simplify_network':And(dict) , 
        'cluster_network':And(dict), 
        'alternative_clustering':And(bool), 
        'distribute_cluster':And(list), 
        'out_logging':And(bool), 
        'aggregation_strategies':Schema({str:{str:str}})
    }

build_shape_schema = Schema(
    {

    }
)

clean_osm_data_schema = Schema({

    })

scenario_schema = {
        "clusters": [int],
        "simpl": [str],
        "ll": [str],
        "opts": [str]
    }

snapshot_schema = {
        "start": datetime_schema,
        "end": datetime_schema,
        "inclusive": Or("left", "right") # not sure what this arg is
    }

crs_schema = {
       k:Or(SUPPORTED_PROJECTIONS) for k in ["geo_crs", "distance_crs", "area_crs"]
    }

augmented_line_connection_schema = {

    }

lines_schema = {
    "ac_types": Or(SUPPORTED_AC_LINE_TYPES),
    "dc_types": Or(SUPPORED_DC_LINE_TYPES),
    "s_max_pu": float, 
    "s_nom_max": float,
    "length_factor": float,
    "under_construction": Or("zero", "remove", "keep")
    }

links_schema = {
    "p_max_pu": float, 
    "p_nom_max": float,
    "under_construction": Or("zero", "remove", "keep")
    }

transformer_schema = {
    "x": float, 
    "s_nom": int, 
    "type": str
}

atlite_schema = {
    "nprocesses": int, 
    "cutouts":{
        str:{
            "module": str,
            "dx": float,
            "dy": float
        }
    }
}

RENEWABLE_RESOURCES = ["onwind","offwind-ac","offwind-dc","solar","hydro","csp"]

renewable_schema = {Optional(k):renewable_subschema for k in RENEWABLE_RESOURCES}


# it would be nice to organize these between required for all renewables and only for some types
renewable_subschema = {
    "cutout": lambda s: open(f"../cutouts/{s}.nc"),
    "resource": {
        "method": str,
        Optional("turbine"): str,
        Optional("panel"): str, 
        Optional("orientation"): str,
        Optional("slope"): float,
        Optional("azimuth"): float,
        Optional("hydrobasins"): open(str)
        Optional("flowspeed"): float,
        Optional("weight_with_height"): bool,
        Optional("show_progress"): bool
    }
    "carriers": [str]
    Optional("PHS_max_hours"): int,
    Optional("hydro_max_hours"): str,
    Optional("hydro_max_hours_default"): float
    Optional("clip_min_inflow"): float,
    "capacity_per_sqkm": lambda n: 1.7 <= n <= 4.6,
    "correction_factor": float,
    "copernicus": {
        "grid_codes": [int],
        "distance": int,
        "distance_grid_codes": [int]
    }
    "natura": bool,
    "max_depth": int,
    Optional("max_shore_distance"): int,
    "potential": Or("simple", "conservative"),
    "clip_p_max_pu": float,
    "extendable": true,
    "normalization":{
        "method": Or("hydro_capacities", "eia", False),
        Optional("year"): int # probably some bounds on the year
    }
}

costs_schema = {
    "year": int, 
    "version": str,
    "rooftop_share": float,
    "USD2013_to_EUR2013": float,
    "fill_values":{
        "FOM": float,
        "VOM": float,
        "efficiency": float,
        "fuel": float,
        "investments": float,
        "lifetime": int,
        "CO2 intensity": float,
        "discount_rate": float
    }
    "marginal_costs":{ # should these align with the renewable resources above?
        "solar": float,
        "onwind": float,
        "offwind": float,
        "hydro": float,
        "H2": float,
        "electrolysis": float,
        "fuel cell": float,
        "battery": float,
        "battery inverter": float
    }
    "emission_prices": { # only used with option ep
        "co2": float
    }
}



def validate_config_file(file="../configs/config.usa_PE.yaml": str):
    with open(file) as f:
        config = yaml.safe_load(f)

    try:
        schema.validate(config)
    except SchemaEror as e:
        print(e)
