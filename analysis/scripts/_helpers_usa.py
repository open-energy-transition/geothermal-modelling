# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import geopandas as gpd
import numpy as np
import pandas as pd
import re
import yaml
import pathlib
import os
import snakemake as sm
from pypsa.descriptors import Dict
from snakemake.script import Snakemake
import random

generators_aggregation_strategies_dict = {
    "p_nom": "sum",
    "p_nom_max": "sum",
    "p_nom_min": "sum",
    "p_min_pu": "mean",
    "marginal_cost": "mean",
    "committable": "any",
    "ramp_limit_up": "max",
    "ramp_limit_down": "max",
    "efficiency": "mean",
}


def get_colors(n):
    return ["#%06x" % random.randint(0, 0xFFFFFF) for _ in range(n)]


def get_gadm_mapping(gadm_shapes_path):
    gadm_gdp_usa = gpd.read_file(gadm_shapes_path)
    gadm_gdp_usa_state = pd.DataFrame()
    gadm_gdp_usa_state["GID_1_new"] = gadm_gdp_usa["GID_1"].str.replace("USA", "US")
    gadm_gdp_usa_state["state"] = gadm_gdp_usa["ISO_1"].str[-2:]
    usa_state_dict = dict(gadm_gdp_usa_state.values)
    return usa_state_dict


def get_state_node(gadm_shapes_path, state):
    usa_state_dict = {v: k for k, v in get_gadm_mapping(gadm_shapes_path).items()}
    return usa_state_dict[state]


def extract_time_res(network_path):
    # Search for the pattern in the filename
    match = re.search(r"(\d+)(H|SEG)", str(network_path.name))

    if match:
        number = np.float64(match.group(1))
        unit = match.group(2)

        if unit == "H":
            return number
        elif unit == "SEG":
            return 8760 / number
    else:
        return None


def config(config_path):
    path_config = pathlib.Path(config_path)
    with open(path_config) as file:
        config_dict = yaml.safe_load(file)
    return config_dict


def rename_carrier(x):
    if x == "ccgt":
        return "CCGT"
    elif x == "phs":
        return "PHS"
    else:
        return x


def eia_to_pypsa_terminology():
    eia_to_pypsa_dict = {
        "Nuclear": "nuclear",
        "Onshore Wind Turbine": "onwind",
        "Conventional Hydroelectric": "hydro",
        "Conventional Steam Coal": "coal",
        "Hydroelectric Pumped Storage": "PHS",
        "Solar Photovoltaic": "solar",
        "Geothermal": "geothermal",
        "Offshore Wind Turbine": "offwind-ac",  # check if it can be differentiated between ac & dc from EIA data
        "Wood/Wood Waste Biomass": "biomass",
        "Natural Gas Fired Combined Cycle": "CCGT",
        "Natural Gas Fired Combustion Turbine": "OCGT",
        "Natural Gas Steam Turbine": "OCGT",
        "Natural Gas Internal Combustion Engine": "OCGT",
        "Petroleum Liquids": "oil",
        "Petroleum Coke": "oil",
        "Batteries": "battery",
        "other": "All Other",
    }

    return eia_to_pypsa_dict


def mock_snakemake(
    rulename, root_dir=None, submodule_dir=None, configfile=None, **wildcards
):
    """
    This function is expected to be executed from the "scripts"-directory of "
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards**.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    configfile: str
        path to config file to be used in mock_snakemake
    wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    script_dir = pathlib.Path(__file__).parent.resolve()
    if root_dir is None:
        root_dir = script_dir.parent.parent
    else:
        root_dir = pathlib.Path(root_dir).resolve()

    user_in_script_dir = pathlib.Path.cwd().resolve() == script_dir
    if str(submodule_dir) in __file__:
        # the submodule_dir path is only need to locate the project dir
        os.chdir(pathlib.Path(__file__[: __file__.find(str(submodule_dir))]))
    elif user_in_script_dir:
        os.chdir(root_dir)
    elif pathlib.Path.cwd().resolve() != root_dir:
        raise RuntimeError(
            "mock_snakemake has to be run from the repository root"
            f" {root_dir} or scripts directory {script_dir}"
        )
    try:
        for p in sm.SNAKEFILE_CHOICES:
            if os.path.exists(p):
                snakefile = p
                break

        if isinstance(configfile, str):
            with open(configfile, "r") as file:
                configfile = yaml.safe_load(file)

        workflow = sm.Workflow(
            snakefile,
            overwrite_configfiles=[],
            rerun_triggers=[],
            overwrite_config=configfile,
        )
        workflow.include(snakefile)
        workflow.global_resources = {}
        try:
            rule = workflow.get_rule(rulename)
        except Exception as exception:
            print(
                exception,
                f"The {rulename} might be a conditional rule in the Snakefile.\n"
                f"Did you enable {rulename} in the config?",
            )
            raise
        dag = sm.dag.DAG(workflow, rules=[rule])
        wc = Dict(wildcards)
        job = sm.jobs.Job(rule, dag, wc)

        def make_accessable(*ios):
            for io in ios:
                for i in range(len(io)):
                    io[i] = os.path.abspath(io[i])

        make_accessable(job.input, job.output, job.log)
        snakemake = Snakemake(
            job.input,
            job.output,
            job.params,
            job.wildcards,
            job.threads,
            job.resources,
            job.log,
            job.dag.workflow.config,
            job.rule.name,
            None,
        )
        snakemake.benchmark = job.benchmark

        # create log and output dir if not existent
        for path in list(snakemake.log) + list(snakemake.output):
            pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)

    finally:
        if user_in_script_dir:
            os.chdir(script_dir)
    return snakemake


def get_component(network, component):
    """
    Returns PyPSA model component

    Parameters
    ----------
    network: PyPSA model object
        pypsa model
    component: str
        PyPSA component to be returned

    Returns
    -------
    network[component]: pandas dataframe
        component of the PyPSA model
    """

    if component == "generators":
        return network.generators
    elif component == "links":
        return network.links
    elif component == "stores":
        return network.stores
    elif component == "storage_units":
        return network.storage_units
    elif component == "generators_t":
        return network.generators_t.p
    elif component == "links_t":
        return network.links_t
    elif component == "stores_t":
        return network.stores_t.p
    elif component == "storage_units_t":
        return network.storage_units_t.p


def get_component_list(energy_carrier):
    """
    Relevant PyPSA model components for each energy carrier

    Parameters
    ----------
    energy_carrier: str

    Returns
    -------
    list of components relevant to that particular energy carrier
    """

    energy_carrier_component_dict = {
        "electricity": ["generators", "links", "storage_units", "stores"],
        "heat": ["links", "stores"],
        "H2": ["links", "stores"],
        "water": ["links", "stores"],
        "oil": ["links", "stores"],
        "biomass": ["links", "stores"],
        "gas": ["links", "stores"],
        "battery": ["links", "stores"],
        "cooling": ["links"],
    }
    return energy_carrier_component_dict[energy_carrier]


def get_energy_carriers_key(energy_carrier):
    """
    List of PyPSA carriers for each energy carrier

    Parameters
    ----------
    energy_carrier: str

    Returns
    -------
    key to filter relevant PyPSA carriers
    """
    if energy_carrier == "electricity":
        return "AC|low voltage"
    else:
        return energy_carrier


def drop_carriers(df, key):
    """
    To drop irrelevant technologies
    """
    drop_carriers_list = [
        "electricity distribution grid",
        "H2 pipeline",
        "H2 pipeline repurposed",
        "DC",
        # "B2B",
        # "urban central gas CHP",
    ]
    for car in drop_carriers_list:
        if car in df.index.tolist() and key == "rows":
            df = df.drop(index=car)
        elif car in df.columns.tolist() and key == "columns":
            df = df.drop(car, axis=1)

    return df
