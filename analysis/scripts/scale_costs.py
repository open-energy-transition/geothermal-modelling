"""
scripts to scale cost assumptions for renewable tech

Relevant Settings
-----------------

.. code:: yaml


.. seealso::
    Documentation of the configuration file ``config.s6_2050_high_cost_ren.yaml``


Inputs
------

- PyPSA Model



Outputs
-------

- PyPSA Model



Description
-----------
"""

import pathlib
import pandas as pd
import datetime as dt
import pypsa

def parse_inputs(default_path):
    """
    Load all input data required for preprocessing and computing utility level demand data

    Parameters
    ----------
    default_path: str
        current total directory path

    Returns
    -------
    pypsa_network: PyPSA Model prenetwork file
    """
    pypsa_network_path = pathlib.Path(
        default_path, snakemake.input.pypsa_network_path[0]
    )

    pypsa_network = pypsa.Network(pypsa_network_path)
    return pypsa_network

def scale_renewable_costs(pypsa_network, scaling_factor, technology_list):
    """
    Renewable cost scaling

    Parameters
    ----------
    pypsa_network: pypsa model
    scaling_factor: factor by which to scale cost assumptions for renewable
    technology_list: tech list for which costs need to be scaled

    Returns
    -------
    pypsa_network: PyPSA Model modified prenetwork file
    """
    for tech in technology_list:
        pypsa_network.generators.loc[pypsa_network.generators.carrier == tech,"marginal_cost"] *= scaling_factor
        pypsa_network.generators.loc[pypsa_network.generators.carrier == tech,"capital_cost"] *= scaling_factor

    return pypsa_network

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("scale_cost_assumptions")

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "scale_cost_assumptions")

    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(
        log_path, f"output_scale_cost_assumptions_{today_date[:10]}.txt"
    )
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    scaling_factor = snakemake.params.cost_scaling_factor

    pypsa_network = parse_inputs(default_path)

    technology_list = ['solar','csp','onwind','offwind-ac','offwind-dc']

    scaled_pypsa_network = scale_renewable_costs(pypsa_network, scaling_factor, technology_list)

    pypsa_output_network_path = pathlib.Path(
        default_path, snakemake.output.pypsa_network_modified_path[0]
    )
    scaled_pypsa_network.export_to_netcdf(pypsa_output_network_path)
