"""
Plotting scripts to visualize model results

Relevant Settings
-----------------

.. code:: yaml


.. seealso::
    Documentation of the configuration file ``config.usa_baseline.yaml``


Inputs
------

- ``workflow/pypsa-earth/resources/{run_name}/demand_profiles.csv``



Outputs
-------

- ``analysis/plots/summary_plots``
- ``analysis/outputs/summary_outputs``


Description
-----------
"""

import pathlib
import pandas as pd
import datetime as dt


def parse_inputs(default_path):
    """
    Load all input data required for preprocessing and computing utility level demand data

    Parameters
    ----------
    default_path: str
        current total directory path

    Returns
    -------
    df_demand: pandas dataframe
        Electricity demand profiles
    energy_totals: pandas dataframe
        Total energy requirement for each carrier
    industrial_demand: pandas dataframe
        Various components of industrial demand
    """
    demand_profiles_path = pathlib.Path(
        default_path, snakemake.input.demand_profile_path
    )
    energy_totals_path = pathlib.Path(
        default_path, snakemake.input.energy_totals_path[0]
    )
    industrial_demand_path = pathlib.Path(
        default_path, snakemake.input.industrial_demand_path[0]
    )

    df_demand = pd.read_csv(demand_profiles_path, index_col="time")
    energy_totals = pd.read_csv(energy_totals_path, index_col="Unnamed: 0.1")
    industrial_demand = pd.read_csv(industrial_demand_path)

    return df_demand, energy_totals, industrial_demand


# Currently handles only a single country
# Replacing time varying electricity load totals
# The services electricity and electricity residential values are replaced to match demand projections based
# on NREL EFS study after deduction of constant loads
def modify_electricity_totals(df_demand, energy_totals, industry_demand, country):
    """
    To modify energy_totals_{demand}_{planning_horizons}.csv

    Parameters
    ----------
    df_demand: pandas dataframe
        Electricity demand profiles from EIA (scaled by NREL EFS projections for future years)
    energy_totals: pandas dataframe
        Annual energy demand for each country across various technologies and sectors
    industry_demand: pandas dataframe
        Industrial demands
    country: str
        Country code

    Returns
    -------
    energy_totals: pandas dataframe
        Modified annual energy demands
    """
    elec_cols = [
        x
        for x in energy_totals.columns
        if "electricity" in x
        and x not in ["electricity residential", "services electricity"]
    ]

    elec_services = energy_totals.loc[country, "services electricity"]
    elec_residential = energy_totals.loc[country, "electricity residential"]
    service_elec_ratio = elec_services / (elec_services + elec_residential)
    elec_residential_ratio = elec_residential / (elec_services + elec_residential)

    # industry_electricity_demand = industry_demand["electricity"].sum() / 1e6  # in TWh
    industry_electricity_demand = 0  # already included in the demand from EIA

    # demand profiles have one hour granularity
    total_electricity_demand = df_demand.sum().sum() / 1e6  # in TWh
    replace_demand = (
        total_electricity_demand
        - energy_totals.loc[country, elec_cols].sum()
        - industry_electricity_demand
    )

    energy_totals.loc[country, "electricity residential"] = (
        replace_demand * elec_residential_ratio
    )
    energy_totals.loc[country, "services electricity"] = (
        replace_demand * service_elec_ratio
    )

    return energy_totals


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("modify_energy_totals")

    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "modify_energy_totals")

    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(
        log_path, f"output_modify_energy_totals_{today_date[:10]}.txt"
    )
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    output_path = pathlib.Path(default_path, snakemake.output.energy_totals_path[0])
    pathlib.Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Only one country - needs to be changed if multiple countries are being accounted for
    country = snakemake.params.country[0]
    df_demand, energy_totals, industrial_demand = parse_inputs(default_path)

    modified_energy_totals = modify_electricity_totals(
        df_demand, energy_totals, industrial_demand, country
    )
    modified_energy_totals.to_csv(output_path)
