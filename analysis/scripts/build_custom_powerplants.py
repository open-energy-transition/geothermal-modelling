import pathlib
import datetime as dt
import pandas as pd
from _helpers_usa import eia_to_pypsa_terminology
import numpy as np


def parse_inputs(base_path, log_file):
    """
    The function parses the necessary inputs for the analysis
    """
    log_file.write("        \n")
    log_file.write("        \n")
    log_file.write("Parse inputs for build powerplants \n")
    eia_generators_path = snakemake.input.eia_generators_data_path
    eia_plants_path = snakemake.input.eia_plants_data_path
    ror_custom_path = snakemake.input.ror_custom_powerplants_path

    eia_generators_df = pd.read_excel(
        eia_generators_path, sheet_name="Operable", skiprows=1
    )
    eia_plants_df = pd.read_excel(eia_plants_path, skiprows=1)
    ror_custom_ppl_df = pd.read_csv(ror_custom_path)

    return eia_generators_df, eia_plants_df, ror_custom_ppl_df


def build_custom_pp(eia_generators_df, eia_plants_df, log_file):
    log_file.write("        \n")
    log_file.write("        \n")
    log_file.write("Build custom powerplants file from EIA 860 database \n")
    eia_to_pypsa_dict = eia_to_pypsa_terminology()
    df_eia_reqd = eia_generators_df.loc[
        eia_generators_df.Technology.isin(eia_to_pypsa_dict.keys())
    ]
    df_eia_reqd = df_eia_reqd[
        [
            "Plant Name",
            "Technology",
            "Nameplate Capacity (MW)",
            "Operating Year",
            "Planned Retirement Year",
        ]
    ]
    df_eia_reqd["Technology"] = df_eia_reqd["Technology"].map(eia_to_pypsa_dict)
    # df_eia_reqd = df_eia_reqd.rename(columns={'Plant Name':'Name','Technology':'Fueltype','Op'})

    df_ppl = pd.DataFrame(
        columns=[
            "Name",
            "Fueltype",
            "Technology",
            "Set",
            "Country",
            "Capacity",
            "Efficiency",
            "Duration",
            "Volume_Mm3",
            "DamHeight_m",
            "StorageCapacity_MWh",
            "DateIn",
            "DateRetrofit",
            "DateMothball",
            "DateOut",
            "lat",
            "lon",
            "EIC",
            "projectID",
            "bus",
        ]
    )
    df_ppl["Name"] = df_eia_reqd["Plant Name"].tolist()
    df_ppl["Fueltype"] = df_eia_reqd["Technology"].tolist()
    df_ppl["Capacity"] = df_eia_reqd["Nameplate Capacity (MW)"].tolist()
    df_ppl["DateIn"] = df_eia_reqd["Operating Year"].replace(" ", np.nan).tolist()
    df_ppl["DateOut"] = (
        df_eia_reqd["Planned Retirement Year"].replace(" ", np.nan).tolist()
    )
    df_ppl["Country"] = "US"
    df_ppl["Set"] = "PP"

    df_ppl.loc[df_ppl.Fueltype == "battery", "Set"] = "Store"
    df_ppl.loc[df_ppl.Fueltype == "hydro", "Technology"] = "Reservoir"
    df_ppl.loc[df_ppl.Fueltype == "OCGT", "Technology"] = "OCGT"
    df_ppl.loc[df_ppl.Fueltype == "CCGT", "Technology"] = "CCGT"
    df_ppl.loc[df_ppl.Fueltype == "PHS", "Technology"] = "Pumped Storage"
    df_ppl.loc[df_ppl.Fueltype == "PHS", "Set"] = "Store"
    df_ppl.loc[df_ppl.Fueltype == "PHS", "Fueltype"] = "hydro"

    df_ppl.loc[:, ["lat", "lon"]] = df_ppl.swifter.apply(
        lambda x: eia_plants_df.query('`Plant Name` == @x["Name"]').iloc[0][
            ["Latitude", "Longitude"]
        ],
        axis=1,
    ).rename(columns={"Latitude": "lat", "Longitude": "lon"})
    df_ppl["cumcount"] = df_ppl.groupby("Name").cumcount().add(1)
    df_ppl["Name"] = df_ppl.apply(
        lambda x: x["Name"] + " " + str(x["cumcount"]), axis=1
    )

    return df_ppl


def identify_ror_plants(custom_ppl_eia_df, custom_ppl_ror_df, log_file):
    log_file.write("        \n")
    log_file.write("        \n")
    log_file.write("Identify and reassign Run-Of-River powerplants \n")
    ror_plant_names = custom_ppl_ror_df.Name.tolist()
    custom_ppl_eia_df.loc[
        (custom_ppl_eia_df["Name"].isin(ror_plant_names))
        & (custom_ppl_eia_df["Fueltype"] == "hydro"),
        "Technology",
    ] = "Run-Of-River"

    return custom_ppl_eia_df


if __name__ == "__main__":
    # set relevant paths
    default_path = pathlib.Path(__file__).parent.parent.parent
    log_path = pathlib.Path(default_path, "analysis", "logs", "custom_powerplants")
    plot_path = pathlib.Path(default_path, "analysis", "plots", "custom_powerplants")
    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(plot_path).mkdir(parents=True, exist_ok=True)
    today_date = str(dt.datetime.now())
    log_output_file_path = pathlib.Path(
        log_path, f"output_build_custom_powerplants_{today_date[:10]}.txt"
    )
    log_output_file_path.touch(exist_ok=True)
    log_output_file = open(log_output_file_path, "w")

    eia_generators_df, eia_plants_df, ror_custom_ppl = parse_inputs(
        default_path, log_output_file
    )
    df_custom_ppl_eia = build_custom_pp(
        eia_generators_df, eia_plants_df, log_output_file
    )
    df_custom_ppl = identify_ror_plants(
        df_custom_ppl_eia, ror_custom_ppl, log_output_file
    )

    df_custom_ppl.to_csv(snakemake.output.output_filepath)

    log_output_file.write("        \n")
    log_output_file.write("        \n")
    log_output_file.write("Complete build powerplants \n")

    log_output_file.close()
