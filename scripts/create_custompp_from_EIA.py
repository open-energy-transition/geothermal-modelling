import pandas as pd
import numpy as np

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
df_eia_gen = pd.read_excel(
    "../Analysis/data/EIA generators/eia8602021/3_1_Generator_Y2021.xlsx",
    sheet_name="Operable",
    skiprows=1,
)
df_eia_plant = pd.read_excel(
    "../Analysis/data/EIA generators/eia8602021/2___Plant_Y2021.xlsx", skiprows=1
)

pypsa_tech = {
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

df_eia_reqd = df_eia_gen.loc[df_eia_gen.Technology.isin(pypsa_tech.keys())]
df_eia_reqd = df_eia_reqd[
    [
        "Plant Name",
        "Technology",
        "Nameplate Capacity (MW)",
        "Operating Year",
        "Planned Retirement Year",
    ]
]
df_eia_reqd["Technology"] = df_eia_reqd["Technology"].map(pypsa_tech)
# df_eia_reqd = df_eia_reqd.rename(columns={'Plant Name':'Name','Technology':'Fueltype','Op'})

df_ppl["Name"] = df_eia_reqd["Plant Name"].tolist()
df_ppl["Fueltype"] = df_eia_reqd["Technology"].tolist()
df_ppl["Capacity"] = df_eia_reqd["Nameplate Capacity (MW)"].tolist()
df_ppl["DateIn"] = df_eia_reqd["Operating Year"].replace(" ", np.nan).tolist()
df_ppl["DateOut"] = df_eia_reqd["Planned Retirement Year"].replace(" ", np.nan).tolist()
df_ppl["Country"] = "US"
df_ppl["Set"] = "PP"

df_ppl.loc[df_ppl.Fueltype == "battery", "Set"] = "Store"
df_ppl.loc[df_ppl.Fueltype == "hydro", "Technology"] = "Reservoir"
df_ppl.loc[df_ppl.Fueltype == "OCGT", "Technology"] = "OCGT"
df_ppl.loc[df_ppl.Fueltype == "CCGT", "Technology"] = "CCGT"
df_ppl.loc[df_ppl.Fueltype == "PHS", "Technology"] = "Pumped Storage"
df_ppl.loc[df_ppl.Fueltype == "PHS", "Set"] = "Store"
df_ppl.loc[df_ppl.Fueltype == "PHS", "Fueltype"] = "hydro"

# df_ppl[['lat','lon']] = df_ppl.swifter.apply(lambda x: df_eia_plant.query('`Plant Name` == @x["Name"]')[['Latitude','Longitude']],axis=1)
df_ppl.loc[:, ["lat", "lon"]] = df_ppl.swifter.apply(
    lambda x: df_eia_plant.query('`Plant Name` == @x["Name"]').iloc[0][
        ["Latitude", "Longitude"]
    ],
    axis=1,
).rename(columns={"Latitude": "lat", "Longitude": "lon"})
df_ppl["cumcount"] = df_ppl.groupby("Name").cumcount().add(1)
df_ppl["Name"] = df_ppl.apply(lambda x: x["Name"] + " " + str(x["cumcount"]), axis=1)
df_ppl.to_csv("../Data/custom_powerplants_eia.csv")
