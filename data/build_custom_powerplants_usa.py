#!/usr/bin/env python
# coding: utf-8

# In[199]:


import argparse
import datetime as dt
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import numpy as np
import pypsa


# In[200]:


base_path = pathlib.Path(__file__).parent.parent.parent

# Read reference data from plants_merged.csv (PyPSA-USA)
df_plants_usa = pd.read_csv('plants_merged.csv', low_memory=False)


# In[201]:


# Prepare plants_merged.csv for the custompowerplants.csv format
rename_map = {"plant_name": "Name", "fuel_type": "Fueltype", "nameplate_capacity_mw": "Capacity", "generator_id": "projectID", "operating_year": "DateIn", "planned_retirement_year": "DateOut", "latitude": "lat", "longitude": "lon"}

df_plants_usa = df_plants_usa.rename(columns = rename_map)
df_plants_usa


# In[202]:


# Generate custom_powerplants.csv
column_names = [
    "Name", "Fueltype", "Technology", "Set", "Country", "Capacity", 
    "Efficiency", "Duration", "Volume_Mm3", "DamHeight_m", 
    "StorageCapacity_MWh", "DateIn", "DateRetrofit", 
    "DateMothball", "DateOut", "lat", "lon", "EIC", "projectID", "bus"
]

df_custom_ppl_usa = pd.DataFrame(columns=column_names)


# In[203]:


# Copy the columns in plants_merged.csv into custom_powerplants.csv
df_plants_usa_aligned = df_plants_usa.reindex(columns=df_custom_ppl_usa.columns)
df_custom_ppl_usa = pd.concat([df_custom_ppl_usa, df_plants_usa_aligned], ignore_index=True).fillna('')


# In[204]:


# Set missing values and distinguish Hydro power plants according to the technology type
df_custom_ppl_usa['Country'] = 'US'
df_custom_ppl_usa['Fueltype'] = df_custom_ppl_usa['Fueltype'].replace('wind', 'onwind')
df_custom_ppl_usa['Set'] = df_custom_ppl_usa['Fueltype'].replace('wind', 'onwind')
df_custom_ppl_usa['Technology'] = df_custom_ppl_usa['Fueltype']
df_custom_ppl_usa['Technology'] = df_custom_ppl_usa['Technology'].replace('hydro', 'Reservoir')

phs = ['Horse Mesa',
'Mormon Flat',
'Castaic',
'Thermalito',
'Cabin Creek',
'Flatiron',
'Rocky River (CT)',
'Northfield Mountain',
'Ludington',
'Taum Sauk',
'Blenheim Gilboa',
'Lewiston Niagara',
'Hiwassee Dam',
'Salina',
'Muddy Run',
'Jocassee',
'Smith Mountain',
'Wallace Dam',
'Helms Pumped Storage',
'Fairfield Pumped Storage',
'Carters',
'Richard B Russell',
'Clarence Cannon',
'Harry Truman',
'Raccoon Mountain',
'Bath County',
'Rocky Mountain Hydroelectric Plant',
'Mount Elbert',
'Yards Creek Energy',
'Bad Creek',
'Bear Swamp',
'Seneca Generation LLC',
'Lake Hodges Hydroelectric Facility',
'Youngs Creek Hydroelectric Project']

# List of PHS plants from https://clui.org/projects/offstream/pumped-storage-facilities-usa
df_custom_ppl_usa.loc[df_custom_ppl_usa['Name'].isin(phs), 'Technology'] = 'Pumped Storage'

df_custom_ppl_usa['Set'] = np.where(df_custom_ppl_usa['Technology'] == 'Pumped Storage', 'Store', 'PP')
df_custom_ppl_usa


# In[205]:


# Generate custom_powerplants.csv
df_custom_ppl_usa.to_csv('custom_powerplants.csv', index=False)


# In[ ]:




