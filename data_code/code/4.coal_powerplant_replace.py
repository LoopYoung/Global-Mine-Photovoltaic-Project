# Import packages
import pandas as pd
import geopandas as gpd
import os

# Set the script directory as the working directory
working_dir = os.getcwd()
os.chdir(working_dir)

# Plant level power plants data
PowerPlant = pd.read_csv('data/coal_fired_powerplants/final_power_database_Feb2021_update.csv')
CoalPowerPlants = PowerPlant[(PowerPlant['status_adj'] == 'Operating') & (PowerPlant['fuel_class'] == 'COAL')]

# Change the dtype of data from str to float
CoalPowerPlants['mw'] = CoalPowerPlants['mw'].astype(float)
CoalPowerPlants['load_factor'] = CoalPowerPlants['load_factor'].astype(float)

# Calculate the electricity generated from coal fired power plants
CoalPowerPlants['electricity_gen(GWh)'] = CoalPowerPlants['mw'] * 24 * 365 * CoalPowerPlants['load_factor']/1000

# Calculate sum of electricity generated grouped by Country
Electricity_CoalPowerPlants = CoalPowerPlants.groupby('country_iso3')['electricity_gen(GWh)'].sum()

# Screen out coal mine
df = pd.read_excel('data/coal_fired_powerplants/MineClassify.xlsx')
# List of distance columns
distance_columns = ['NearDis_coal', 'NearDisFac', 'NearDisChi', 'NearDisCan', 'NearDisDep', 'NearDismrd']
# List of corresponding mine type columns
mine_type_columns = ['Cola Mine Accuracy', 'primary_co', '矿种', 'CommoCanad', 'Dep_commod', 'mrds_commo']

# Calculate the smallest distance and its source column
df['Smallest_Distance'] = df[distance_columns].min(axis=1)
df['Smallest_DisSource'] = df[distance_columns].idxmin(axis=1)

# Extract mine type values based on the source column
df['Mine_Type'] = df.apply(lambda row: row[mine_type_columns[distance_columns.index(row['Smallest_DisSource'])]], axis=1)

df['CoalMine'] = (
    (df['Mine_Type'].str.contains('Coal|煤|coal|COA|Approximate') & (df['Smallest_Distance'] <= 20000)) |
    (df['Mine_Type'].str.contains('Exact') & (df['Smallest_Distance'] <= 10000))
).astype(int)

# Correlate coal mine（FID is 1）with "numbering" in shp file
fid_to_keep = df[df['CoalMine'] == 1]['FID'].unique()
gdf = gpd.read_file('data/landuse_occupation/MineArea_Centroid.shp')
gdf['Numbering'] = gdf['Numbering'].astype(type(df['FID'].iloc[0]))
CoalMine_gdf = gdf[gdf['Numbering'].isin(fid_to_keep)]

# Import energy output data
Energy = pd.read_csv('data/coal_fired_powerplants/Final_version.csv')
Energy['Electricity_Gen_Gwh'] = Energy['energy_v4'] * Energy['Area_3d_v1']*1000000/1e6

# Merge two datasets
CoalMine_gdf_Energy = pd.merge(CoalMine_gdf, Energy, on='Numbering', how='left')

# Electrictity generation from coal fired power plants by country
Electricity_CoalMine = CoalMine_gdf_Energy.groupby('Country_code')['Electricity_Gen_Gwh'].sum()

# Compare CMPV and coal-fired powerplants
CoalMine_CoalFirePowerPlant = pd.concat([Electricity_CoalPowerPlants, Electricity_CoalMine], axis=1)
CoalMine_CoalFirePowerPlant.columns = ['Electricity_CoalPowerPlants(GWh)', 'Electricity_CoalMinePV(GWh)']
CoalMine_CoalFirePowerPlant['Ratio'] = (CoalMine_CoalFirePowerPlant['Electricity_CoalMinePV(GWh)'] / CoalMine_CoalFirePowerPlant['Electricity_CoalPowerPlants(GWh)'])*100
CoalMine_CoalFirePowerPlant['Ratio'] = CoalMine_CoalFirePowerPlant['Ratio'].clip(upper=100)
CoalMine_CoalFirePowerPlant_sorted = CoalMine_CoalFirePowerPlant.sort_values(by='Ratio', ascending=False)
CoalMine_CoalFirePowerPlant_sorted.reset_index(inplace=True)
CoalMine_CoalFirePowerPlant_sorted.rename(columns={'index': 'ISO_A3'}, inplace=True)

# Import world country polygons
world = gpd.read_file('data/coal_fired_powerplants/WB_countries_Admin0_10m.shp')

# Electricity from CMPV by continent
Electricity_CoalMine_Continent = CoalMine_gdf_Energy.groupby('Continent').agg({
    'Electricity_Gen_Gwh': 'sum',
    'geom_Area_x': 'sum'
}).reset_index()
Electricity_CoalMine_Continent.columns = ['Continent', 'Total_Electricity_Gen_Gwh', 'Total_geom_Area_x']

# Electricity from CMPV by country
Electricity_CoalMine_Country = CoalMine_gdf_Energy.groupby('Country_code').agg({
    'Electricity_Gen_Gwh': 'sum',
    'geom_Area_x': 'sum'
}).reset_index()
Electricity_CoalMine_Country.columns = ['Country_code', 'Total_Electricity_Gen_Gwh', 'Total_geom_Area_x']

# Combine CMPV data and coal-fired powerplants data together
merged_df = pd.merge(world, CoalMine_CoalFirePowerPlant_sorted, on='ISO_A3', how='left')
Country_mpv_coalfire = merged_df[merged_df['TYPE'].isin(['Sovereign country', 'Country'])]
Country_mpv_coalfire = Country_mpv_coalfire.dropna(subset=['Electricity_CoalPowerPlants(GWh)', 'Electricity_CoalMinePV(GWh)'], how='all')
Country_mpv_coalfire[['Electricity_CoalPowerPlants(GWh)', 'Electricity_CoalMinePV(GWh)']] = Country_mpv_coalfire[['Electricity_CoalPowerPlants(GWh)', 'Electricity_CoalMinePV(GWh)']].fillna(0)
Country_mpv_coalfire = Country_mpv_coalfire[['FORMAL_EN','ISO_A3', 'CONTINENT', 'Electricity_CoalPowerPlants(GWh)', 'Electricity_CoalMinePV(GWh)', 'Ratio']]
Country_mpv_coalfire = Country_mpv_coalfire[Country_mpv_coalfire['CONTINENT'] != 'Seven seas (open ocean)']
# NLD and NZL has duplicates, remove them
Country_mpv_coalfire = Country_mpv_coalfire.drop_duplicates(subset='ISO_A3')
