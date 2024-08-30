# -*- coding: utf-8 -*-
import geopandas as gpd
import rasterio
from rasterio.features import geometry_mask
from shapely.geometry import mapping
import numpy as np
import pandas as pd
from shapely.geometry import Point
from scipy.spatial import distance_matrix
import os

# Set the script directory as the working directory
working_dir = os.getcwd()
os.chdir(working_dir)

#################################### Landuse of existing pv #######################################
# installed pv data
installed_pv = gpd.read_file('mpv_shp_file_with_former_landuse')
installed_pv["ISO3"] = installed_pv["ISO3"].replace('-99', "FRA") # replace "-99" with "FRA"
area_sum_by_lc_vis = installed_pv.groupby('landuse_ty')['area'].sum()
area_sum_by_country = installed_pv.groupby('ISO3')['area'].sum()
area_sum_by_country_lc = installed_pv.groupby(['ISO3', 'landuse_ty'])['area'].sum( )

# pv area by country
total_area_by_ISO3 = installed_pv.groupby('ISO3')['area'].sum()

# calculate ratio by landuse type
area_ratio = installed_pv.groupby(['ISO3', 'landuse_ty'])['area'].sum().div(total_area_by_ISO3, level='ISO3')

pvlanduse_by_country = pd.merge(area_sum_by_country_lc, area_ratio, on=['ISO3', 'landuse_ty'])
pvlanduse_by_country.rename(columns={"area_x": "m2", "area_y": "ratio"}, inplace=True)
# pv land by country and landuse type
pvlanduse_by_country = pvlanduse_by_country[pvlanduse_by_country['ratio'] != 1]
pvlanduse_by_country = pd.merge(area_sum_by_country_lc, area_ratio, on=['ISO3', 'landuse_ty'])
pvlanduse_by_country.rename(columns={"area_x": "m2", "area_y": "ratio"}, inplace=True)
# remove countries only install pv panel on only one type of land
pvlanduse_by_country = pvlanduse_by_country[pvlanduse_by_country['ratio'] != 1]

######################################################### Mine area #################################
# mine area polygon file
gdf = gpd.read_file('data/landuse_occupation/MineArea_Centroid.shp')
# excel with energy output
df = pd.read_csv('data/landuse_occupation/Final_version.csv')
# calculate energy output
df['Electricity_Gen_Gwh'] = df['energy_v4'] * df['Area_3d_v1']*1000000/1e6
# world background map
world = gpd.read_file('data/landuse_occupation/continent.shp')
merged_df = pd.merge(gdf, df, on='Numbering', how='left')

# calculate suitable MPV area for each country
suitable_mine_area_by_country = merged_df.groupby('Country_code')['Area_3d_v1'].sum()
suitable_mine_area_by_country.sort_values(inplace=True,ascending=False)
suitable_mine_area_by_country = pd.DataFrame(suitable_mine_area_by_country)
suitable_mine_area_by_country.rename_axis(index={'Country_code':'ISO3'}, inplace=True)

#### Estimate pv land occupation for each country (countries without data use surrouding-3-countries's average value ###
# Get the lists of ISO3 values from both DataFrames
mine_iso3_list = suitable_mine_area_by_country.index.tolist()
pv_iso3_list = pvlanduse_by_country.index.get_level_values('ISO3').tolist()

# Convert the lists to sets for easy comparison
mine_iso3 = set(mine_iso3_list)
pv_iso3 = set(pv_iso3_list)

# Find common and unique ISO3 values
common_iso3 = mine_iso3.intersection(pv_iso3)
unique_mine_iso3 = mine_iso3.difference(pv_iso3)

# Create a list of tuples with ISO3 and the corresponding value (0 for common, 1 for unique to suitable_mine_area)
iso3_values = [(iso3, 0) for iso3 in common_iso3] + [(iso3, 1) for iso3 in unique_mine_iso3]

# Create a new DataFrame from the list of tuples
country_with_label = pd.DataFrame(iso3_values, columns=['ISO3', 'label']).set_index('ISO3')

# Load the shapefile
shapefile_path = 'global_shp_map_at_country_level'
country_code = gpd.read_file(shapefile_path)
country_code = country_code.rename(columns={'Alpha_3_co': 'ISO3'})

# Merge the GeoDataFrame with the result DataFrame on 'ISO3'
country_code = country_code.merge(country_with_label, on='ISO3')

# Filter the countries with a value of 1 and 0
countries_with_value_1 = country_code[country_code['label'] == 1]
countries_with_value_0 = country_code[country_code['label'] == 0]

# Create lists to store results
nearest_countries = []

# Calculate distances from each country with a value of 1 to all countries with a value of 0
for idx, row in countries_with_value_1.iterrows():
    country_iso3 = row['ISO3']
    country_geometry = row['geometry']

    # Calculate distances to all countries with a value of 0
    countries_with_value_0['distance'] = countries_with_value_0['geometry'].apply(
        lambda x: country_geometry.distance(x))

    # Sort by distance and get the nearest 3 countries
    nearest = countries_with_value_0.nsmallest(3, 'distance')

    # Store the results
    nearest_countries.append({
        'ISO3': country_iso3,
        'nearest_countries': nearest['ISO3'].tolist(),
        'distances': nearest['distance'].tolist()
    })

# Convert the results to a DataFrame
nearest_countries_df = pd.DataFrame(nearest_countries)

# Calculate the average "m2" value for each combination of "ISO3" and "landuse_ty"
avg_m2 = pvlanduse_by_country.groupby(['ISO3', 'landuse_ty'])['m2'].mean().reset_index()

mid_results = []

# Iterate over the nearest_countries_df to calculate the average "m2" values
for idx, row in nearest_countries_df.iterrows():
    iso3 = row['ISO3']
    nearest_countries = row['nearest_countries']

    # Filter avg_m2 for the nearest countries
    filtered_avg_m2 = avg_m2[avg_m2['ISO3'].isin(nearest_countries)]

    # Calculate the average "m2" for each combination of "ISO3" and "landuse_ty"
    for landuse_ty in filtered_avg_m2['landuse_ty'].unique():
        landuse_avg = filtered_avg_m2[filtered_avg_m2['landuse_ty'] == landuse_ty]['m2'].mean()
        mid_results.append({
            'ISO3': iso3,
            'landuse_ty': landuse_ty,
            'm2': landuse_avg
        })

# Convert the results to a DataFrame
country_1_label = pd.DataFrame(mid_results)

# sum area by country
total_area_by_ISO3_1label = country_1_label.groupby('ISO3')['m2'].sum()

# ratio by country and landuse type
area_ratio_1_label = country_1_label.groupby(['ISO3', 'landuse_ty'])['m2'].sum().div(total_area_by_ISO3_1label,
                                                                                     level='ISO3')

country_1_label_group = country_1_label.groupby(['ISO3', 'landuse_ty'])['m2'].sum().reset_index()
country_1_label_group = country_1_label_group.set_index(['ISO3', 'landuse_ty'])

pvlanduse_by_country_1label = pd.merge(country_1_label_group, area_ratio_1_label, on=['ISO3', 'landuse_ty'])
pvlanduse_by_country_1label.rename(columns={"m2_x": "m2", "m2_y": "ratio"}, inplace=True)

# merge countries with pv installation data and those without pv installation
all_country_landuse = pd.concat([pvlanduse_by_country, pvlanduse_by_country_1label])

################################ Estimate MPV ability to avoid land occupation  ###################
# load country ISO3 data
country_table = pd.read_excel('data/landuse_occupation/Country_table.xlsx',usecols=['ISO_A3'])
country_table = country_table.drop_duplicates(subset=['ISO_A3'])
country_table.rename(columns={'ISO_A3':'ISO3'}, inplace=True)

# combine country code and suitable mine area together
suitable_mine_area_by_country1 = pd.merge(suitable_mine_area_by_country, country_table, how='left', on='ISO3')
suitable_mine_area_by_country1.set_index(['ISO3'], inplace=True)
suitable_mine_area_by_country1

# combine landuse occupation data by existing pv and mine area data together
pvlanduse_iso3_values = all_country_landuse.index.get_level_values('ISO3')
pvland_country_mine = all_country_landuse.join(suitable_mine_area_by_country1, how='inner')
pvland_country_mine.drop(columns='m2', inplace=True)
pvland_country_mine_reset = pvland_country_mine.reset_index()
pvland_country_mine_restructured = pvland_country_mine_reset.set_index(['ISO3', 'landuse_ty'])
landuse_avoid = pvland_country_mine_restructured
landuse_avoid['area(km2)'] = landuse_avoid['ratio']*landuse_avoid['Area_3d_v1']
landuse_avoid_sorted = landuse_avoid.sort_values(by='Area_3d_v1', ascending=False)
landuse_avoid_sorted.drop(columns='Area_3d_v1', inplace=True)

######################### calculate biomass loss due to land occupation by pv ####################
# Projected coordinate system to geographic coordinate system
country_code = country_code.to_crs(epsg=4326)


# Calculate the centroid of each polygon
country_code['centroid'] = country_code.geometry.centroid

# Extract the latitude of each centroid
country_code['latitude'] = country_code['centroid'].apply(lambda x: x.y)

# Classify based on the latitude
def classify_latitude(lat):
    if abs(lat) <= 23.5:
        return 'Tropical'
    elif 23.5 < abs(lat) <= 40:
        return 'Subtropical'
    elif 40 < abs(lat) <= 60:
        return 'Temperate'
    elif abs(lat) > 60:
        return 'Boreal'
    else:
        return 'Unknown'

country_code['climate_class'] = country_code['latitude'].apply(classify_latitude)


# climate type of each country
country_climate = pd.read_csv('data/landuse_occupation/biomass_calculation.csv')
landuse_avoid_sorted_reindex = landuse_avoid_sorted.reset_index(level='landuse_ty')
country_climate_area = pd.merge(landuse_avoid_sorted_reindex, country_climate, on='ISO3', how='left')

# select green land
values_to_keep = ['grasslands', 'shrub/herbaceous/sparse', 'treecover']
country_climate_area_greenland = country_climate_area[country_climate_area['landuse_ty'].isin(values_to_keep)]

# get the biomass
conditions = [
    (country_climate_area_greenland['landuse_ty'] == 'grasslands') & (country_climate_area_greenland['climate_class'] == 'Boreal'),
    (country_climate_area_greenland['landuse_ty'] == 'grasslands') & (country_climate_area_greenland['climate_class'] == 'Temperate'),
    (country_climate_area_greenland['landuse_ty'] == 'grasslands') & (country_climate_area_greenland['climate_class'] == 'Subtropical'),
    (country_climate_area_greenland['landuse_ty'] == 'grasslands') & (country_climate_area_greenland['climate_class'] == 'Tropical'),
    (country_climate_area_greenland['landuse_ty'] == 'shrub/herbaceous/sparse') & (country_climate_area_greenland['region'] == 'Asia'),
    (country_climate_area_greenland['landuse_ty'] == 'shrub/herbaceous/sparse') & (country_climate_area_greenland['region'] == 'Africa'),
    (country_climate_area_greenland['landuse_ty'] == 'shrub/herbaceous/sparse') & (country_climate_area_greenland['region'] == 'Americas'),
    (country_climate_area_greenland['landuse_ty'] == 'shrub/herbaceous/sparse') & (country_climate_area_greenland['region'] == 'Europe'),
    (country_climate_area_greenland['landuse_ty'] == 'shrub/herbaceous/sparse') & (country_climate_area_greenland['region'] == 'Oceania'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Asia') & (country_climate_area_greenland['climate_class'] == 'Tropical'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Americas') & (country_climate_area_greenland['climate_class'] == 'Tropical'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Africa') & (country_climate_area_greenland['climate_class'] == 'Tropical'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Oceania') & (country_climate_area_greenland['climate_class'] == 'Tropical'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Europe') & (country_climate_area_greenland['climate_class'] == 'Tropical'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Asia') & (country_climate_area_greenland['climate_class'] == 'Subtropical'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Americas') & (country_climate_area_greenland['climate_class'] == 'Subtropical'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Africa') & (country_climate_area_greenland['climate_class'] == 'Subtropical'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Oceania') & (country_climate_area_greenland['climate_class'] == 'Subtropical'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Europe') & (country_climate_area_greenland['climate_class'] == 'Subtropical'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Asia') & (country_climate_area_greenland['climate_class'] == 'Temperate'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Americas') & (country_climate_area_greenland['climate_class'] == 'Temperate'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Africa') & (country_climate_area_greenland['climate_class'] == 'Temperate'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Oceania') & (country_climate_area_greenland['climate_class'] == 'Temperate'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Europe') & (country_climate_area_greenland['climate_class'] == 'Temperate'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Asia') & (country_climate_area_greenland['climate_class'] == 'Boreal'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Americas') & (country_climate_area_greenland['climate_class'] == 'Boreal'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Africa') & (country_climate_area_greenland['climate_class'] == 'Boreal'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Oceania') & (country_climate_area_greenland['climate_class'] == 'Boreal'),
    (country_climate_area_greenland['landuse_ty'] == 'treecover') & (country_climate_area_greenland['region'] == 'Europe') & (country_climate_area_greenland['climate_class'] == 'Boreal')
]

choices = [
    1.7, 1.7, 1.6, 2.3,
    70, 70, 80, 60, 60,
    130, 210, 120, 130, 'NoData',
    130, 210, 140, 130, 130,
    20, 60, 'NoData', 20, 20,
    10, 10, 'NoData', 'NoData', 10,
]

country_climate_area_greenland['abg_biomass'] = np.select(conditions, choices, default=np.nan)



# get the carbon fraction
country_climate_area_greenland['carbon_fraction'] = np.where(
    country_climate_area_greenland['landuse_ty'] == 'grasslands',
    0.5,
    np.where(
        country_climate_area_greenland['landuse_ty'].isin(['shrub/herbaceous/sparse', 'treecover']),
        0.47,
        np.nan  # Default value if none of the conditions match
    )
)

# calculate CO2 release, the  unit of abg_biomass is (toones d.m./ha);
country_climate_area_greenland['area(km2)'] = pd.to_numeric(country_climate_area_greenland['area(km2)'], errors='coerce')
country_climate_area_greenland['abg_biomass'] = pd.to_numeric(country_climate_area_greenland['abg_biomass'], errors='coerce')
country_climate_area_greenland['carbon_fraction'] = pd.to_numeric(country_climate_area_greenland['carbon_fraction'], errors='coerce')
country_climate_area_greenland['CO2'] = (
    country_climate_area_greenland['area(km2)'] *
    country_climate_area_greenland['abg_biomass'] *
    100 *
    country_climate_area_greenland['carbon_fraction'] *
    44 / 12
)

# CO2 captured by each type of land
CO2_landuse = country_climate_area_greenland.groupby('landuse_ty')['CO2'].sum().reset_index()

# country level analysis
Cropland_country = country_climate_area[country_climate_area['landuse_ty'].isin(['cropland'])].sort_values(by='area(km2)', ascending=False)


