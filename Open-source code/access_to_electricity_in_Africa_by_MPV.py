# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import geopandas as gpd
import fiona
import rasterio
from rasterio.mask import mask
from shapely.geometry import shape
import rasterio
from rasterio.features import geometry_mask
from shapely.geometry import box
from shapely.ops import unary_union
from rasterstats import zonal_stats
import os

# Set the script directory as the working directory
working_dir = os.getcwd()
os.chdir(working_dir)

# africa part
world = gpd.read_file('world_map_at_country_level.shp')
africa = world[world["CONTINENT"] == "Africa"]
gdf = gpd.read_file('mine_area_map.shp')
df = pd.read_csv('the_energy_potential_of_each_mpv.csv')
df['Electricity_Gen_Gwh'] = df['energy_v4'] * df['Area_3d_v1'] * 1000000 / 1e6
df_africa = df[df["Continent"] == "Africa"]
merged_df = gdf.merge(df_africa, on='Numbering', how='inner')

################## installed PV in Africa ###################
installed_pv = gpd.read_file('global_soalr_PV_installation_map.shp')
installed_pv_africa = installed_pv[installed_pv['Continent_'] == 'Africa']
installed_pv_africa_area = installed_pv_africa.groupby('ISO3')['area'].sum() / 1000000  # convert to km2
mine_pv_africa_area = merged_df.groupby('Country_code')['Area_3d_v1'].sum()

# Convert Series to DataFrames
installed_pv_africa_area_df = installed_pv_africa_area.reset_index()
mine_pv_africa_area_df = mine_pv_africa_area.reset_index()

# Merge the DataFrames on 'Country_code' and 'ISO3'
minepv_installed_comparison = pd.merge(installed_pv_africa_area_df, mine_pv_africa_area_df, left_on='ISO3',
                                       right_on='Country_code')

# Drop redundant columns (one of 'ISO3' or 'Country_code')
minepv_installed_comparison.drop(columns=['Country_code'], inplace=True)
minepv_installed_comparison = minepv_installed_comparison.rename(
    columns={'area': 'Current installed areas', 'Area_3d_v1': 'Potential mine areas for installation'})

# Sort the DataFrame by the 'Potential mine areas for installation' column
minepv_installed_comparison_sorted = minepv_installed_comparison.sort_values(by='Potential mine areas for installation',
                                                                             ascending=True)
minepv_installed_comparison_sorted.to_csv('comparing_existing_pv_and_MPV.csv')

######################## determine main grid diffusion distance #######################
elec_countries = pd.read_csv('results_from_Step_"maingrid_difussion_distance_estimation".csv')


# Define the function to find the nearest value and its corresponding distance
def find_nearest(row):
    access_value = row['access_to_electricity']
    # Extract relevant columns and their corresponding distances
    distances = {
        500: row['access_ratio_500'],
        1000: row['access_ratio_1000'],
        2000: row['access_ratio_2000'],
        3000: row['access_ratio_3000'],
        4000: row['access_ratio_4000'],
        5000: row['access_ratio_5000'],
        6000: row['access_ratio_6000'],
        7000: row['access_ratio_7000'],
        8000: row['access_ratio_8000']
    }
    # Find the key (distance) with the nearest value
    nearest_distance = min(distances, key=lambda k: abs(distances[k] - access_value))
    return nearest_distance


# Apply the function to each row
elec_countries['suit_distance'] = elec_countries.apply(find_nearest, axis=1)

elec_countries.to_csv('output.csv')

################### clip africa population map by country ##########################
# Load the shapefile
shapefile_path = 'Africa_map.shp'
shapefile = gpd.read_file(shapefile_path)

# List of country ISO3 codes
iso3_list = ['ZAF', 'GHA', 'NAM', 'BWA', 'MAR', 'ZWE', 'AGO', 'MLI', 'BFA', 'NER', 'MRT', 'SEN', 'EGY', 'SDN', 'TUN'
    , 'MDG', 'KEN', 'NGA', 'DZA', 'ERI', 'COG', 'UGA', 'SOM', 'RWA', 'BDI']

# Filter the shapefile to include only the polygons with ISO3 codes in the list
filtered_shapefile = shapefile[shapefile['ISO_A3'].isin(iso3_list)]

raster_path = 'africa_population2022_map.tif'
raster = rasterio.open(raster_path)

output_directory = 'data/MPV_africa'
os.makedirs(output_directory, exist_ok=True)

# Iterate over each country in the filtered shapefile
for _, row in filtered_shapefile.iterrows():
    country_iso3 = row['ISO_A3']
    country_name = f"{country_iso3}_2022_population.tif"

    # Get the geometry of the country
    geom = [row['geometry']]

    # Clip the raster with the geometry
    out_image, out_transform = mask(raster, geom, crop=True)

    out_meta = raster.meta.copy()

    # Update the metadata to match the shape of the clipped raster
    out_meta.update({
        "driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform
    })

    output_path = os.path.join(output_directory, country_name)
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(out_image)

raster.close()

############## using diffusion distance for each country to calculate access to electricity rate ########
iso3_list = [
    'ZAF', 'GHA', 'NAM', 'BWA', 'MAR', 'ZWE', 'AGO', 'MLI', 'BFA', 'NER', 'MRT', 'SEN', 'EGY', 'SDN', 'TUN',
    'MDG', 'KEN', 'NGA', 'DZA', 'ERI', 'COG', 'UGA', 'SOM', 'RWA', 'BDI'
]

# Load the dataframe with suit_distance values
suit_distance_df = pd.read_csv('data/MPV_africa/Africa_electricity_v3.csv')

# Directory where the country-specific raster files are stored
country_raster_dir = 'data/MPV_africa'

# Directory where the buffered line shapefiles are stored
buffered_lines_dir = 'data/MPV_africa/'

# Directory to save the masked rasters
output_directory = 'data/MPV_africa'
os.makedirs(output_directory, exist_ok=True)


# Function to apply buffer and mask operations
def process_country_raster(country_iso3, suit_distance):
    # Load the corresponding country raster
    raster_path = os.path.join(country_raster_dir, f"{country_iso3}_2022_population.tif")
    with rasterio.open(raster_path) as raster:
        raster_data = raster.read(1)  # Read the first band
        raster_profile = raster.profile

        # Load the appropriate buffered line shapefile
        buffered_line_path = os.path.join(buffered_lines_dir, f"africa_maingrid_{suit_distance}mbuffer.shp")
        buffered_lines = gpd.read_file(buffered_line_path)

        # Clip the buffered lines by the raster bounds
        lines_clipped = gpd.clip(buffered_lines, box(*raster.bounds))

        if not lines_clipped.empty:
            # Create a merged geometry from buffered lines
            merged_geometry = lines_clipped.unary_union

            # Create a mask from the merged geometry
            mask = geometry_mask([merged_geometry], out_shape=raster_data.shape,
                                 transform=raster_profile['transform'], invert=True)

            # Apply the mask to the raster data
            raster_data[mask] = 0

            # Update the profile for the output raster
            raster_profile.update(dtype=rasterio.float32, count=1)

            # Save the masked raster
            output_path = os.path.join(output_directory, f"{country_iso3}_population2022_without_elec.tif")
            with rasterio.open(output_path, "w", **raster_profile) as dest:
                dest.write(raster_data, 1)


# Iterate over each country in the suit_distance_df
for _, row in suit_distance_df.iterrows():
    country_iso3 = row['ISO_A3']
    if country_iso3 in iso3_list:
        suit_distance = row['suit_distance']
        process_country_raster(country_iso3, suit_distance)

############## calculate number of people get access to electricity by MPV ####################
# Define paths
polygon_shapefile_path = 'mpv_with_1000m_buffer.shp'
raster_dir = 'directory_where_the rater_layer_are'
output_shapefile_path = 'mpv_in_africa_with_people_electrified_by_mpv.shp'

# Load the polygon shapefile
polygons = gpd.read_file(polygon_shapefile_path)

# Initialize a new column for population sum
polygons['Pop_wo_elec'] = 0.0


# Function to calculate population sum within a polygon
def calculate_population_sum(polygon, raster_path):
    with rasterio.open(raster_path) as src:
        # Mask the raster with the polygon
        out_image, out_transform = mask(src, [polygon], crop=True)
        out_image = out_image[0]  # First band

        # Calculate the sum of the population values
        population_sum = out_image[out_image > 0].sum()
        return population_sum


# Iterate over each polygon in the shapefile
for idx, row in polygons.iterrows():
    country_code = row['CountryCod']
    polygon = row['geometry']

    # Construct the raster file path
    raster_path = os.path.join(raster_dir, f"{country_code}_population2022_without_elec.tif")

    if os.path.exists(raster_path):
        # Calculate the population sum for the current polygon
        pop_sum = calculate_population_sum(polygon, raster_path)
        polygons.at[idx, 'Pop_wo_elec'] = pop_sum
    else:
        print(f"Raster file for {country_code} not found.")

# Save the updated GeoDataFrame to a new shapefile
polygons.to_file(output_shapefile_path)

print(f"Updated shapefile saved to {output_shapefile_path}")

# number of people get access to electricity by MPV for each country
pop_elec_from_MPV = gpd.read_file('population_get_access_to_electricity_from_MPV.shp')
selected_columns_pop_elec_from_MPV = pop_elec_from_MPV[['Electrici', 'CountryCod', 'Pop_wo_ele']]
elec_bycountry = selected_columns_pop_elec_from_MPV.groupby('CountryCod')[['Pop_wo_ele', 'Electrici']].sum()
