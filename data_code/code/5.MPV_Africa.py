# Import packages
import pandas as pd
import os
from rasterio.features import geometry_mask
import geopandas as gpd
from shapely.geometry import box
from shapely.ops import unary_union
from rasterstats import zonal_stats
import rasterio
from rasterio.mask import mask

# Set the script directory as the working directory
working_dir = os.getcwd()
os.chdir(working_dir)

#  Create raster map outside maingride buffer
with rasterio.open('data/MPV_africa/africa_population2022_1.tif') as src:
    raster_profile = src.profile
    raster_data = src.read(1)  # Assuming it's a single band raster

# Loop through shapefiles with different buffer distances
for buffer_distance in [500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]:
    ## Open the shapefile
    shp_path = 'data/MPV_africa/africa_maingrid_{buffer_distance}mbuffer.shp'
    gdf = gpd.read_file(shp_path)

    ## Merge all geometries into a single polygon
    merged_geometry = unary_union(gdf.geometry)

    ## Create a mask from the merged geometry
    mask = geometry_mask([merged_geometry], out_shape=raster_data.shape, transform=raster_profile['transform'], invert=True)

    ## Apply the mask to the raster data
    raster_data[mask] = 0

    ## Write the modified raster to a new file
    output_raster_path = 'data/MPV_africa/africa_access_elec_{buffer_distance}m.tif'
    with rasterio.open(output_raster_path, 'w', **raster_profile) as dst:
        dst.write(raster_data, 1)
# Calculate access to electricity rate under different scenarioes
world = gpd.read_file('data/MPV_africa/WB_countries_Admin0_10m.shp')
africa_countries = world[world["CONTINENT"] == "Africa"]

raster_500= 'data/MPV_africa/africa_access_elec_500m.tif'
raster_1000= 'data/MPV_africa/africa_access_elec_1000m.tif'
raster_2000= 'data/MPV_africa/africa_access_elec_2000m.tif'
raster_3000= 'data/MPV_africa/africa_access_elec_3000m.tif'
raster_4000= 'data/MPV_africa/africa_access_elec_4000m.tif'
raster_5000= 'data/MPV_africa/africa_access_elec_5000m.tif'
raster_6000= 'data/MPV_africa/africa_access_elec_6000m.tif'
raster_7000= 'data/MPV_africa/africa_access_elec_7000m.tif'
raster_8000= 'data/MPV_africa/africa_access_elec_8000m.tif'
raster_All = 'data/MPV_africa/africa_population2022_1.tif'

# Using rasterstats to calculate the sum of raster values for each polygon
stats500 = zonal_stats(africa_countries, raster_500, stats="sum", all_touched=True)
# Adding the statistics back to the original GeoDataFrame
africa_countries['population_500'] = [stat['sum'] for stat in stats500]
# Using rasterstats to calculate the sum of raster values for each polygon
stats1000 = zonal_stats(africa_countries, raster_1000, stats="sum", all_touched=True)
# Adding the statistics back to the original GeoDataFrame
africa_countries['population_1000'] = [stat['sum'] for stat in stats1000]
# Using rasterstats to calculate the sum of raster values for each polygon
stats2000 = zonal_stats(africa_countries, raster_2000, stats="sum", all_touched=True)
# Adding the statistics back to the original GeoDataFrame
africa_countries['population_2000'] = [stat['sum'] for stat in stats2000]
# Using rasterstats to calculate the sum of raster values for each polygon
stats3000 = zonal_stats(africa_countries, raster_3000, stats="sum", all_touched=True)
# Adding the statistics back to the original GeoDataFrame
africa_countries['population_3000'] = [stat['sum'] for stat in stats3000]
# Using rasterstats to calculate the sum of raster values for each polygon
stats4000 = zonal_stats(africa_countries, raster_4000, stats="sum", all_touched=True)
# Adding the statistics back to the original GeoDataFrame
africa_countries['population_4000'] = [stat['sum'] for stat in stats4000]
# Using rasterstats to calculate the sum of raster values for each polygon
stats5000 = zonal_stats(africa_countries, raster_5000, stats="sum", all_touched=True)
# Adding the statistics back to the original GeoDataFrame
africa_countries['population_5000'] = [stat['sum'] for stat in stats5000]
# Using rasterstats to calculate the sum of raster values for each polygon
stats6000 = zonal_stats(africa_countries, raster_6000, stats="sum", all_touched=True)
# Adding the statistics back to the original GeoDataFrame
africa_countries['population_6000'] = [stat['sum'] for stat in stats6000]
# Using rasterstats to calculate the sum of raster values for each polygon
stats7000 = zonal_stats(africa_countries, raster_7000, stats="sum", all_touched=True)
# Adding the statistics back to the original GeoDataFrame
africa_countries['population_7000'] = [stat['sum'] for stat in stats7000]
# Using rasterstats to calculate the sum of raster values for each polygon
stats8000 = zonal_stats(africa_countries, raster_8000, stats="sum", all_touched=True)
# Adding the statistics back to the original GeoDataFrame
africa_countries['population_8000'] = [stat['sum'] for stat in stats8000]
# Using rasterstats to calculate the sum of raster values for each polygon
statsAll = zonal_stats(africa_countries, raster_All, stats="sum", all_touched=True)
# Adding the statistics back to the original GeoDataFrame
africa_countries['population_all'] = [stat['sum'] for stat in statsAll]
column_names = africa_countries.columns[-10:]
column_names =  column_names.insert(0, 'ISO_A3')
df_last_10 = africa_countries[column_names]
df_last_10['access_ratio_500']= (df_last_10['population_all']- df_last_10['population_500'])/df_last_10['population_all']
df_last_10['access_ratio_1000']= (df_last_10['population_all']- df_last_10['population_1000'])/df_last_10['population_all']
df_last_10['access_ratio_2000']= (df_last_10['population_all']- df_last_10['population_2000'])/df_last_10['population_all']
df_last_10['access_ratio_3000']= (df_last_10['population_all']- df_last_10['population_3000'])/df_last_10['population_all']
df_last_10['access_ratio_4000']= (df_last_10['population_all']- df_last_10['population_4000'])/df_last_10['population_all']
df_last_10['access_ratio_5000']= (df_last_10['population_all']- df_last_10['population_5000'])/df_last_10['population_all']
df_last_10['access_ratio_6000']= (df_last_10['population_all']- df_last_10['population_6000'])/df_last_10['population_all']
df_last_10['access_ratio_7000']= (df_last_10['population_all']- df_last_10['population_7000'])/df_last_10['population_all']
df_last_10['access_ratio_8000']= (df_last_10['population_all']- df_last_10['population_8000'])/df_last_10['population_all']
real = pd.read_csv('data/MPV_africa/africa_access_elc.csv')
elec_countries = pd.merge(real, df_last_10, on = 'ISO_A3', how='left')
elec_countries.to_csv('data/MPV_africa/africa_access_electricity.csv')

# Find the diffusion distance for each country
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

elec_countries['suit_distance'] = elec_countries.apply(find_nearest, axis=1)
elec_countries.to_csv('data/MPV_africa/africa_access_electricity_v2.csv')

# Clip population raster data for each country
# Load the shapefile
shapefile_path = 'data/MPV_africa/Africa_continent.shp'
shapefile = gpd.read_file(shapefile_path)

# List of country ISO3 codes
iso3_list = ['ZAF', 'GHA', 'NAM', 'BWA', 'MAR', 'ZWE', 'AGO', 'MLI', 'BFA', 'NER', 'MRT', 'SEN', 'EGY', 'SDN', 'TUN'
    , 'MDG', 'KEN', 'NGA', 'DZA', 'ERI', 'COG', 'UGA', 'SOM', 'RWA', 'BDI']

# Filter the shapefile to include only the polygons with ISO3 codes in the list
filtered_shapefile = shapefile[shapefile['ISO_A3'].isin(iso3_list)]

# Load the raster layer
raster_path = 'data/MPV_africa/africa_population2022.tif'
raster = rasterio.open(raster_path)

# Directory to save the output rasters
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

    # Copy the metadata of the original raster
    out_meta = raster.meta.copy()

    # Update the metadata to match the shape of the clipped raster
    out_meta.update({
        "driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform
    })

    # Save the clipped raster
    output_path = os.path.join(output_directory, country_name)
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(out_image)

# Close the raster file
raster.close()

# Determine population without access to electricity
iso3_list = [
    'ZAF', 'GHA', 'NAM', 'BWA', 'MAR', 'ZWE', 'AGO', 'MLI', 'BFA', 'NER', 'MRT', 'SEN', 'EGY', 'SDN', 'TUN',
    'MDG', 'KEN', 'NGA', 'DZA', 'ERI', 'COG', 'UGA', 'SOM', 'RWA', 'BDI'
]

# Load the dataframe with suit_distance values
suit_distance_df = pd.read_csv(
    'data/MPV_africa/Africa_electricity_v3.csv')

# Directory where the country-specific raster files are stored
country_raster_dir = 'data/MPV_africa'

# Directory where the buffered line shapefiles are stored
buffered_lines_dir = 'data/MPV_africa/'

# Directory to save the masked rasters
output_directory = 'data/MPV_africa'



# Function to apply buffer and mask operations
def process_country_raster(country_iso3, suit_distance):
    ## Load the corresponding country raster
    raster_path = os.path.join(country_raster_dir, f"{country_iso3}_2022_population.tif")
    with rasterio.open(raster_path) as raster:
        raster_data = raster.read(1)  # Read the first band
        raster_profile = raster.profile

        ### Load the appropriate buffered line shapefile
        buffered_line_path = os.path.join(buffered_lines_dir, f"africa_maingrid_{suit_distance}mbuffer.shp")
        buffered_lines = gpd.read_file(buffered_line_path)

        #### Clip the buffered lines by the raster bounds
        lines_clipped = gpd.clip(buffered_lines, box(*raster.bounds))

        if not lines_clipped.empty:
            ##### Create a merged geometry from buffered lines
            merged_geometry = lines_clipped.unary_union

            ##### Create a mask from the merged geometry
            mask = geometry_mask([merged_geometry], out_shape=raster_data.shape,
                                 transform=raster_profile['transform'], invert=True)

            ##### Apply the mask to the raster data
            raster_data[mask] = 0

            ##### Update the profile for the output raster
            raster_profile.update(dtype=rasterio.float32, count=1)

            ##### Save the masked raster
            output_path = os.path.join(output_directory, f"{country_iso3}_population2022_without_elec.tif")
            with rasterio.open(output_path, "w", **raster_profile) as dest:
                dest.write(raster_data, 1)


# Iterate over each country in the suit_distance_df
for _, row in suit_distance_df.iterrows():
    country_iso3 = row['ISO_A3']
    if country_iso3 in iso3_list:
        suit_distance = row['suit_distance']
        process_country_raster(country_iso3, suit_distance)

# Estimate number of people who get access from MPV
# Define paths
polygon_shapefile_path = 'data/MPV_africa/Africa_minePV_1000mbuffer_v1.shp'
raster_dir = 'data/MPV_africa'
output_shapefile_path = 'data/MPV_africa/Africa_minePV_1000mbuffer_pop_wo_elec.shp'

# Load the polygon shapefile
polygons = gpd.read_file(polygon_shapefile_path)

# Initialize a new column for population sum
polygons['Pop_wo_elec'] = 0.0

# Function to calculate population sum within a polygon
def calculate_population_sum(polygon, raster_path):
    with rasterio.open(raster_path) as src:
        ## Mask the raster with the polygon
        out_image, out_transform = mask(src, [polygon], crop=True)
        out_image = out_image[0]  # First band

        ## Calculate the sum of the population values
        population_sum = out_image[out_image > 0].sum()
        return population_sum

# Iterate over each polygon in the shapefile
for idx, row in polygons.iterrows():
    country_code = row['CountryCod']
    polygon = row['geometry']

    ## Construct the raster file path
    raster_path = os.path.join(raster_dir, f"{country_code}_population2022_without_elec.tif")

    if os.path.exists(raster_path):
        ### Calculate the population sum for the current polygon
        pop_sum = calculate_population_sum(polygon, raster_path)
        polygons.at[idx, 'Pop_wo_elec'] = pop_sum
    else:
        print(f"Raster file for {country_code} not found.")

# Save the updated GeoDataFrame to a new shapefile
polygons.to_file(output_shapefile_path)

# Number of people get access to electricity by MPV for each country
pop_elec_from_MPV = gpd.read_file('data/MPV_africa/Africa_minePV_1000mbuffer_pop_wo_elec.shp')
lines = gpd.read_file('data/MPV_africa/africa_maingrid.shp')
# Check whether mine polygons overlap with the main grid
polygons['ol_grid'] = polygons.apply(lambda poly: 1 if lines.intersects(poly.geometry).any() else 0, axis=1)
# Save the updated polygons to a new shapefile
polygons.to_file('data/MPV_africa/mine_pop_elec_v2.shp')
polygons[['Elec', 'CountryCod','Pop_wo_ele','population','ol_grid']].to_csv('data/MPV_africa/mine_pop_elec_v2.csv')
Access_to_elec_df = polygons[['Elec', 'CountryCod','Pop_wo_ele','ol_grid']]
Access_to_elec_df_1 = Access_to_elec_df[Access_to_elec_df['ol_grid'] == 0]
Access_to_elec = Access_to_elec_df_1.groupby('CountryCod')[['Pop_wo_ele','Elec']].sum()
Access_to_elec_1 =  Access_to_elec[Access_to_elec['Pop_wo_ele'] != 0]
Access_to_elec_2 =  Access_to_elec_1[Access_to_elec_1['Elec'] != 0]
# Delete this number (1991) from DZA, EGY, and MAR, as they have achieved 99% access to electricity rate
indexes_to_remove = ['DZA', 'EGY',  'MAR']
Access_to_elec_3 = Access_to_elec_2.drop(indexes_to_remove)
# Final results
Access_to_elec_3
