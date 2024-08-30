# -*- coding: utf-8 -*-

import rasterio
from rasterio.features import geometry_mask
import geopandas as gpd
from shapely.geometry import box
from shapely.ops import unary_union
from rasterstats import zonal_stats
import pandas as pd
import numpy as np
import fiona
import rasterio
from rasterio.mask import mask
from shapely.geometry import shape
import os

# Set the script directory as the working directory
working_dir = os.getcwd()
os.chdir(working_dir)

################  Create raster map outside maingride buffer ####################
# Open the raster file
with rasterio.open('the_population_raster_layer.tif') as src:
    raster_profile = src.profile
    raster_data = src.read(1)  # Assuming it's a single band raster

# Loop through shapefiles with different buffer distances
for buffer_distance in [500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]:
    # Open the shapefile
    shp_path = 'africa_maingrid_{buffer_distance}mbuffer.shp'
    gdf = gpd.read_file(shp_path)

    # Merge all geometries into a single polygon
    merged_geometry = unary_union(gdf.geometry)

    # Create a mask from the merged geometry
    mask = geometry_mask([merged_geometry], out_shape=raster_data.shape, transform=raster_profile['transform'], invert=True)

    # Apply the mask to the raster data
    raster_data[mask] = 0

    # Write the modified raster to a new file
    output_raster_path = 'output_{buffer_distance}m.tif'
    with rasterio.open(output_raster_path, 'w', **raster_profile) as dst:
        dst.write(raster_data, 1)
########################### calculate access to electricity rate under different scenarioes ###################
world = gpd.read_file('wordmap_at_country_level.shp')
africa_countries = world[world["CONTINENT"] == "Africa"]

raster_500= 'africa_access_elec_500m.tif'
raster_1000= 'africa_access_elec_1000m.tif'
raster_2000= 'africa_access_elec_2000m.tif'
raster_3000= 'africa_access_elec_3000m.tif'
raster_4000= 'africa_access_elec_4000m.tif'
raster_5000= 'africa_access_elec_5000m.tif'
raster_6000= 'africa_access_elec_6000m.tif'
raster_7000= 'africa_access_elec_7000m.tif'
raster_8000= 'africa_access_elec_8000m.tif'
raster_All = 'africa_population2022_1.tif'

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
real = pd.read_csv('final_data.csv')# same as africa_access_el_v2.csv