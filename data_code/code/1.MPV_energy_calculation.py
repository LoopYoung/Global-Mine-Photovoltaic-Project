#Import packages
import pvlib
from pvlib.modelchain import ModelChain
from pvlib.location import Location
from pvlib.pvsystem import PVSystem
from pvlib.temperature import TEMPERATURE_MODEL_PARAMETERS
import pandas as pd
from tqdm import tqdm
from pvlib.solarposition import get_solarposition
from timezonefinder import TimezoneFinder
from datetime import datetime
import pytz
import math
import os
# Function to calculate the optimal tilt
def calculate_tilt(latitude):
    tilt = (1.3793 + latitude * (1.2011 + latitude * (-0.014404 + 0.000080509 * latitude))) if latitude > 0 else (0 if latitude == 0 else (-(-0.41657 + latitude * (1.4216 + latitude * (0.024051 + latitude * 0.00021828)))))
    return tilt # based on equation1
#Function to calculate solar altitude and azimuth at 3:00PM at winter solstice
def calculate_solar_position_solstice(latitude, longitude):
    ## Use TimezoneFinder to determine the timezone
    tf = TimezoneFinder()
    timezone_str = tf.timezone_at(lat=latitude, lng=longitude)
    timezone = pytz.timezone(timezone_str)
    
    ## Determine the year for calculation
    year = datetime.now().year

    ## Determine winter solstice date based on hemisphere
    if latitude >= 0:  # Northern Hemisphere
        winter_solstice_date = f"{year}-12-21"
    else:  # Southern Hemisphere
        winter_solstice_date = f"{year}-06-21"

    ## Create a naive datetime object for 3:00 PM on the winter solstice
    naive_time = datetime.strptime(winter_solstice_date + ' 15:00:00', '%Y-%m-%d %H:%M:%S')

    ## Localize the naive datetime to the found timezone
    local_time = timezone.localize(naive_time)

    ## Calculate the solar position
    solar_position = get_solarposition(local_time, latitude, longitude)

    ## Extract solar altitude and azimuth
    altitude = solar_position['apparent_elevation'].values[0]
    azimuth = solar_position['azimuth'].values[0]

    return altitude, azimuth
# Function to calculate packing factor
def calculate_packing_factor(tilt, altitude, azimuth):
    pf = 1/(math.cos(tilt/180*math.pi)+(math.sin(tilt/180*math.pi)/math.tan(altitude/180*math.pi))*math.cos(azimuth/180*math.pi))
    return pf

# Set the script directory as the working directory
working_dir = os.getcwd()
os.chdir(working_dir)
MineLocation = pd.read_csv("data/mine_energy/MineAreaCentroid.csv")
for i in tqdm(range(len(MineLocation)), desc='Processing locations', position=0, mininterval=300):
    latitude = MineLocation.at[i,"latitude"]
    longitude = MineLocation.at[i,"longitude"]
    try:
        altitude = pvlib.location.lookup_altitude(latitude=latitude, longitude=longitude)
        location = Location(latitude =latitude, longitude =longitude, altitude =altitude)
        #PVmodule
        sandia_modules  = pvlib.pvsystem.retrieve_sam("SandiaMod")
        #Inverters
        cec_inverters = pvlib.pvsystem.retrieve_sam("CECInverter")
        module = sandia_modules["Panasonic_VBHN235SA06B__2013_"]
        inverter = cec_inverters["Trina_Energy_Storage_Solutions__Jiangsu__Co___Ltd__TB6000SHU__240V_"]
        #pv panel
        tilt = calculate_tilt(latitude)
        surface_azimuth = 180 if latitude >=0 else 0
        
        #module temperature
        temperature_parameters = TEMPERATURE_MODEL_PARAMETERS["sapm"]["open_rack_glass_glass"]
        
        system = PVSystem(surface_tilt = tilt, surface_azimuth = surface_azimuth,
                          module_parameters = module, inverter_parameters = inverter,
                          temperature_model_parameters = temperature_parameters,
                          modules_per_string = 8, strings_per_inverter = 1)
        modelchain = ModelChain(system, location)
        #Solar data
        tmy_data, months_selected, inputs, metadata = pvlib.iotools.get_pvgis_tmy(
                                       latitude=latitude, longitude=longitude,
                                       outputformat='json', usehorizon=True, 
                                       userhorizon=None, url='https://re.jrc.ec.europa.eu/api/v5_2/', 
                                       map_variables=True, timeout=30)
        selected_columns = ["temp_air", "ghi", "dni","dhi", "wind_speed"]
        tmy_data = tmy_data[selected_columns]
        tmy_data.index = pd.to_datetime(tmy_data.index,format="%Y%m%d:%H%M")
        
        modelchain.run_model(tmy_data)
        # Replace negative values with 0
        modelchain.results.ac = modelchain.results.ac.apply(lambda x: max(0, x))
        
        #Calculate solar altitude and azimuth at winter solstice 3:00pm and packing factor
        solar_altitude, solar_azimuth = calculate_solar_position_solstice(latitude, longitude)
        azimuth_adjust = (180-solar_azimuth) if latitude>=0 else (solar_azimuth-360) # pvlib与公式中的azimuth不一样
        pf = calculate_packing_factor(tilt, solar_altitude, azimuth_adjust)
        #Calculate energy generation(kwh)/m^2/year(original unit is XX)
        energy = modelchain.results.ac.sum() / ((module.loc["Area"]*8)/pf) /(10**3)
        
        # Store the calculated energy value in the "energy" column for the current row
        MineLocation.at[i, "energy"] = energy
    except Exception as e:
        MineLocation.at[i, "energy"] = pd.NA
    # Write the current row to the CSV file after each iteration
    MineLocation.iloc[i:i + 1].to_csv(os.path.join("data/mine_energy/mine_energy.csv"), mode='a', header=False, index=False)
