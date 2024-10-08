# Import packages
import pandas as pd
import os
from pvlib.solarposition import get_solarposition
from timezonefinder import TimezoneFinder
from datetime import datetime
import pytz
import math

# Set the script directory as the working directory
working_dir = os.getcwd()
os.chdir(working_dir)

# Calculate cost elements
    ## Install cost
install_cost = pd.read_excel("data/lcoe_estimation/LCOE_estimation.xlsx", sheet_name='Install')
    ## MPV data
base = pd.read_excel("data/lcoe_estimation/LCOE_estimation.xlsx", sheet_name='base')
    ##OECD country list
OECD_list = pd.read_excel("data/lcoe_estimation/LCOE_estimation.xlsx", sheet_name='OECD_list')
    ##WACC for different countries
WACC = pd.read_excel("data/lcoe_estimation/LCOE_estimation.xlsx", sheet_name='WACC')

    ## Install cost (for countries without install cost data, use continental average value)
continent_mean_install = install_cost.groupby('Continent')['Install_cost(2022 USD/kW)'].mean()
base['Install_cost'] = base.apply(lambda row: install_cost.loc[install_cost['Country_code'] == row['Country_code'], 'Install_cost(2022 USD/kW)'].values[0]
                                   if not pd.isnull(row['Country_code']) and not pd.isnull(install_cost.loc[install_cost['Country_code'] == row['Country_code'], 'Install_cost(2022 USD/kW)'].values).all()
                                   else continent_mean_install[row['Continent']] if not pd.isnull(row['Continent'])
                                   else pd.NA, axis=1)
    ## O&M cost; numbers come from International Renewable Cost Database
base['O&M_cost'] = base.apply(lambda row: 17.8 if row['Country_code'] in OECD_list['Country_code'].values else 9.2, axis=1)

    ## Discount rate
continent_mean_WACC = WACC.groupby('Continent')['real after-tax WACC'].mean()
base['Discount_rate'] = base.apply(lambda row: WACC.loc[WACC['Country_code'] == row['Country_code'], 'real after-tax WACC'].values[0]
                                   if not pd.isnull(row['Country_code']) and not pd.isnull(WACC.loc[WACC['Country_code'] == row['Country_code'], 'real after-tax WACC'].values).all()
                                   else continent_mean_WACC[row['Continent']] if not pd.isnull(row['Continent'])
                                   else pd.NA, axis=1)

# Function to calculate optimal tilt for PV panels
def calculate_tilt(latitude):
    tilt = (1.3793 + latitude * (1.2011 + latitude * (-0.014404 + 0.000080509 * latitude))) if latitude > 0 else (
        0 if latitude == 0 else (-(-0.41657 + latitude * (1.4216 + latitude * (0.024051 + latitude * 0.00021828)))))
    return tilt

# Function to calculate the solar altitude and azimuth at 3:00PM at winter solstice
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

# Calculate packing factor(panel area versus land occupation)
def calculate_packing_factor1(latitude, longitude):
    tilt = calculate_tilt(latitude)
    altitude, azimuth = calculate_solar_position_solstice(latitude, longitude)
    azimuth_v2 = azimuth if azimuth <= 90 else (180 - azimuth if azimuth <= 270 else azimuth - 360)
    pf = 1 / (math.cos(tilt / 180 * math.pi) + (
                math.sin(tilt / 180 * math.pi) / math.tan(altitude / 180 * math.pi)) * math.cos(
        azimuth_v2 / 180 * math.pi))
    return pf

base['packing factor'] = base.apply(lambda row: calculate_packing_factor1(row['latitude'], row['longitude']), axis=1)

# Calculate lcoe
    ## module parameters: Area(m2),Imp(A),Vmp(V)
Module_Area = 1.26
Imp = 5.5383
Vmp = 43.1204
RatedPower = Imp*Vmp

    ## Calculate LCOE (units:$/MWh)
n = 25 # assume lifetime span 25 years
base['LCOE_v2'] = 1000*((base['Install_cost'] * base['Area_3d_v1']* base['packing factor'] * 1000000 / Module_Area * RatedPower / 1000
          + sum(base['O&M_cost'] * base['Area_3d_v1'] * base['packing factor'] * 1000000 / Module_Area * RatedPower / 1000 / pow(1 + base['Discount_rate'], i) for i in range(1, n+1)))
          /sum(base['energy_v4'] * base['Area_3d_v1'] * 1000000  / pow(1 + base['Discount_rate'], i) for i in range(1, n + 1)))

# Remove values smaller than 0 (caused by too little energy output)
base.loc[(base['energy_v4'] < 0) | (base['packing factor'] < 0), 'LCOE_v2'] = ''
base.loc[base['packing factor'] < 0, 'energy_v4'] = ''

base.to_csv("data/lcoe_estimation/LCOE_estimation_result.csv")
