# Download ERA5 Met Data
# Bharat Sharma
# sharma.bha@ornl.gov

# https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview
# ======================

import sys

import cdsapi

varname = [
    "surface_pressure",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "surface_solar_radiation_downwards",
    "surface_thermal_radiation_downwards",
]
var = ["sp", "u10", "v10", "ssrd", "strd"]
save_path = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/ERA5_Duke_Met"

### Lat - Lon
# Source: https://www.ornl.gov/content/duke-forest-face-site-characteristics
# The Blackwood Division of the Duke Forest is near Chapel Hill, in Orange County, North Carolina (35° 58' 41.430"N, 79° 05' 39.087" W, 163 m asl)
# i.e. Latitude: 35.978175° N Longitude: -79.09419083333333° W


def download_era5land_var(varname, year, month):
    c = cdsapi.Client()

    c.retrieve(
        "reanalysis-era5-land",
        {
            "format": "netcdf",
            "variable": varname,
            "year": str(year),
            "month": str(month),
            "day": [
                "01",
                "02",
                "03",
                "04",
                "05",
                "06",
                "07",
                "08",
                "09",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "17",
                "18",
                "19",
                "20",
                "21",
                "22",
                "23",
                "24",
                "25",
                "26",
                "27",
                "28",
                "29",
                "30",
                "31",
            ],
            "time": [
                "00:00",
                "01:00",
                "02:00",
                "03:00",
                "04:00",
                "05:00",
                "06:00",
                "07:00",
                "08:00",
                "09:00",
                "10:00",
                "11:00",
                "12:00",
                "13:00",
                "14:00",
                "15:00",
                "16:00",
                "17:00",
                "18:00",
                "19:00",
                "20:00",
                "21:00",
                "22:00",
                "23:00",
            ],
            "area": [
                36.1,
                -79.4,
                35.7,
                -78.7,
            ],
        },
        f"{save_path}/{year}_{month}_daily_duke_met_ERA5.nc",
    )


year_list = [f"{i:d}" for i in range(1996, 2013)]
month_list = [f"{i:02d}" for i in range(1, 13)]

# Downloading
for year in year_list:
    for month in month_list:
        download_era5land_var(varname, year, month)
        download_era5land_var(varname, year, month)
