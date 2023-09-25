Duke Forest FACE, North Carolina, USA
Site ID: DUKE
Site coordinate: latitude = 35.9782, longitude = -79.0942
File Origin: this file was created at Oak Ridge National Laboratory for FACE modeling
synthesis
Date: Sep 25, 2023
Contact: Anthony Walker (walkerap@ornl.gov), Bharat Sharma (sharmabd@ornl.gov)
Missing value: -6999

I. Half-hourly file - DUKE_forcing_h.txt

'YEAR':'Year of measurement',
'DTIME':'Fractional day of year',
'DOY':'Day of year',
'HRMIN':'Hour:minute, marked at the middle of measurement interval with last two digits as minute',
'Rainf':'Total Precipitation over a time step of measurement',
'Rainf_f': 'gap-filling flag, 0 = measured, 1 = derived from other variables, 2 = filled by \
exising FACEMDS data, 3 = filled by data from nearby weather station, 4 = \
filled by using ERA5 data',
'Tair':'Mean air temperature over a time step of measurement',
'Tair_f':'gap-filling flag',
'RH':'Mean relative humidity over a time step of measurement',
'RH_f':'gap-filling flag',
'VPD':'Vapor pressure deficit, kPa',
'VPD_f':'gap-filling flag',
'PAR':'Incident or downward photosynthetically active radiation',   
'PAR_f':'gap-filling flag',    
'SM':'Soil Moisture integrates measurements from 0 to 30cm depth',
'SM_f':'gap-filling flag', 
'SWP':' Soil Water Potential',
'SWP_f':'gap-filling flag', 
'SVP':'Saturated Vapor Pressure',
'SVP_f':'gap-filling flag', 
'Rn':'Net Radiation',
'Rn_f':'gap-filling flag', 
'SLT':'Soil Temperature',   
'SLT_f':'gap-filling flag', 
'Wind':'Mean wind speed over a time step of measurement, m/s',
'Wind_f':'gap-filling flag',
'PSurf': 'Surface barometric pressure, Pa',
'PSurf_f':'gap-filling flag',
'aCO2': 'Daily mean ambient CO2 concentration in daytime (solar angle > 15), ppmv',
'eCO2': 'Daily mean elevated treatment CO2 concentration in daytime (solar angle >15), ppmv',
'Ndep': 'Total N deposition over a time step of measurement (30 minutes), g/m2/(30-minute)',
'SolarElevation': 'Solar elevation angle, degree'


II. Daily file – DUKE_forcing_d.txt


'YEAR':'Year of measurement',
'DOY':'Day of year',
'Rainf':'Daily Total Precipitation',
'Tair':'Daily Mean air temperature',
'RH':'Daily mean relative humidity, %',
'VPD':'Daily mean vapor pressure deficit, Pa',
'Wind':'Daily mean wind speed, m/s',
'SWdown':'Daily total incident or downward short-wave radiation,W/m2',
'PAR':'Daily total incident or downward photosynthetically active radiation,mol/m2/day',
'LWdown:'Daily total ncident or downward short-wave radiation,W/m2',
'Psurf':'Daily mean surface barometric pressure, Pa',
'aCO2':'Daily mean ambient CO2 concentration in daytime (solar angle > 15), ppmv',
'eCO2':'Daily mean elevated treatment CO2 concentration in daytime (solar angle >15), ppmv',
'Ndep':'Total N deposition over a day, g/m2/day'


III. Annual file - DUKE_forcing_y.txt

'Rainf':'Yearly Total Precipitation',
'Tair':'Yearly Mean air temperature',
'PSurf':'Yearly mean surface barometric pressure, Pa',
'aCO2':'Yearly mean ambient CO2 concentration in daytime (solar angle > 15), ppmv',
'eCO2':'Yearly mean elevated treatment CO2 concentration in daytime (solar angle >15), ppmv',
'Ndep':'Total N deposition over a year, g/m2/day'


IV. Corrections, repair and manipulations to the original data

Please see https://github.com/sharma-bharat/PrepMetData/blob/main/README.md
Issued addressed:
- Fixing Duplicate Index
- Fixing Incorrect Time Data
- Formula Used to derive Variables like SWdown and LWdown
- Codes to download pressure and Wind Data from ERA5