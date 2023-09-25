# MET Data Processing to FACEMDS Format
#=============================
#by Bharat Sharma and Anthony Walker

## Requisites
# ============
#The processed files from the plot data are saved at `/Users/ud4/Documents/FACEMDS/MET_Data_Processing/Oren_2022_Met_Data_processed/` <br>
#There are two files: <br>
#1. `Processed_Duke_Met_Data_All_Vars_30m_FV.csv`: FillValue/Missing = -6999.0
#2. `Processed_Duke_Met_Data_All_Vars_30m.csv`: FillValue/Missing = NaN or empty

# importing libraries
import xarray as xr
import glob
from datetime import datetime
import cftime
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr

# paths

paths = {}
paths ["ELM-DUKE"] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/ELM_Data/data/atm/datm7/CLM1PT_data/1x1pt_US-DUK/"
paths ["FACEMDS_Walker2018"] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/Walker_2018_FATES_MDS/data/"
paths ["DukeFACE_Oren2022"] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/Oren_2022_DUKE_Met/data/"
paths ["Save_Processed"] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/Oren_2022_Met_Data_processed/"

# Read the Processed Data
df_all_vars_30m_FV = pd.read_csv(f"{paths['Save_Processed']}Processed_Duke_Met_Data_All_Vars_30m_FV.csv",index_col=0)
df_all_vars_30m = pd.read_csv(f"{paths['Save_Processed']}Processed_Duke_Met_Data_All_Vars_30m.csv",index_col=0)

col_names = [
    'Year',
    'DOY',
    'Time',
    'Rainf',
    'Tair',
    'RH',
    'VPD',
    'PAR',
    'SM',
    'SWP',
    'SVP',
    'Rn',
    'SLT',
    'Wind',
    'PSurf',
    'aCO2',
    'eCO2',
    'Ndep',
    'SolarElevation',
]

## Creating DUKE_forcing_h.txt
# ----------------------------

# Dictionary of columns and their descriptions
dict_cols = {
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
'SWdown':'Incident or downward short-wave radiation, W/m2',
'SWdown_f':'gap-filling flag', 
'LWdown':'Incident or downward long-wave radiation, W/m2',
'LWdown_f':'gap-filling flag', 
'Wind':'Mean wind speed over a time step of measurement, m/s',
'Wind_f':'gap-filling flag',
'PSurf': 'Surface barometric pressure, Pa',
'PSurf_f':'gap-filling flag',
'aCO2': 'Daily mean ambient CO2 concentration in daytime (solar angle > 15), ppmv',
'eCO2': 'Daily mean elevated treatment CO2 concentration in daytime (solar angle >15), ppmv',
'Ndep': 'Total N deposition over a time step of measurement (30 minutes), g/m2/(30-minute)',
'SolarElevation': 'Solar elevation angle, degree',
}

# Units of Saving File Format
# Change: Rainf, mm to kg/m2/s; factor_add = 0,  factor_multiple = 1;  
# Change: Tair,  deg C to K   ; factor_add = 273.15,  factor_multiple = 1;
# Change: RH,    '' to %      ; factor_add = 0,  factor_multiple = 100;
# Change: VPD,   kPa to Pa    ; factor_add = 0,  factor_multiple = 1000;
# Change: SLT,    deg C to K   ; factor_add = 273.15,  factor_multiple = 1;
dict_units = {
    'Rainf':'kg/m2/s', # 'mm' = 'kg/m2/s'
    'Tair':'K',
    'RH':'%',
    'VPD':'Pa',
    'PAR':'umol/m2/s',
    'SM':'',
    'SWP':'',
    'SVP':'kPa',
    'Rn':'umol/m2/s',
    'SLT':'K',
    'SWdown':'W/m2',
    'LWdown':'W/m2',
    'Wind':'m/s',
    'PSurf': 'Pa',
    'aCO2': 'ppmv',
    'eCO2': 'ppmv',
    'Ndep': 'g/m2/(30-minute)',
    'SolarElevation':'degree',
}



df_h = pd.DataFrame(columns=dict_cols.keys())
df_h['YEAR'] = df_all_vars_30m['Year'].astype(int)
df_h['DOY'] = df_all_vars_30m['DOY'].astype(int)
df_h['HRMIN'] = df_all_vars_30m['Time'].astype(int)
df_h['DTIME'] = (df_all_vars_30m['Time']/(24*60)+df_all_vars_30m['DOY']).astype(float)
df_h[list(dict_units.keys())] = df_all_vars_30m[list(dict_units.keys())]
gap_fill_cols = [col + '_f' for col in list(dict_units.keys())]
df_h[gap_fill_cols] = 0

for k_var in ['Rainf']:
    factor_add = 0
    factor_multiple = 1/(60*30) # From per day to per second for 48 timesteps in a day i.e. 24*60*60/48
    df_h[k_var] = df_h[k_var]*factor_multiple+factor_add
for k_var in ['Tair','SLT']:
    factor_add = 273.15
    factor_multiple = 1
    df_h[k_var] = df_h[k_var]*factor_multiple+factor_add
for k_var in ['RH']:
    factor_add = 0
    factor_multiple = 100
    df_h[k_var] = df_h[k_var]*factor_multiple+factor_add
for k_var in ['VPD']:
    factor_add = 0
    factor_multiple = 1000
    df_h[k_var] = df_h[k_var]*factor_multiple+factor_add
# Adding the unit row below first row
unit_row = pd.Series([dict_units.get(col, '') for col in df_h.columns], index=df_h.columns)
df_h_save = pd.concat([pd.DataFrame([unit_row]), df_h.iloc[:]]).reset_index(drop=True)

### Correcting the Flags
'''
**'gap-filling flag,** <br>
0 = measured, <br>
1 = derived from other variables, <br>
2 = filled by exising FACEMDS data, <br>
3 = filled by data from nearby weather station, <br>
4 = filled by using ERA5 data'<br>
'''

import pickle
fname_start_dates = f"{paths['Save_Processed']}start_dates_Oren.pkl"
# Open the pickle file for reading in binary mode
with open(fname_start_dates, 'rb') as file:
    dict_start_dates = pickle.load(file)

# Copied from FACE MDS existing till the start date of dataset
for col_name in ['PAR','RH','Tair','Rainf','VPD'] :
    flag_col = f'{col_name}_f'
    year_tmp = dict_start_dates [col_name].Year
    doy_tmp = dict_start_dates [col_name].DOY
    time_tmp = dict_start_dates [col_name].Time
    filter_tmp = df_h_save[(df_h_save['YEAR'] == year_tmp) & (df_h_save['DOY'] == doy_tmp) & (df_h_save['HRMIN'] == time_tmp)]
    df_h_save[flag_col].iloc[1:filter_tmp.index[0]] = 2
    
# ERA5 Source
for col_name in ['Wind','PSurf']:
    flag_col = f'{col_name}_f'
    df_h_save[flag_col][1:] = 4
    
# Copied from FACE MDS existing
for col_name in ['aCO2', 'eCO2', 'Ndep', 'SolarElevation']:
    flag_col = f'{col_name}_f'
    df_h_save[flag_col][1:] = 2
    
# Derived based on a formula and other variabls
for col_name in ['SWdown', 'LWdown']:
    flag_col = f'{col_name}_f'
    df_h_save[flag_col][1:] = 1

# Saving the output
df_h_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_h.csv")
df_h_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_h.txt")
fill_value = -6999.
df_h_fv_save = df_h_save.fillna(fill_value)
df_h_fv_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_h_fv.csv")
df_h_fv_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_h_fv.txt")

# with time 
df_h_wTime = df_h.copy(deep=True)
df_h_wTime['Time'] = df_all_vars_30m['Date']
df_h_wTime['Time'] = pd.to_datetime(df_h_wTime['Time'])
df_h_wTime = df_h_wTime.set_index('Time')
df_h_wTime
# ----------------------------

## Creating DUKE_forcing_d.txt
# by talking average across all variables except 
# PAR, Rn, and Rainf where we computed the sum, resulting in Total daily metrics
# ----------------------------

# Units of Saving File Format
# Change: Rainf, kg/m2/s to 'kg/m2/d'; factor_add = 0,  factor_multiple = 30*60 (b/c 24*60*60/48);  
# Change: PAR/Rn, umol/m2/s to 'mol/m2/d'; factor_add = 0,  factor_multiple = 30*60 *10(--6);  
dict_units_d = {
    'Rainf':'kg/m2/d', # 'kg/m2/s' to 'kg/m2/d'
    'Tair':'K',
    'RH':'%',
    'VPD':'Pa',
    'PAR':'mol/m2/d',
    'SM':'',
    'SWP':'',
    'SVP':'kPa',
    'Rn':'mol/m2/d',
    'SLT':'K',
    'SWdown':'W/m2',
    'LWdown':'W/m2',
    'Wind':'m/s',
    'PSurf': 'Pa',
    'aCO2': 'ppmv',
    'eCO2': 'ppmv',
    'Ndep': 'g/m2/d',
    'SolarElevation':'degree',
}

# Using the datetime index to calculate means
df_d_wTime = df_h_wTime.resample('D').mean()
df_d_wTime= df_d_wTime.drop(['DTIME','HRMIN'], axis=1)
df_d_wTime['DOY'] = round(df_d_wTime['DOY']).astype('int')

# For Variables to convert from per second to per day if the original Tstep is of 30 mins
factor_multiple_s2d = 30*60

# For Rainf we need to take sum
df_d_wTime['Rainf'] = df_h_wTime[['YEAR', 'DOY', 'Rainf']].resample('D').sum()['Rainf'] * factor_multiple_s2d
# For PAR we need to take sum
df_d_wTime['PAR'] = df_h_wTime[['YEAR', 'DOY', 'PAR']].resample('D').sum()['PAR'] * factor_multiple_s2d *10**(-6)
# For Rn we need to take sum
df_d_wTime['Rn'] = df_h_wTime[['YEAR', 'DOY', 'Rn']].resample('D').sum()['Rn'] * factor_multiple_s2d *10**(-6)
# For Ndep we need to take sum
df_d_wTime['Ndep'] = df_h_wTime[['YEAR', 'DOY', 'Ndep']].resample('D').sum()['Ndep'] * factor_multiple_s2d 

df_d_wTime['YEAR'] =  df_d_wTime['YEAR'].astype('int')

# Reset index to columns
df_d = df_d_wTime.reset_index()
df_d = df_d.drop('Time', axis =1)

# Adding the unit row below first row
unit_row = pd.Series([dict_units_d.get(col, '') for col in df_d.columns], index=df_d.columns)
df_d_save = pd.concat([pd.DataFrame([unit_row]), df_d.iloc[:]]).reset_index(drop=True)

#Saving the dataframes
df_d_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_d.csv")
df_d_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_d.txt")
fill_value = -6999.
df_d_fv_save = df_d_save.fillna(fill_value)
df_d_fv_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_d_fv.csv")
df_d_fv_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_d_fv.txt")

# ----------------------------

## Creating DUKE_forcing_y.txt
# ----------------------------
# Units of Saving File Format
# Change: Rainf, kg/m2/d to 'kg/m2/y'; ;  
# Units of Saving File Format
# Change: Rainf, kg/m2/d to 'kg/m2/y'; ;  
dict_units_y = {
    'Rainf':'kg/m2/y', # 'kg/m2/d' to 'kg/m2/y'
    'Tair':'K',
    'PSurf': 'Pa',
    'aCO2': 'ppmv',
    'eCO2': 'ppmv',
    'Ndep': 'g/m2/y',
}
# Using the datetime index to calculate means
df_y_wTime = df_d_wTime.resample('Y').mean()
df_y_wTime= df_y_wTime.drop('DOY', axis=1)
df_y_wTime['YEAR'] = round(df_y_wTime['YEAR']).astype('int')

# For Rainf we need to take sum
df_y_wTime['Rainf'] = df_d_wTime[['YEAR', 'Rainf']].resample('Y').sum()['Rainf']
df_y_wTime['YEAR'] =  df_y_wTime['YEAR'].astype('int')

# For Ndep we need to take sum
df_y_wTime['Ndep'] = df_d_wTime[['YEAR', 'Ndep']].resample('Y').sum()['Ndep']
df_y_wTime['YEAR'] =  df_y_wTime['YEAR'].astype('int')


# Reset index to columns
df_y = df_y_wTime.reset_index()
df_y = df_y.drop('Time', axis =1)

# Adding the unit row below first row
unit_row = pd.Series([dict_units_y.get(col, '') for col in df_y.columns], index=df_y.columns)
df_y_save = pd.concat([pd.DataFrame([unit_row]), df_y.iloc[:]]).reset_index(drop=True)

#Saving the dataframes
df_y_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_y.csv")
df_y_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_y.txt")
fill_value = -6999.
df_y_fv_save = df_y_save.fillna(fill_value)
df_y_fv_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_y_fv.csv")
df_y_fv_save.to_csv(f"{paths['Save_Processed']}DUKE_forcing_y_fv.txt")

# Making netcdf files

## Creating DUKE_forcing_h.nc
cols_netcdf_h = ['YEAR',
 'DOY',
 'HRMIN',
 'Rainf',
  'Tair',
  'RH',
  'VPD',
  'PAR',
  'SM',
  'SWP',
  'SVP',
  'Rn',
  'SLT',
  'SWdown',
  'LWdown',
  'Wind',
  'PSurf',
  'aCO2',
  'eCO2',
  'Ndep',
  'SolarElevation'
                ]

# Convert the DataFrame to an xarray Dataset

df_h_wTime_fv = df_h_wTime.fillna(fill_value)
ds = xr.Dataset.from_dataframe(df_h_wTime_fv[cols_netcdf_h])
ds = ds.expand_dims({'lat': 1}).assign_coords({'lat':[35.9782] })
ds = ds.expand_dims({'lon': 1}).assign_coords({'lon':[-79.0942] })
ds['lon'].attrs['units'] = "degrees_east"
ds['lon'].attrs['long_name'] = "Longitude"
ds['lat'].attrs['units'] = "degrees_north"
ds['lat'].attrs['long_name'] = "Latitude"
ds['YEAR'].attrs['units'] = ""
ds['YEAR'].attrs['long_name'] = "Year of Measurement"
ds['DOY'].attrs['units'] = ""
ds['DOY'].attrs['long_name'] = "Day of Measurement"
ds['HRMIN'].attrs['units'] = ""
ds['HRMIN'].attrs['long_name'] = "Hour:minute - marked at the middle of measurement interval with last two digits as minute"
# Adding attributes to variables
for k_var in cols_netcdf_h: # looping of columns
    if k_var in ['YEAR','DOY','HRMIN']:
        pass
    else:
        ds[k_var].attrs['units'] = dict_units[k_var]
        ds[k_var].attrs['missing_value'] = fill_value
        ds[k_var].attrs['long_name'] = dict_cols[k_var]
        ds[k_var].attrs['associate'] = "Time lat lon"
        ds[k_var].attrs['axis'] = "TYX"

    
    
# Add global attributes to the Dataset
ds.attrs['site_id'] = "DUKE"
ds.attrs['title'] = "Half-hourly forcing data from DUKE Forest FACE, North Carolina, USA"
ds.attrs['history'] = "File Origin - This file was created at Oak Ridge National Laboratory for FACE model data synthesis"
ds.attrs['creation_date'] = "Sep 25, 2023" ;
ds.attrs['contact'] = 'Bharat Sharma (sharmabd@ornl.gov), Anthony Walker (walkerp@ornl.gov)'
ds.to_netcdf(f"{paths['Save_Processed']}DUKE_forcing_h.nc")

print(f"Dataset saved to {paths['Save_Processed']}DUKE_forcing_h.nc")

## Creating DUKE_forcing_d.nc

cols_netcdf_d = ['YEAR',
 'DOY',
 'Rainf',
  'Tair',
  'RH',
  'VPD',
  'PAR',
  'SM',
  'SWP',
  'SVP',
  'Rn',
  'SLT',
  'SWdown',
  'LWdown',
  'Wind',
  'PSurf',
  'aCO2',
  'eCO2',
  'Ndep',
  'SolarElevation'
                ]

# Convert the DataFrame to an xarray Dataset
df_d_wTime_fv = df_d_wTime.fillna(fill_value)
ds = xr.Dataset.from_dataframe(df_d_wTime_fv[cols_netcdf_d])
ds = ds.expand_dims({'lat': 1}).assign_coords({'lat':[35.9782] })
ds = ds.expand_dims({'lon': 1}).assign_coords({'lon':[-79.0942] })
ds['lon'].attrs['units'] = "degrees_east"
ds['lon'].attrs['long_name'] = "Longitude"
ds['lat'].attrs['units'] = "degrees_north"
ds['lat'].attrs['long_name'] = "Latitude"
ds['YEAR'].attrs['units'] = ""
ds['YEAR'].attrs['long_name'] = "Year of Measurement"
ds['DOY'].attrs['units'] = ""
ds['DOY'].attrs['long_name'] = "Day of Measurement"
# Adding attributes to variables
for k_var in cols_netcdf_d: # looping of columns
    if k_var in ['YEAR','DOY','HRMIN']:
        pass
    else:
        ds[k_var].attrs['units'] = dict_units_d[k_var]
        ds[k_var].attrs['missing_value'] = fill_value
        ds[k_var].attrs['long_name'] = dict_cols[k_var]
        ds[k_var].attrs['associate'] = "Time lat lon"
        ds[k_var].attrs['axis'] = "TYX"

# Add global attributes to the Dataset
ds.attrs['site_id'] = "DUKE"
ds.attrs['title'] = "Daily forcing data from DUKE Forest FACE, North Carolina, USA"
ds.attrs['history'] = "File Origin - This file was created at Oak Ridge National Laboratory for FACE model data synthesis"
ds.attrs['creation_date'] = "Sep 25, 2023" ;
ds.attrs['contact'] = 'Bharat Sharma (sharmabd@ornl.gov), Anthony Walker (walkerp@ornl.gov)'
ds.to_netcdf(f"{paths['Save_Processed']}DUKE_forcing_d.nc")

print(f"Dataset saved to {paths['Save_Processed']}DUKE_forcing_d.nc")

## Creating DUKE_forcing_y.nc

cols_netcdf_y = ['YEAR',
 'Rainf',
  'Tair',
  'PSurf',
  'aCO2',
  'eCO2',
  'Ndep',
                ]

# Convert the DataFrame to an xarray Dataset
df_y_wTime_fv = df_y_wTime.fillna(fill_value)
ds = xr.Dataset.from_dataframe(df_y_wTime_fv[cols_netcdf_y])
ds = ds.expand_dims({'lat': 1}).assign_coords({'lat':[35.9782] })
ds = ds.expand_dims({'lon': 1}).assign_coords({'lon':[-79.0942] })
ds['lon'].attrs['units'] = "degrees_east"
ds['lon'].attrs['long_name'] = "Longitude"
ds['lat'].attrs['units'] = "degrees_north"
ds['lat'].attrs['long_name'] = "Latitude"
ds['YEAR'].attrs['units'] = ""
ds['YEAR'].attrs['long_name'] = "Year of Measurement"
# Adding attributes to variables
for k_var in cols_netcdf_y: # looping of last 10 columns
    if k_var in ['YEAR','']:
        pass
    elif k_var in ['PAR','Rn']:
        ds[k_var].attrs['units'] = dict_units_y[k_var]
        ds[k_var].attrs['missing_value'] = fill_value
        ds[k_var].attrs['long_name'] = f"Mean {dict_cols[k_var]}"
        ds[k_var].attrs['associate'] = "Time lat lon"
        ds[k_var].attrs['axis'] = "TYX"        
    else:
        ds[k_var].attrs['units'] = dict_units_y[k_var]
        ds[k_var].attrs['missing_value'] = fill_value
        ds[k_var].attrs['long_name'] = dict_cols[k_var]
        ds[k_var].attrs['associate'] = "Time lat lon"
        ds[k_var].attrs['axis'] = "TYX"

# Add global attributes to the Dataset
ds.attrs['site_id'] = "DUKE"
ds.attrs['title'] = "Yearly forcing data from DUKE Forest FACE, North Carolina, USA"
ds.attrs['history'] = "File Origin - This file was created at Oak Ridge National Laboratory for FACE model data synthesis"
ds.attrs['creation_date'] = "Sep 25, 2023" ;
ds.attrs['contact'] = 'Bharat Sharma (sharmabd@ornl.gov), Anthony Walker (walkerp@ornl.gov)'
ds.to_netcdf(f"{paths['Save_Processed']}DUKE_forcing_y.nc")

print(f"Dataset saved to {paths['Save_Processed']}DUKE_forcing_y.nc")

