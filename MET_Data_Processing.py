# Met Data Analysis
# by Bharat Sharma and Anthony Walker <br>
# sharmabd@ornl.gov <br>
# Site: US-DUKE
# Features:
## It will also print any duplicate values in the original data
## We have a fix that issue (read Readme_MetData.md)


# importing libraries
import xarray as xr
import glob
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

## Files
# - **ELM-DUKE** : has the nc files that we use to run the current version of ELM-FATES <br>
# - **[FACEMDS_Walker2018](https://data.ess-dive.lbl.gov/view/ess-dive-7807cf86f1dd42a-20181127T173047368940)** : Processed Data (My focus: DUKE) <br>
# - **[DukeFACE_Oren2022](https://data.ess-dive.lbl.gov/view/doi:10.15485/1895465)**: Updated DUKE met Data

# Make plots yes/on
make_plots = 'no'

# paths

paths = {}
paths[
    "ELM-DUKE"
] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/ELM_Data/data/atm/datm7/CLM1PT_data/1x1pt_US-DUK/"
paths[
    "FACEMDS_Walker2018"
] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/Walker_2018_FATES_MDS/data/"
paths[
    "DukeFACE_Oren2022"
] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/Oren_2022_DUKE_Met/data/"
paths[
    "Save_Processed"
] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/Oren_2022_Met_Data_processed/"

## ELM Duke Data

key = "ELM-DUKE"
ds_elm_all = xr.open_mfdataset(glob.glob(f"{paths[key]}*.nc"))

# cftime to datetime


def cftime_to_dtime(cftime_object):
    # Convert cftime.DatetimeNoLeap to datetime
    datetime_object = datetime(
        year=cftime_object.year,
        month=cftime_object.month,
        day=cftime_object.day,
        hour=cftime_object.hour,
        minute=cftime_object.minute,
        second=cftime_object.second,
        microsecond=cftime_object.microsecond,
    )
    return datetime_object

if make_plots in ['y','yes']:
    time_objects = [cftime_to_dtime(t) for t in ds_elm_all.time.values]
    fig1 = plt.figure(figsize=(20, 8))
    plt.scatter(
        x=time_objects,
        y=ds_elm_all.PRECTmms,
        c=ds_elm_all.PRECTmms,
        cmap="Reds",
        s=22,
        marker="o",
        label=ds_elm_all.PRECTmms.units,
    )
    plt.title(f"{key}")
    plt.legend()
    plt.show()

## Investigating FACEMDS_Walker2018

key = "FACEMDS_Walker2018"
ds_Walker_h = xr.open_dataset(f"{paths[key]}DUKE_forcing_h.nc", decode_times=False)

from datetime import datetime, timedelta


def seconds_to_datetime(seconds, reference_date):
    return reference_date + timedelta(seconds=seconds)


reference_date = datetime(1996, 1, 1, 0, 0)

list_seconds = ds_Walker_h.TIME[...].values
list_seconds = np.asarray(list_seconds, dtype=float)
resulting_datetime = [
    seconds_to_datetime(seconds, reference_date) for seconds in list_seconds
]
if make_plots in ['y','yes']:
    fig2 = plt.figure(figsize=(20, 8))
    plt.scatter(
        x=resulting_datetime,
        y=ds_Walker_h["Rainf"],
        c=ds_Walker_h["Rainf"],
        cmap="Reds",
        s=22,
        marker="o",
        label=ds_Walker_h["Rainf"].units,
    )
    plt.title(f"{key}")
    plt.legend()
    plt.show()

### FACEMDS_Walker CSV files
key = "FACEMDS_Walker2018"
#### Reading Hourly files
df_FACEMDS = {}
df_FACEMDS["h"] = pd.read_csv(glob.glob(f"{paths[key]}DUKE*_h.txt")[0], delimiter="\t")
print(len(df_FACEMDS["h"].columns))
df_FACEMDS["h"].columns


## Oren New Met Data
# Dictionary that will contain the dataframes of all the variables from Oren New Met Data
dict_dfs_common = {}

# Specify the file path for storing any duplicate index information
DuplicateDataFilename = f'{paths["Save_Processed"]}DuplicateDukeData.txt'
# This file path is to document the duplicate Data in Duke Data

key = "DukeFACE_Oren2022"

### AT : Tair
var_key = "AT"
face_var_key = "Tair"

plots_cols = ["R1uat", "R2uat", "R3uat", "R4uat", "R5uat", "R6uat", "R7uat", "R8uat"]

"""
Oren Data has 3 dirs for AT for different plots, time periods, and sensors. <br>
I intend to use `*_gl.csv` files; i believe these are gap filled. <br>
I will save the mean of the plots for FACEMDS <br>
The common sensor data will be used among all the files.
"""
# gap filled files
files = sorted(glob.glob(f"{paths[key]}DukeFACE_{var_key}*/*_gf.csv"))


# Define a custom sorting key function to extract the year from the file path
def extract_year(filepath):
    return int(filepath.split(f"{var_key}")[-1][:4])  # Year


# Sort the list of file paths based on the Year
sorted_filepaths = sorted(files, key=extract_year)

# Create a common dataframe
# List of column names
common_columns = ["Year", "JDT", "DOY", "Time", f"{face_var_key}"]

# Create an empty DataFrame with the specified columns
df_tmp_common_gf = pd.DataFrame(columns=common_columns)

selected_columns = plots_cols
# Open the file in append mode and write text
with open(DuplicateDataFilename, "a") as duplicate_data:
    for file in sorted_filepaths:
        i_dup_file_name = 0  # to save file name once for multiple dupliate entries
        df_tmp = pd.read_csv(file)
        for lbl, gr in df_tmp.groupby(["Year", "DOY", "Time"]):
            if len(gr) > 1:
                if i_dup_file_name == 0:
                    duplicate_data.write(f"{'/'.join(file.split('/')[-3:])}\n")
                duplicate_data.write(f"{gr.iloc[:,:4]}\n")
                i_dup_file_name += 1

        # Calculate the mean of selected columns
        df_tmp[f"{face_var_key}"] = round(df_tmp[selected_columns].mean(axis=1), 2)
        # only saving the common columns
        df_tmp = df_tmp[common_columns]
        # Appending all the common columns to the common dataframe
        df_tmp_common_gf = df_tmp_common_gf.append(df_tmp)
    duplicate_data.write(f"\n")
dict_dfs_common[f"{face_var_key}"] = df_tmp_common_gf

### Precip : Rainf
# units: mm

key = "DukeFACE_Oren2022"
var_key = "Precip"
face_var_key = "Rainf"

plots_cols = ["FACE.PO"]

# gap filled files
files = sorted(glob.glob(f"{paths[key]}DukeFACE_{var_key}*/*_gf.csv"))

# Sort the list of file paths based on the Year
sorted_filepaths = sorted(files, key=extract_year)

# Create a common dataframe
# List of column names
common_columns = ["Year", "JDT", "DOY", "Time", f"{face_var_key}"]

# Create an empty DataFrame with the specified columns
df_tmp_common_gf = pd.DataFrame(columns=common_columns)


selected_columns = plots_cols
# Open the file in append mode and write text
with open(DuplicateDataFilename, "a") as duplicate_data:
    for file in sorted_filepaths:
        i_dup_file_name = 0  # to save file name once for multiple dupliate entries
        df_tmp = pd.read_csv(file)
        for lbl, gr in df_tmp.groupby(["Year", "DOY", "Time"]):
            if len(gr) > 1:
                if i_dup_file_name == 0:
                    duplicate_data.write(f"{'/'.join(file.split('/')[-3:])}\n")
                duplicate_data.write(f"{gr.iloc[:,:4]}\n")
                i_dup_file_name += 1

        # Calculate the mean of selected columns
        df_tmp[f"{face_var_key}"] = round(df_tmp[selected_columns].mean(axis=1), 2)
        # only saving the common columns
        df_tmp = df_tmp[common_columns]
        # Appending all the common columns to the common dataframe
        df_tmp_common_gf = df_tmp_common_gf.append(df_tmp)
    duplicate_data.write(f"\n")
dict_dfs_common[f"{face_var_key}"] = df_tmp_common_gf

### RH : RH
# Relative Humidity
key = "DukeFACE_Oren2022"
var_key = "RH"
face_var_key = "RH"
plots_cols = ["R1urh", "R2urh", "R3urh", "R4urh", "R5urh", "R6urh", "R7urh", "R8urh"]
# gap filled files
files = sorted(glob.glob(f"{paths[key]}DukeFACE_{var_key}*/*_gf.csv"))

# Sort the list of file paths based on the Year
sorted_filepaths = sorted(files, key=extract_year)

# Create a common dataframe
# List of column names
common_columns = ["Year", "JDT", "DOY", "Time", f"{face_var_key}"]

# Create an empty DataFrame with the specified columns
df_tmp_common_gf = pd.DataFrame(columns=common_columns)


selected_columns = plots_cols
# Open the file in append mode and write text
with open(DuplicateDataFilename, "a") as duplicate_data:
    for file in sorted_filepaths:
        i_dup_file_name = 0  # to save file name once for multiple dupliate entries
        df_tmp = pd.read_csv(file)
        for lbl, gr in df_tmp.groupby(["Year", "DOY", "Time"]):
            if len(gr) > 1:
                if i_dup_file_name == 0:
                    duplicate_data.write(f"{'/'.join(file.split('/')[-3:])}\n")
                duplicate_data.write(f"{gr.iloc[:,:4]}\n")
                i_dup_file_name += 1

        # Calculate the mean of selected columns
        df_tmp[f"{face_var_key}"] = round(df_tmp[selected_columns].mean(axis=1), 2)
        # only saving the common columns
        df_tmp = df_tmp[common_columns]
        # Appending all the common columns to the common dataframe
        df_tmp_common_gf = df_tmp_common_gf.append(df_tmp)
    duplicate_data.write(f"\n")
dict_dfs_common[f"{face_var_key}"] = df_tmp_common_gf

### SM:SM
# (not in existing WalkerFACEMDS Data)<br>
# Soil moisture integrates measurements from 0 to 30cm depth

key = "DukeFACE_Oren2022"
var_key = "SM"
face_var_key = "SM"
plots_cols = ["R1tdr", "R2tdr", "R3tdr", "R4tdr", "R5tdr", "R6tdr", "R7tdr", "R8tdr"]

# gap filled files
files = sorted(glob.glob(f"{paths[key]}DukeFACE_{var_key}*/*_gf.csv"))

# Sort the list of file paths based on the Year
sorted_filepaths = sorted(files, key=extract_year)

# Create a common dataframe
# List of column names
common_columns = ["Year", "JDT", "DOY", "Time", f"{face_var_key}"]

# Create an empty DataFrame with the specified columns
df_tmp_common_gf = pd.DataFrame(columns=common_columns)


selected_columns = plots_cols
# Open the file in append mode and write text
with open(DuplicateDataFilename, "a") as duplicate_data:
    for file in sorted_filepaths:
        i_dup_file_name = 0  # to save file name once for multiple dupliate entries
        df_tmp = pd.read_csv(file)
        for lbl, gr in df_tmp.groupby(["Year", "DOY", "Time"]):
            if len(gr) > 1:
                if i_dup_file_name == 0:
                    duplicate_data.write(f"{'/'.join(file.split('/')[-3:])}\n")
                duplicate_data.write(f"{gr.iloc[:,:4]}\n")
                i_dup_file_name += 1

        # Calculate the mean of selected columns
        df_tmp[f"{face_var_key}"] = round(df_tmp[selected_columns].mean(axis=1), 2)
        # only saving the common columns
        df_tmp = df_tmp[common_columns]
        # Appending all the common columns to the common dataframe
        df_tmp_common_gf = df_tmp_common_gf.append(df_tmp)
    duplicate_data.write(f"\n")
dict_dfs_common[f"{face_var_key}"] = df_tmp_common_gf

### SWP: SWP
# not in existing Walker 2018 <br>
# (Soil water potential) <br>
# from 2007 to 2012

key = "DukeFACE_Oren2022"
var_key = "SWP"
face_var_key = "SWP"
plots_cols = ["R1swp", "R2swp", "R3swp", "R4swp", "R5swp", "R6swp"]

# gap filled files
files = sorted(glob.glob(f"{paths[key]}DukeFACE_{var_key}*/*_gf.csv"))

# Sort the list of file paths based on the Year
sorted_filepaths = sorted(files, key=extract_year)

# Create a common dataframe
# List of column names
common_columns = ["Year", "JDT", "DOY", "Time", f"{face_var_key}"]

# Create an empty DataFrame with the specified columns
df_tmp_common_gf = pd.DataFrame(columns=common_columns)


selected_columns = plots_cols
# Open the file in append mode and write text
with open(DuplicateDataFilename, "a") as duplicate_data:
    for file in sorted_filepaths:
        i_dup_file_name = 0  # to save file name once for multiple dupliate entries
        df_tmp = pd.read_csv(file)
        for lbl, gr in df_tmp.groupby(["Year", "DOY", "Time"]):
            if len(gr) > 1:
                if i_dup_file_name == 0:
                    duplicate_data.write(f"{'/'.join(file.split('/')[-3:])}\n")
                duplicate_data.write(f"{gr.iloc[:,:4]}\n")
                i_dup_file_name += 1

        # Calculate the mean of selected columns
        df_tmp[f"{face_var_key}"] = round(df_tmp[selected_columns].mean(axis=1), 2)
        # only saving the common columns
        df_tmp = df_tmp[common_columns]
        # Appending all the common columns to the common dataframe
        df_tmp_common_gf = df_tmp_common_gf.append(df_tmp)
    duplicate_data.write(f"\n")
dict_dfs_common[f"{face_var_key}"] = df_tmp_common_gf


### SVP: SVP

# Saturated Vapor Pressure; Not in Walker 2018<br>
# units kPa

key = "DukeFACE_Oren2022"
var_key = "SVP"
face_var_key = "SVP"
plots_cols = [
    "R1usvp",
    "R2usvp",
    "R3usvp",
    "R4usvp",
    "R5usvp",
    "R6usvp",
    "R7usvp",
    "R8usvp",
]
# gap filled files
files = sorted(glob.glob(f"{paths[key]}DukeFACE_{var_key}*/*_gf.csv"))

# Sort the list of file paths based on the Year
sorted_filepaths = sorted(files, key=extract_year)

# Create a common dataframe
# List of column names
common_columns = ["Year", "JDT", "DOY", "Time", f"{face_var_key}"]

# Create an empty DataFrame with the specified columns
df_tmp_common_gf = pd.DataFrame(columns=common_columns)


selected_columns = plots_cols
# Open the file in append mode and write text
with open(DuplicateDataFilename, "a") as duplicate_data:
    for file in sorted_filepaths:
        i_dup_file_name = 0  # to save file name once for multiple dupliate entries
        df_tmp = pd.read_csv(file)
        for lbl, gr in df_tmp.groupby(["Year", "DOY", "Time"]):
            if len(gr) > 1:
                if i_dup_file_name == 0:
                    duplicate_data.write(f"{'/'.join(file.split('/')[-3:])}\n")
                duplicate_data.write(f"{gr.iloc[:,:4]}\n")
                i_dup_file_name += 1

        # Calculate the mean of selected columns
        df_tmp[f"{face_var_key}"] = round(df_tmp[selected_columns].mean(axis=1), 2)
        # only saving the common columns
        df_tmp = df_tmp[common_columns]
        # Appending all the common columns to the common dataframe
        df_tmp_common_gf = df_tmp_common_gf.append(df_tmp)
    duplicate_data.write(f"\n")
dict_dfs_common[f"{face_var_key}"] = df_tmp_common_gf


### VPD: VPD
# Vapor pressure deficit <br>
# Units kPa <br>

key = "DukeFACE_Oren2022"
var_key = "VPD"
face_var_key = "VPD"
plots_cols = [
    "R1uvpd",
    "R2uvpd",
    "R3uvpd",
    "R4uvpd",
    "R5uvpd",
    "R6uvpd",
    "R7uvpd",
    "R8uvpd",
]

# gap filled files
files = sorted(glob.glob(f"{paths[key]}DukeFACE_{var_key}*/*_gf.csv"))

# Sort the list of file paths based on the Year
sorted_filepaths = sorted(files, key=extract_year)

# Create a common dataframe
# List of column names
common_columns = ["Year", "JDT", "DOY", "Time", f"{face_var_key}"]

# Create an empty DataFrame with the specified columns
df_tmp_common_gf = pd.DataFrame(columns=common_columns)


selected_columns = plots_cols
# Open the file in append mode and write text
with open(DuplicateDataFilename, "a") as duplicate_data:
    for file in sorted_filepaths:
        i_dup_file_name = 0  # to save file name once for multiple dupliate entries
        df_tmp = pd.read_csv(file)
        for lbl, gr in df_tmp.groupby(["Year", "DOY", "Time"]):
            if len(gr) > 1:
                if i_dup_file_name == 0:
                    duplicate_data.write(f"{'/'.join(file.split('/')[-3:])}\n")
                duplicate_data.write(f"{gr.iloc[:,:4]}\n")
                i_dup_file_name += 1

        # Calculate the mean of selected columns
        df_tmp[f"{face_var_key}"] = round(df_tmp[selected_columns].mean(axis=1), 2)
        # only saving the common columns
        df_tmp = df_tmp[common_columns]
        # Appending all the common columns to the common dataframe
        df_tmp_common_gf = df_tmp_common_gf.append(df_tmp)
    duplicate_data.write(f"\n")
dict_dfs_common[f"{face_var_key}"] = df_tmp_common_gf

### SLT: SLT

# Walker 2018 does not have this data <br>
# soil temperature from plots 1-6, one sensor per plot at 15 cm depth, measured in Degree Celsius. <br>
# I will be averaging over plots 2-6 because 1 was upgraded over time and 15 cm depth was not avaiable later

key = "DukeFACE_Oren2022"
var_key = "SLT"
face_var_key = "SLT"
plots_cols = ["R2slt", "R3slt", "R4slt", "R5slt", "R6slt"]

# gap filled files
files = sorted(glob.glob(f"{paths[key]}DukeFACE_{var_key}*/*_gf.csv"))

# Sort the list of file paths based on the Year
sorted_filepaths = sorted(files, key=extract_year)

# Create a common dataframe
# List of column names
common_columns = ["Year", "JDT", "DOY", "Time", f"{face_var_key}"]

# Create an empty DataFrame with the specified columns
df_tmp_common_gf = pd.DataFrame(columns=common_columns)


selected_columns = plots_cols
# Open the file in append mode and write text
with open(DuplicateDataFilename, "a") as duplicate_data:
    for file in sorted_filepaths:
        i_dup_file_name = 0  # to save file name once for multiple dupliate entries
        df_tmp = pd.read_csv(file)
        for lbl, gr in df_tmp.groupby(["Year", "DOY", "Time"]):
            if len(gr) > 1:
                if i_dup_file_name == 0:
                    duplicate_data.write(f"{'/'.join(file.split('/')[-3:])}\n")
                duplicate_data.write(f"{gr.iloc[:,:4]}\n")
                i_dup_file_name += 1

        # Calculate the mean of selected columns
        df_tmp[f"{face_var_key}"] = round(df_tmp[selected_columns].mean(axis=1), 2)
        # only saving the common columns
        df_tmp = df_tmp[common_columns]
        # Appending all the common columns to the common dataframe
        df_tmp_common_gf = df_tmp_common_gf.append(df_tmp)
    duplicate_data.write(f"\n")
dict_dfs_common[f"{face_var_key}"] = df_tmp_common_gf

### PAR: PAR

# The data from 2008-12 was averaged to get PAN and Rn <br>
# PAR - Photosynthetically active radiation, umol/m^2\*s <br>
# PAR had values from only one plot for 1997-2007 and from two plots for 2008-2012, which were averaged <br>

key = "DukeFACE_Oren2022"
var_key = "PAR"
face_var_key = "PAR"
plots_cols = ["PAR"]

exp_filename = "Rad"  # Exception in the filename

# gap filled files
files = sorted(glob.glob(f"{paths[key]}DukeFACE_{exp_filename}*/*_gf.csv"))


# Define a custom sorting key function to extract the year from the file path
def extract_year(filepath):
    return int(filepath.split(f"{exp_filename}")[-1][:4])  # Year


# Sort the list of file paths based on the Year
sorted_filepaths = sorted(files, key=extract_year)

# Create a common dataframe
# List of column names
common_columns = ["Year", "JDT", "DOY", "Time", f"{face_var_key}"]

# Create an empty DataFrame with the specified columns
df_tmp_common_gf = pd.DataFrame(columns=common_columns)


selected_columns = plots_cols
# Open the file in append mode and write text
with open(DuplicateDataFilename, "a") as duplicate_data:
    for file in sorted_filepaths:
        i_dup_file_name = 0  # to save file name once for multiple dupliate entries
        df_tmp = pd.read_csv(file)
        for lbl, gr in df_tmp.groupby(["Year", "DOY", "Time"]):
            if len(gr) > 1:
                if i_dup_file_name == 0:
                    duplicate_data.write(f"{'/'.join(file.split('/')[-3:])}\n")
                duplicate_data.write(f"{gr.iloc[:,:4]}\n")
                i_dup_file_name += 1

        # Calculate the mean of selected columns
        df_tmp[f"{face_var_key}"] = round(df_tmp[selected_columns].mean(axis=1), 2)
        # only saving the common columns
        df_tmp = df_tmp[common_columns]
        # Appending all the common columns to the common dataframe
        df_tmp_common_gf = df_tmp_common_gf.append(df_tmp)
    duplicate_data.write(f"\n")
dict_dfs_common[f"{face_var_key}"] = df_tmp_common_gf


### Rn: Rn

# Not in Walker 2018 <br>
# The data from 2008-12 was averaged to get PAN and Rn <br>
# Rn  - Net radiation. Q7 sensor before 2004, CNR1 thereafter, W/m^2  <br>
# Rn had values from only one plot for 1997-2007 and from two plots for 2008-2012, which were averaged <br>

key = "DukeFACE_Oren2022"
var_key = "Rn"
face_var_key = "Rn"
plots_cols = ["Rn"]

exp_filename = "Rad"  # Exception in the filename

# gap filled files
files = sorted(glob.glob(f"{paths[key]}DukeFACE_{exp_filename}*/*_gf.csv"))


# Define a custom sorting key function to extract the year from the file path
def extract_year(filepath):
    return int(filepath.split(f"{exp_filename}")[-1][:4])  # Year


# Sort the list of file paths based on the Year
sorted_filepaths = sorted(files, key=extract_year)

# Create a common dataframe
# List of column names
common_columns = ["Year", "JDT", "DOY", "Time", f"{face_var_key}"]

# Create an empty DataFrame with the specified columns
df_tmp_common_gf = pd.DataFrame(columns=common_columns)


selected_columns = plots_cols
# Open the file in append mode and write text
with open(DuplicateDataFilename, "a") as duplicate_data:
    for file in sorted_filepaths:
        i_dup_file_name = 0  # to save file name once for multiple dupliate entries
        df_tmp = pd.read_csv(file)
        for lbl, gr in df_tmp.groupby(["Year", "DOY", "Time"]):
            if len(gr) > 1:
                if i_dup_file_name == 0:
                    duplicate_data.write(f"{'/'.join(file.split('/')[-3:])}\n")
                duplicate_data.write(f"{gr.iloc[:,:4]}\n")
                i_dup_file_name += 1

        # Calculate the mean of selected columns
        df_tmp[f"{face_var_key}"] = round(df_tmp[selected_columns].mean(axis=1), 2)
        # only saving the common columns
        df_tmp = df_tmp[common_columns]
        # Appending all the common columns to the common dataframe
        df_tmp_common_gf = df_tmp_common_gf.append(df_tmp)
    duplicate_data.write(f"\n")
dict_dfs_common[f"{face_var_key}"] = df_tmp_common_gf

## Modification of the variables
### PAR: PAR 
#Saving the existing PAR in a new variable 'PAR_ori'. <br>
#Correcting existing values of PAR that reduce over time because of the degrading sensor <br>

dict_dfs_common['PAR_ori'] = dict_dfs_common['PAR']
dict_dfs_common['PAR_ori']['PAR_ori'] = dict_dfs_common['PAR_ori']['PAR']
dict_dfs_common['PAR_ori'] = dict_dfs_common['PAR_ori'].drop(['PAR'], axis=1)


for lbl,gr in dict_dfs_common['PAR'].groupby(['Year']):
    percentiles = np.arange(10, 101, 5)
    percentile_values = np.percentile(gr['PAR'], percentiles)

# Normalizing based on the 95th percent value
# Assuming that the sensor was good in the first 3 years
# Taking the average for first 3 years.
df_par_per = pd.DataFrame(data= np.array(per_95)[:,1], columns=['PAR_95p'])
df_par_per.index = np.asarray(np.array(per_95)[:,0], dtype=int)
df_par_per.index.name = 'Year'
base_95th = df_par_per.iloc[0:3].mean()

# Factor Reduction in amplitude
df_par_div_factor = df_par_per/base_95th

df_par_div_factor.columns = ['PAR_divisible_factor']
# we need to divide the original data with the factor corresponding the year of above array.

# Merge the original DataFrame with the factors DataFrame based on the 'Year' column
merged_df = dict_dfs_common['PAR'].merge(df_par_div_factor, on='Year', how='left')

# Calculate the modified 'PAR' column by dividing 'PAR' with 'PAR_95p'
merged_df['Modified_PAR'] = merged_df['PAR'] / merged_df['PAR_divisible_factor']

# Replacing the original PAR with new PAR values
merged_df['PAR'] = merged_df['Modified_PAR']
merged_df = merged_df.drop(['Modified_PAR', 'PAR_divisible_factor'], axis=1)
dict_dfs_common['PAR'] = merged_df

##  Calculating new variables
### SWdown: SWdown 
#Calculating it from PAR using the following function: <br>
#Using the corrected par (umolm-2s-1) we calculate SWDown in (W/m2) <br>

def PAR2SWdown (data, out_units='umolm-2s-1'):
    """
    by BS
    
    input units: umolm-2s-1
    -----------------------
    
    Returns:
    --------
    PAR in desired out_units.
    
    Source: https://www.researchgate.net/post/Can-I-convert-PAR-photo-active-radiation-value-of-micro-mole-M2-S-to-Solar-radiation-in-Watt-m2/59ca6422217e201e2b23415f/citation/download
    It says PAR is 45% of total Solar radiation. However, most models use a factor of 0.5 or 50%. We will go with 0.5. 
    And 1 W/m2 ≈ 4.6 μmole/m2/s 
    
    """
    if out_units == 'umolm-2s-1':
        conversion_factor = 2 # (1/0.5)
    if out_units == 'W/m2':
        conversion_factor = 2/4.6 # (1/0.5/4.6)
        
    return data*conversion_factor
    



### Making a common Dataframe of 30 min data
"""
Largest timeseries is from "PAR" <br>
1997 - 2558.00 - 1 - 30 to 2012 - 8401.98 - 366 - 2400 <br>
Using these Index to fill in the rest <br>
FillValue = -6999.0
"""

# list of variables
keys_vars = list(dict_dfs_common.keys())

# Making a copy of the dataframe with most data
df_all_vars_30m = dict_dfs_common["PAR"].copy(deep=True)

# Dropping JDT column since it is not same across vars based on Dates and time, due to which I get NaNs during merge
df_all_vars_30m = df_all_vars_30m.drop("JDT", axis=1)

for k in keys_vars:
    df_all_vars_30m = pd.merge(
        df_all_vars_30m, dict_dfs_common[k].drop("JDT", axis=1), how="left"
    )


fill_value = -6999.0
df_all_vars_30m["Rn"] = df_all_vars_30m["Rn"].replace(fill_value, np.nan)

# Fill NaNs with -6999.0
df_all_vars_30m_FV = df_all_vars_30m.fillna(fill_value)

# Adding a datetime column
# Convert DOY and Time to timedelta

# Convert 'Time' values to HH:MM format
time_str = df_all_vars_30m["Time"].astype(str)
time_str = time_str.str.zfill(4)  # Ensure all times are 4 digits

# Making Sure time_str in has integer values only
# time_int = time_str.astype(float).astype(int)
time_int = np.asarray(np.asarray(time_str, dtype=float), dtype=int)

# Extract hours and minutes
hours = np.asarray(time_int // 100, dtype=int)
minutes = np.asarray(time_int % 100, dtype=int)

# Calculate the total minutes
total_minutes = hours * 60 + minutes

df_all_vars_30m["Date"] = pd.to_timedelta(
    df_all_vars_30m["DOY"] - 1, unit="D"
) + pd.to_timedelta(total_minutes, unit="m")

# Add Year to the Date
df_all_vars_30m["Date"] = (
    pd.to_datetime(df_all_vars_30m["Year"].astype(float).astype(int).astype(str))
    + df_all_vars_30m["Date"]
)
df_all_vars_30m_FV["Date"] = df_all_vars_30m["Date"]

df_all_vars_30m_FV.to_csv(
    f"{paths['Save_Processed']}Processed_Duke_Met_Data_All_Vars_30m_FV.csv"
)
df_all_vars_30m.to_csv(
    f"{paths['Save_Processed']}Processed_Duke_Met_Data_All_Vars_30m.csv"
)
