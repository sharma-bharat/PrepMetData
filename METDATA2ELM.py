## Preparing the Met data for ELM
# by Bharat Sharma

# Source : https://github.com/dmricciuto/OLMT/blob/master/metdata_tools/site/data_to_elmbypass.py
# importing libraries
import xarray as xr
import glob
from datetime import datetime
import cftime
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
paths["Raleigh_Airport_Met"] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/"
paths["ERA5_Met"] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/ERA5_Duke_Met/"
paths["NARR_Met"] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/NARR_Met/"
paths["AmeriFlux"] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/AmeriFlux/"
paths[
    "Save_Processed"
] = "/Users/ud4/Documents/FACEMDS/MET_Data_Processing/Oren_2022_Met_Data_processed/"


df_h = pd.read_csv(f"{paths['Save_Processed']}DUKE_forcing_h.csv", index_col=0)
# Drop the first row (index 0) i.e. of units
df_h = df_h.drop(0)

cols_data = ["YEAR", "DOY", "HRMIN", "Tair", "RH", "Wind", "PSurf", "SWdown", "Rainf"]

df_h_elm = df_h[cols_data]

cols_4_elm = ["YEAR", "DOY", "HRMIN", "TA", "RH", "WS", "PA", "PPFD_OUT", "H2O"]

# Rename the columns
df_h_elm.columns = cols_4_elm
df_h_elm
# df_h


def retun_n_days(yr):
    if int(yr) % 4 == 0:
        days = 366
    else:
        days = 365

    return days


# saving dataevery year
fname_template = "SITE-MetData-YYYY.csv"
year_list = [f"{i:d}" for i in range(1996, 2013)]
site = "US-DUK"

for yr in year_list:
    yr_filter = df_h_elm.YEAR == float(yr)
    print(f"Total count of rows {yr} = {sum(yr_filter)} == {retun_n_days(yr)*24*2}")
    df_h_elm_y = df_h_elm[yr_filter]
    # Removing the Feb 29th of leap years
    if yr in ["1996", "2000", "2004", "2008", "2012"]:
        df_h_elm_y = df_h_elm_y.drop(df_h_elm_y[(df_h_elm_y.DOY == 60)].index, axis=0)
    df_h_elm_y = df_h_elm_y.reset_index(drop=True)
    fname = fname_template.replace("SITE", site).replace("YYYY", str(yr))
    df_h_elm_y[["YEAR", "DOY", "HRMIN"]] = df_h_elm_y[["YEAR", "DOY", "HRMIN"]].astype(
        int
    )

    print(df_h_elm_y.shape)
    df_h_elm_y.to_csv(f'{paths ["Save_Processed"]}ELM_MET/raw_{fname}', index=False)


# Running the Code to convert Data to ELM Met data¶

from netCDF4 import Dataset
import os, sys

sys.path.append("/Users/ud4/repos/GitHub/PrepMetData/in_data")
import gapfill
import write_elm_met

# ------- user input -------------
site = "US-DUK"
start_year = 1996
end_year = 2012
time_offset = -5  # Standard time offset from UTC (e.g. EST is -5)
npd = 48  # number of time steps per day (48 = half hourly)
mylon = 280.9058  # site longitude (0 to 360)
mylat = 35.9782  # site latitude
measurement_height = 2  # tower height (m)
fname_template = f"{paths ['Save_Processed']}ELM_MET/raw_US-DUK-MetData-YYYY.csv"
calc_flds = True  # use T and RH to comput FLDS (use if data missing or sparse)
leapdays = False  # input data has leap days (to be removed for ELM)
outdir = (
    f"{paths ['Save_Processed']}/1x1pt_" + site + "/"
)  # Desired directory for ELM met inputs


metdata = {}
# outvars   - met variables used as ELM inputs
# invars    - corresponding variables to be read from input file
# conv_add  - offset for converting units (e.g. C to K)
# conv_mult - multiplier for converting units (e.g. hPa to Pa, PAR to FSDS)
# valid_min - minimum acceptable value for this variable (set as NaN outside range)
# valid_max - maximum acceptable value for this variable (set as NaN outside range)

# Note - FLDS not included here (calculated)
outvars = ["TBOT", "RH", "WIND", "PSRF", "FSDS", "PRECTmms"]
invars = ["TA", "RH", "WS", "PA", "PPFD_OUT", "H2O\n"]  # matching header of input file
conv_add = [0, 0, 0, 0, 0, 0]
conv_mult = [1, 1, 1, 1, 1, 1]
valid_min = [180.00, 0, 0, 8e4, 0, 0]
valid_max = [350.00, 100.0, 80, 1.5e5, 2500, 15]

# ELM Variable names and units
# TBOT:     Air temperature at measurement (tower) height (K)
# RH:       Relative humidity at measurment height (%)
# WIND:     Wind speeed at measurement height (m/s)
# PSRF:     air pressure at surface  (Pa)
# FSDS:     Incoming Shortwave radiation  (W/m2)
# FLDS:     Incoming Longwave radiation   (W/m2)
# PRECTmms: Precipitation       (kg/m2/s)

os.system("mkdir -p " + outdir)
for v in outvars:
    metdata[v] = []

# Load the data
for y in range(start_year, end_year + 1):
    if (y % 4) == 0:
        isleapyear = True
    filename = fname_template.replace("SITE", site).replace("YYYY", str(y))
    lnum = 0
    myfile = open(filename, "r")
    for s in myfile:
        if lnum == 0:
            header = s.split(",")
        else:
            # skip leap days
            if not leapdays or (
                not isleapyear or (isleapyear and (lnum - 1) / npd != 59)
            ):
                data = s.split(",")
                for v in range(0, len(invars)):
                    for h in range(0, len(header)):
                        if header[h] == invars[v]:
                            # print(header)
                            try:
                                val = float(data[h]) * conv_mult[v] + conv_add[v]
                                if val >= valid_min[v] and val <= valid_max[v]:
                                    metdata[outvars[v]].append(val)
                                else:
                                    metdata[outvars[v]].append(np.NaN)
                            except:
                                metdata[outvars[v]].append(np.NaN)
        lnum = lnum + 1

# Fill missing values with diurnal mean
for key in metdata:
    print(key)
    gapfill.diurnal_mean(metdata[key], npd=npd)


out_fname = outdir + f"/all_hourly{start_year}_{end_year}.nc"
write_elm_met.bypass_format(
    out_fname,
    metdata,
    mylat,
    mylon,
    start_year,
    end_year,
    edge=0.1,
    time_offset=time_offset,
    calc_qbot=False,
    calc_lw=calc_flds,
    zbot=measurement_height,
)


print(f"file saved at {out_fname}")
