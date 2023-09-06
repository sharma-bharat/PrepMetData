# PrepMetData
This repository is created to help one with tools to process meterological data into continuous timeseries in formats that are in land model readable format.

# Processing of Duke Met Data
by Bharat Sharma and Anthony Walker 
sharmabd@ornl.gov
=============================

## Fixing the duplicate time index
There are some duplicate/incorrect time information in the original data in the year 2007 and 2009 across all variables. <br>
We wrote a code to fix the files : `FixingDupicateDuke.py` <br>
Usage example:
`python FixingDupicateDuke.py -file_ci /Users/ud4/repos/GitHub/FATESFACE/Jupyter_Notebooks/DuplicateDukeDataCorrectIndexOnly.txt -path_data /Users/ud4/Documents/FACEMDS/MET_Data_Processing/Oren_2022_DUKE_Met/data/ -replace_file yes` <br>
This will replace tge gap filled files with correct values <br>
You will also need to download the file `DuplicateDukeDataCorrectIndexOnly.txt` <br>

## Sub-hourly

### Variable "AT" 
- Air Temperature (Degree Celsius)
- Renamed to 'Tair'
- Calculated the mean of the following plots every time step for the from 1997 to 2012
  - 'R1uat', 'R2uat', 'R3uat', 'R4uat','R5uat', 'R6uat', 'R7uat', 'R8uat'

### Variable "Precip"
- Precipitation (mm)
- Renamed to 'Rainf'
- Calculated the mean of the following plots every time step for the from 1997 to 2012
  - 'FACE.PO'

### Variable "RH"
- Relative Humidity
- Renamed to 'RH'
- Calculated the mean of the following plots every time step for the from 1997 to 2012
  - 'R1urh', 'R2urh', 'R3urh', 'R4urh','R5urh', 'R6urh', 'R7urh', 'R8urh'

### SM:SM 
- Soil moisture integrates measurements from 0 to 30cm depth 
- Renamed to 'SM'
- Calculated the mean of the following plots every time step for the from 1997 to 2012
  - 'R1urh', 'R2urh', 'R3urh', 'R4urh','R5urh', 'R6urh', 'R7urh', 'R8urh'

### SWP:SWP 
- Soil water potential
- Renamed to 'SWP'
- Calculated the mean of the following plots every time step for the from 2007 to 2012
  - 'R1swp', 'R2swp', 'R3swp', 'R4swp','R5swp', 'R6swp'

### Variable "SVP"
- Saturated Vapor Pressure (kPa)
- Renamed to 'SVP'
- Calculated the mean of the following plots every time step for the from 1997 to 2012
  - 'R1usvp', 'R2usvp', 'R3usvp', 'R4usvp','R5usvp', 'R6usvp', 'R7usvp', 'R8usvp'

### Variable "VPD"
- Vapor pressure deficit (kPa)
- Renamed to 'VPD'
- Calculated the mean of the following plots every time step for the from 1997 to 2012
  - 'R1uvpd', 'R2uvpd', 'R3uvpd', 'R4uvpd','R5uvpd', 'R6uvpd', 'R7uvpd', 'R8uvpd'


### Variable "SLT"
- Soil Temperature (Degree Celsius) at 15 cm Depth
- Renamed to 'SLT'
- Calculated the mean of the following plots every time step for the from 1997 to 2012
  - 'R2slt', 'R3slt', 'R4slt','R5slt', 'R6slt'
  - Other plots had varied depths over time hence, ignored them for this calculation


### Variable "PAR"
- Photosynthetically active radiation (umol/m^2\*s) 
- Renamed to 'PAR'
- Calculated the mean of the following plots every time step for the from 1997 to 2012
  - 'PAR'
  - 'PAR' had values from only one plot for 1997-2007 and from two plots for 2008-2012, which were averaged 

### Variable "Rn"
- Net radiation (W/m^2) 
- Renamed to 'Rn'
- Calculated the mean of the following plots every time step for the from 1997 to 2012
  - 'Rn'
  - 'Rn' had values from only one plot for 1997-2007 and from two plots for 2008-2012, which were averaged 
