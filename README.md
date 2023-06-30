# CMCC Fire Detection
by AdaptiveMeteo S.r.l.

## 1. Install dependencies

> conda env create -f env/cmcc_env_linux.yml
or
> conda env create -f env/cmcc_env_win.yml

## 2. Download data

> python read_database.py

**Options:**
- --lon_min/--lon_max : Minimum/Maximum Longitude [Default: -30°/70°]
- --lat_min/--lat_max : Minimum/Maximum Latitude [Default: -30°/81°]
- --start/--end       : Temporal extention of the query [Format: YYYY-mm-ddTHH:MM:SS]
- --domain_tag        : Tag for custom domains
- -p/--products       : List of fire products to be downloaded [Default: AFM, AFM_RS, FDEM, FRP-PIX]

## 3. Visualize Data

[Grafana Client](http://elena.hopto.org:50022)


**User:**  CMCC_Guest 
**Psswd:** W3lc0m3
