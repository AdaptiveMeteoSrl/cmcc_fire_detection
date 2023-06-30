#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 15:42:56 2023

@author: avalletti
"""
import os, sys
import pygrib
import h5py
import shutil
import pandas as pd
import numpy  as np
from datetime import datetime, timedelta
from influxdb import DataFrameClient
from glob import glob
from satpy import Scene
import warnings
warnings.simplefilter("ignore") 
import matplotlib.pyplot as plt
import pyresample
import xarray as xr
from bs4 import BeautifulSoup
from scipy.spatial.distance import cdist
os.environ['XRIT_DECOMPRESS_PATH'] = "/home/cmcc/anaconda3/envs/cmcc/bin/xRITDecompress"

name_conventions_HDF = {"FDeM"          : "FDEM",
                        "FRP-PIXEL"     : "FRP-PIX",
                        "FRP-GRID"      : "FRP-GRID",
                        "FRM-"          : "FRM"}

name_conventions_grb = {"4.bin"    : "AFM",
                        "3.bin"    : "AFM_RS"}

name_conventions = {"S-LSA_-HDF5_LSASAF_MSG_"   : name_conventions_HDF , 
                    ".bin"                      : name_conventions_grb}

ZONES = ["EURO"] #, "NAFR", "SAFR", "SAME", "DISK"]
       

def parser():
    import argparse
    
    parser=argparse.ArgumentParser()
    parser.add_argument('-fd',
                    '--filedir',
                    type=str,
                    required = True,
                    help='Input files directory')
    parser.add_argument('--filelist',
                    type=str,
                    required=False,
                    default="./sent_filelist.txt",
                    help='Sent File List')
    parser.add_argument('-gpd',
                        '--geoprojdir',
                    type=str,
                    required=False,
                    default="/home/avalletti/CMCC/GEO Projection in HDF5/",
                    help='LSA-SAF Geo Projection in HDF5')
    return parser.parse_args()

def decompress_bzip2(path, filelist):
    import subprocess
    print("Decompressing bz2 files...")
    for i,f in enumerate(filelist):
        if ".bz2" in f:
            if not os.path.isfile(path+'/'+f.replace('.bz2','')):
                subprocess.run(["bzip2","-d","{}/{}".format(path,f)])
    listdir = os.listdir(path)
    return listdir

def set_HDF5_coords(geoprojdir):
    listdir = os.listdir(geoprojdir)
    
    if any(".bz2" in f for f in listdir):
        listdir = decompress_bzip2(geoprojdir, listdir)
    COORDS = {}
    for f in listdir:
        filename = geoprojdir+f
        try:
           coord = f.split("_")[3].upper()
           area  = f.split("_")[4].upper()
        
           x = h5py.File(filename, "r")
           COORDS[area+"_"+coord] = x[coord][:]
        except:
           print('Skipped HDF5 coord file ',f)
    return COORDS

def read_grb(filename,data_type,coords):
    # extreme coords:
    # Latitudes
    N = +81
    S = +33
    # # Longitudes
    E = -30
    W = +70

    grbs = pygrib.open(filename)
    grbs.seek(0)    
    grb = grbs[1]
    time = grb.analDate

    lons = -grb.longitudes
    #lons[lons>180] -= 360
    #ds = xr.open_dataset(filename, engine="cfgrib")
    #values = ds[ [ k for k in ds.keys() if 'p' in k ][0]].values.flatten()
    #lons   = ds['longitude'].values.flatten()
    #lats   = ds['latitude'].values.flatten()

    #lons = coords['lon'].values.flatten()
    #lons = np.mod(lons + 180.0, 360.0) - 180.0
    lats   = grb.latitudes
    values = grb.values.flatten()    
    #lats = coords['lat'].values.flatten()

    #a = np.unique(values)
    #if 1 in a or 2 in a:
        #filt = np.logical_and(values>0, values <3)
        # print(lats[filt],lons[filt],values[filt])
    
    lldf = pd.DataFrame({"time":time,"lat":lats,"lon":lons,"category":values})
    #lldf.drop(lldf[np.logical_or(lldf["lat"]>N,lldf["lat"]<S)].index, inplace=True)
    #lldf.drop(lldf[np.logical_or(lldf["lon"]>W,lldf["lon"]<E)].index,inplace=True)
    
    #filt = np.logical_and(lldf["category"]>0, lldf["category"] <3)
    #lldf.drop(lldf[np.logical_or(lldf["category"]==0,lldf["category"]>2)].index,inplace=True)
    
    indices = np.array(list(lldf.index))
    ncol  = grb.values.shape[1]
    icols = indices% ncol
    irows = indices//ncol
    #lons = lldf["lon"]
    #lons[lons<0] += 360
    #lldf["lon"] = lons
    #lldf["row_idx"] = irows
    #lldf["col_idx"] = icols
    lldf = lldf[ lldf['category'].isin([1.,2.])]
    if len(lldf)==0:
        print(datetime.now(),"No fire detected within the target domain")
        return None
    
    #df = pd.concat([df,lldf],ignore_index=True)
    print(datetime.now(),"Data added to output dataframe")

    return lldf

def read_HDF5(filename,COORDS):

    area      = filename.split("/")[-1].split("_")[5].upper()
    data_type = filename.split("/")[-1].split("_")[4].upper()
    
    if area not in ZONES and "GRID" not in data_type:
        return None
    
    time = filename.split("/")[-1].split("_")[-1]
    time = datetime.strptime(time,"%Y%m%d%H%M")
    
    x = h5py.File(filename, "r")

    
    if data_type == "FDEM":
        
        """
        
        https://nextcloud.lsasvcs.ipma.pt/s/zXJoeTf6HByE6RP?dir=undefined&path=%2FFDeM%2FATBD&openfile=27139
        [pag. 27] 
        
        FD&M product is provided in the HDF5 format as requested by the LSA-SAF system.
        This format allows defining a set of attributes that provide the relevant information.
        Table 8 shows description of dataset (named CF) that composes FD&M file (Table 8).
        
        TABLE 8:
       
          Class |  Description
        _________________________
            0   | Water
            1   | Land
            2   | Land with fire
            
        """
        prev_date = time - timedelta(seconds=15)
        key = "CF"
        categories = [2]
        lats = COORDS[area+"_LAT"]/100
        lons = COORDS[area+"_LON"]/10000
        data = x[key][:]
        
        lldf = pd.DataFrame({"lat":lats.flatten(),"lon":lons.flatten(),"category":data.flatten()})
        indices = np.array(list(lldf.index))
        ncol  = lats.shape[1]
        icols = indices% ncol
        irows = indices//ncol

        #lldf["row_idx"] = irows
        #lldf["col_idx"] = icols
        lldf = lldf[ lldf['category'].isin(categories)]

    elif "FRP-PIX" in data_type:

        if len(x)==0:
            return None
        """
        
        https://nextcloud.lsasvcs.ipma.pt/s/zXJoeTf6HByE6RP?dir=undefined&path=%2FFRPPIXEL_FRPGRID%2FPUM&openfile=27166
        [pag. 25-29] 
        
        from 'ProductList' file: FRP (range > 0, Scale factor: 10)
        
        from 'QualityProduct':  QUALITY FLAG 
       
        EXAMPLE QUALITY FLAG CATEGORY:
        ###########
        # QF == 1 #
        ###########
        NAME    : FRP_APL_FRP  
        VALUE   : 1 
        STATUS  : FRP Estimated 
        REASON  : Successful fire detection and FRP estimation.
            
        """
        lats = x["LATITUDE"][:]/100
        lons = x["LONGITUDE"][:]/100
        data = x["FRP"][:]/10
        
        lldf = pd.DataFrame({"time":time,"lat":lats.flatten(),"lon":lons.flatten(),"FRP":data.flatten()})
        
        
        KEYS = ['ACQTIME', 'ERR_ATM_TRANS', 'ERR_BACKGROUND', 
                'ERR_FRP_COEFF', 'ERR_RADIOMETRIC',
                'ERR_VERT_COMP', 'FIRE_CONFIDENCE', 
                'FRP_UNCERTAINTY', 'RAD_BCK', 'RAD_PIX', 
                'STD_BCK']
        
            
        path        = ("/").join(filename.split("/")[:-1])+"/"
        FRP_PIX_QC  = [qc for qc in os.listdir(path) if "Quality" in qc and "FRP-PIX" in qc and '.bz2' not in qc]
        qc = []
        for f in FRP_PIX_QC:
            if time.strftime("%Y%m%d%H%M") in f and area in f.upper():
                qc = h5py.File(path+'/'+f, "r")
                
                
        lldf["row_idx"] = x["REL_LINE"][:].flatten()
        lldf["col_idx"] = x["REL_PIXEL"][:].flatten()
        
        for k in KEYS:
            lldf[k] = x[k][:].flatten()

        prev_date = time - timedelta(minutes=15)
   
        
    elif "FRP-GRID" in data_type:
        
        if len(x)==0:
            return None
        
        """
        
        https://nextcloud.lsasvcs.ipma.pt/s/zXJoeTf6HByE6RP?dir=undefined&path=%2FFRPPIXEL_FRPGRID%2FPUM&openfile=27166
        [pag. 30-...] 
        Ref. SAF/LAND/KCL/PUM_FRP/2.2

        Every hour the FRP product processing chain generates one external FRP-GRID product output file, 
        according to the following name convention:
            
                    HDF5_LSASAF_MSG_FRP_Frp_Grid_Global_YYYYMMDDSSEE
                    
        where YYYY, MM, DD, SS and EE respectively denote the year, the month, the day,
        the Starting and End hour of the analysed period.
            
        """
        
        lats = x["LATITUDE"][:]/100
        lons = x["LONGITUDE"][:]/100
        data = x["GFRP"][:]
        
        KEYS = ['ATMTRANS', 'BURNTSURF', 'GFRP_CLOUD_CORR', 
                'GFRP_ERROR', 'GFRP_ERR_FRP', 'GFRP_QI', 
                'GFRP_RANGE', 'GRIDPIX', 'NUMFIRES', 
                'NUMIMG', 'NUMSATFIRES']
        
        lldf = pd.DataFrame({"time":time,"lat":lats.flatten(),"lon":lons.flatten(),"category":data.flatten()})
        for k in KEYS:
            lldf[k] = x[k][:].flatten()
        # lldf["row_idx"] = x["REL_LINE"].flatten()
        # lldf["col_idx"] = x["REL_PIXEL"].flatten()
        lldf.drop(lldf[np.logical_or(lldf["category"]==0,lldf["category"]==32767)].index,inplace=True)

        prev_date = time - timedelta(minutes=60)
        
    elif "FRM" in data_type:
        if len(x)==0:
            return None
        
        """
        
        https://nextcloud.lsasvcs.ipma.pt/s/zXJoeTf6HByE6RP?dir=undefined&path=%2FFRM%2FPUM&openfile=27152
        pag. 31
        
        Output data are coded in HDF5 format. The HDF5 files in LSA SAF system have the following structure:
        A common set of attributes for all kind of data, containing general information
        about the data (including metadata compliant with U-MARF requirements);
        A dataset for the parameter values;
        Additional datasets for metadata (e.g., quality flags).
        
        The FRM product, estimated once a day, is available in an HDF5 file containing ten
        datasets:
             FWI;
             DSR;
             BUI;
             ISI;
             FFMC;
             DMC;
             DC;
             Risk;
             Table_Ref;
             Quality Control information;
        Namely, the three fuel moisture codes (FFMC, DMC, DC), the three fire behaviour
        indices (ISI, BUI, FWI), the Daily Severity Rating (DSR), the Risk (low, moderate and
        high, 10, 20 and 30, respectively), the Table_Ref (10, 20, 30, 41, 42, 43 and 50 for
        shrub, tree cover broad-leaved, tree cover needle-leaved, cultivated and managed areas
        for Iberian Peninsula and France, cultivated and managed areas for the remaining
        regions and all other types of vegetation, respectively) and a flag indicating processed
        pixels (pixels in land) and non-processed pixels (water pixels or pixels out of Europe or
        pixels out of MSG disk).            
        """        
        lats = COORDS[area+"_LAT"]/100
        lons = COORDS[area+"_LON"]/10000
        
        lldf = pd.DataFrame({"time":time,"lat":lats.flatten(),"lon":lons.flatten()})
        lldf['Risk'] = x['Risk'][:].flatten()
        #lldf = lldf[ lldf['Risk'].isin(categories)]
        lldf.drop(lldf[lldf["Risk"]<0].index,inplace=True)
        prev_date = time - timedelta(days=1)


    # Filter out points outside the target domain (defualt: EU)
    lldf = lonlat_filter(lldf)
 
    n_points = lldf.index.size
    if n_points !=0 :

         # Assign a tag time to each pixel 
         n_sec  = (time - prev_date).total_seconds()/n_points 
         lldf['time'] = [ time - timedelta(seconds=i*n_sec) for i in range(n_points) ]
         lldf = lldf.set_index('time')

         print(datetime.now(),"Data added to output dataframe")

         return lldf

    else:
        print(datetime.now(),"No fire detected within the target domain")
        return None


def debug_plot(df,title):
    import cartopy.crs as ccrs
    plt.clf()
    print('DEBUG')
    data = df.to_xarray()
    globe = ccrs.Globe(semimajor_axis=6378137, flattening=1/298.257223563)
    myccrs=ccrs.PlateCarree( globe=globe)           
    plt.figure(figsize=(15, 15))
    ax = plt.axes(projection=myccrs)
    gc=ax.coastlines(resolution='10m',linewidth=.75,color='orange')
    #gl=ax.gridlines(linewidth=.5,color='black')
    gl=ax.gridlines(draw_labels=True, linewidth=.75, color='black')
    gl.bottom_labels = False
    gl.left_labels = False
    lon = data['lon']
    lat = data['lat']
    xmin, xmax = lon.min(), lon.max()
    ymin, ymax = lat.min(), lat.max()

    ax.set_extent([ -90, 90, -60,60],crs=ccrs.PlateCarree())
    ax.scatter(lon,lat,marker='o',c='red')
    file = './'+title+'.png'
    plt.savefig(file)
    print("Saved ",file)
   
    return

def HRIT_img(scn,composite,date,save = True,outdir = "/."):
    from satpy.writers import to_image

    scn.load([composite])
    img = to_image(scn[composite])

    img.stretch("linear") 
    img.gamma(2)

    outfile = f'{outdir}/{date[:-4]}/{composite}_{date}.png'
    img.save(outfile)
    print(datetime.now(),"Saved satellite image ",outfile)

    return

def read_HRIT(file, path):
    outdir = '/'.join(path.split('/')[:-2]) + '/processed'
    if not os.path.isdir(outdir):
       os.system(f"mkdir -p {outdir}")

    date = file.split('-')[-2]
    header  = file[:17]

    channels  = ['HRV','VIS006','VIS008','IR_016','IR_108','IR_039','IR_120','IR_087','IR_134','WV_062','WV_073']
    epilogue  = file.replace('-PRO__','-EPI__')
    cloudmask = f"W_XX-EUMETSAT-Darmstadt,SING+LEV+SAT,MET10+CLM_C_EUMG_{date}00_3.bin"

    if not os.path.isfile(path + '/' + epilogue):
           print(f"Missing EPILOGUE segment for file {file}")
           return False
    if not os.path.isfile(path + '/' + cloudmask):
           print(f"Missing Cloud Mask for file {file}")
           return False
    with pygrib.open(path + '/' + cloudmask) as grbfile:
          grbfile.seek(0)
          grb = grbfile[1]
          num_cm = grb.values
          cm = np.zeros_like(num_cm,dtype=bool)
          cm[ num_cm >= 2 ] = True
          cm = xr.DataArray(cm,dims=('lon','lat'))
    os.system(f"mkdir -p {outdir}/{date[:-4]}")
    related_segments = [  f for f in glob(path + '/' + header + '*' + date +'*') if f.split('-')[-4].strip('_') in channels]
    related_segments += [ path + '/' + file, path + '/' + epilogue]   
    scn = Scene(filenames = related_segments, reader = 'seviri_l1b_hrit')
    scn.load(channels)
    scn.compute()
    width, height = scn['VIS006'].shape
    area = pyresample.geometry.create_area_def( area_id='mygrid',projection={'proj': 'latlong', 'lon_0': 0},width=width, height=height)
    scn = scn.resample(area)
    
    # SHOW AND SAVE IMAGES
    """
    HRIT_img(scn,'natural_color',date,True,outdir)
    HRIT_img(scn,'cloudtop',date,True,outdir)
    HRIT_img(scn,'convection',date,True,outdir)
    HRIT_img(scn,'night_fog',date,True,outdir)
    HRIT_img(scn,'snow',date,True,outdir)
    HRIT_img(scn,'day_microphysics',date,True,outdir)
    HRIT_img(scn,'fog',date,True,outdir)
    HRIT_img(scn,'realistic_colors',date,True,outdir)
    HRIT_img(scn,'overview',date,True,outdir)
    HRIT_img(scn,"natural_enh",date,True,outdir)
    """
    outfile = f'{outdir}/{date[:-4]}/ndvi_{date}.png'    
    scn['ndvi'] = (scn['VIS008'] - scn['VIS006'])/(scn['VIS006'] + scn['VIS008'])
    if not os.path.isfile(outfile):

        plt.imshow(scn['ndvi'][:-1:,:])
        plt.title('NDVI '+date)
        plt.colorbar()
        plt.xticks([])
        plt.yticks([])
        plt.savefig(outfile,dpi=300)
        print(datetime.now(),f'Saved NDVI image {outfile}')
        plt.clf()

    scn['ndvi'] = scn['ndvi'].rename( {'y': 'lat', 'x': 'lon'})
    scn['ndvi'] = scn['ndvi'].where( ~ cm, drop=False)
    scn['ndvi'] = scn['ndvi'].where( scn['ndvi'] >= 0.0, drop=False)

    #scn['ndvi'] = scn['ndvi'].to_masked_array(copy=True)
    scn['cm'] = cm

    ncfile =  f'{outdir}/{date[:-4]}/ndvi_{date}.nc'
    if not os.path.isfile(ncfile):
        scn.save_datasets(datasets=['ndvi','cm'],
                       filename=ncfile,
                       flatten_attrs=True,
                       exclude_attrs=['raw_metadata'],
                       include_lonlats=False)
        print(datetime.now(),f'Saved NetCDF file {ncfile}')
 
    return True

def read_CAP(file, path, client):
    filetype = 'AFM_RS' if '3' in os.path.basename(file).split('_')[-1] else 'AFM' 

    date = os.path.basename(file).split('_')[-2]
    date = datetime.strptime(date,'%Y%m%d%H%M%S')
    prev_date = date - timedelta(minutes=5) if filetype == 'AFM_RS' else date + timedelta(minutes=15)
    with open(file, 'r') as f:
         data = f.read()
    bs_data = BeautifulSoup(data, 'xml')
    categories  = []
    status_list = []
    lats        = []
    lons        = []
    sizes       = []

    cases = bs_data.find_all('info')
    for case in cases:
       
       status = str(case.find('certainty').text)
       category = 2.0 if 'Like' in status else 1.0
       circles = case.find_all('circle')
       for circle in circles:
             lat,lon = circle.text.split()[0].split(',')
             size = circle.text.split()[-1]
             lats.append(float(lat))
             lons.append(float(lon))
             sizes.append(float(size))
             categories.append(float(category))
             status_list.append(status)
    n_points = len(lats)
    n_sec = (date - prev_date).total_seconds()/n_points if n_points != 0 else 0
    times = np.array([ date - timedelta(seconds=i*n_sec) for i in range(n_points)  ])

    """
    print('times: ',len(times))
    print('lon: ',len(lons))
    print('lat: ',len(lats))
    print('size: ',len(sizes))
    print('certainty: ',len(status_list))
    print(file, date, next_date)
    print('last_time: ',times[-1])
    """
    df = pd.DataFrame({'lat'       : np.array(lats,dtype=float),
                       'lon'       : np.array(lons,dtype=float),
                       'category'  : np.array(categories,dtype=float),
                       'size'      : np.array(sizes,dtype=float),
                       'certainty' : np.array(status_list,dtype=str),
                       'time'      : times })
    df = df.set_index('time')
    df.title = filetype
    df = lonlat_filter(df)

    if df.index.size == 0:
        print(datetime.now(), "No AFM detected within the target domain in ",os.path.basename(file))
        return None
    else:
       print(datetime.now(), "Data added to output dataframe")
       return df 


def move_file_to_sent_dir(file,path,newpath):
    
    shutil.move(path+'/'+file, newpath)    
    
def lonlat_filter(df, zone = 'EU'):
    if zone == 'EU':
         # Latitudes
         N = +81
         S = +33

         # # Longitudes
         E = -30
         W = +70
    else:
         sys.exit(f"Zone {zone} not implemented")
    df.drop(df[np.logical_or(df["lat"]>N,df["lat"]<S)].index, inplace=True)
    df.drop(df[np.logical_or(df["lon"]>W,df["lon"]<E)].index,inplace=True)

    return df

def retrieve_FRM_and_NDVI(df, FRM, NDVI):
    
    # Read FRM
    target_date = df.index[0].date().strftime('%Y%m%d')
    if target_date in FRM:
        FRM_df = FRM[target_date]
        distances = cdist(df[['lon','lat']].values,FRM_df[['lon','lat']].values)
        closest_points = distances.argmin(axis=1)
        for var in FRM_df.columns:
            if var not in df.columns:
                df[var] = FRM_df.iloc[closest_points][var].values
    else:
        print(datetime.now(),"Associated Fire Risk Map (FRM) file not found")
        FRM_df = FRM[list(FRM.keys())[0]]
        for var in FRM_df.columns:
            if var not in df.columns:
                df[var] = [np.nan for i in df.index]

    # Read NDVI
    target_date = (df.index[0] - timedelta(days=1)).date().strftime('%Y%m%d')
    if target_date in NDVI:
        NDVI_df = NDVI[target_date]
        distances = cdist(df[['lon','lat']].values,NDVI_df[['lon','lat']].values)
        closest_points = distances.argmin(axis=1)
        df['NDVI'] = NDVI_df.iloc[closest_points]['ndvi'].values
    else:
        print(datetime.now(),"Daily NDVI file not found for ",target_date)
        df['NDVI'] = [ np.nan for i in df.index]

    return df

        
    
    return

def send_data(path,listdir, log_file, client, FRM, NDVI, HDF5COORDS = None):
    listdir.sort()
    filelist = np.loadtxt(log_file,dtype=str)
    count = 0
    if len(listdir)==0:
        return "EMPTY FILES LIST"
    # df = pd.DataFrame(columns = ["time","lat","lon","category","row_idx","col_idx"])
    listdir.sort()
    
    # Read FRP-PIXEL or FDEM
    if listdir[0][:23] in name_conventions.keys():
        
        data_type = listdir[0].split("_")[4]
        if "FRP-PIXEL" in data_type:
            data_type = "FRP-PIXEL"  
        elif "FDeM" in data_type:
            data_type = "FDeM"  
        else:
            sys.exit("Wrong HDF5 product!")  
            
        if HDF5COORDS is None:
            print(datetime.now(),f"DATA TYPE: {data_type} \nTHERE AREN'T COORDS MATRICES")
            return "NO HDF COORDINATES DIRECTORY"
        
        for f in listdir:
            if "QualityProduct" in f:
                continue
            if f not in filelist:
                print(datetime.now(),"Reading file:",f)       
                #try:
                # df = pd.DataFrame(columns = ["time","lat","lon","category","row_idx","col_idx"])
                df = read_HDF5(path+'/'+f,HDF5COORDS)

                if not (df is None) and df.size > 0:

                        # Find corresponding FRM and NDVI
                        df = retrieve_FRM_and_NDVI(df,FRM,NDVI)

                        # Updte InfluxDB
                        client.write_points( df,  name_conventions_HDF[data_type], 
                                                  batch_size = 1000,
                                                 tags = { 'type' : name_conventions_HDF[data_type]})
                        print(datetime.now(),"Updated Database")

                # Update list with read files
                with open(log_file,'a') as outfile:
                           outfile.write(os.path.basename(f)+'\n')

                #except Exception as error:
                #   print(f"Skipped file {f}: {error}")
            # move_file_to_sent_dir(f, path, newpath)
            
    # Read AFM/AFM_RS
    elif '.txt' in listdir[0] and 'FIRC' in listdir[0]:
        data_type = listdir[0].split("_")[-1]

        for f in listdir:
            if f not in filelist:
                print(datetime.now(),"Reading file:",f)
                #try:
                df = read_CAP(path+'/'+f, path, client)
                if not (df is None) and df.size > 0:
                    
                    # Find corresponding FRM and NDVI
                    df = retrieve_FRM_and_NDVI(df,FRM,NDVI)

                    # Update InfluxDB
                    client.write_points( df,   df.title, batch_size = 1000, tags = { 'type' : df.title})
                    print(datetime.now(),"Updated Database")

                # Update list with read files
                with open(log_file,'a') as outfile:
                          outfile.write(os.path.basename(f)+'\n')

                #except Exception as error:
                #    print(f"Skipped file {f}: {error}")

    elif listdir[0][:5] == 'H-000':
        print(datetime.now(),"Reading HRIT...")
        for f in listdir:
          if f not in filelist:
             print(datetime.now(),"Reading file:",f)
             try:
                flag = read_HRIT(f,path)
                
                # Update list with read files
                with open(log_file,'a') as outfile:
                       outfile.write(os.path.basename(f)+'\n')

             except Exception as error:
                 print(datetime.now(),f'Skipped file {f}: {error}')
    else:
        sys.exit(f"Wrong File Format in {listdir[0]}")
    
    
    del(filelist)        
    return f"SENT {count} FILES"
          
def read_all_FRM(FRM_list, path, HDF5COORDS):
    FRM_list.sort()
    print(datetime.now(),"Reading Fire Risk Map files...")

    out_dict = {}
    for f in FRM_list:
        if "QualityProduct" in f:
            continue
        print(datetime.now(),"Reading file:",f)       
        #try:
        # df = pd.DataFrame(columns = ["time","lat","lon","category","row_idx","col_idx"])
        df = read_HDF5(path+'/'+f,HDF5COORDS)

        if not (df is None):
                out_dict[df.index[0].date().strftime('%Y%m%d')] = df
    
    return out_dict
    
def read_all_NDVI(processed_days, path ):
    processed_days.sort()
    print(datetime.now(),"Reading Daily NDVI files...")

    out_dict = {}
    for day in processed_days:
          print(datetime.now(),"Scanning processed folder ",day)
          ndvi_file = day + '/daily_ndvi.nc'
          if os.path.isfile(ndvi_file):
               try:
                   data = xr.open_dataset(ndvi_file,engine='netcdf4')
                   day_str = os.path.basename(day) if day[-1] != '/' else os.path.basename(day[:-1])

                   lon, lat = np.meshgrid( data['lon'].values, data['lat'].values)
                   df = pd.DataFrame({'lon':lon.flatten(), 'lat':lat.flatten(),'ndvi':data['ndvi'].values.flatten()})
                   out_dict[day_str] = df                   
                   del(data)
               except Exception as err:
                   print(datetime.now(),f"Skipped file {ndvi_file}: ",err)
          else:
               print(datetime.now(), 'Daily NDVI file not found')
               continue
    return out_dict

if __name__ == "__main__" :

    args = parser()
    path            = args.filedir
    if path[-1] != '/':
        path += '/'
    sentfilelist    = args.filelist
    geoprojHDF5dir  = args.geoprojdir
    
    listdir = os.listdir(path)
    listdir.sort()

    listbz2 = [f for f in listdir if ".bz2" in f]
    if len(listbz2)!=0:
        listdir = decompress_bzip2(path, listdir)
    del(listbz2)

    # Read File list for each product
    FRM       = [f for f in listdir if "FRM" in f and '.bz2' not in f]
    HRIT      = [f for f in listdir if "-PRO______" in f]
    AFM       = [    f for f in listdir if "FIRC" in f and f.split('.')[-1] == 'txt' ]
    FDEM      = [f for f in listdir if "FDeM" in f and ( "Quality" not in f and '.bz2' not in f)]
    FRP_PIX   = [f for f in listdir if ("FRP-PIXEL" in f) and ('.bz2' not in f)]
    PROC_DAYS = glob(path.replace('default','processed/*'))

    # Read auxiliary data for HDF5 products
    COORDS = set_HDF5_coords(geoprojHDF5dir)

    # Read all available Fire Risk Map and Daily NDVI
    FRM_data = read_all_FRM(FRM, path, COORDS)
    NDVI_data = read_all_NDVI(PROC_DAYS, path)

    client = DataFrameClient('elena.hopto.org',port=40022, database = "Fire_Detection")
    
    log_AFM = send_data(path, AFM,sentfilelist, client, FRM_data, NDVI_data) 
    log_FDEM = send_data(path, FDEM, sentfilelist, client, FRM_data, NDVI_data, COORDS)
    log_FRP = send_data(path, FRP_PIX, sentfilelist, client, FRM_data, NDVI_data, COORDS)
    log_HRIT = send_data(path, HRIT, sentfilelist, client, FRM_data, NDVI_data)
    
    client.close()
    del(FRM_data, NDVI_data)
    
    
