#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 17:39:05 2023

@author: paoloscaccia
"""
import pandas as pd
import numpy  as np
from datetime import datetime
from influxdb import DataFrameClient
import sys, re
import warnings
warnings.simplefilter("ignore")

DATE_MASK = r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}'


def parser():
    import argparse
    
    parser=argparse.ArgumentParser()
    parser.add_argument('-o','--outdir',type=str,default='./',
                        help='Output Directory')  
    parser.add_argument('--start',type=str,default='1970-01-01T00:00:01',
                        help='Starting date for the query')  
    parser.add_argument('--end',type=str,default='2100-01-01T00:00:01',
                        help='Ending date for the query')  
    parser.add_argument('--lon_min',type=float,default=-30,
                        help='Minimum Longitude (default: -30°)')  
    parser.add_argument('--lon_max',type=float,default=70,
                        help='Maximum Longitude (default: 70°)')  
    parser.add_argument('--lat_min',type=float,default=-30,
                        help='Minimum Latitude (default: 33°)')  
    parser.add_argument('--lat_max',type=float,default=81,
                        help='Maximum Latitude (default: 81°)')  
    parser.add_argument('--domain',type=str,default="XXXX",
                        help='Custom Domain Tag')  
    parser.add_argument('-p','--products',
                        type=str,nargs='+',
                        required=False, default = ['FDEM','AFM','AFM_RS','FRP-PIX'],
                        help='Product to be downloaded from the remote DB')
    
    return parser.parse_args()

def print_config(args):

    print()    
    print("\n"+"#"*16)
    print("#"+" "*4+"QUERY"+" "*5+"#")
    print("#"*16)
    print(f"  Products: {', '.join(args.products)}")
    print(f"  Start:    {args.start}")
    print(f"  End:      {args.end}")
    print(f"  Lon Min:  {args.lon_min}°")
    print(f"  Lon Max:  {args.lon_max}°")
    print(f"  Lat Min:  {args.lat_min}°")
    print(f"  Lat Max:  {args.lat_max}°")
    print()

    return


def read_database(args):
    import xarray as xr
    list_of_products = '"' + '\", \"'.join(args.products) + '"'
    
    query_str = f'SELECT * FROM {list_of_products} WHERE '\
                f"( time<=\'{args.end}\' AND time>=\'{args.start}\' ) AND "\
                f'( \"lon\" >= {args.lon_min} AND \"lon\" <= {args.lon_max} ) AND '\
                f'( \"lat\" >= {args.lat_min} AND \"lat\" <= {args.lat_max} ) '
    dataframe_dict = client.query(query_str, epoch = 'ns', raise_errors = False)
    print('\nFire Products found in database: ')
    for product in dataframe_dict:
        print('     ',product)
        
    
    # Build ncfile header
    domain = 'EURO' if [args.lat_min, args.lat_max, args.lon_min, args.lon_max] == [-30, 81,-30, 70] else args.domain
    if len(domain) < 4:
        domain += '_'*(4 - len(domain) )
    else:
        domain = domain[:4]
    ncfile_header = f'_{domain}_{args.start.replace(" ","T").replace(":","").replace("-","")}_{args.end.replace(" ","T").replace(":","").replace("-","")}'

    print("\nCreating netcdf files...")    
    for product in dataframe_dict:
        
        # Read DB as Pandas DataFrame
        df = dataframe_dict[product]
        df['time'] = pd.to_datetime(df.index, unit='s')
        df = df.set_index(['time'])
        columns = [x for x in df.columns if x != 'type']
        df = df[columns]

        # Convert to Xarray
        if 'row_idx' in df:
            df = df.drop(columns=['row_idx','col_idx'])

        xr_data = df.to_xarray()
        xr_data['time'] = pd.DatetimeIndex(xr_data['time'].values)
        xr_data['time'].attrs['long_name'] = 'Time [UTC]'
        
        if 'FRP' in product:
            good_vars = ['time','FRP','NDVI','Risk','lat','lon']
            bad_vars = [x for x in xr_data.keys() if x not in good_vars]
            xr_data = xr_data.drop(labels=bad_vars)
            
        xr_data['lon'].attrs['long_name'] = 'longitude'
        xr_data['lon'].attrs['units'] = 'degree_east'
        xr_data['lat'].attrs['long_name'] = 'latitude'
        xr_data['lat'].attrs['units'] = 'degree_north'
        xr_data['NDVI'] = xr_data['NDVI'].astype(np.int8)
        xr_data['Risk'] = xr_data['Risk'].astype(np.int8)
        try:
            xr_data['category'].attrs['description'] = '1: Probable 2: Likely'
            xr_data['category'].attrs['long_name'] = 'Fire Certainty Class'
            xr_data['certainty'].attrs['long_name'] = 'Fire Certainty'
        except:
            pass
        try:
            xr_data['FRP'].attrs['long_name'] = 'Fire Radiative Power'
            xr_data['FRP'].attrs['units'] = 'MW'
        except:
            pass
        
        try:
            xr_data['size'].attrs['long_name'] = 'Pixel Size'
        except:
            pass
        xr_data['NDVI'].attrs['long_name'] = 'Normalized Difference Vegetation Index'
        xr_data['Risk'].attrs['long_name'] = 'Fire Risk Map (LSA-SAF)'
        
        # Save Netcdf4 (CF compliant)
        ncfile = args.outdir + '/' + str(product)+'_'*(6 - len(product)) + ncfile_header + '.nc'
        
        xr_data.to_netcdf(ncfile)
        print(f"Saved {ncfile}")
        
    del(df,xr_data)
    print()
    
    
    return

if __name__ == '__main__':
    
    # Parse Shell Argument
    args = parser()
    
    # Print Selected Configuration
    print_config(args)
    
    # Check Date Format
    if not re.match(DATE_MASK,args.start) or not re.match(DATE_MASK,args.end):
        sys.exit("Wrong time format for START/END date:  %Y-%m-%dT%H:%M:%S")
    else:
        args.start = args.start.replace('T',' ')
        args.end   = args.end.replace('T',' ')
        
    # Init InfluxDB client
    client = DataFrameClient('elena.hopto.org',
                             port     = 40022, 
                             database = "Fire_Detection")

    # Send Query to DB
    read_database(args)    
    
    client.close()
    
