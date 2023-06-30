#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 16:11:25 2023

@author: paoloscaccia
"""
import numpy as np
from   glob import glob
import xarray as xr
from   datetime import datetime
import os


def parser():
    import argparse
    
    parser=argparse.ArgumentParser()
    parser.add_argument('-d',
                    '--directory',
                    type=str,
                    required = True,
                    help='Directory with ndvi')
    return parser.parse_args()

def aggregate_netcdf(files):
    files.sort()
    for i, file in enumerate(files):
        if i == 0:
             agg_ndvi = xr.open_dataset(file)
             agg_ndvi['times'] = [ os.path.basename(file).split('_')[-1].replace('.nc','') ]
             agg_ndvi['cm']    = agg_ndvi['cm'].astype(int)
             agg_ndvi['ndvi']  = agg_ndvi['ndvi'].fillna(0)
        else:
             with xr.open_dataset(file) as tmp_dataset:
                 tmp_dataset['ndvi'] = tmp_dataset['ndvi'].fillna(0)
                 agg_ndvi['ndvi'] += tmp_dataset['ndvi']
                 agg_ndvi['cm']    += tmp_dataset['cm'].astype(int)
                 agg_ndvi['times'] = np.append( agg_ndvi['times'], os.path.basename(file).split('_')[-1].replace('.nc','') )
 
    denom = len(files) - agg_ndvi['cm']
    print(agg_ndvi['ndvi'])
    agg_ndvi['ndvi'] /= denom

    agg_ndvi['n_good_quarters'] = denom
    agg_ndvi['cm'].attrs['longname']   = "Daily Aggregated Cloud Mask"
    agg_ndvi['n_good_quarters'].attrs['longname']   = "Number of masked segments (max 96)"
    agg_ndvi['ndvi'].attrs['longname'] = "Daily Average NDVI"
             
    return agg_ndvi


if __name__ == '__main__':
    args = parser()
    
    ndvi_files = glob(args.directory + '/ndvi_*.nc')
    if len(ndvi_files) != 96:
        print("Warning: NDVI data do not cover the enitire day: " + repr(len(ndvi_files)) +' out of 96 files')

    agg_ndvi = aggregate_netcdf(ndvi_files)

    outfile = args.directory + '/daily_ndvi.nc'
    agg_ndvi.to_netcdf(outfile)
    
    print(datetime.now(), "Saved aggregated NDVI in ",outfile)

        
    
