#!/usr/env python
"""
A program to download the requested HRRR forecasts


.. codeauthor:: Tony Lorenzo <alorenzo175@gmail.com>
"""


import sys
import os
import logging
import traceback
import datetime as dt
import time
try:
    import ConfigParser as configparser
except ImportError:
    import configparser
import ftplib
import shutil


import numpy as np
import pandas as pd
import requests
import pygrib
import netCDF4 as nc4

from prepare_template import make_netcdf


class Grib2NC(object):
    def __init__(self, init_time, level):
        self.logger = logging.getLogger()
        self.config = configparser.ConfigParser()
        self.config.read(os.path.join(
            os.path.dirname(os.path.realpath('__file__')),
            'settings.txt'))
        self.download_dict = dict(self.config.items('download_settings'))
        self.hrrr_type_dict = dict(self.config.items('output_types'))
        if level not in self.hrrr_type_dict:
            raise Exception
        else:
            self.level = self.hrrr_type_dict[level]
        self.init_time = init_time
        self.grib_vars = self.config.items(
            '{level}_settings'.format(level=level))

        self.base_path = os.path.join(self.download_dict['folder'], 
                                      init_time.strftime('%Y'),
                                      init_time.strftime('%m'),
                                      init_time.strftime('%d'),
                                      init_time.strftime('%Hz'))
        self.grib_path = os.path.join(self.base_path, 'grib')
        self.netcdf_path = os.path.join(self.base_path, 'netcdf')
        if not os.path.isdir(self.grib_path):
            os.makedirs(self.grib_path)
        if not os.path.isdir(self.netcdf_path):
            os.makedirs(self.netcdf_path)
        self.ncfilename = 'hrrr.{init_time}.{level}.nc'.format(
            init_time=init_time.strftime('%Y%m%d%H'), level=level)

    def read_index(self):
        idx_files = [afile for afile in os.listdir(self.grib_path) 
                     if afile.endswith('.idx') and self.level in afile]
        index_df = None
        for idx_file in sorted(idx_files):
            grib_file_desc =  pd.read_table(os.path.join(self.grib_path,
                                                         idx_file), 
                                            sep=':', engine='c', 
                                            lineterminator='\n', header=None, 
                                            names=['grib_level','start_byte', 
                                                   'init_time', 'field', 
                                                   'vertical_layer', 
                                                   'forecast_hour', 'end_byte'],
                                            index_col=False)
            grib_file_desc['end_byte'] = grib_file_desc.start_byte.shift(-1) - 1
            filename_series = pd.Series([idx_file[:-4] for i in range(len(grib_file_desc['grib_level']))])
            grib_file_desc['filename'] = filename_series
            
            if index_df is None:
                index_df = grib_file_desc
            else:
                index_df = index_df.append(grib_file_desc, ignore_index=True)

        self.index_df = index_df
        return index_df

    def setup_netcdf(self):
        # first check the grid
        if not hasattr(self, 'index_df'):
            self.read_index()

        for filename in sorted(self.index_df['filename'].unique()):
            try:
                grbs = pygrib.open(os.path.join(self.grib_path, filename))
                grb = grbs[1]
            except IOError:
                continue
            lats, lons = grb.latlons()

            lat1 = self.config.getint('subdomain', 'lat1')
            lat2 = self.config.getint('subdomain', 'lat2')
            lon1 = self.config.getint('subdomain', 'lon1')
            lon2 = self.config.getint('subdomain', 'lon2')

            clons = lons[int(np.where((lats >= lat1) 
                                      & (lats <= lat2)
                                  )[0].mean()), :]
            clats = lats[:,int(np.where((lons >= lon1 ) 
                                        & (lons <= lon2))[0].mean())]
            latlims = np.where((clats >= lat1) & (clats <= lat2))[0]
            lonlims = np.where((clons >= lon1) & (clons <= lon2))[0]
            grbs.close()
            break
        
        self.max_lat = latlims.max() + 1
        self.min_lat = latlims.min()
        self.max_lon = lonlims.max() + 1
        self.min_lon = lonlims.min()

        path = os.path.join(self.netcdf_path, self.ncfilename)
        make_netcdf(path, self.max_lat-self.min_lat,
                    self.max_lon-self.min_lon)
        self.ncfile = nc4.Dataset(path, 'a', format='NETCDF4')
        nclats = self.ncfile.variables['XLAT']
        nclons = self.ncfile.variables['XLONG']
        nclats[0] = lats[self.min_lat:self.max_lat,
                         self.min_lon:self.max_lon]
        nclons[0] = lons[self.min_lat:self.max_lat,
                         self.min_lon:self.max_lon]
        
    def load_into_netcdf(self):
        if not hasattr(self, 'ncfile'):
            self.setup_netcdf()
        
        field_dict = {}
        relevant_df = None
        for nc_field, grib_f in self.grib_vars:
            field, vertical_layer = grib_f.split(',')
            if vertical_layer is not '':
                field_dict[(field, vertical_layer)] = nc_field
                if relevant_df is None:
                    relevant_df = self.index_df[
                        (self.index_df['field'] == field) & 
                        (self.index_df['vertical_layer'] == vertical_layer)]
                else:
                    relevant_df = relevant_df.append(self.index_df[
                        (self.index_df['field'] == field) & 
                        (self.index_df['vertical_layer'] == vertical_layer)], 
                                                     ignore_index=True)
            else:
                field_dict[field] = nc_field
                if relevant_df is None:
                    relevant_df = self.index_df[
                        (self.index_df['field'] == field)]
                else:
                    relevant_df = relevant_df.append(self.index_df[
                        (self.index_df['field'] == field)]
                                                     ,ignore_index=True)

        relevant_df.set_index('filename', inplace=True)
        times = []
        for filename in sorted(relevant_df.index.unique()):
            try:
                grbs = pygrib.open(os.path.join(self.grib_path, filename))
            except IOError:
                continue

            for filename, series in relevant_df.loc[filename].iterrows():
                try:
                    grb = grbs[series.grib_level]
                except IOError:
                    continue

                thetime = grb.validDate
                if thetime not in times:
                    timed = len(times)
                    times.append(thetime)
                else:
                    timed = times.index(thetime)

                self.ncfile.variables['Times'][timed] = nc4.stringtoarr(thetime.strftime('%Y-%m-%d_%H:%M:%S'), 19)

                if (series.field, series.vertical_layer) in field_dict:
                    nc_field = field_dict[(series.field, 
                                           series.vertical_layer)].upper()
                elif series.field in field_dict:
                    nc_field = field_dict[series.field].upper()
                else:
                    self.logger.warning("Why is this file open?")
                    continue

                if nc_field not in self.ncfile.variables:
                    self.ncfile.createVariable(nc_field, 'f4', ('Time', 
                        'south_north', 'west_east'), zlib=True, complevel=1)
                    var = self.ncfile.variables[nc_field]
                    var.description = grb.name
                    var.units = grb.units
                    var.MemoryOrder = 'XY'
                    var.FieldType = 104

                var = self.ncfile.variables[nc_field]
                var[timed] = grb.values[self.min_lat:self.max_lat,
                                        self.min_lon:self.max_lon]
        

def main():
    logging.basicConfig(level=logging.DEBUG)
    g2nc = Grib2NC(dt.datetime(2014,10,7,1), 'subhourly')
    g2nc.read_index()
    g2nc.setup_netcdf()
    g2nc.load_into_netcdf()
    g2nc.ncfile.close()
    return g2nc.index_df


def exceptionlogging(*exc_info):
    text = "".join(traceback.format_exception(*exc_info))
    logger = logging.getLogger()
    logger.exception("Unhandled exception: %s" % text)


if __name__ == '__main__':
    sys.excepthook = exceptionlogging
    main()
