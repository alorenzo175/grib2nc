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


import numpy as np
import pandas as pd
import requests
import pygrib
import netCDF4 as nc4


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
            self.level = level
        self.init_time = init_time

        self.base_path = os.path.join(self.download_dict['folder'], 
                                      init_time.strftime('%Y'),
                                      init_time.strftime('%m'),
                                      init_time.strftime('%d'),
                                      init_time.strftime('%Hz'),
                                      'grib')
        if not os.path.isdir(self.base_path):
            os.makedirs(self.base_path)
        self.ncfilename = 'hrrr.{init_time}.{level}.nc'.format(
            init_time=init_time.strftime('%Y%m%d%H'), level=level)

    def read_index(self):
        idx_files = [afile for afile in os.listdir(self.base_path) 
                     if afile.endswith('.idx')]
        index_df = None
        for idx_file in sorted(idx_files):
            grib_file_desc =  pd.read_table(os.path.join(self.base_path, 
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
        self.ncfile = nc4.Dataset(os.path.join(self.base_path, self.ncfilename), 
                              'w', format='NETCDF4')
        self.ncfile.createDimension('Time', None)
        self.ncfile.createDimension('time', 2)
        self.ncfile.createDimension('DateStrLen', 19)
        self.ncfile.createDimension('west_east', None)
        self.ncfile.createDimension('south_north', None)
        self.ncfile.createVariable('Times', 'S1', ('Time', 'DateStrLen'), zlib=True)
        self.ncfile.createVariable('XLAT', 'f4', ('time', 'south_north', 'west_east'), zlib=True)
        self.ncfile.createVariable('XLON', 'f4', ('time', 'south_north', 'west_east'), zlib=True)
        for ncfield in self.config.options('surface_settings'):
            self.ncfile.createVariable(ncfield, 'f4', ('Time', 'south_north', 'west_east'), zlib=True)

    def convert(self):
        times = self.ncfile.variables['Times']
        fhs= self.index_df['forecast_hour'].unique()
        for fh in fhs:
            if isinstance(fh, float):
                continue
            if 'hour' in fh and '-' not in fh:
                thetime = self.init_time + dt.timedelta(hours=int(fh[:2]))
            elif 'anl' in fh:
                thetime = self.init_time
            else:
                thetime = None
            if thetime is not None:
                row_time = nc4.stringtoarr(
                    thetime.strftime('%Y-%m-%d_%H:%M:%S'), 19)
                times[len(times)] = row_time

        lats = self.ncfile.variables['XLAT']
        lons = self.ncfile.variables['XLON']

        for nc_field, grib_f in self.config.items('surface_settings'):
            field, vertical_layer = grib_f.split(',')
            relevant_df = self.index_df[(self.index_df['field'] == field) & 
                (self.index_df['vertical_layer'] == vertical_layer)]
            relevant_df.set_index('filename', inplace=True)
            for filename, series in relevant_df.iterrows():
                grbs = pygrib.open(os.path.join(self.base_path, filename))

                try:
                    grb = grbs[series.grib_level]
                except IOError:
                    continue
                print grb.latlons()[0]
                print grb.latlons()[1]
                print grb.values.shape

                data, lat, lon = grb.data(
                    lat1=self.config.getfloat('subdomain', 'lat1'),
                    lat2=self.config.getfloat('subdomain', 'lat2'),
                    lon1=self.config.getfloat('subdomain', 'lon1'),
                    lon2=self.config.getfloat('subdomain', 'lon2'))


                print data.shape
                lats[:] = lat
                lons[:] = lon
                fvar = self.ncfile.variables[field]
                fvar[len(fvar),:,:] = data
                grbs.close()

            


def main():
    logging.basicConfig(level=logging.DEBUG)
    g2nc = Grib2NC(dt.datetime(2014,10,5,1), 'surface')
    g2nc.read_index()
    g2nc.setup_netcdf()
    g2nc.convert()
    g2nc.ncfile.close()
    return g2nc.index_df


def exceptionlogging(*exc_info):
    text = "".join(traceback.format_exception(*exc_info))
    logger = logging.getLogger()
    logger.exception("Unhandled exception: %s" % text)


if __name__ == '__main__':
    sys.excepthook = exceptionlogging
    main()
