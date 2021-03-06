#!/usr/env python
"""
Convert grib forecasts to netCDF


.. codeauthor:: Tony Lorenzo <alorenzo175@gmail.com>
"""


import sys
import os
import logging
import traceback
import datetime as dt
try:
    import ConfigParser as configparser
except ImportError:
    import configparser
import subprocess


from pkg_resources import resource_filename, Requirement
import numpy as np
import pandas as pd
import pygrib
import netCDF4 as nc4


WRGRIB2 = '/usr/bin/wgrib2'


class EmptyError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return str(self.value)


class Grib2NC(object):
    """Grib2NC converts many grib files into a single netCDF
    forecast file.

    Parameters
    ----------
    init_time : datetime
        The initilization time of the forecast to convert
    level : str
        The HRRR forecast "level" i.e. surface, native, pressures, or subhourly
    """

    def __init__(self, init_time, level, config_path=None, grib_path=None,
                 netcdf_path=None, ncfilename=None):
        self.logger = logging.getLogger()
        self.config = configparser.ConfigParser()
        if config_path is not None:
            self.config_path = config_path
        else:
            self.config_path = resource_filename(
                Requirement.parse('grib2nc'), 'settings.txt')
        self.config.read(self.config_path)
        self.download_dict = dict(self.config.items('download_settings'))
        self.hrrr_type_dict = dict(self.config.items('output_types'))
        if level.lower() not in self.hrrr_type_dict:
            raise Exception
        else:
            self.level = self.hrrr_type_dict[level.lower()]
            if level in ['native', 'pressure']:
                self.vertical = True
            else:
                self.vertical = False
        self.init_time = init_time
        self.grib_vars = self.config.items(
            '{level}_settings'.format(level=level))

        self.grib_path = (grib_path or
                          self.download_dict['grib_base_folder'].format(
                              year=init_time.strftime('%Y'),
                              month=init_time.strftime('%m'),
                              day=init_time.strftime('%d'),
                              hour=init_time.strftime('%H'),
                              level=level))

        self.netcdf_path = (netcdf_path or
                            self.download_dict['netcdf_base_folder'].format(
                                year=init_time.strftime('%Y'),
                                month=init_time.strftime('%m'),
                                day=init_time.strftime('%d'),
                                hour=init_time.strftime('%H')))

        if not os.path.isdir(self.grib_path):
            old_umask = os.umask(0)
            os.makedirs(self.grib_path, mode=0o775)
            os.umask(old_umask)
        if not os.path.isdir(self.netcdf_path):
            old_umask = os.umask(0)
            os.makedirs(self.netcdf_path, mode=0o775)
            os.umask(old_umask)
        self.ncfilename = (ncfilename or
                           self.download_dict['netcdf_filename'].format(
                               init_time=init_time.strftime('%Y%m%d%H'),
                               level=level))

    def make_index_files(self):
        grib2_files = [afile for afile in os.listdir(self.grib_path)
                       if afile.endswith('.grib2')]
        if not grib2_files:
            raise EmptyError('No .grib2 files in path {}'.format(
                self.grib_path))
        for afile in grib2_files:
            command = '{wgrib} -set center 7 -s {file} > {new_file}'.format(
                wgrib=WRGRIB2, file=os.path.join(self.grib_path, afile),
                new_file=os.path.join(self.grib_path, afile[:-6]+'.idx'))
            try:
                subprocess.check_output([command, ], shell=True,
                                        stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                self.logger.error(e.output)
                raise e

    def read_index(self):
        """Read the .idx files

        """
        idx_files = [afile for afile in os.listdir(self.grib_path)
                     if afile.endswith('.idx') and self.level in afile]

        if not idx_files:
            self.make_index_files()
            idx_files = [afile for afile in os.listdir(self.grib_path)
                         if afile.endswith('.idx')]
            if not idx_files:
                raise EmptyError('No index files found')

        index_df = None
        for idx_file in sorted(idx_files):
            grib_file_desc = pd.read_table(
                os.path.join(self.grib_path, idx_file), sep=':', engine='c',
                lineterminator='\n', header=None,
                names=['grib_level', 'start_byte', 'init_time', 'field',
                       'vertical_layer', 'forecast_hour', 'end_byte'],
                index_col=False)
            grib_file_desc['end_byte'] = grib_file_desc.start_byte.shift(-1) - 1  # NOQA
            filename_series = pd.Series(
                [idx_file[:-4] for i in range(
                    len(grib_file_desc['grib_level']))])
            grib_file_desc['filename'] = filename_series

            if index_df is None:
                index_df = grib_file_desc
            else:
                index_df = index_df.append(grib_file_desc, ignore_index=True)

        if index_df['grib_level'].dtype == float:
            mapping = {}
            for i, floatnum in enumerate(sorted(
                    index_df['grib_level'].unique())):
                mapping[floatnum] = i+1

            index_df['grib_level'] = index_df['grib_level'].map(mapping)

        self.index_df = index_df
        return index_df

    def setup_netcdf(self):
        """Setup the netCDF file and add lat/lon from the grib files

        """
        # first check the grid
        if not hasattr(self, 'index_df'):
            self.read_index()

        for filename in sorted(self.index_df['filename'].unique()):
            try:
                grbs = pygrib.open(os.path.join(self.grib_path, filename))
                grb = grbs[1]
            except IOError:
                try:
                    grbs = pygrib.open(os.path.join(self.grib_path,
                                                    filename+'.grib2'))
                    grb = grbs[1]
                except IOError:
                    continue

            if 'ELON' in self.index_df['field'] and 'NLAT' in self.index_df['field']:  # NOQA
                lats_ind = self.index_df[
                    (self.index_df['field'] == 'NLAT') &
                    (self.index_df['filename'] == filename)].grib_level
                lats = grbs[lats_ind].values

                lons_ind = self.index_df[
                    (self.index_df['field'] == 'ELON') &
                    (self.index_df['filename'] == filename)].grib_level
                lons = grbs[lons_ind].values

                self.index_df = self.index_df[
                    (self.index_df['field'] != 'ELON') &
                    (self.index_df['field'] != 'NLAT')]
            else:
                lats, lons = grb.latlons()

            lat1 = self.config.getint('subdomain', 'lat1')
            lat2 = self.config.getint('subdomain', 'lat2')
            lon1 = self.config.getint('subdomain', 'lon1')
            lon2 = self.config.getint('subdomain', 'lon2')

            clons = lons[int(np.where((lats >= lat1) &
                                      (lats <= lat2))[0].mean()), :]
            clats = lats[:, int(np.where((lons >= lon1) &
                                         (lons <= lon2))[0].mean())]
            latlims = np.where((clats >= lat1) & (clats <= lat2))[0]
            lonlims = np.where((clons >= lon1) & (clons <= lon2))[0]
            start_date = grb.analDate
            dx = grb.Dx*1.0/1000
            dy = grb.Dy*1.0/1000
            grbs.close()
            break

        self.max_lat = latlims.max() + 1
        self.min_lat = latlims.min()
        self.max_lon = lonlims.max() + 1
        self.min_lon = lonlims.min()

        # move vvel somewhere in config
        if self.vertical:
            nvert = self.index_df.groupby(['filename', 'field']
                                          ).size().loc[:, 'VVEL'].iloc[0]
        else:
            nvert = None

        path = os.path.join(self.netcdf_path, self.ncfilename)
        # check here if it already exists and just append
        self.ncwriter = NCWriter()
        self.ncwriter.make_netcdf(path, self.max_lat-self.min_lat,
                                  self.max_lon-self.min_lon,
                                  self.vertical, nvert)

        self.ncwriter.set_variable('XLAT', lats[self.min_lat:self.max_lat,
                                                self.min_lon:self.max_lon], 0)
        self.ncwriter.set_variable('XLONG', lons[self.min_lat:self.max_lat,
                                                 self.min_lon:self.max_lon], 0)
        self.ncwriter.set_attribute('START_DATE',
                                    start_date.strftime('%Y-%m-%d_%H:%M:%S'))
        self.ncwriter.set_attribute('DX', dx)
        self.ncwriter.set_attribute('DY', dy)
        self.ncwriter.set_attribute('WEST-EAST_GRID_DIMENSION',
                                    self.max_lon - self.min_lon + 1)
        self.ncwriter.set_attribute('SOUTH-NORTH_GRID_DIMENSION',
                                    self.max_lat - self.min_lat + 1)
        self.ncwriter.set_attribute('BOTTOM-TOP_GRID_DIMENSION', nvert or 0)

    def load_into_netcdf(self):
        """Load the grib files into the netCDF file that has been setup

        """

        if not hasattr(self, 'ncwriter'):
            self.setup_netcdf()

        self.logger.info('Making netCDF file %s' % self.ncfilename)
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
                        (self.index_df['field'] == field)], ignore_index=True)

        relevant_df.set_index('filename', inplace=True)
        times = []
        levels = []
        for filename in sorted(relevant_df.index.unique()):
            try:
                grbs = pygrib.open(os.path.join(self.grib_path, filename))
            except IOError:
                try:
                    grbs = pygrib.open(os.path.join(self.grib_path,
                                                    filename+'.grib2'))
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
                    self.ncwriter.set_variable('Times', nc4.stringtoarr(
                        thetime.strftime('%Y-%m-%d_%H:%M:%S'), 19), timed)
                else:
                    timed = times.index(thetime)

                if self.vertical:
                    thelevel = grb.level
                    if thelevel not in levels:
                        leveld = len(levels)
                        levels.append(thelevel)
                    else:
                        leveld = levels.index(thelevel)

                    nc_field = field_dict[series.field].upper()
                else:
                    nc_field = field_dict[(series.field,
                                           series.vertical_layer)].upper()

                if not self.ncwriter.check_variable(nc_field):
                    self.ncwriter.add_variable(nc_field, description=grb.name,
                                               units=grb.units,
                                               vertical=self.vertical)

                if not self.vertical:
                    self.ncwriter.set_variable(
                        nc_field, grb.values[self.min_lat:self.max_lat,
                                             self.min_lon:self.max_lon],
                        timed)
                else:
                    ivals = grb.values[self.min_lat:self.max_lat,
                                       self.min_lon:self.max_lon]
                    try:
                        self.ncwriter.set_variable(nc_field, ivals,
                                                   [timed, leveld])

                    except RuntimeError:
                        self.logger.debug(
                            'Error with leveld = %s and var = %s' %
                            (leveld, grb.cfName))
                    else:
                        if self.level == 'wrfprs':
                            if not self.ncwriter.check_variable('P'):
                                self.nc_writer.add_variable('P', units='mb',
                                                            vertical=True)

                            self.ncwriter.set_variable(
                                'P', (np.ones(ivals.shape) * grb.level),
                                [timed, leveld])
        self.ncwriter.close()
        return os.path.join(self.netcdf_path, self.ncfilename)


class NCWriter(object):
    """NCWriter handles all the writing and setup of the netCDF file

    """

    def __init__(self):
        self.logger = logging.getLogger('NCWriter')

    def make_netcdf(self, ncfilename='base.nc', nlats=1059, nlons=1799,
                    vertical=False, nvert=40):
        """Create the netCDF file and initilize dimensions

        Parameters
        ----------
        ncfilename : str
        nlats : int
            length of lat (south_north) dimension
        nlons : int
            length of lon (west_east) dimension
        vertical : boolean
            whether the bottom_top dimension should be used
        nvert : int
            length of vertical (bottom_top) dimension
        """
        XLONG_DICT = {'FieldType': 104, 'MemoryOrder': "XY",
                      'stagger': "", 'units': 'degree_east',
                      'description': "LONGITUDE, WEST IS NEGATIVE"}
        XLAT_DICT = {'FieldType': 104, 'MemoryOrder': "XY",
                     'stagger': "", 'units': 'degree_north',
                     'description': "LATITUDE, SOUTH IS NEGATIVE"}

        self.netc = nc4.Dataset(ncfilename, 'w', format='NETCDF4')
        self.netc.createDimension('Time', None)
        self.netc.createDimension('time', 2)
        self.netc.createDimension('south_north', nlats)
        self.netc.createDimension('west_east', nlons)
        self.netc.createDimension('DateStrLen', 19)
        if vertical:
            self.netc.createDimension('bottom_top', nvert)

        self.netc.createVariable('Times', 'S1', ('Time', 'DateStrLen'),
                                 zlib=True, complevel=1)
        lons = self.netc.createVariable('XLONG', 'f4',
                                        ('time', 'south_north',
                                         'west_east'), zlib=True, complevel=1)
        lats = self.netc.createVariable('XLAT', 'f4', ('time', 'south_north',
                                                       'west_east'),
                                        zlib=True, complevel=1)
        lons.setncatts(XLONG_DICT)
        lats.setncatts(XLAT_DICT)
        self.netc.sync()

    def add_variable(self, var_name, field_type=104, mem_order="XY",
                     stagger='', units='', description='', dformat='f4',
                     vertical=False):
        """Create a new variable in the netCDF file

        Parameters
        ----------
        var_name : str
            name of the new variable
        field_type : int or str
            value for netCDF FieldType attribute
        mem_order : str
            value for netCDF MemoryOrder attribute
        stagger : str
            value for netCDF stagger attribute
        units : str
            units of the variable
        description : str
            description of the variable
        dformat : str
            data type for the variable
        vertical : boolean
            whether the new variable should be index by the bottom_top dim
        """
        if vertical:
            var = self.netc.createVariable(var_name, dformat,
                                           ('Time', 'bottom_top',
                                            'south_north', 'west_east'),
                                           zlib=True, complevel=1)
        else:
            var = self.netc.createVariable(
                var_name, dformat, ('Time', 'south_north', 'west_east'),
                zlib=True, complevel=1)
        var.description = description
        var.units = units
        var.MemoryOrder = mem_order
        var.FieldType = field_type
        var.stagger = stagger
        self.netc.sync()

    def set_variable(self, var_name, data, *position):
        """Inserts data into the variable

        Parameters
        ----------
        var_name : str
            variable name
        data : numpy array
            data to insert
        *position : int or tuple or list
            numpy slicing index for position to place data in variable
        """
        self.logger.debug('Putting data into %s at %s' % (var_name, position))
        var = self.netc.variables[var_name]
        var[position] = data
        return var

    def set_attribute(self, attr_name, value):
        """Set file attributes

        Parameters
        ----------
        attr_name: str
            attribute name
        value : str or int or float
            value of the attribute
        """
        self.logger.debug('Setting attribute {} to {}'.format(attr_name,
                                                              value))
        self.netc.setncattr(attr_name, value)

    def check_variable(self, var_name):
        """Check if the variable has been added to the file yet
        """
        return var_name in self.netc.variables

    def close(self):
        """Close the netCDF file
        """
        if hasattr(self, 'netc'):
            self.netc.close()


def main():
    logging.basicConfig(level=logging.DEBUG)
    g2nc = Grib2NC(dt.datetime(2014, 10, 9, 0), 'subhourly')
    g2nc.load_into_netcdf()
    return g2nc.index_df


def exceptionlogging(*exc_info):
    text = "".join(traceback.format_exception(*exc_info))
    logger = logging.getLogger()
    logger.exception("Unhandled exception: %s" % text)


if __name__ == '__main__':
    sys.excepthook = exceptionlogging
    main()
