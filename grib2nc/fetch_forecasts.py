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


import pandas as pd
import requests


class RequestError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class Fetcher(object):
    def __init__(self, init_time, level, method='ftp'):
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

    def connect_ftp(self):
        self.ftp = ftplib.FTP(self.download_dict['ftp_host'])
        self.ftp.login()
        self.ftp.cwd(self.download_dict['ftp_dir'])
        stored_days = self.ftp.nlst()
        it_date = self.init_time.strftime('hrrr.%Y%m%d')
        if it_date not in stored_days:
            raise RequestError('Requested forecast initilization day not on NCEP server')
        self.ftp.cwd(it_date)

    def read_index_ftp(self):
        self.fetch_ftp('grib2.idx', False)
        idx_files = [afile for afile in os.listdir(self.base_path) 
                     if afile.endswith('.idx')]
        
        self.index_panel = None
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
            grib_file_desc.set_index('field', inplace=True)
            if self.index_panel is None:
                self.index_panel = pd.Panel({idx_file[:-4]:grib_file_desc})
            else:
                self.index_panel[idx_file[:-4]] = grib_file_desc

        return self.index_panel
        
    def fetch_ftp(self, extension='grib2', overwrite=True):
        if not hasattr(self, 'ftp'):
            self.connect_ftp()
        base_name = 'hrrr.t{init_hour}z.{level}f*.{extension}'.format(
            init_hour=self.init_time.strftime('%H'), 
            level=self.hrrr_type_dict[self.level],
            extension=extension)
        files = self.ftp.nlst(base_name)
        total_size = 0
        for afile in files:
            total_size += self.ftp.size(afile)
        self.logger.info('Attempting to retrieve %s bytes of data' % total_size)
        start = time.time()
        for afile in files:
            local_path = os.path.join(self.base_path,  afile)
            if os.path.isfile(local_path) and not overwrite:
                continue
            with open(local_path, 'wb') as f:
                def callback(data):
                    f.write(data)
                self.ftp.retrbinary('RETR %s' % afile, callback)
        end = time.time()
        self.logger.info('Retrieved %s files in %s seconds' % (len(files), 
                                                               end-start))
        self.logger.debug('Files retrieved are %s' % os.listdir(self.base_path))

    def partial_fields_http(self, field, chunk_size=4098):
        to_retrieve = self.index_panel.iloc[0].loc[field]
        afile = self.index_panel.items[0]

        start = time.time()
        for i in range(len(to_retrieve)):
            df = to_retrieve.iloc[i]
            grib_level = df.grib_level
            local_path = os.path.join(self.base_path, afile + '.%s.%s' % 
                                      (field, grib_level))
            url = (self.download_dict['html_site'] + 
                   self.init_time.strftime('hrrr.%Y%m%d') + 
                   '/' + afile)
            headers = {'range': 'bytes=%s-%s' % (df.start_byte, df.end_byte)}
            r = requests.get(url, headers=headers, stream=True)
            with open(local_path, 'wb') as fd:
                for chunk in r.iter_content(chunk_size):
                    fd.write(chunk)
            break
        end = time.time()
        self.logger.info('Retrieved files in %s seconds' % (end-start))
        

def main():
    logging.basicConfig(level=logging.DEBUG)
    f = Fetcher(dt.datetime(2014,10,7,1), 'subhourly')
    #f.read_index_ftp()
    f.fetch_ftp()


def exceptionlogging(*exc_info):
    text = "".join(traceback.format_exception(*exc_info))
    logger = logging.getLogger()
    logger.exception("Unhandled exception: %s" % text)


if __name__ == '__main__':
    sys.excepthook = exceptionlogging
    main()
