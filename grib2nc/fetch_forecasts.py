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
import argparse
try:
    import ConfigParser as configparser
except ImportError:
    import configparser


class RequestError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class HRRRFetcher(object):
    def __init__(self, init_time, level, method='ftp', config_path=None):
        self.logger = logging.getLogger('HRRRFetcher')
        self.config = configparser.ConfigParser()
        config_path = config_path or os.path.join(
            os.path.dirname(os.path.realpath('__file__')),
            'settings.txt')
        self.config.read(config_path)
        self.download_dict = dict(self.config.items('download_settings'))
        self.hrrr_type_dict = dict(self.config.items('output_types'))
        if level not in self.hrrr_type_dict:
            raise Exception
        else:
            self.level = level
        
        if isinstance(init_time, str):
            try:
                import dateutil.parser as dparser
            except ImportError:
                raise ImportError(
                    'dateutil is required for string datetime input')
            self.init_time = dparser.parse(init_time)
        elif isinstance(init_time, dt.datetime):
            self.init_time = init_time
        else:
            raise TypeError('init_time must be a datetime string or datetime object')

        self.default_method = method
        self.base_path = os.path.join(self.download_dict['folder'], 
                                      self.init_time.strftime('%Y'),
                                      self.init_time.strftime('%m'),
                                      self.init_time.strftime('%d'),
                                      self.init_time.strftime('%Hz'),
                                      'grib')
            
        if not os.path.isdir(self.base_path):
            os.makedirs(self.base_path)

    def fetch(self, init_time=None, level=None, overwrite=True):
        if self.default_method == 'ftp':
            return self.fetch_ftp(init_time, level, overwrite)
        elif self.default_method == 'http':
            return self.fetch_http(init_time, level, overwrite)
        else:
            raise AttributeError('Method must be ftp or http')

    def connect_ftp(self, init_time=None):
        try:
            import ftplib
        except ImportError:
            raise ImportError('FTP requires the ftplib module')
        self.ftp = ftplib.FTP(self.download_dict['ftp_host'])
        self.ftp.login()
        self.ftp.cwd(self.download_dict['ftp_dir'])
        stored_days = self.ftp.nlst()
        init_time = init_time or self.init_time
        it_date = init_time.strftime('hrrr.%Y%m%d')
        if it_date not in stored_days:
            raise RequestError('Requested forecast initilization day '+
                               'not on NCEP server')
        self.ftp.cwd(it_date)

    def fetch_ftp(self, init_time=None, level=None, overwrite=True):
        level = level or self.level
        init_time = init_time or self.init_time
        if not hasattr(self, 'ftp'):
            self.connect_ftp(init_time)

        forecast_name ='hrrr.t{init_hour}z.{level}f*'.format(
            init_hour=init_time.strftime('%H'),
            level=self.hrrr_type_dict[level])

        files = self.ftp.nlst(forecast_name)
        total_size = []
        expected_size = 0
        nfiles = 0
        for afile in files:
            expected_size += self.ftp.size(afile)
        self.logger.info(('Attempting to retrieve %0.2f MB of data in %s'+
                         ' files over FTP')
                         % ((1.0*expected_size/1024**2), len(files)))
        start = time.time()
        for afile in files:
            self.logger.debug('Retrieving %s' % afile)
            local_path = os.path.join(self.base_path,  afile)
            if os.path.isfile(local_path) and not overwrite:
                self.logger.warning('%s already exists!' % afile)
                continue
            with open(local_path, 'wb') as f:
                def callback(data):
                    f.write(data)
                    total_size.append(len(data))
                try:
                    self.ftp.retrbinary('RETR %s' % afile, callback)
                except KeyboardInterrupt:
                    break
            nfiles += 1
        end = time.time()
        total_size = sum(total_size)
        self.logger.info('Retrieved %s files in %s seconds' % (nfiles, 
                                                               end-start))
        self.logger.info('Download rate: %0.2f MB/s' % 
                         (1.0*total_size/1024**2/(end-start)))
        self.logger.debug('Files retrieved are %s' % os.listdir(self.base_path))

    def fetch_http(self, init_time=None, level=None, overwrite=True):
        try:
            import requests
            import pandas as pd
        except ImportError:
            raise ImportError('HTTP requires the requests and BeautifulSoup4')

        init_time = init_time or self.init_time
        level = level or self.level

        forecast_name ='hrrr.t{init_hour}z.{level}f'.format(
            init_hour=init_time.strftime('%H'),
            level=self.hrrr_type_dict[level])
        
        html_folder = (self.download_dict['html_site'] + 
            'hrrr.{dtime}'.format(dtime=init_time.strftime('%Y%m%d')))

        if not hasattr(self, 'session'):
            self.session = requests.Session()

        forecast_table = self.session.get(html_folder)

        if forecast_table.status_code != 200:
            raise RequestError('Invalid URL/date')
            
        table = pd.io.html.read_html(forecast_table.content,
                                     header=0, parse_dates=True,
                                     infer_types=False)[0]
        
        files_df= table[table['Name'].str.contains(
            forecast_name)]
        
        def process_size(size):
            if size.endswith('K'):
                out = float(size[:-1])*1024
            elif size.endswith('M'):
                out = float(size[:-1])*1024**2
            else:
                out = float(size)
            return out

        expected_size = files_df['Size'].map(process_size).sum()
        self.logger.info('Downloading %0.2f MB of data' 
                         % (1.0*expected_size/1024**2))
        start = time.time()
        total_size = 0
        nfiles = 0
        for filename in files_df['Name']:
            self.logger.debug('Retrieving %s' % filename)
            local_path = os.path.join(self.base_path,  filename)
            if os.path.isfile(local_path) and not overwrite:
                self.logger.warning('Failed to retrieve %s' % filename)
                continue
            try:
                r = self.session.get(html_folder + '/' + filename, stream=True)

                with open(local_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1024):
                        if chunk:
                            f.write(chunk)
                            f.flush()
                            total_size += len(chunk)
                nfiles += 1
            except KeyboardInterrupt:
                break

        end = time.time()
        self.logger.info('Retrieved %s files in %s seconds' % 
                         (nfiles, end-start))
        self.logger.info('Download rate: %0.2f MB/s' % 
                         (1.0*total_size/1024**2/(end-start)))
        self.logger.debug('Files retrieved are %s' % os.listdir(self.base_path))

    def close(self):
        if hasattr(self, 'ftp'):
            self.ftp.close()
        if hasattr(self, 'session'):
            self.session.close()

def main():
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s')
    
    argparser = argparse.ArgumentParser(
        description='Download HRRR grib files and store select fields in netCDF format')
    argparser.add_argument('-v', '--verbose', help='Increase logging verbosity',
                           action='count')
    argparser.add_argument('-p', '--protocol', help='Download protocol',
                           choices=['ftp', 'http'], default='ftp')
    argparser.add_argument('-c', '--config', help='Config file path')
    argparser.add_argument('-l', '--level', help='HRRR level subtype',
                           choices=['surface', 'subhourly','pressure', 
                                    'native'],
                           default='subhourly')
    argparser.add_argument('INITDT', help='Initilization datetime')
    args = argparser.parse_args()
    
    if args.verbose == 1:
        logging.getLogger().setLevel(logging.INFO)
    elif args.verbose > 1:
        logging.getLogger().setLevel(logging.DEBUG)
        
    f = HRRRFetcher(args.INITDT, args.level, args.protocol, 
                    args.config)
    f.fetch()
    f.close()


def exceptionlogging(*exc_info):
    text = "".join(traceback.format_exception(*exc_info))
    logger = logging.getLogger()
    logger.exception("Unhandled exception: %s" % text)


if __name__ == '__main__':
    sys.excepthook = exceptionlogging
    main()
