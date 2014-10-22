#!/usr/bin/env python
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
import fcntl
import argparse
import warnings
try:
    import ConfigParser as configparser
except ImportError:
    import configparser


import pandas as pd
import numpy as np
try:
    import concurrent.futures as cf
except ImportError:
    import futures as cf


class RequestError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class HRRRFetcher(object):
    """Fetch the specified HRRR forecasts from the NCEP servers
    
    Parameters
    ----------
    init_time : str or datetime
        Initilization time of the forecast you want
    level : str
        "Level" of the HRRR forecast i.e. surface, subhourly, native, or pressure
    method : str
        Method/protocol to download files. Currently ftp or http
    config_path : str
        Path to the config file
    """
    def __init__(self, init_time, level, method='ftp', config_path=None, 
                 threads=None):
        self.logger = logging.getLogger('HRRRFetcher')
        self.config = configparser.ConfigParser()
        self.config_path = config_path or os.path.join(
            os.path.dirname(os.path.realpath('__file__')),'..'
            'settings.txt')
        self.config.read(self.config_path)
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

        self.base_path = self.download_dict['grib_base_folder'].format( 
            year=self.init_time.strftime('%Y'), month=self.init_time.strftime('%m'),
            day=self.init_time.strftime('%d'), hour=self.init_time.strftime('%H'))
            
        if not os.path.isdir(self.base_path):
            os.makedirs(self.base_path)
        self.downloaded_files = []

        self.nthreads = threads or self.config.getint(
            'download_settings', 'workers')

    def _http_setup(self, init_time, level, forecast_name):
        try:
            global requests
            import requests
        except ImportError:
            raise ImportError('HTTP requires requests')
        else:
            logging.getLogger('requests.packages.urllib3'
                              ).setLevel(logging.WARNING)

        html_folder = (self.download_dict['html_site'] + 
            'hrrr.{dtime}'.format(dtime=init_time.strftime('%Y%m%d')))

        forecast_table = requests.get(html_folder)

        if forecast_table.status_code != 200:
            raise RequestError('Invalid URL/date')
            
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
        files = files_df['Name']
        return (files, expected_size)

    def _ftp_setup(self, init_time, level, forecast_name):
        try:
            global ftplib
            import ftplib
        except ImportError:
            raise ImportError('FTP fetching requires ftplib')

        self.logger.info('Connecting to FTP site')
        ftp = self.connect_ftp(init_time)
        files = ftp.nlst(forecast_name + '*')
        expected_size = 0
        for afile in files:
            expected_size += ftp.size(afile)
        return (files, expected_size)

    def connect_ftp(self, init_time=None):
        """Connect to the FTP site and change directories
        """        
        ftp = ftplib.FTP(self.download_dict['ftp_host'], 
                         timeout=float(self.download_dict['ftp_timeout']))
        ftp.login()
        ftp.cwd(self.download_dict['ftp_dir'])
        stored_days = ftp.nlst()
        init_time = init_time or self.init_time
        it_date = init_time.strftime('hrrr.%Y%m%d')
        if it_date not in stored_days:
            raise RequestError('Requested forecast initilization day '+
                               'not on NCEP server')
        ftp.cwd(it_date)
        return ftp

    def fetch(self, init_time=None, level=None, overwrite=True):
        """Fetch the forecasts useing a threaded fetcher
        """

        level = level or self.level
        init_time = init_time or self.init_time
        forecast_name ='hrrr.t{init_hour}z.{level}f'.format(
            init_hour=init_time.strftime('%H'),
            level=self.hrrr_type_dict[level])

        html_folder = (self.download_dict['html_site'] + 
            'hrrr.{dtime}'.format(dtime=init_time.strftime('%Y%m%d')))

        def _ftp_fetch(filename, local_path):
            size = []
            ftp = self.connect_ftp(init_time)
            self.logger.debug('Retrieving %s' % filename)
            with open(local_path, 'wb') as f:
                def callback(data):
                    f.write(data)
                    size.append(len(data))
                ftp.retrbinary('RETR %s' % filename, callback)
            return sum(size)

        def _http_fetch(filename, local_path):
            size = 0
            self.logger.debug('Retrieving %s' % filename)
            r = requests.get(html_folder + '/' + filename, stream=True)
            with open(local_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=1024):
                    if chunk:
                        f.write(chunk)
                        f.flush()
                        size += len(chunk)
            return size

        if self.default_method == 'ftp':
            files, expected_size = self._ftp_setup(
                init_time, level, forecast_name)
            fetcher = _ftp_fetch
        elif self.default_method == 'http':
            files, expected_size = self._http_setup(
                init_time, level, forecast_name)
            fetcher = _http_fetch
        else:
            raise AttributeError('method must be ftp or http not %s' % 
                                 self.default_method)
        
        self.logger.info('Downloading %0.2f MB of data for forecast %s' 
                         % ((1.0*expected_size/1024**2), forecast_name))

        total_size = 0
        nfiles = 0
        # threading loop
        with cf.ThreadPoolExecutor(max_workers=self.nthreads) as executor:
            worker_dict = {}
            for filename in files:
                local_path = os.path.join(self.base_path,  filename)
                if os.path.isfile(local_path) and not overwrite:
                    self.logger.warning('%s already exists!' % filename)
                    continue
                wk = executor.submit(fetcher, filename, local_path)
                worker_dict[wk] = (filename, local_path)

            requeue = []
            start = time.time()
            for future in cf.as_completed(worker_dict):
                thefile, local_path = worker_dict[future]

                try:
                    size = future.result()
                except Exception as e:
                    self.logger.exception('Exception when fetching %s' % thefile)
                    requeue.append((thefile, local_path))
                else:
                    self.downloaded_files.append(local_path)
                    total_size += size
                    nfiles += 1

        if requeue:
            self.logger.warning('Must try to refetch some files')
            for filename, local_path in requeue:
                self.logger.info('Trying to refetch %s' % filename)
                try:
                    size = fetcher(filename, local_path)
                except Exception as e:
                    self.logger.exception('Failed to retrieve %s' % filename)
                else:
                    self.downloaded_files.append(local_path)
                    total_size += size
                    nfiles += 1
                    
        end = time.time()
        self.logger.info('Retrieved %s files in %0.1f seconds' % 
                         (nfiles, end-start))
        self.logger.info('Approx. download rate: %0.2f MB/s' % 
                         (1.0*total_size/1024**2/(end-start)))
        self.logger.debug('Files retrieved are %s' % self.downloaded_files)


def check_availability(config_path, levels):
    """Check the HRRR availability website and compare
    to already archived forecasts to fetch new forecasts
    """
    
    config = configparser.ConfigParser()
    config.read(config_path)
    avail_sites = config.get('download_settings', 'availability_sites')
    dfs = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for site in avail_sites.split(','):
            dfs.append(pd.io.html.read_html(site, header=0,
                                            flavor='bs4')[0])
    avail_df = pd.concat(dfs, ignore_index=True)
    utcdate = dt.datetime.utcnow().date()
    newi = [dt.datetime.combine(utcdate, 
                                dt.datetime.strptime(e, '%H UTC HRRR').time()) 
            for e in avail_df['EVENT']]
    newi[-1] = newi[-1] - dt.timedelta(days=1)
    avail_df.index = pd.Index(newi)

    invalid = [i for i, val in enumerate(avail_df['STATUS']) if not "COMPLETE" in val]
    avail_df.iloc[invalid] = np.nan
    avail_df = avail_df.dropna()

    netcdf_folder = config.get('download_settings', 'netcdf_base_folder')
    netcdf_filename = config.get('download_settings', 'netcdf_filename')

    lsers = {}
    for level in levels:
        lsers[level] = pd.Series(np.zeros(len(avail_df.index)), 
                                 index=avail_df.index)
    
    for atime in avail_df.index:
        netcdf_path = netcdf_folder.format(year=atime.strftime('%Y'),
                                           month=atime.strftime('%m'),
                                           day=atime.strftime('%d'),
                                           hour=atime.strftime('%H'))
        if os.path.isdir(netcdf_path):
            for level in levels:
                ncfilename = netcdf_filename.format(level=level, 
                    init_time=atime.strftime('%Y%m%d%H'))
                if os.path.isfile(os.path.join(netcdf_path, ncfilename)):
                    lsers[level].loc[atime] = 1
                

    avail_df = avail_df.join(pd.DataFrame(lsers))
        
    to_get = {}
    for level in levels:
        to_get[level] = sorted(avail_df[(avail_df[level] == 0)
                                    ].index.to_pydatetime())
        
    return to_get


class Locker(object):
    """Places a file lock so only one instance of the program
    runs at a time
    """
    def __init__(self, lock_path='/tmp/fetch_forecasts.lock'):
        self.lock_path = lock_path
        try:
            self.locked = open(self.lock_path, 'w+')
            fcntl.lockf(self.locked, fcntl.LOCK_EX|fcntl.LOCK_NB)
        except IOError:
            logging.debug('File (%s) locked. Program probably running.' % 
                          self.lock_path)
            sys.exit(0)
    def unlock(self):
        fcntl.lockf(self.locked, fcntl.LOCK_UN)
        os.remove(self.lock_path)


def main():
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s')

    argparser = argparse.ArgumentParser(
        description='Download HRRR grib files and store select fields in netCDF format')
    
    argparser.add_argument('-v', '--verbose', help='Increase logging verbosity',
                           action='count')
    argparser.add_argument('-p', '--protocol', help='Download protocol',
                           choices=['ftp', 'http'], default='ftp')
    argparser.add_argument('-c', '--config', help='Config file path',
                           default='../settings.txt')
    argparser.add_argument('-l', '--levels', default='subhourly',
                           help='HRRR level subtype; split multiple with commas')
    argparser.add_argument('-r', '--remove', action='store_true',
                           help='Remove grib files after conversion')
    argparser.add_argument('-t', '--threads', type=int,
                           help='Number of threads to use to fetch data')
    argparser.add_argument('--no-convert', action='store_true',
                           help="Don't convert the grib files to netCDF")
    argparser.add_argument('--no-overwrite', action='store_false',
                           help="Don't overwrite any files already downloaded")
    argparser.add_argument('--all', action='store_true',
                           help="Get all the HRRR files not already archived")
    argparser.add_argument('-i', '--init-times', 
        help='Initilization datetimes; split multiple with commas')
    args = argparser.parse_args()
    
    if args.verbose == 1:
        logging.getLogger().setLevel(logging.INFO)
    elif args.verbose > 1:
        logging.getLogger().setLevel(logging.DEBUG)

    global lock
    lock = Locker()

    def loop(atime, alevel, args):
        f = HRRRFetcher(atime, alevel, args.protocol, 
                        args.config, args.threads)
        f.fetch(overwrite=args.no_overwrite)

        if not args.no_convert:
            from grib2nc import Grib2NC
            g2nc = Grib2NC(f.init_time, f.level, f.config_path)
            g2nc.load_into_netcdf()
            del g2nc

        if args.remove and not args.no_convert:
            for apath in f.downloaded_files:
                os.remove(apath)
                if not os.listdir(os.path.dirname(apath)):
                    os.rmdir(os.path.dirname(apath))


    if args.all:
        files_to_get = check_availability(args.config, args.levels.split(','))
        for alevel, thetimes in files_to_get.items():
            for atime in thetimes:
                loop(atime, alevel, args)
    elif args.init_times is not None:
        for atime in args.init_times.split(','):
            for alevel in args.levels.split(','):
                loop(atime,alevel, args)
    else:
        argparser.error('Either --all or -i/--init-times must be specified')
    lock.unlock()


def exceptionlogging(*exc_info):
    text = "".join(traceback.format_exception(*exc_info))
    logger = logging.getLogger()
    logger.exception("Unhandled exception: %s" % text)
    lock.unlock()


if __name__ == '__main__':
    sys.settrace
    sys.excepthook = exceptionlogging
    main()
