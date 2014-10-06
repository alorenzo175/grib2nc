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


class Fetcher(object):
    def __init__(self, method='ftp'):
        self.logger = logging.getLogger()
        self.config = configparser.ConfigParser()
        self.config.read(os.path.join(
            os.path.dirname(os.path.realpath('__file__')),
            'settings.txt'))
        self.ftp_dict = dict(self.config.items('download_settings'))
        self.hrrr_type_dict = dict(self.config.items('output_types'))

    def fetch_ftp(self, init_time, level):
        if level not in self.hrrr_type_dict:
            raise Exception
        ftp = ftplib.FTP(self.ftp_dict['ftp_host'])
        ftp.login()
        ftp.cwd(self.ftp_dict['ftp_dir'])
        stored_days = ftp.nlst()
        it_date = init_time.strftime('hrrr.%Y%m%d')
        if it_date not in stored_days:
            raise Exception('l')
        ftp.cwd(it_date)
        base_name = 'hrrr.t{init_hour}z.{level}f*.grib2*'.format(
            init_hour=init_time.strftime('%H'), 
            level=self.hrrr_type_dict[level])
        files = ftp.nlst(base_name)
        total_size = 0
        for afile in files:
            total_size += ftp.size(afile)
        
        self.logger.info('Attempting to retrieve %s data' % total_size)
        base_path = os.path.join(self.ftp_dict['folder'], 
                                 init_time.strftime('%Y'),
                                 init_time.strftime('%m'),
                                 init_time.strftime('%d'),
                                 init_time.strftime('%Hz'),
                                 'grib')
        if not os.path.isdir(base_path):
            os.makedirs(base_path)

        start = time.time()
        for afile in files:
            local_path = os.path.join(base_path,  afile)
            with open(local_path, 'wb') as f:
                def callback(data):
                    f.write(data)
                ftp.retrbinary('RETR %s' % afile, callback)
            break
        end = time.time()
        self.logger.info('Retrieved %s files in %s seconds' % (len(files), 
                                                               end-start))
        

def main():
    logging.basicConfig(level=logging.INFO)
    f = Fetcher()
    f.fetch_ftp(dt.datetime(2014,10,2,1), 'surface')


def exceptionlogging(*exc_info):
    text = "".join(traceback.format_exception(*exc_info))
    logger = logging.getLogger()
    logger.exception("Unhandled exception: %s" % text)


if __name__ == '__main__':
    sys.excepthook = exceptionlogging
    main()
