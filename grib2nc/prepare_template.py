
import os
import sys
import logging
import argparse


import netCDF4 as nc4


XLONG_DICT = {'FieldType':104, 'MemoryOrder':"XY", 
              'stagger':"", 'units':'degree_east',
              'description':"LONGITUDE, WEST IS NEGATIVE"}
XLAT_DICT = {'FieldType':104, 'MemoryOrder':"XY", 
              'stagger':"", 'units':'degree_north',
              'description':"LATITUDE, SOUTH IS NEGATIVE"}

def main():
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', 
                        level=logging.WARNING)
    argparser = argparse.ArgumentParser(
        description='Program to make an empty netCDF file for grib conversion')
    argparser.add_argument('-v', '--verbose', help='Verbose logging', 
                           action='count')
    argparser.add_argument('--ncfilename', help='NetCDF filename', 
                           default='base.nc')
    argparser.add_argument('--nlats', help='Number of north/south coords',
                           type=int, default=1059)
    argparser.add_argument('--nlons', help='Number of west/east coords',
                           type=int, default=1799)
    

    args = argparser.parse_args()
    if args.verbose == 1:
        logging.getLogger().setLevel(logging.INFO)
    elif args.verbose > 1:
        logging.getLogger().setLevel(logging.DEBUG)

    make_netcdf(args.ncfilename, args.nlats, args.nlons)
        
        
def make_netcdf(ncfilename='base.nc', nlats=1059, nlons=1799):    
    netc = nc4.Dataset(ncfilename, 'w', format='NETCDF4')
    netc.createDimension('Time', None)
    netc.createDimension('time', 2)
    netc.createDimension('south_north', nlats)
    netc.createDimension('west_east', nlons)
    netc.createDimension('DateStrLen', 19)

    netc.createVariable('Times', 'S1', ('Time', 'DateStrLen'), 
                        zlib=True, complevel=1)
    lons = netc.createVariable('XLONG', 'f4', ('time', 'south_north', 'west_east'), 
                               zlib=True, complevel=1)
    lats = netc.createVariable('XLAT', 'f4', ('time', 'south_north', 'west_east'),
                        zlib=True, complevel=1)    
    lons.setncatts(XLONG_DICT)
    lats.setncatts(XLAT_DICT)
    netc.sync()
    netc.close()


if __name__ == '__main__':
    main()
