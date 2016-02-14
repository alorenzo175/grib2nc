import sys
from distutils.core import setup


import versioneer


PACKAGE = 'grib2nc'


if sys.version_info[:2] < (2, 7):
    sys.exit('%s requires Python 2.7 or higher.' % PACKAGE)


SHORT_DESC = 'Convert grib2 files to netCDF format'
AUTHOR = 'Tony Lorenzo'
MAINTAINER_EMAIL = 'alorenzo175@gmail.com'
URL = 'https://github.com/alorenzo175/grib2nc'


setup(
    name=PACKAGE,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=SHORT_DESC,
    author=AUTHOR,
    maintainer_email=MAINTAINER_EMAIL,
    url=URL,
    packages={'grib2nc': 'grib2nc'})
