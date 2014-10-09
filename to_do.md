HRRR [website](http://www.nco.ncep.noaa.gov/pmb/products/hrrr/)

1) Only partial grib file retrieval *DOESN'T seem to actually work
 a) load .idx files to determine ranges; see [here](http://nomads.ncep.noaa.gov/txt_descriptions/fast_downloading_grib.shtml)
  i) read into pandas dataframe for manipulation
  ii) ~~perhaps not to limit dependencies~~
 b) HTTP 'Range' header
 c) FTP REST and abort
 d) OPeNDAP
2) Conversion/combination
 a) looks like [regional subsetting is still required](http://www.cpc.ncep.noaa.gov/products/wesley/fast_downloading_grib.html)
 b) netCDF4-classic probably safer? 
 c) test compression

