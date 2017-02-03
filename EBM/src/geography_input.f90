! Read land masks from the preprocess into the EBM main code

      subroutine geography_input(geography)
      implicit none
      include 'netcdf.inc'

!     This is the name of the data file we will create.
      character*(*) FILE_NAME
      parameter (FILE_NAME = '../input/geography.nc')
      integer ncid
      real error 
      parameter(error = 1.0e-10)

!     We are reading a geography map with spatial resolutions of 2.8125 degrees both longitudinal and latitudinal
!     So we will have a grid of 128 longitude X 65 latitude
      integer NDIMS
      parameter (NDIMS = 2)
      integer NLATS, NLONS
      parameter (NLATS = 65, NLONS = 128)
      character*(*) LAT_NAME, LON_NAME
      parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
      integer lon_dimid, lat_dimid

      integer start(NDIMS), count(NDIMS)
      data start /1,1/
      data count /NLONS,NLATS/

      real lats(NLATS), lons(NLONS)
      integer lon_varid, lat_varid

      character*(*) GEOG_NAME
      parameter (geog_NAME='landmask')
      integer GEOG_varid
      integer dimids(NDIMS)

!     We define each variable a "units" attribute.
      character*(*) UNITS
      parameter (UNITS = 'units')
      character*(*) GEOG_UNITS, LAT_UNITS, LON_UNITS
      parameter (GEOG_UNITS = 'MASK: 1. land;  2. sea ice; 3. land ice; 5 seawater. ')
      parameter (LAT_UNITS = 'degrees_north')
      parameter (LON_UNITS = 'degrees_east')

      integer geography(NLONS, NLATS)

!     Our grid starts from uppler right corner at longitude 0 degree and latitude 90 degree (0,90)
      integer START_LAT, START_LON
      parameter (START_LAT = 90.0, START_LON = -180.0)

!     Loop indices.
      integer lat, lon, i

!     Error handling.
      integer retval


!     Create grid in real longitude and latitude.
      do lat = 1, NLATS
         lats(lat) = START_LAT - (lat - 1) * 2.8125
      end do
      do lon = 1, NLONS
         lons(lon) = START_LON + lon * 2.8125
      end do

!     Open the file. 
      retval = nf_open(FILE_NAME, nf_nowrite, ncid)

!     Get the varids of the latitude and longitude coordinate variables.
      retval = nf_inq_varid(ncid, LAT_NAME, lat_varid)
      retval = nf_inq_varid(ncid, LON_NAME, lon_varid)
  
!     Read the latitude and longitude data.
      retval = nf_get_var_real(ncid, lat_varid, lats)
      retval = nf_get_var_real(ncid, lon_varid, lons)

!     Assign units attributes to coordinate variables.
      retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
      retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
  
!     Check to make sure we got what we expected 
!     The input netcdf should STRICTLY follow predefined lon-lat grid 
      do lat = 1, NLATS
         if (abs(lats(lat)-(START_LAT - (lat - 1) * 2.8125)) .gt. error) STOP 'Error: Wrong latitude'
      end do
      do lon = 1, NLONS
         if (abs(lons(lon)-(START_LON + lon  * 2.8125)) .gt. error) STOP 'Error: Wrong longitude'
      end do

!     Get the varids of the geography netCDF variables.
      retval = nf_inq_varid(ncid, GEOG_NAME, GEOG_varid)
      
!     Read the surface ALBEDO data from the file
      retval = nf_get_vara_int(ncid, GEOG_varid, start, count, geography)
 
!     Close the dat file and make sure your data are really written to disk.
      retval = nf_close(ncid)
      
      print *,'SUCCESSLY reading geography from [', FILE_NAME,'] into simulation!'
      write(2,*)'SUCCESSLY reading geography from [', FILE_NAME,'] into simulation!'      
      end

