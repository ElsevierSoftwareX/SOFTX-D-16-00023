! Read albedo values from the preprocessed into EBM main program

      subroutine albedo_input(Pcoalbedo)

      implicit none
      include 'netcdf.inc'

!     This is the name of the data file we will read.
      logical Albedo_2D
      character*(*) FILE_NAME
      parameter (FILE_NAME='../input/albedo.nc')
      integer ncid
      real error 
      parameter(error = 1.0e-10)

!     We are reading 4D data, a 128 x 65 x 1 x 48 lon-lat-lvl-timestep grid

      integer NDIMS 
      parameter (NDIMS = 4)
      integer NRECS, NLVLS, NLATS, NLONS
      parameter (NRECS = 48, NLVLS = 1, NLATS = 65, NLONS = 128)
      character*(*) LVL_NAME, LAT_NAME, LON_NAME, REC_NAME
      parameter (LVL_NAME = 'level')
      parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
      parameter (REC_NAME = 'time')
      integer lvl_dimid, lon_dimid, lat_dimid, rec_dimid

!     The start and count arrays to read netCDF library 
      integer start(NDIMS), count(NDIMS)
      data start /1,1,1,1/
      data count /NLONS,NLATS,NLVLS,NRECS/

!     latitude and longitude variables to check the input netcdf correct
      real lats(NLATS), lons(NLONS)
      integer lon_varid, lat_varid

!     We will read surface ALBEDO from the preprocessed 
      character*(*) ALBEDO_NAME
      parameter (ALBEDO_NAME='ALBEDO')
      integer ALBEDO_varid
      integer dimids(NDIMS)

!     variable attribute of units
      character*(*) UNITS
      parameter (UNITS = 'units')
      character*(*) ALBEDO_UNITS, LAT_UNITS, LON_UNITS
      parameter (ALBEDO_UNITS = 'dimensionless')
      parameter (LAT_UNITS = 'degrees_north')
      parameter (LON_UNITS = 'degrees_east')

!     ALBEDO to be read from the input netcdf
      real ALBEDO_in(NLONS, NLATS, NLVLS, NRECS),Pcoalbedo(NLONS, NLATS, NRECS)

!     Use these predefined grids to calculate the values we expect to find
      integer START_LAT, START_LON
      parameter (START_LAT = 90.0, START_LON = -180.0)

!     Loop indices.
      integer lvl, lat, lon, ts, rec, i

!     Error handling.
      integer retval

!     Open the file. 
      retval = nf_open(FILE_NAME, nf_nowrite, ncid)
      
!     Get the varids of the latitude and longitude coordinate variables.
      retval = nf_inq_varid(ncid, LAT_NAME, lat_varid)
      retval = nf_inq_varid(ncid, LON_NAME, lon_varid)

!     Read the latitude and longitude data.
      retval = nf_get_var_real(ncid, lat_varid, lats)
      retval = nf_get_var_real(ncid, lon_varid, lons)
      
!     Check to make sure we got what we expected 
!     The input netcdf should STRICTLY follow predefined lon-lat grid 
      do lat = 1, NLATS
         if (abs(lats(lat)-(START_LAT - (lat - 1) * 2.8125)) .gt. error) STOP 'Error: Wrong latitude'
      end do
      do lon = 1, NLONS
         if (abs(lons(lon)-(START_LON + lon * 2.8125)) .gt. error) STOP 'Error: Wrong longitude'
      end do

!     Get the varids of the pressure and ALBEDO netCDF variables.
      retval = nf_inq_varid(ncid, ALBEDO_NAME, ALBEDO_varid)

!     Read the surface ALBEDO data from the file
      retval = nf_get_vara_real(ncid, ALBEDO_varid, start, count, ALBEDO_in)
      do ts=1,NRECS
        do lat=1,NLATS
         do lon=1,NLONS
            Pcoalbedo(lon,lat,ts) = 1.0-ALBEDO_IN(lon,lat,1,ts)
         end do
        end do
      end do

!     Close the file
      retval = nf_close(ncid)

!     Now everything worked as expected
      print *,'SUCCESSLY reading albedo from [', FILE_NAME,'] into simulation!'
      write(2,*)'SUCCESSLY reading albedo from [', FILE_NAME,'] into simulation!'
      end


