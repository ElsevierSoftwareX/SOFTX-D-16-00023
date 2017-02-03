! Output the simulated results into netcdf

      subroutine monthly_output

      implicit none
      include 'netcdf.inc'

!     This is the name of the data file that will be created.
      character*(*) FILE_NAME
      parameter (FILE_NAME = '../output/monthly-output.nc')
      integer ncid

      integer NDIMS 
      parameter (NDIMS = 4)
      integer NRECS, NLVS, NLATS, NLONS
      parameter (NRECS = 12, NLVS=1, NLATS = 65, NLONS = 128)
      character*(*) LVL_NAME, LAT_NAME, LON_NAME, REC_NAME
      parameter (LVL_NAME='level', LAT_NAME = 'latitude', LON_NAME = 'longitude')
      parameter (REC_NAME = 'time')
      integer lvl_dimid, lon_dimid, lat_dimid, rec_dimid
      
!     The start and count arrays will tell the netCDF library where to
!     write our data.
      integer start(NDIMS), count(NDIMS)
      data start /1,1,1,1/
      data count /NLONS,NLATS,NLVS,NRECS/

      real lats(NLATS), lons(NLONS)
      integer lon_varid, lat_varid

!     We will create the temperature field.
      character*(*) TEMP_NAME
      parameter (TEMP_NAME='temperature')
      integer temp_varid
      integer dimids(NDIMS)

!     We recommend that each variable carry a "units" attribute.
      character*(*) UNITS
      parameter (UNITS = 'units')
      character*(*) TEMP_UNITS, LAT_UNITS, LON_UNITS
      parameter (TEMP_UNITS = 'celsius')
      parameter (LAT_UNITS = 'degrees_north')
      parameter (LON_UNITS = 'degrees_east')

      real temp_out(NLONS,NLATS,NLVS,NRECS),temp(NLONS,NLATS)

!     Use these to construct some latitude and longitude data for this
!     example.
      integer START_LAT, START_LON
      parameter (START_LAT = 90.0, START_LON = -180.0)

!     Loop indices.
      integer lvl, lat, lon, rec, i

!     Error handling.
      integer retval
      character(len=3):: months(12) = (/'jan','feb','mar','apr','may','jun',  &
                                  'jul','aug','sep','oct','nov','dec'/)
      integer ts
      character::  tsnum*2, BINfile*30
      do ts = 1, NRECS
         BINfile='../output/'//months(ts)//'.bin'
         open (unit=1, file=BINfile, status='old')
         read (1,*) temp
         close (1,status='delete')
         temp_out(:,:,NLVS,ts) = temp(:,:)
      end do

!     Create real data to write.
      do lat = 1, NLATS
         lats(lat) = START_LAT - (lat - 1) * 2.8125
      end do
      do lon = 1, NLONS
         lons(lon) = START_LON + lon * 2.8125
      end do

!     Create the file. 
      retval = nf_create(FILE_NAME, nf_clobber, ncid)

!     Define the dimensions
      retval = nf_def_dim(ncid, LVL_NAME, NLVS, lvl_dimid)
      retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
      retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
      retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
  
!     Define the coordinate variables
      retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, lat_varid)
      retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, lon_varid)

!     Assign units attributes to coordinate variables.
      retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
      retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)

!     The dimids array 
      dimids(1) = lon_dimid
      dimids(2) = lat_dimid
      dimids(3) = lvl_dimid
      dimids(4) = rec_dimid
 
!     Define the netCDF variables for the temperature data.
      retval = nf_def_var(ncid, TEMP_NAME, NF_REAL, NDIMS, dimids, temp_varid)
      
!     Assign units attributes to the netCDF variables.
      retval = nf_put_att_text(ncid, temp_varid, UNITS, len(TEMP_UNITS), TEMP_UNITS)
 
!     End define mode.
      retval = nf_enddef(ncid)
      
!     Write the coordinate variable data. This will put the latitudes
!     and longitudes of our data grid into the netCDF file.
      retval = nf_put_var_real(ncid, lat_varid, lats)
      retval = nf_put_var_real(ncid, lon_varid, lons)
 
!     write temperature output
  
      retval = nf_put_vara_real(ncid, temp_varid, start, count, temp_out)

!     Close the file. This causes netCDF to flush all buffers and make
!     sure your data are really written to disk.
      retval = nf_close(ncid)
      print *,'Monthly temperature results stored successfully in ', FILE_NAME, '!'
      write(2,*)'Monthly temperature results stored successfully in ', FILE_NAME, '!'      
      
      end

