! program name: prepare_albedo.f90
! purpose: turn ascii albedo into NetCDF
! input: albedo.dat
! output: albedo.nc

! Be Aware!
! This is only to the annual mean albedo
! for T48 albedo need modify a little bit

! Compile
!ifort prepare_albedo.f90 -o prepare_albedo -I/home/kelin/netcdf/include -L/home/kelin/netcdf/lib -lnetcdff -lnetcdf

      program monthly_output
      implicit none
      include 'netcdf.inc'

!     This is the name of the data file we will create.
      character*(*) FILE_NAME
      parameter (FILE_NAME = 'albedo.nc')
      integer ncid

      integer NDIMS 
      parameter (NDIMS = 4)
      integer NRECS, NLVS, NLATS, NLONS
      parameter (NRECS = 48, NLVS=1, NLATS = 65, NLONS = 128)
      character*(*) LVL_NAME, LAT_NAME, LON_NAME, REC_NAME
      parameter (LVL_NAME='level', LAT_NAME = 'latitude', LON_NAME = 'longitude')
      parameter (REC_NAME = 'time')
      integer lvl_dimid, lon_dimid, lat_dimid, rec_dimid
      
!     The start and count arrays will tell the netCDF library where to
!     write our data.
      integer start(NDIMS), count(NDIMS)
      data start /1,1,1,1/
      data count /NLONS,NLATS,NLVS,NRECS/

!     These program variables hold the latitudes and longitudes.
      real lats(NLATS), lons(NLONS)
      integer lon_varid, lat_varid

!     We will create the albedo field.
      character*(*) ALBEDO_NAME
      parameter (ALBEDO_NAME='ALBEDO')
      integer ALBEDO_varid
      integer dimids(NDIMS)

!     We recommend that each variable carry a "units" attribute.
      character*(*) UNITS
      parameter (UNITS = 'units')
      character*(*) ALBEDO_UNITS, LAT_UNITS, LON_UNITS
      parameter (ALBEDO_UNITS = 'dimensionless')
      parameter (LAT_UNITS = 'degrees_north')
      parameter (LON_UNITS = 'degrees_east')

!     Program variables to hold the data we will write out. We will only
!     need enough space to hold one timestep of data; one record.
!      real ALBEDO_out(NLONS, NLATS,NRECS)
      real ALBEDO_out(NLONS,NLATS,NLVS,NRECS),ALBEDO(NLONS,NLATS,NRECS)

!     Use these to construct some latitude and longitude data for this
!     example.
!      integer START_LAT, START_LON
      real START_LAT, START_LON,longitude
      parameter (START_LAT = 90.0, START_LON = -180.0)

!     Loop indices.
      integer lvl, lat, lon, rec, i

!     Error handling.
      integer retval
      
      integer ts
      character::  tsnum*2, BINfile*30

      open(10,file='albedo.dat',status='old')
      
! Annual mean albedo    
      do lat=1,65
           Read(10,101) (albedo(lon,lat,1),lon=1,128)
      end do      

      do ts=2,48
        do lat=1,65
        do lon=1,128
             albedo(lon,lat,ts)=albedo(lon,lat,1)
        end do
        end do
      end do
     
101 format(128F10.2)
    close(10)

      do ts=1,NRECS
        do lat=1,NLATS
         do lon=1,NLONS
            ALBEDO_OUT(lon,lat,1,ts) = albedo(lon,lat,ts)
         end do
        end do
      end do

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
 
!     Define the netCDF variables for the albedo data.
      retval = nf_def_var(ncid, ALBEDO_NAME, NF_REAL, NDIMS, dimids, ALBEDO_varid)
      
!     Assign units attributes to the netCDF variables.
      retval = nf_put_att_text(ncid, ALBEDO_varid, UNITS, len(ALBEDO_UNITS), ALBEDO_UNITS)
 
!     End define mode.
      retval = nf_enddef(ncid)
      
!     Write the coordinate variable data. This will put the latitudes
!     and longitudes of our data grid into the netCDF file.
      retval = nf_put_var_real(ncid, lat_varid, lats)
      retval = nf_put_var_real(ncid, lon_varid, lons)
 
!     write albedo output
  
      retval = nf_put_vara_real(ncid, ALBEDO_varid, start, count, ALBEDO_out)

!     Close the file. This causes netCDF to flush all buffers and make
!     sure your data are really written to disk.
      retval = nf_close(ncid)
      print *,'SUCCESS in writing albedo to ', FILE_NAME, '!'
      print *,'Preprocess DONE! '
      print *, '-------------------------------------------'
      print *     
      end

