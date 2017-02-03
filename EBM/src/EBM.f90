!*************************************************************************
! This is the main program
! Final version: May 29, 2014 at Dickinson College

! A netCDF version of the Energy Balance Model 
!         using the full multigrid solver

!  EBM_netCDF is a free software.
!  You can redistribute it and/or modify it 
!  under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
! 
!  EBM_netCDF is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

!  For compilation, please see Makefile in the src subdirectory.
!  This code has been tested under Ubuntu 12.04LTS, CentOS6.5 and MacOSX
!  using intel fortran compiler and gfortran with netcdf3.6.3 and netcdf4.1. 

!  Users need to modify the netcdf library and include paths in Makefile 
!  to your own installed directories. 
!  For netcdf3.6.3 and below, the -lnetcdff flag is not needed and should be
!  removed from Makefile.

!***************************************************************************
PROGRAM EBM

IMPLICIT NONE
INCLUDE 'memory.inc'
INCLUDE 'ebm.inc'

! ======= variable explanations               ============================

! There are 48 time steps in the year. The astronomical calculations of the 
! Earth's orbit begin at the vernal equinox. This is used to calculate the TOA
! radiation. The first time step is at the vernal equinox. One month is then
! 4 time steps. Model years use the astronomical year, NOT the calendar year.
!
!       Vernal Equinox     1
!       Summer Solstice   13
!       Autumnal Equinox  25
!       Winter Solstice   37
!
!       Month     Time Steps
!        Jan    38, 39, 40, 41
!        Feb    42, 43, 44, 45
!        Mar    46, 47, 48   1
!        Apr     2,  3,  4,  5
!        May     6,  7,  8,  9
!        Jun    10, 11, 12, 13
!        Jul    14, 15, 16, 17
!        Aug    18, 19, 20, 21
!        Sep    22, 23, 24, 25
!        Oct    26, 27, 28, 29
!        Nov    30, 31, 32, 33
!        Dec    34, 35, 36, 37
! 
!  All model output is in binary files. The longitude-latitude grid is
!  128 x 65, with points at both poles. There grid interval is 2.8125
!  degrees. The arrangement of the data is:
!  
!  data(1,1) = 180W, 90N ...... data(128,1) = 177.1875E, 90N
!  ...                        ....
!  data(1,65) = 180W, 90S.......data(128,65) = 177.1875E, 90S

! NX6:  number of grid points in longitude on the finest grid level 
! NY6:  number of grid points in latitude on the finest grid level
! NT:  number of time steps per model year
! rhs(NX6,NY6):  right hand side of equation
! Temp(NX6,NY6):  surface temperatures 
!  Lastrhs(NX6,NY6):  last values of the RHS 
!  F(NX6,NY6,NT):  total radiative forcing 
!  SF(NX6,NY6,NT):  TOA solar radiative forcing 
!  Pcoalbedo(NX6,NY6,NT):  Planetary Coalbedos
!  HeatCap(NX6,NY6):  heat capacities
!  sum(NX6,NY6):  sum over time steps
!  latd(NY6):  latitude (degrees)
!  lond(NX6+1):  longitude (degrees)
!  cost(NY6):  cosine of colatitude
!  Global:  function names
!  GTemp:   global mean temperature
!  NHTemp:   NH mean temperature
!  SHTemp:   SH mean temperature
!  AnnTemp:   global annual mean temp
!  NHAnnTemp:   NH annual mean temp
!  SHAnnTemp:   SH annual mean temp
!  Last_AnnTemp:   last annual temperature
!  tau_land:   relaxation time for land
!  tau_snow:   relaxation time for snow
!  tau_sea_ice:   relaxation time for sea ice
!  tau_mixed_layer:   relaxation time for mixed layer
!  A:   OLR coefficient, A+BT
!  spinup:  model spin-up time
!  geography(NX6,NY6):  masks for the whole globe      
!  RelErr = 2.0e-5:  relative error between years
!  nx(NG), ny(NG):  grid points on each level
!  h(NG):  dx,dy on each grid level
!  geom, GCnp, GCsp:  geometry and terms at poles
!  ecc:  orbital elements, eccentricity
!  ob:                     obliquity 
!  per:                   long of perihelion
!  Solar_Constant:  mean value of solar constant S0 (W m^-2)
!  dt = 1.0/real(NT):  time step as fraction of a year (1/48)
!  dy = pi/real(NY6-1):  increment of latitude in radians
!  B = 2.15:  radiation damping coefficient, A+BT
!  KlandSP:  land diffusion coefficient at South Pole
!  KlandNP:  land diffusion coefficient at North Pole
!  Kland:  land diffusion coefficient at equator
!  Keq:  ocean diffusion coefficient at equator
!  Kocean:  ocean diffusion coefficient at poles 


!--------------------- general EBM variables -----------------------------
real::  rhs(NX6,NY6)                        
real::  Temp(NX6,NY6)                       
real::  Lastrhs(NX6,NY6)                    
real::  F(NX6,NY6,NT)                       
real::  SF(NX6,NY6,NT)                      
real::  Pcoalbedo(NX6,NY6,NT)               
real::  HeatCap(NX6,NY6)                     
real::  sum(NX6,NY6)                        
real::  latd(NY6)                            
real::  lond(NX6+1)                          
real::  cost(NY6)                            
real::  Global                              
real::  GTemp                                
                           
real::  AnnTemp                              
real::  Last_AnnTemp                        
real::  tau_land                             
real::  tau_snow                             
real::  tau_sea_ice                          
real::  tau_mixed_layer                      
real::  A                                    
real::  Monthly, ryear    
real::  tscount
real::  latitude
real::  solar0(NY6,NT)
real::  CO2ppm
integer:: initial_year

real:: SECNDS, secs, elapsed_time           
integer:: geography(NX6,NY6)                
integer:: i, j, l, n, yr, tstep, year, ts
integer::  Maxyrs, mcount,nf
logical:: Equilibrium 
character:: datafile*40, step*2, YearChar*4, RunType*5
character(len=3):: months(12) = (/'jan','feb','mar','apr','may','jun',  &
                                  'jul','aug','sep','oct','nov','dec'/)
                                  
!------------------ Set Some Parameters ----------------------------------
real,parameter::  RelErr = 2.0e-5   

!------------------  variables for FMG solver  ---------------------------
integer,dimension(NG):: nx, ny               
real:: h(NG)                                 
real,dimension(NG):: geom, GCnp, GCsp        
logical:: Converged 

!------------------  variables for orbital forcing  ----------------------
           
real::  ecc, ob, per                         

!------------------  variables for solar forcing  ------------------------
real,parameter:: Solar_Constant =1371.685 !100%  
!real,parameter:: Solar_Constant =1331.685 !3% less   
!real,parameter:: Solar_Constant = 1316.685  !4%   
!real,parameter:: Solar_Constant = 1344.685 !2%   
!real,parameter:: Solar_Constant =1357.685 !1%   

real::  S0(NT)                               
!-----------------  variables for mixed-layer heatflux  ------------------ 

!-----------------  Old Input Namelist for EBM  -------------------------- 
integer:: FirstYr 
logical:: Albedo_2D
                
Albedo_2D = .true.

! output model parameters and global temperature of each step         
open (unit=2, file='../output/Briefing.out', status='replace') 

!*************************************************************************    
!                        BEGIN EXECUTABLE STATEMENTS
!************************************************************************* 

 secs = SECNDS(0.0)                 

!------------------  Initialize common block parameters  -----------------    
pi = acos(-1.0)
dt = 1.0/real(NT)                
dy = pi/real(NY6-1)               
B = 2.15                         
KlandSP = 0.20                    
KlandNP = 0.28                    
Kland   = 0.65                   
Keq     = 0.65                    
Kocean  = 0.40 

!------------------------  Initialize some variables  --------------------         

CO2ppm=315.0 !AD1950
initial_year=1950
!CO2ppm=315.0 !9kaBP
!initial_year=-9000
!CO2ppm=315.0 !21kaBP
!initial_year=-21000

Maxyrs = 100
FirstYr = Maxyrs

!------------------------  Initialize some arrays  ----------------------- 
do ts = 1, NT
  S0(ts) = Solar_Constant                    
end do
do j = 1, NY6
  latd(j) = (pi/2.0 - dy*real(j-1))*180./pi    
  cost(j) = cos(dy*real(j-1))                   
end do
do i = 1, NX6+1
  lond(i) = dy*real(i-1)*180./pi               
end do  

!------------------ Read geography and coalbedos from input --------------
CALL geography_input(geography)
CALL albedo_input(Pcoalbedo)  

!-------------------  Calculate the heat capacities  --------------------- 
CALL HeatCapacities (HeatCap, geography, tau_land, tau_snow, tau_sea_ice, &
                     tau_mixed_layer, .FALSE.)

!---------------------- Setup for the FMG Solver  ------------------------
CALL FMG_Setup (nx, ny, h, geom, Heatcap, geography, GCnp, GCsp) 

!------------  Calculate the background planetary coalbedos  -------------
!CALL Coalbedos (Albedo_2D, Pcoalbedo)
 
!--------  ALLOCATE AND INITIALIZE ARRAYS FOR THE CLIMATE FORCINGS  ------ 
!---------------  Calculate the background planetary coalbedos  ---------- 
 
!--------  ALLOCATE AND INITIALIZE ARRAYS FOR THE CLIMATE FORCINGS  ------ 
    
!------  change this value according to greenhouse gas concentrations ---- 
!------  Depending upon Myrhe et al. (1998)                           ----      

call A_value(CO2ppm, A)

!-------------------------- Write input parameters  ---------------------- 
     
CALL WriteParameters (tau_land, tau_snow, tau_sea_ice, tau_mixed_layer, A, &
           Albedo_2D)  
!-----------------------  Get the orbital elements  ----------------------
call orbital_params(initial_year,ecc,ob,per)

! HERE THERE IS AN OPTION TO DIRECTLY ASSIGN ORBITAL PARAMETERS
! BY COMMENTING THE ABOVE LINE "call orbital_params(initial_year, ecc,ob,per)"
! AND DIRECTLY ASSIGN THE ORBITS BELOW

!ecc = 0.019456 !18k            
!ob  = 0.409153 
!per = 2.853565 

!ecc = 0.016740 !1950AD            
!ob  = 0.409253
!per = 1.783037  

!------------------  Initailize the TOA solar forcing  ------------------- 
CALL Solar_Forcing (0, .false., S0, .false., Pcoalbedo, A, ecc, ob, per, SF)

!---- Initialize the temperature field from Legates and Willmott data ---- 
CALL Initial_Temp (Temp)            

year = 1

sum = 0.0
rhs = 0.0
LastRhs = 0.0
GTemp = 0.0                   

write (*,23)
write (2,23)
23 format(/,' Start Energy Balance Model simulation:')

write (*,24)
write (2,24)
24 format(/,'  Year   Global Temperature')
!---------------------- START LOOP OVER MODEL YEARS  --------------------- 
Equilibrium = .false. 
DO yr = 1, Maxyrs            

  F = SF                             
   
!--------------------  START LOOP OVER MODEL TIME STEPS  ----------------- 
  DO tstep = 1, NT             
     
    Converged = .false.
    if (tstep == 38 .and. Equilibrium) then      
      write (*,58) yr-1, AnnTemp                            
      write (2,58) yr-1, AnnTemp
58    format (/,'EQUILIBRIUM REACHED AFTER ',i3,' YEARS. GLOBAL TEMP= ',f7.4)
      goto 80                       
    else if (tstep == 38) then     
      tscount = 1.0
      mcount = 1
      AnnTemp = GTemp/real(NT)               
         
      if (abs(AnnTemp-Last_AnnTemp) <= RelErr) then
        Equilibrium = .true.
      end if
      if (year < FirstYr) then
        write (*,60) yr, AnnTemp
      else
        if (year==FirstYr) print *,'Begin writing output data'
        if (year>=FirstYr) then
 
          write (*,60) year-FirstYr+1, Anntemp
        else
          write (*,60) year, Anntemp
        end if
      end if
      write (2,60) year, AnnTemp
60    format (1x,i4,2x,f10.5)
      year = year + 1
      Last_AnnTemp = AnnTemp
      GTemp  = 0.0                
         
    end if
   
    CALL UpdateRHS (tstep, HeatCap, Temp, F, rhs, LastRhs)
      
!    CALL FMG_Solver (nx, ny, h, geom, GCnp, GCsp, Converged, rhs, Temp, .FALSE.)
      CALL FMG_Solver (nx, ny, h, geom, GCnp, GCsp, Converged, rhs, Temp, .FALSE.)

    if (.NOT.Converged) then
        write (*,61) 
        write (*,59) yr, tstep 
59    format ('Year =',i7,'Time Step = ',i2)
61    format ('NO Convergence within max number of V cycles')
      STOP
    end if

!---------------------    run time step data  ----------------------------
    if (yr==Maxyrs.or.Equilibrium) then     
      write (step, '(i2)') tstep 
      if (tstep < 10) step(1:1) = '0'
      datafile = '../output/t'//step//'.bin'      
      open (unit=7, file=datafile, status='replace')            
      write (7,*) Temp
      close(7)

      sum = sum + 0.25*Temp         

      if (mod(tscount,4.0) == 0.0) then
        datafile = '../output/'//months(mcount)//'.bin'
        open (unit=7, file=datafile, status='replace')            
        write (7,*) sum
        close(7)
        Monthly = Global (sum, z(iaw(NG)), 0)
        ryear = real(year) + 0.041667 + 0.083333*real(mcount-1)
        mcount = mcount + 1
        sum = 0.0            
      end if
    end if
    GTemp =  GTemp  + Global (Temp, z(iaw(NG)), 0)     
    tscount = tscount + 1.0           
  
  END DO                       
END DO                        


 80 elapsed_time = SECNDS(secs)      

if (.NOT.Equilibrium) then
  write (*,82) Maxyrs, AnnTemp
  write (2,82) Maxyrs, AnnTemp
82 format (/,' WARNING: EQUILIBRIUM NOT REACHED AFTER ',i3,' YEARS!!!!', &
          /,'         GLOBAL TEMPERATURE = ',f8.4)
end if

write (*,90) elapsed_time/60.
write (2,90) elapsed_time/60. 
90 format(/,' Elapsed time:',f10.2,' minutes') 

call monthly_output
call timesteps_output

close (2)

END PROGRAM EBM


