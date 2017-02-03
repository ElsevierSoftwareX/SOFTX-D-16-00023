!-----------> EBM.inc - Include file for the EBM and FMG Solver <-------------

integer,parameter:: NX6=128, NY6=65, NT=48 

! NX6 is the number of grid points in longitude on the finest grid level 
! NY6 is the number of grid points in latitude on the finest grid level
! NT is the number of time steps per model year
  
real:: pi               ! pi = circumference/(2*radius)
real:: dt               ! time step interval
real:: dy               ! latitude increment
real:: Keq              ! coefficinet for diffusion at equator 
real:: Kocean           ! coefficient for ocean diffusion 
real:: Kland            ! coefficinet for land diffusion 
real:: KlandNP          ! coefficinet for land diffusion 
real:: KlandSP          ! coefficinet for land diffusion 
real:: B                ! radiation damping 

common /constants/pi, dt, dy, Kocean, Kland, KlandNP, KlandSP, Keq, B
!------------------------------------------------------------------------------
