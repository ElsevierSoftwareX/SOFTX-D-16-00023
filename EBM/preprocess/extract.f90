! This program calculates albedo values based on land masks and Legende polynomials
! See reference:
! Gerald R. North and James A. Coakley, Jr., 1979,
! Differences between Seasonal and Mean Annual Energy Balance Model Calculations
! of Climate and Climate Sensistivity.
! Journal of the Atmospheric Sciences
! Vol 36, pp 1189-1204.
! 

program extract
implicit none

integer,parameter:: nx=128, ny=65
real:: albedo(nx,ny), legende
integer:: i, j, alb
integer:: geography(nx,ny)

open(20,file='The_World.dat')

  do j = ny, 1, -1
    read (20,100) (geography(i,j),i=1,nx)
100  format (128I1)
  end do

! ================  newly added ==================================
do j=1,ny
   legende=0.5*(3*sin((90.0-(j-1)*2.8125)*3.1415926/180.0)*sin((90.0-(j-1)*2.8125)*3.1415926/180.0)-1)
do i=1,nx
   alb=geography(i,j)
  if(alb.eq.1)  albedo(i,j)=0.30 +0.09*legende
  if(alb.eq.2)  albedo(i,j)=0.60
  if(alb.eq.3)  albedo(i,j)=0.70 !0.68 
  if(alb.eq.5)  albedo(i,j)=0.29 +0.09*legende
end do
end do

open(21,file='albedo.dat')
do j=ny, 1, -1
  write(21,101) (albedo(i,j),i=1,nx)
101 format(128F10.2)
end do

close(21)
close(20)
end program
