! Here list procedures to materialize the FMG in EBM
! See references:
!Huang, J., Bowman, K.P., 1992. The small ice cap instability in seasonal energy balance models. 
!Climate Dynamics 7, 205-215. 
!Bowman, K. and J. Huang, 1991: A multigrid solver for the Helmholtz equation 
!on a semiregular grid on the sphere. Monthly Weather Review, 119, 769-775.
!Stevens, M.J., and G.R. North, 1996: Detection of the climate response 
!to the solar cycle. Journal of Atmospheric Sciences, 53, 2594-2608.

! addcor.f90
      
subroutine AddCorrection (nxc, nyc, uc, area, nxf, nyf, res, uf, Debug)
                               

implicit none
include 'ebm.inc'

integer nxc, nyc, nxf, nyf, i, j
real    uc(nxc,nyc), uf(nxf,nyf), res(nxf,nyf), area(nyc)
logical Debug
!-----------------------------------------------------------------------
! Do interpolation (prolongation) of the correction from the coarse grid
! to the fine grid and then add the result to the solution. 

      call Interpolation (nxc, nyc, area, uc, nxf, nyf, res, debug)

      do j = 1, nyf
        do i = 1, nxf
          uf(i,j) = uf(i,j) + res(i,j)        
        end do
      end do


      if (Debug) then
        write (4,2)
 2      format (/,'ADD CORRECTION to solution')
      end if

      end

! part 3
! coarse.f90
subroutine CoarseSol (nx, ny, h, geom, GCnp, GCsp, rhs, u, Debug)

implicit none
include 'memory.inc'
include 'ebm.inc'

integer nx, ny
real    rhs(nx,ny), u(nx,ny)
real    GCnp, GCsp, h, geom
integer itmax, n
parameter (itmax = 3)
logical  Debug
!-------------------------------------------------------------------------
! Calculate the initial solution on the coarsest grid (4 x 3).

      if (Debug) then
!        write (4,2) nx, ny
 2      format (/,'COARSE SOLUTION for nx = ',i1,' ny = ',i1)
      end if

! Initialize solution to zero
      u = 0.0

      do n = 1, itmax
        call XRelaxation (nx, ny, h, u, rhs, geom, z(ic0(1)), &
        z(ic1(1)), z(ic2(1)), z(ic3(1)), z(ic4(1)), GCnp, GCsp, z(iDnp(1)), &
        z(iDsp(1)), Debug)

      end do 

      end

! part 4
! constb.f90
!****************************************************************************

  subroutine ConstB (nx, ny, h, Heatcap, geom, MDnp, MDsp, &
    			 GCnp, GCsp, DiffC, csc2, c0) 

      implicit none
      include 'ebm.inc'

      integer  nx, ny, i, j
      real   Heatcap(nx,ny), DiffC(nx,ny)
      real   c0(nx,ny), csc2(ny), h, hh
      real   geom, GCnp, GCsp, MDnp, MDsp
!-------------------------------------------------------------------------
! Initialize for the case of time-independent B
      hh = h*h
      
      GCnp = geom * MDnp + B + 2.*Heatcap(1,1)/dt     ! north pole 
      GCsp = geom * MDsp + B + 2.*Heatcap(1,ny)/dt    ! south pole

! Calculate differential operator which depends on B     
      do j = 2, ny-1
        do i = 1, nx
          c0(i,j) = 2.*DiffC(i,j)*(1.+csc2(j)) + hh*(2.*Heatcap(i,j)/dt + B)
        end do
      end do

      end

! part 5
! constdif.f90
!******************************************************************************

      subroutine ConstDiffOper (nx, ny, h, DiffC, csc2, cot, c1, c2, c3, c4)

      implicit none

      integer nx, ny, i, j 
      real    csc2(ny), cot(ny)
      real    DiffC(nx,ny), c1(nx,ny), c2(nx,ny), c3(nx,ny), c4(nx,ny)
      real    h, h2, dTheta, dPhi
!--------------------------------------------------------------------------
! Calculate the coefficients of the diffusion operator at each colatitude and
! longitude on the sphere for the current grid level.

      h2 = 2.0*h
      
      do j = 2, ny-1
        do i = 1, nx
          if (i .eq. 1) then
            dPhi = (DiffC(2,j) - DiffC(nx,j))/h2
          else if (i .eq. nx) then
            dPhi = (DiffC(1,j) - DiffC(nx-1,j))/h2
          else
            dPhi = (DiffC(i+1,j) - DiffC(i-1,j))/h2
          end if

          c1(i,j) = csc2(j)*(DiffC(i,j) - h/2.0 * dPhi)
          c3(i,j) = csc2(j)*(DiffC(i,j) + h/2.0 * dPhi)

          dTheta= (DiffC(i,j+1) - DiffC(i,j-1))/h2
          c2(i,j) = DiffC(i,j) - h/2.0*(DiffC(i,j)*cot(j) + dTheta)  
          c4(i,j) = DiffC(i,j) + h/2.0*(DiffC(i,j)*cot(j) + dTheta)  
        end do
      end do

      end

! part 6
! converge.f90
!*****************************************************************************

      subroutine Convergence (u, lastu, tolmax, Converged)

      implicit none

      include 'ebm.inc'

      logical Converged
      real   u(NX6,NY6), lastu(NX6,NY6), diff(NX6,NY6)
      real   tolmax, relerr, maxd, maxu
      integer i, j
!-----------------------------------------------------------------------------
! Calculate the maximum relative error


      do j = 1, NY6
        do i = 1, NX6
          diff(i,j) = abs(u(i,j)-lastu(i,j))
        end do
      end do

      maxd = 0.0
      maxu = 0.0
      do j = 1, NY6
        do i = 1, NX6
          if (diff(i,j) .gt. maxd) then
            maxd = diff(i,j)
          end if
          if (abs(u(i,j)) .gt. maxu) then
            maxu = abs(u(i,j))
          end if 
        end do
      end do

	if (maxu .ne. 0.0) then
        relerr = maxd/maxu
!       write (1,10) relerr
!10     format ('Relative Error = ',f11.8)
      end if

      if (relerr .le. tolmax) then
        Converged = .true.
      end if

      end

! part 7
! diffcoef.f90
!**************************************************************************

  subroutine DiffCoeffs (nx, ny, DiffC, Geography, area, DT2np, DT2sp, &
                         MDnp, MDsp)  

      implicit none
      include 'ebm.inc'

      integer nx, ny, i, j
      real   DiffC(nx,ny), Geography(nx,ny)                      
      real   area(ny), DT2np(nx), DT2sp(nx)
      real   Tarea, MDnp, MDsp, Geo
      real   theta, colat_term
!-------------------------------------------------------------------------
! Calculate the diffusion coefficients at each grid level

      do j = 1, ny
       theta = pi*real(j-1)/real(ny-1)
       colat_term = sin(theta)**5                    ! centered on equator
       do i = 1, nx
         Geo = Geography(i,j)
         if (Geo>=5.0 .and. Geo<=7.0) then           ! oceans 
           DiffC(i,j) = (Keq-Kocean)*colat_term + Kocean 
         else                                        ! land, sea ice, etc
           if (j <= 33) then                         ! northern hemisphere
             DiffC(i,j) = (Kland-KlandNP)*colat_term + KlandNP  
           else                                      ! southern hemisphere
             DiffC(i,j) = (Kland-KlandSP)*colat_term + KlandSP  
	     end if
         end if
       end do
      end do

! Calculate the area weighted mean of the diffusion at the mid-point between
! the pole and the first ring of grid points.
      MDnp = 0.0
      MDsp = 0.0
      Tarea = area(1) + area(2)
      do i = 1, nx
        DT2np(i) = (area(1)*DiffC(1,1)  + area(2)*DiffC(i,2))/Tarea 
        DT2sp(i) = (area(1)*DiffC(1,ny) + area(2)*DiffC(i,ny-1))/Tarea 
        MDnp = MDnp + DT2np(i)
        MDsp = MDsp + DT2sp(i)
      end do

      end

! part 8
! fmgsetup.f90
!*****************************************************************************

SUBROUTINE FMG_Setup (nx, ny, h, geom, Heatcap, Geography, GCnp, GCsp)

IMPLICIT NONE
INCLUDE 'memory.inc'
INCLUDE 'ebm.inc'

integer,dimension(NG):: nx, ny
integer:: Geography(NX6,NY6)
real,dimension(NG):: h, geom, GCnp, GCsp, MDnp, MDsp
real:: Heatcap(NX6,NY6)
real,allocatable:: Geog(:,:)
integer:: grid

!-------------- Calculate the time-independent stuff for FMG Solver ----------
  DO grid = 1, NG
    nx(grid) = 2**(grid+1)
    ny(grid) = 2**grid+1
    h(grid) = pi/real(2**grid)
  END DO

! Initialize memory storage
  mindex = 0
  CALL GetMemory 
      
! Save heatcap at finest grid level into memory storage
  CALL CopytoMemory (NX6, NY6, heatcap, z(ihc(NG)))
      
! Save geography at finest grid level into memory storage
  allocate (Geog(NX6,NY6)) 
  Geog = real(geography)                   ! convert to real
  CALL CopytoMemory (NX6, NY6, Geog, z(iGeo(NG)))
  deallocate (Geog)

! Calculate grid point areas and trig functions on the grids. Resample the geography
! and heat capacity onto lower grid levels.

  DO grid = NG, 1, -1
    CALL TrigFuncs (ny(grid), h(grid), z(icsc(grid)), z(icsc2(grid)), &
                    z(icot(grid))) 
    CALL Geometry (nx(grid), ny(grid), h(grid), z(iaw(grid)), geom(grid))     
    if (grid .gt. 1) then
      CALL ResampleGrid (nx(grid), ny(grid), z(ihc(grid)), nx(grid-1), &
                         ny(grid-1), z(ihc(grid-1)))
      CALL ResampleGrid (nx(grid), ny(grid), z(iGeo(grid)), nx(grid-1), &
                         ny(grid-1), z(iGeo(grid-1)))
    end if
  END DO

!Calculate the diffusion coefficients, time-independant differential operators, 
!and B on each grid level
  DO grid = NG, 1, -1
    CALL DiffCoeffs (nx(grid), ny(grid), z(iD(grid)), z(iGeo(grid)),  &
         z(iaw(grid)), z(iDnp(grid)), z(iDsp(grid)), MDnp(grid), MDsp(grid))
     
    CALL ConstDiffOper(nx(grid),ny(grid),h(grid),z(iD(grid)),z(icsc2(grid)), &
         z(icot(grid)), z(ic1(grid)), z(ic2(grid)), z(ic3(grid)), z(ic4(grid)))
     
    CALL ConstB (nx(grid), ny(grid), h(grid), z(ihc(grid)), geom(grid), &
                 MDnp(grid), MDsp(grid), GCnp(grid), GCsp(grid), &
                 z(iD(grid)), z(icsc2(grid)), z(ic0(grid)))
  END DO
      
END SUBROUTINE

! part 9
! fmgsolve.f90
!******************************************************************************

subroutine FMG_Solver (nx, ny, h, geom, GCnp, GCsp, Converged, rhs, u, Debug)

! Solve the elliptic PDE for the energy balance model on the sphere using the
! Full MultiGrid method. The grid levels and the number of grid points in each
! direction are as shown:
!
!   grid  nx   ny      
!    6   128   65   'finest'     
!    5    64   33    
!    4    32   17   
!    3    16    9  
!    2     8    5 
!    1     4    3   'coarsest'
!
! The finest grid will be considered to be at the 'top' and the coarsest
! grid at the 'bottom' of a stack of grids. Examples for the bottom two grid
! levels are:
!
!   Grid 2 is 8 x 5
!
!   colatitude  180   *   *   *   *   *   *   *   *  
!               135   *   *   *   *   *   *   *   *
!                90   *   *   *   *   *   *   *   *
!                45   *   *   *   *   *   *   *   *
!                 0   *   *   *   *   *   *   *   *
!                     -180                       180
!                                 latitude
!
!   Grid 1 is 4 x 3
!
!                180  *       *       *       *
!                 90  *       *       *       *
!                  0  *       *       *       *
!
!                     -180                       180
!-------------------------------------------------------------------------------
      implicit none
      include 'memory.inc'
      include 'ebm.inc'

      integer  nx(NG)                   ! x dimensions at all grid levels
      integer  ny(NG)                   ! y dimensions at all grid levels
      real     h(NG)                    ! grid increment at each level
      real     geom(NG)                 ! geometry at poles for each grid level
      real     GCnp(NG)
      real     GCsp(NG)
      real     u(NX6,NY6)               ! solution on finest grid level 
      real     rhs(NX6,NY6)             ! rhs on the finest grid level 
      real     LastSol(NX6,NY6)         ! saved last solution
      real     lasterr
      real     tolmax                   ! maximum relative error tolerance
	parameter (tolmax = 1.0e-3)
      integer  MaxCycles                ! maximum number of V cycles
	parameter (MaxCycles = 15)
      integer  NPRE                     ! number of relaxation sweeps before 
      parameter (NPRE = 3)              ! restriction of residual
      integer  NPOST                    ! number of relaxation sweeps after  
      parameter (NPOST = 2)             ! restriction of residual
      integer  grid                     ! current grid level
      integer  Vcycles                  ! number of V cycles in FMG
      logical  Converged                ! has solution converged or not?
      logical  Relax_XY                 ! relax in both x and y directions?
      logical  Debug
      integer  j, n 
!----------------------------------------------------------------------------
      Relax_XY = .true.
      grid = NG

! Copy the temperature into the top grid memory storage
      call CopytoMemory (nx(grid), ny(grid), u, z(iu(grid)))

! Restrict temperature from top grid to the next grid down.
      call Restriction (nx(grid), ny(grid), z(iaw(grid)), z(iu(grid)), &
                       nx(grid-1), ny(grid-1), z(iu(grid-1)), Debug)

! Continue to recursively restrict from current grid to the 
! next grid down until the bottom (coarsest) grid is reached.
      do grid = NG-1, 2, -1
        call Restriction (nx(grid), ny(grid), z(iaw(grid)), z(iu(grid)), &
                         nx(grid-1), ny(grid-1), z(iu(grid-1)), Debug)
      end do

! Calculate the diffusion coefficients at each grid level and then the 
! coefficients for the diffusion operator at each grid level. Only do
! this once for the constant diffusion case. For the time-dependent diffusion 
! coefficients you must discretize the operators at each time step. 

      grid = NG

! Copy the input rhs into memory storage at finest grid level
      call CopytoMemory (nx(grid), ny(grid), rhs, z(irhs(grid)))

! Restrict rhs from top grid to the next grid down.
      call Restriction (nx(grid), ny(grid), z(iaw(grid)), z(irhs(grid)), &
                       nx(grid-1), ny(grid-1), z(irho(grid-1)), Debug)

! Continue to recursively restrict from current grid to the 
! next grid down until the bottom (coarsest) grid is reached.
      do grid = NG-1, 2, -1
        call Restriction (nx(grid), ny(grid), z(iaw(grid)), z(irho(grid)), &
                         nx(grid-1), ny(grid-1), z(irho(grid-1)), Debug)
      end do

! Get initial solution on the coarsest grid
      call CoarseSol (nx(1), ny(1), h(1), geom(1), GCnp(1), GCsp(1), &
                     z(irho(1)), z(iu(1)), Debug)
 
!----------------------- FULL MULTIGRID ITERATION ------------------------------
  if (Debug) then
!    write (4,1) grid
1   format (/,'>>>>>>>>>> BEGIN FULL MULTIGRID CYCLE on grid: ',i2,'<<<<<<<<<')
  end if

      lasterr = 0.0
      Vcycles = 1 

      do grid = 2, NG 

! Interpolate solution from coarse grid to next finer grid
      call Interpolation(nx(grid-1),ny(grid-1), z(iaw(grid-1)), z(iu(grid-1)), &
                         nx(grid), ny(grid), z(iu(grid)), Debug)

! set up the right hand side
        if (grid .lt. NG) then
          call CopytoMemory (nx(grid), ny(grid), z(irho(grid)), z(irhs(grid))) 
        end if
   
 100    CONTINUE          ! loop for V cycles

        if (Debug) then
          write (4,2)
 2        format (/,'------->> BEGIN DOWN STROKE OF V CYCLE')
        end if

        do n = grid, 2, -1 

          if (Debug) then
            write (4,3) n
 3          format (/,'DOWN_STROKE n = ',i2)
            write (4,4)
 4          format (/,'PRE-RESTRICTION SMOOTHING')
          end if

!      Pre-smoothing
          do j = 1, NPRE
            if (Relax_XY) then
              call YRelaxation (nx(n), ny(n), h(n), z(iu(n)), z(irhs(n)), &
              geom(n), z(ic0(n)), z(ic1(n)), z(ic2(n)), z(ic3(n)), &
              z(ic4(n)), GCnp(n), GCsp(n), z(iDnp(n)), z(iDsp(n)), Debug) 
            end if
            call XRelaxation (nx(n), ny(n), h(n), z(iu(n)), z(irhs(n)), &
            geom(n), z(ic0(n)), z(ic1(n)), z(ic2(n)), z(ic3(n)), & 
            z(ic4(n)), GCnp(n), GCsp(n), z(iDnp(n)), z(iDsp(n)), Debug) 
          end do

          if (Vcycles .eq. 1 .and. grid .eq. NG) then
            call CopytoMemory (nx(grid), ny(grid), z(iu(grid)), LastSol) 
          end if

          call Residual (nx(n), ny(n), h(n), geom(n), GCnp(n), &
           GCsp(n), z(iDnp(n)), z(iDsp(n)), z(ic0(n)), z(ic1(n)), z(ic2(n)), & 
           z(ic3(n)), z(ic4(n)), z(iu(n)), z(irhs(n)), z(ires(n)), Debug)

!      Restrict the residual down to the next coarser grid. This will be the
!      next right hand side.
          call Restriction (nx(n), ny(n), z(iaw(n)), z(ires(n)), nx(n-1),  &
                           ny(n-1), z(irhs(n-1)), Debug)

!     For the next relaxation start with zero initial guess
          call SettoZero (nx(n-1), ny(n-1), z(iu(n-1)))

        end do 

        if (Debug) then
          write (4,7)
 7        format (/,'------->> END DOWN STROKE OF V CYCLE')
        end if

!    At bottom of V solve residual equation on coarsest grid
        call CoarseSol (nx(1), ny(1), h(1), geom(1), GCnp(1), GCsp(1), &
                       z(irhs(1)), z(iu(1)), Debug)

        if (Debug) then
          write (4,8)
 8        format (/,'------>> BEGIN UP STROKE OF V CYCLE')
        end if

        do n = 2, grid
       
          if (Debug) then
            write (4,9) n
 9          format (/,'UP_STROKE n = ',i2)
          end if

          call AddCorrection (nx(n-1), ny(n-1), z(iu(n-1)), z(iaw(n-1)), & 
                             nx(n), ny(n), z(ires(n)), z(iu(n)), Debug)
          if (Debug) then
            write (4,10)
 10         format (/,'POST-PROLONGATION SMOOTHING')
          end if

!      Post-smoothing
          do j = 1, NPOST
            if (Relax_XY) then
              call YRelaxation (nx(n), ny(n), h(n), z(iu(n)), z(irhs(n)), &
              geom(n), z(ic0(n)), z(ic1(n)), z(ic2(n)), z(ic3(n)), &
              z(ic4(n)), GCnp(n), GCsp(n), z(iDnp(n)), z(iDsp(n)), Debug) 
            end if
            call XRelaxation (nx(n), ny(n), h(n), z(iu(n)), z(irhs(n)),  &
            geom(n), z(ic0(n)), z(ic1(n)), z(ic2(n)), z(ic3(n)), & 
            z(ic4(n)), GCnp(n), GCsp(n), z(iDnp(n)), z(iDsp(n)), Debug) 
          end do

        end do 

        if (Debug) then
          write (4,11)
 11       format (/,'------>> END UP STROKE OF V CYCLE')
        end if

! Determine if convergence criteria are met on finest grid
        if (grid .eq. NG) then
          call Convergence (z(iu(NG)), LastSol, tolmax, Converged)

          if (.NOT.Converged .and. (Vcycles .lt. MaxCycles)) then
            call CopytoMemory (nx(NG), ny(NG), z(iu(NG)), LastSol)
            Vcycles = Vcycles + 1
            goto 100 
          else if (.NOT.Converged .and. (Vcycles .eq. MaxCycles)) then
            goto 200 
          else
!           write (2,16) Vcycles 

 16         format ('CONVERGENCE AFTER ',i2,' V CYCLES')
          end if
        end if
         
      end do 

!-----------------------------------------------------------------------------
! Save final solution
      call CopytoMemory (nx(NG), ny(NG), z(iu(NG)), u)

200   continue
      end

! part 10
! geometry.f90
      subroutine Geometry (nx, ny, h, area, geom)

      implicit none

      integer ny, nx, j
      real  geom, h, area(ny)
!-----------------------------------------------------------------------------
! Geometry and needed constants at each grid level.

! Grid point fractional area for interior points
      do j = 2, ny-1
        area(j) =  sin(h/2.0)*sin(h*float(j-1))/float(nx)
      end do

! Fractional area for the poles
      area(1) = 0.50*(1.0 - cos(h/2.0))
      area(ny) = area(1)

      geom = sin(h/2.0)/area(1)

      end

! part 11
! interp.f90
!**************************************************************************

      subroutine Interpolation (nxc, nyc, area, uc, nxf, nyf, uf, debug)

      implicit none

      include 'ebm.inc'

      integer  nxc, nyc, nxf, nyf
      real     uc(nxc,nyc), uf(nxf,nyf), area(nyc)
      integer  ci, cj, fi, fj
	logical  Debug
!----------------------------------------------------------------------------
      if (Debug) then
        write (4,2) nxc, nxf
 2      format (/,'PROLONGATION of correction from grid nxc = ',i2, &
                ' to grid nxf = ',i2)
      end if

! Copy points which coincide, which are the odd rows and odd columns on 
! the fine grid.
      do cj = 1, nyc
        if (cj .eq. 1) then                ! north pole
          do fi = 1, nxf
            uf(fi,1) = uc(1,1)
          end do
        else if (cj .eq. nyc) then         ! south pole
          do fi = 1, nxf
            uf(fi,nyf) = uc(1,nyc)
          end do
        else                             ! interior
          fj = 2*cj - 1
          do ci = 1, nxc
            fi = 2*ci - 1
            uf(fi,fj) = uc(ci,cj)
          end do
        end if
      end do

! interpolate horizontally at same latitude on the coarse grid 
      do cj = 2, nyc-1
        fj = 2*cj-1                                     ! odd rows on fine grid
        do ci = 1, nxc
          fi = 2*ci                                     ! even columns on fine grid
          if (ci .eq. nxc) then                           ! right periodic boundary
            uf(fi,fj) = 0.50 * (uc(ci,cj) + uc(1,cj)) 
          else
            uf(fi,fj) = 0.50 * (uc(ci,cj) + uc(ci+1,cj)) 
          end if
        end do
      end do

! interpolate vertically at the same longitude on the coarse grid 
      do ci = 1, nxc
        fi = 2*ci-1                                    ! odd columns on fine grid
        do cj = 1, nyc-1 
          fj = 2*cj                                   ! even rows on fine grid
          uf(fi,fj) = (area(cj)*uc(ci,cj) + area(cj+1)*uc(ci,cj+1))/  &
                    (area(cj) + area(cj+1)) 
        end do
      end do

! Fill in the fine grid holes at even columns, even rows
      do fj = 2, nyf-1, 2
        do fi = 2, nxf, 2
         if (fi .eq. nxf) then
          uf(fi,fj) = 0.5*((area(fj/2+1)*uf(fi,fj+1) + area(fj/2)*uf(fi,fj-1))/ &
                      (area(fj/2+1) + area(fj/2)) + 0.5*(uf(fi-1,fj) + uf(1,fj)))
         else
          uf(fi,fj) = 0.5*((area(fj/2+1)*uf(fi,fj+1) + area(fj/2)*uf(fi,fj-1))/ &
                     (area(fj/2+1) + area(fj/2)) + 0.5*(uf(fi-1,fj) + uf(fi+1,fj)))
          end if
        end do
      end do

      if (Debug) then
        write (4,3)
 3      format (/,'INTERPOLATED RESULT on finer grid:')
        do fj = 1, nyf
          write (4,10) (uf(fi,fj), fi = 1, nxf)
        end do
 10     format (5(f13.8,1x))
      end if

      end

! part 13
! memory.f90
!************************************************************************
      subroutine GetMemory 

      implicit none
      include 'memory.inc'

      integer  nx, ny, grid, Memalloc
!----------------------------------------------------------------------------
! Get the indicies for the memory array z for each grid level

      do grid = 1, NG
        nx = 2**(grid+1)
        ny = 2**grid + 1
        iaw(grid)   = Memalloc (ny)
        icsc(grid)  = Memalloc (ny)
        icsc2(grid) = Memalloc (ny)
        icot(grid)  = Memalloc (ny)
        iDnp(grid)  = Memalloc (nx)
        iDsp(grid)  = Memalloc (nx)
        ihc(grid)   = Memalloc (nx*ny)
        iD(grid)    = Memalloc (nx*ny)
        iGeo(grid)  = Memalloc (nx*ny)
        ic0(grid)   = Memalloc (nx*ny)
        ic1(grid)   = Memalloc (nx*ny)
        ic2(grid)   = Memalloc (nx*ny)
        ic3(grid)   = Memalloc (nx*ny)
        ic4(grid)   = Memalloc (nx*ny)
        if (grid .lt. NG) then
        irho(grid)  = Memalloc (nx*ny)
        end if
        irhs(grid)  = Memalloc (nx*ny)
        ires(grid)  = Memalloc (nx*ny)
        iu(grid)    = Memalloc (nx*ny)
      end do

      end
!****************************************************************************

      integer function Memalloc (length)
 
! This function mimics dynamical storage allocation of memory. It returns
! the index to the starting element of 'length' array elements in the array
! z. The preceeding array element is filled with the value of 'length', and
! the variable mindex is updated to point to the last element of z that has
! been used.

      implicit none
      include 'memory.inc'

      integer length                 ! length of array to be allocated memory
!----------------------------------------------------------------------------
     
      if ((mindex + 1 + length) .gt. MEMLEN) then
        write (*,10) mindex+1, length
 10     format (/,'mindex = ',i7,' length = ',i6)
        STOP 'Memory allocation is insufficient'
      else
        z(mindex+1) = length         ! store length of memory to be allocated
        Memalloc = mindex + 2        ! return index of first element of memory
        mindex = mindex + 1 + length ! point to last element of memory used
      end if
    
      end

! part 14
! memutils.f90
!**********************************************

SUBROUTINE CopytoMemory (nx, ny, ain, aout)

IMPLICIT NONE
integer::  nx, ny
real,dimension(nx,ny):: ain, aout 

aout = ain

END SUBROUTINE

!---------------------------------------------

SUBROUTINE SettoZero (nx, ny, a)

IMPLICIT NONE
integer::  nx, ny
real,dimension(nx,ny):: a

a = 0.0

END SUBROUTINE

! part 15
! regrid.f90
!******************************************************************************
      subroutine ResampleGrid (nxf, nyf, fhc, nxc, nyc, chc)

      implicit none

      integer nxf, nyf, nxc, nyc
      real    fhc(nxf,nyf), chc(nxc,nyc)
      integer ci, cj, fi, fj
!-----------------------------------------------------------------------------
! Resample the array values from fine to coarse grid

      do cj = 1, nyc
        fj = 2*cj - 1
        do ci = 1, nxc
          fi = 2*ci - 1
          chc(ci,cj) = fhc(fi,fj)
        end do
      end do

      end

! part 16
! relaxx.f90
!**************************************************************************
   subroutine XRelaxation (nx, ny, h, u, rhs, geom, c0, c1, c2, c3, c4, &
                             GCnp, GCsp, DT2np, DT2sp, Debug)

      implicit none
      include 'ebm.inc'

      integer nx, ny, col, row, first_row
      real    u(nx,ny), rhs(nx,ny) 
      real    c0(nx,ny), c1(nx,ny), c2(nx,ny), c3(nx,ny), c4(nx,ny)
      real    h, geom, NPsum, SPsum, hh 
      real    GCnp, GCsp, DT2np(nx), DT2sp(nx)
      real   Upper(NX6), Diag(NX6), Lower(NX6), R(NX6), X(NX6)
      integer i, j 
	logical Debug
!----------------------------------------------------------------------------

      if (Debug) then
        write (4,2) nx, ny
 2      format (/,'X LINE RELAXATION  nx = ',i3,' ny = ',i2)
      end if

      hh = h*h 

! Do Gauss Seidel Line Relaxation (Zebra relaxation) 

      do first_row = 2, 3           ! relax even rows first

        do row = first_row, ny-1, 2

! Build LHS operator, these can depend upon lat, long, and time
          do col = 1, nx-1
            Lower(col) = -c1(col+1,row)         ! lower diagonal
            Diag(col) =  c0(col,row)            ! main diagonal
            Upper(col) = -c3(col,row)           ! upper diagonal
          end do
          Diag(nx) = c0(nx,row)
          Lower(nx) = -c3(nx,row)               ! periodic BC
          Upper(nx) = -c1(1,row)                ! periodic BC

! Build RHS for this row
          do col = 1, nx
!      rows above and below
            R(col) = -hh*rhs(col,row) + c2(col,row)*u(col,row-1) +  &
                                      c4(col,row)*u(col,row+1)
          end do

! Do the Gaussian elimination to put the array into an upper triangular matrix
! form, then solve by back substitution.
          call SolveSparse (nx, Upper, Lower, Diag, R, X)

! Update solution array for this row
          do col = 1, nx
            u(col,row) = X(col)
          end do

        end do 

        if (first_row .eq. 2) then                ! relax the poles

          NPsum = 0.0
          SPsum = 0.0
          do i = 1, nx
            NPsum = NPsum + DT2np(i) * u(i,2)
            SPsum = SPsum + DT2sp(i) * u(i,ny-1)
          end do

          u(1,1) =  (geom*NPsum - rhs(1,1))/GCnp
          u(1,ny) = (geom*SPsum - rhs(1,ny))/GCsp

          do i = 1, nx
            u(i,1)  = u(1,1) 
            u(i,ny) = u(1,ny)
          end do

          if (Debug) then
            write (4,3) u(1,1) 
 3          format (/,'SOLUTION FOR NORTH POLE',f13.8)
            write (4,4) u(1,ny) 
 4          format ('SOLUTION FOR SOUTH POLE',f13.8)
          end if

        end if

        if (ny .eq. 3) goto 100 

      end do 

100   continue

      if (Debug) then
        write (4,9)
 9      format (/,'SOLUTION :')
        do j = 1, ny
          write (4,10) (u(i,j), i = 1, nx)
        end do
 10     format (5(f13.8,1x))
      end if

      end

! part 17
! relaxy.f90
!**************************************************************************
  subroutine YRelaxation (nx, ny, h, u, rhs, geom, c0, c1, c2, c3, c4, &
                         GCnp, GCsp, DT2np, DT2sp, Debug)

      implicit none

      include 'ebm.inc'

      integer nx, ny, col, row, first_col
      real    u(nx,ny), rhs(nx,ny) 
      real    c0(nx,ny), c1(nx,ny), c2(nx,ny), c3(nx,ny), c4(nx,ny)
      real    h, geom, NPsum, SPsum, hh 
      real    GCnp, GCsp, DT2np(nx), DT2sp(nx)
      real   Upper(NY6-2), Diag(NY6-2), Lower(NY6-2), R(NY6-2)
      integer i, j 
	logical Debug
!----------------------------------------------------------------------------
      if (Debug) then
        write (4,2) nx, ny
2       format (/,'Y LINE RELAXATION  nx = ',i3,' ny = ',i2)
      end if

      hh = h*h 

! Do Gauss Seidel Line Relaxation (Zebra relaxation) 

      do first_col = 1, 2           ! relax odd columns first

        do col = first_col, nx, 2

! Build LHS operator, these can depend upon lat, long, and time
          do row = 2, ny-2
            Lower(row-1) = -c2(col,row+1)         ! lower diagonal
            Diag (row-1) =  c0(col,row)           ! main diagonal
            Upper(row-1) = -c4(col,row)           ! upper diagonal
          end do
          Diag(ny-2) = c0(col,ny-1)

! Build RHS for this column
          do row = 2, ny-1 
!      cols left and right 
            if (col .eq. 1) then            ! periodic boundary
              R(row-1) = -hh*rhs(col,row) + c1(col,row)*u(nx,row) +  &
                                          c3(col,row)*u(col+1,row)
            else if (col .eq. nx) then      ! periodic boundary
              R(row-1) = -hh*rhs(col,row) + c1(col,row)*u(col-1,row) +  & 
                                          c3(col,row)*u(1,row)
            else
              R(row-1) = -hh*rhs(col,row) + c1(col,row)*u(col-1,row) +  &
                                          c3(col,row)*u(col+1,row)
            end if
          end do
! Boundaries at the poles
          R(1) = R(1) + c2(col,2)*u(1,1)
          R(ny-2) = R(ny-2) + c4(col,ny-1)*u(1,ny)

! Solve the tridiagonal matrix
          call TriDiagonal (ny-2, Upper, Lower, Diag, R)

! Update solution array for this column 
          do row = 1, ny-2 
            u(col,row+1) = R(row)
          end do

        end do 

      end do 

! Relax the poles

      NPsum = 0.0
      SPsum = 0.0
      do i = 1, nx
        NPsum = NPsum + DT2np(i) * u(i,2)
        SPsum = SPsum + DT2sp(i) * u(i,ny-1)
      end do

      u(1,1) =  (geom*NPsum - rhs(1,1))/GCnp
      u(1,ny) = (geom*SPsum - rhs(1,ny))/GCsp

      do i = 1, nx
        u(i,1)  = u(1,1) 
        u(i,ny) = u(1,ny)
      end do

      if (Debug) then
        write (4,3) u(1,1) 
 3      format (/,'SOLUTION FOR NORTH POLE',f13.8)
        write (4,4) u(1,ny) 
 4      format ('SOLUTION FOR SOUTH POLE',f13.8)
        write (4,9)
 9      format (/,'SOLUTION :')
        do j = 2, ny-1
          write (4,10) (u(i,j), i = 1, nx)
        end do
 10     format (5(f13.8,1x))
      end if

      end

! part 18
! residual.f90
!******************************************************************************
    subroutine residual (nx, ny, h, geom, GCnp, GCsp, DT2np,  &
                        DT2sp,c0, c1, c2, c3, c4, u, rhs, res, Debug)

! Returns the negative of the residual for the problem (i.e. -res = rhs - Lu).

      implicit none

      include 'ebm.inc'

      integer  nx, ny
      real    c0(nx,ny), c1(nx,ny), c2(nx,ny), c3(nx,ny), c4(nx,ny)
      real    u(nx,ny), rhs(nx,ny), res(nx,ny)
      real    GCnp, GCsp, DT2np(nx), DT2sp(nx)
      real    h, h2, geom, NPsum, SPsum
      integer i, j
	logical Debug
!-------------------------------------------------------------------------------
      if (Debug) then
        write (4,2)
 2      format (/,'Calculate RESIDUAL')
      end if

      h2 = 1.0/(h*h)

! residual at the interior grid points
      do j = 2, ny-1
        do i = 1, nx
          if (i .eq. 1) then
            res(i,j) = rhs(i,j) + h2*(c0(i,j)*u(i,j) - c1(i,j)*u(nx,j) -  &
                      c3(i,j)*u(i+1,j) - c2(i,j)*u(i,j-1) - c4(i,j)*u(i,j+1))
          else if (i .eq. nx) then
            res(i,j) = rhs(i,j) + h2*(c0(i,j)*u(i,j) - c1(i,j)*u(i-1,j) -  &
                      c3(i,j)*u(1,j) - c2(i,j)*u(i,j-1) - c4(i,j)*u(i,j+1))
          else
            res(i,j) = rhs(i,j) + h2*(c0(i,j)*u(i,j) - c1(i,j)*u(i-1,j) -  &
                      c3(i,j)*u(i+1,j) - c2(i,j)*u(i,j-1) - c4(i,j)*u(i,j+1))
          end if
        end do
      end do

! residual at the poles

      NPsum = 0.0
      SPsum = 0.0
      do i = 1, nx
        NPsum = NPsum + DT2np(i) * u(i,2)
        SPsum = SPsum + DT2sp(i) * u(i,ny-1)
      end do

      res(1,1)  =  rhs(1,1)  - geom*NPsum + GCnp*u(1,1)
      res(1,ny) =  rhs(1,ny) - geom*SPsum + GCsp*u(1,ny)

      do i = 1, nx
        res(i,1)  = res(1,1) 
        res(i,ny) = res(1,ny) 
      end do
 
      end

! part 19
! restrict.f90
!***********************************************************************
      subroutine Restriction (nxf, nyf, area, uf, nxc, nyc, uc, Debug)

      implicit none

      include 'ebm.inc'

      integer nxf, nyf, nxc, nyc, fi, fj, ci, cj
      real    area(nyf), uf(nxf,nyf), uc(nxc,nyc), Tarea, MRnp, MRsp
	logical Debug
!------------------------------------------------------------------------
! Do restriction from fine grid to coarse grid for all points excluding the
! poles. Uses half weighting. 

      if (Debug) then
        write (4,2)
 2      format (/,'RESTRICTION from fine to coarser grid')
      end if

      do cj = 2, nyc-1
        fj = 2*cj - 1
        do ci = 1, nxc
          fi = 2*ci - 1
          if (ci .eq. 1) then                       ! periodic boundary on the left
            Tarea = 0.125*(area(fj+1) + 6.0*area(fj) + area(fj-1))
            uc(ci,cj) = 0.125 * (area(fj+1) * uf(fi,fj+1)  &                     
                     + area(fj) * (uf(nxf,fj) + 4.0*uf(fi,fj) + uf(fi+1,fj)) &  
                     + area(fj-1) * uf(fi,fj-1))/Tarea                
          else
            uc(ci,cj) = 0.125 * (area(fj+1) * uf(fi,fj+1)  &                    
                      + area(fj) * (uf(fi-1,fj) + 4.0*uf(fi,fj) + uf(fi+1,fj)) &  
                      + area(fj-1) * uf(fi,fj-1))/Tarea                
          end if
        end do
      end do

! Do restriction at the North and South Poles
      MRnp = 0.0
      MRsp = 0.0
      do fi = 1, nxf
        MRnp = MRnp + uf(fi,2)
        MRsp = MRsp + uf(fi,nyf-1)
      end do
      MRnp = MRnp/float(nxf)           ! mean value of ring of grid points one fine 
      MRsp = MRsp/float(nxf)           ! grid point from the poles
      Tarea = area(1) + area(2)       ! total fractional area of pole and ring points

      uc(1,1) = (area(1)*uf(1,1) + area(2)*MRnp)/Tarea 
      uc(1,nyc)= (area(nyf)*uf(1,nyf) + area(nyf-1)*MRsp)/Tarea

      do ci = 2, nxc
        uc(ci,1) = uc(1,1)
        uc(ci,nyc) = uc(1,nyc)
      end do

      end

! part 20
! sparse.f90
!**********************************************************************
      subroutine SolveSparse (n, U, L, D, R, X)

! Solves a system of equations which are of the form
!  
!    d  c  0  0  0  a     This type of system of equations can arise
!    a  d  c  0  0  0     in solving an elliptical PDE with periodic
!    0  a  d  c  0  0     boundary in the x direction
!    0  0  a  d  c  0
!    0  0  0  a  d  c
!    c  0  0  0  a  d
!
!  The top right corner element must be stored as the element U(n) 
!  The bottom left corner element must be stored as the element L(n)

      implicit none

      include 'ebm.inc'

      integer i, n 
      real    U(n)                     ! Upper diagonal
      real    D(n)                     ! Main diagonal
      real    L(n)                     ! Lower diagonal
      real    R(n)                     ! Right hand side of system
      real    X(n)                     ! Solution vector
      real    mult                     ! multiplier for eliminations
      real    row(NX6-1), col(NX6)     ! temporary storage
!------------------------------------------------------------------
! Initialize the temporary arrays
      do i = 1, NX6-1
        row(i) = 0.0
        col(i) = 0.0
      end do
      col(NX6) = 0.0

! Assign values to temporary arrays
      row(1) =   L(n)                    ! bottom left corner element
      row(n-1) = L(n-1)
      col(1) =   U(n)                    ! top right corner element
      col(n-1) = U(n-1)
      col(n) =   D(n)                    ! last main diagonal element

! First eliminate the element in the last row of the column   
      do i = 1, n-2
        mult = row(i)/D(i)
        row(i+1) = row(i+1) - mult * U(i)
        col(n) = col(n) - mult * col(i)
        R(n) = R(n) - mult * R(i)

! Now work on the sub diagonal in this column
        mult = L(i)/D(i)
        D(i+1) = D(i+1) - mult * U(i)
        col(i+1) = col(i+1) - mult * col(i)
        R(i+1) = R(i+1) - mult * R(i)
      end do

! Now eliminate the last element on the sub diagonal
      mult = row(n-1)/D(n-1)
      D(n) = col(n) - mult * col(n-1)    ! update D array
      R(n) = R(n) - mult * R(n-1) 

      U(n-1) = col(n-1)                  ! update U array

! Next do the back substitution to solve the system
! of equations
      if (D(n) .eq. 0.0) goto 99
      X(n) = R(n)/D(n)

      if (D(n-1) .eq. 0.0) goto 99
      X(n-1) = (R(n-1) - U(n-1)*X(n))/D(n-1)

      do i = 2, n-1
        if (D(n-i) .eq. 0.0) goto 99
        X(n-i) = (R(n-i) - U(n-i)*X(n-i+1) - col(n-i)*X(n))/D(n-i)
      end do

      return
 99   STOP 'SOLVESPARSE: Zero diagonal element'

      end

! part 21
! tridiag.f90
!******************************************************************************
      subroutine TriDiagonal (n, U, L, D, R)

! Solve the tridiagonal system of equations

      implicit none

      include 'ebm.inc'

      integer n                         ! number of equations in system
      real    U(n)                      ! upper diagonal of array
      real    D(n)                      ! main diagonal of array
      real    L(n)                      ! lower diagonal of array
      real    R(n)                      ! right hand side of system
      real    mult                      ! multiplier for eliminations
      integer i
!-----------------------------------------------------------------------------
! Solution is put into vector R

      do i = 2, n
        if (D(i-1) .eq. 0.0) goto 99
        mult = L(i-1)/D(i-1)
        D(i) = D(i) - mult * U(i-1)
        R(i) = R(i) - mult * R(i-1)
      end do

      if (D(n) .eq. 0.0) goto 99
      R(n) = R(n)/D(n)

      do i = 1, n-1
        if (D(n-i) .eq. 0.0) goto 99
        R(n-i) = (R(n-i) - U(n-i) * R(n-i+1))/D(n-i)
      end do

      return
 99   STOP 'TRIDIAGONAL: Zero diagonal element'

      end

! part 22
! trig.f90
!*******************************************************************************
      subroutine TrigFuncs (ny, dy, csc, csc2, cot)

      implicit none

      integer ny, j
      real    csc(ny), csc2(ny), cot(ny)
      real    dy, theta
!------------------------------------------------------------------------------

! Calculate the trig functions as a function of colatitude once 
! at each grid level 
      do j = 2, ny-1
        theta = dy * real(j-1)
        csc(j) = 1.0/sin(theta)
        csc2(j) = csc(j)*csc(j) 
        cot(j) = cos(theta)/sin(theta)
      end do

      end                                                                              
