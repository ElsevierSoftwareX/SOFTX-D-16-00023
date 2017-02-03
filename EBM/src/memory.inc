!---------------------------- memory.inc ---------------------------------------
! Linear memory storage stuff
      integer   NG                           ! number of grid levels
      parameter (NG = 6)
      integer   MEMLEN                       ! total size of memory array
      parameter (MEMLEN = 125*1024)
      integer   mindex                       ! index for memory array
      real      z(MEMLEN)                    ! linear memory array
      common /memory/ z, mindex
!-------------------------------------------------------------------------------
! Indicies into the memory storage
      integer iaw(NG)                   ! index for area weights
      integer icsc(NG)                  ! index for csc(theta)
      integer icsc2(NG)                 ! index for csc(theta)**2
      integer icot(NG)                  ! index for cot(theta)
      integer ihc(NG)                   ! index for heat capacities
      integer irho(NG-1)                ! index for initial rhs (forcing)
      integer iD(NG)                    ! index for diffusion coefficient D
      integer iGeo(NG)                  ! index for geography 
      integer ic0(NG)                   ! index for diff operator coefficient c0
      integer ic1(NG)                   ! index for diff operator coefficient c1
      integer ic2(NG)                   ! index for diff operator coefficient c2
      integer ic3(NG)                   ! index for diff operator coefficient c3
      integer ic4(NG)                   ! index for diff operator coefficient c4
      integer irhs(NG)                  ! index for rhs
      integer ires(NG)                  ! index for residual
      integer iu(NG)                    ! index for solution
      integer iDnp(NG)                  ! constant for north pole relaxation
      integer iDsp(NG)                  ! constant for south pole relaxation

      common /indicies/ iaw, icsc, icot, icsc2, irho, iD, iGeo, iDnp, iDsp, &
                        ihc, ic0, ic1, ic2, ic3, ic4, irhs, ires, iu
!-------------------------------------------------------------------------------
