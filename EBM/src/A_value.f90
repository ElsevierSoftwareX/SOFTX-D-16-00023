subroutine A_value(CO2ppm,A)
!
!  REFERENCE:
!  Myhre, G., E.J. Highwood, K.P. Shine, and F. Stordal, 1998: 
!  New estimates of radiative forcing due to well mixed greenhouse gases. 
!  Geophysical Research Letters, 25, 2715-2718.
!  PURPOSE:
!     CALCULATE A VALUES (A+BT) CAUSED BY CHANGES OF CO2 LEVELS
!  INPUT:
!     CO2ppm
!
!  OUTPUT:
!     A
!
!     Table 3, Myhre et al.(1998)
!     A=210.15 FOR 1950AD WITH CO2ppm=315.0ppmv
      
      real CO2ppm, A
      real CO2_Base, A_Base
      parameter(CO2_Base=315.0, A_Base=210.3)
      
      A=A_Base-5.35*log(CO2ppm/CO2_Base)
      
      return
end subroutine A_value
