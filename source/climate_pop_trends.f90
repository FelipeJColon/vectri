SUBROUTINE climate_pop_trends(iday)
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2012, ICTP
! tompkins@ictp.it
!
! climate-pop subroutine to read in forcing data and setup output files
!
! If simple trends are requested, climate change, or fixed population growth,
! these are implemented here...
!
!---------------------------------------------------------
  USE mo_constants
  USE mo_vectri
  USE mo_climate

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iday

  ! temperature trend and offset
  rtemp=rtemp+rtemperature_trend*REAL(iday)*dt/rdaysinyear+rtemperature_offset

  ! rainfall factor
  rrain=rrain*rrainfall_factor

  ! population 

END SUBROUTINE climate_pop_trends

