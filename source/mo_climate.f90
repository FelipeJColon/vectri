MODULE mo_climate
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2011, ICTP
! tompkins@ictp.it
!
! parameters concerned with the climate input for model
!---------------------------------------------------------

  IMPLICIT NONE

! these are constant rain and temperature amounts
  REAL :: rconst_rain, rconst_temp

! temperature offset and trends
  REAL :: rtemperature_offset=0.0  ! fixed temperature offset
  REAL :: rtemperature_trend=0.0    ! fixed temperature trend K/annum

! rainfall multiplicative parameter
  REAL :: rrainfall_factor=1.0     ! scale rainfall input by this ratio

END MODULE mo_climate
