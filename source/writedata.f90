SUBROUTINE writedata(iday)
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2013, ICTP
! tompkins@ictp.it
!
! writedata subroutine to write to netcdf output
!
!
!---------------------------------------------------------
  USE netcdf
  USE mo_control
  USE mo_climate
  USE mo_constants
  USE mo_vectri
  USE mo_ncdf_tools

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iday

  INTEGER :: i

  DO i=1,ndiag2d
    WHERE(rpopdensity(:,:)<0.0) rdiag2d(:,:,i)=rfillvalue
  ENDDO

  IF (loutput_rain) CALL check( nf90_put_var(ncidout, ncvar_rain(1), rrain, start=(/ 1, 1, iday /) ))
  IF (loutput_t2m) CALL check( nf90_put_var(ncidout, ncvar_t2m(1), rtemp, start=(/ 1, 1, iday /) ))

  IF(loutput_waterfrac) &
  & CALL check(nf90_put_var(ncidout, ncvar_waterfrac(1), rwaterfrac, start=(/ 1, 1, iday /)))
  IF(loutput_vector) &
  & CALL check(nf90_put_var(ncidout, ncvar_vector(1),   rdiag2d(:,:,ncvar_vector(2)), start=(/ 1, 1, iday /)))
  IF(loutput_larvae) &
  & CALL check(nf90_put_var(ncidout, ncvar_larvae(1),   rdiag2d(:,:,ncvar_larvae(2)), start=(/ 1, 1, iday /)))
  IF(loutput_lbiomass) & 
  & CALL check(nf90_put_var(ncidout, ncvar_lbiomass(1),  rdiag2d(:,:,ncvar_lbiomass(2)), start=(/ 1, 1, iday /)))
  IF(loutput_pr) &
  & CALL check(nf90_put_var(ncidout, ncvar_pr(1),  rdiag2d(:,:,ncvar_pr(2)), start=(/ 1, 1, iday /)))
  !     CALL check(nf90_put_var(ncidout, ncvar_vectinfect(1),  rdiag2d(:,:,ncvar_vectinfect(2)), start=(/ 1, 1, iday /)))

  IF(loutput_prd) &
  & CALL check(nf90_put_var(ncidout, ncvar_prd(1),  rdiag2d(:,:,ncvar_prd(2)), start=(/ 1, 1, iday /)))
  IF(loutput_hbr) &
  & CALL check(nf90_put_var(ncidout, ncvar_hbr(1),  rdiag2d(:,:,ncvar_hbr(2)), start=(/ 1, 1, iday /)))
  IF(loutput_cspr) &
  & CALL check(nf90_put_var(ncidout, ncvar_cspr(1),  rdiag2d(:,:,ncvar_cspr(2)), start=(/ 1, 1, iday /)))
  IF(loutput_eir) &
  & CALL check(nf90_put_var(ncidout, ncvar_eir(1),  rdiag2d(:,:,ncvar_eir(2)), start=(/ 1, 1, iday /)))
  IF(loutput_cases) &
  & CALL check(nf90_put_var(ncidout, ncvar_cases(1),  rdiag2d(:,:,ncvar_cases(2)), start=(/ 1, 1, iday /)))
  IF(loutput_immunity) &
  & CALL check(nf90_put_var(ncidout, ncvar_immunity(1),  rdiag2d(:,:,ncvar_immunity(2)), start=(/ 1, 1, iday /)))
  IF(loutput_vecthostratio) &
  & CALL check(nf90_put_var(ncidout, ncvar_vecthostratio(1),  rdiag2d(:,:,ncvar_vecthostratio(2)), start=(/ 1, 1, iday /)))
 
END SUBROUTINE writedata

