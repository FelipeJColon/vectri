MODULE mo_ncdf_tools
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2011, ICTP
! tompkins@ictp.it
!
! tools for ncdf input output
!---------------------------------------------------------

CONTAINS

SUBROUTINE check(status)

  USE netcdf
  USE mo_control

  IMPLICIT NONE

  INTEGER, INTENT(in) :: status

  IF (status /= nf90_noerr) THEN
    WRITE(IOUNIT,*)TRIM(NF90_STRERROR(status))
    STOP 'Bad NETCDF status'
  END IF
END SUBROUTINE check


SUBROUTINE define_ncdf_output(ivarid,name,title,iunits,ndiag,idarray)
  USE netcdf
  USE mo_vectri
  USE mo_control

  IMPLICIT NONE

  INTEGER, INTENT(in) :: idarray(:) ! flexible vector to set coordinate dimensions
  INTEGER, INTENT(inout) :: ivarid(2), ndiag
  CHARACTER (len=*), INTENT(IN) :: name, title, iunits

  CALL check(NF90_DEF_VAR(ncidout, name, nf90_FLOAT, idarray , iVarID(1)))
  CALL check(NF90_PUT_ATT(ncidout, iVarID(1), "title", title) )
  CALL check(NF90_PUT_ATT(ncidout, iVarID(1), "units", iunits) )
  CALL check(NF90_PUT_ATT(ncidout, iVarID(1), "_FillValue", rfillvalue) )
  ndiag=ndiag+1   ! increase the diagnostic counter
  iVarid(2)=ndiag
END SUBROUTINE define_ncdf_output

END MODULE mo_ncdf_tools

