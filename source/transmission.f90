SUBROUTINE transmission(rdailyeir,rprob)
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins A.M. 2011, ICTP
! tompkins@ictp.it
!
! bites per person follow a Poisson distribution
!   - This is important if bite rate low - since it will reduce the mean 
!     tranmission rate (few people get bitten a lot, but can only get infected once!)
! Poisson distribution - used for biting rate
    ! Poisson PDF - this reduces mean transmission, as getting 7 infectious bites not 
    ! not much different probability to getting 4, but this implies 3 other people do 
    ! not get bitten.
!----------------------------------------------------------

  USE mo_constants
  USE mo_vectri

  IMPLICIT NONE

  REAL, INTENT(IN)    :: rdailyeir 
  REAL, INTENT(INOUT) :: rprob     ! mean probability of getting infected, integrated across the PDF

  INTEGER :: i,j

  REAL :: zpk
  REAL :: zpdf(0:nbitepdf) 

  !
  ! probability of N bites distributed by Poisson PDF
  ! with a mean of rdailyeir 
  !
  DO i=0,nbitepdf-1
    zpk=0.0
    DO j=1,i
      zpk=zpk+LOG(REAL(j))
    ENDDO
    zpdf(i)=EXP(REAL(i)*LOG(rdailyeir)-rdailyeir-zpk)
  ENDDO
  zpdf(nbitepdf)=MAX(1.0-SUM(zpdf(0:nbitepdf-1)),0.0)

  !
  ! final transmission rate is the dot product of the bite rate
  ! with the transmission probability array   
  !
  rprob=SUM(zpdf*rpdfvect2host)

  ! bed net use will go in here?
  
  RETURN
  END
