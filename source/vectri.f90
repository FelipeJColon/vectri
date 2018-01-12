PROGRAM VECTRI
  ! 
  ! VECTRI: VECtor borne disease infection model of ICTP TRIeste.
  !
  ! Author: A.M.Tompkins 2011-2016
  !         ICTP, Strada Costiera 11, 34152 Trieste, ITALY
  !         tompkins@ictp.it
  !
  ! --------------------------------------------------------------------------
  ! LICENSE INFORMATION: 
  ! 1) This software is issued under a free-to-use license for research/education purposes
  ! 2) Please do not distribute further without prior permission from the author
  ! 3) Please email a copy of manuscripts using the software to the author 
  ! 4) For the model description cite:
  !      Tompkins and Ermert 2013 for the key model description (v1.2.6)
  !      Tompkins and Di Giuseppe 2015 for the updates to mortality
  !      Asare et al. 2016a GH for revisions to surface hydrology
  !
  !  Documentation for running the model and a short exercise sheet is found in doc/
  !
  ! --------------------------------------------------------------------------
  !
  ! Loosely based on HM2005 and Depinay 62 and Bomblies 2008 
  !
  ! Eventual wishlist of features:
  !   1) DONE: modified survival rate formula (Ermert 2010)
  !   2) advection/diffusion(flight) of vector from cell to cell 
  !   3) modified surface hydrology for egg laying
  !      eventually this will try to represent improved conditions in drought near rivers
  !      and improved conditions in rainy season elsewhere.
  !   4) diffusion of infected human population to account for local migration
  !   5) gravity model/relaxation of proportion p to account for long-distance migration
  !   6) standard deviation of temperature to account for subgrid-scale variability and microclimates
  !   7) DONE:  biting rates account for population density and are Poisson distributed.
  !   8) Immunity
  !
  ! modules explicitly the 
  ! 1) gonotrophic cycle: 
  !   the time taken to prepare a brood (determines biting cycle)
  !
  ! 2) sporogonic cycle (SC) (the process of fertilization
  !   of the macrogametocyte, formation of the oocyst,
  !   ookinete, penetration of the midgut and then the subsequent
  !   development of the sporozoites which dwell in the
  !   salivary glands)
  !
  !   variable convention
  !   -------------------
  !   integer parameter n
  !   integer local     i
  !   real parameter    r
  !   real local        z
  !   
  !   references 
  !
  !   1) Hoshen and Morse (2004)
  !   2) Ermert (2010)
  !   3) Bomblies et al. (2008)
  !   4) garki model (1974)
  !   5) depinay et al. (2004)
  !   6) Tompkins and Ermert (2013) <- model description at v1.2.6
  !   7) Craig et al. 1999 
  !
  ! --- Glossary ---
  !
  !  Entomological data
  !
  !   HBR - host bite rate - number of bites per host per day
  !   EIR - Entomological inoculation rate (n of infectious bites per person)
  !   CSPR - Circumsporozoite protein rate (fraction eir/hbr)
  !
  !  Epidemiological data
  !    
  !
  ! --- order of code ---
  !
  ! 1. write out diagnostics
  ! 2. read in met data for timestep
  ! 3. biting treatment
  !    3.1 - calculation transmission probabilities
  !    3.2 - Vector to human transmission
  !      3.2.1 Evolution of disease in host
  !      3.2.2 clearing rate 
  ! 4.  Sporogonic cycle
  ! 5.  Gonotrophic cycle
  ! 6.  Vector temperature dependent survival
  ! 7.  Vector Oviposition
  ! 8.  Puddle model
  ! 9. Larvae maturation - limitation of resourses
  ! 10. Larvae maturation - larvae progression
  ! 12. Larvae hatching
  ! 15. NCDF OUTPUT
  !
  ! --- release history --- 
  ! v1.0 May 2011 - Basic model consisting of LMM and Bomblies and Depinay.
  !               - with addition parametrization for 
  !                 a) ponds, b) pop density c) zoophily d) bednets
  ! v1.1 Nov 2011 - all params written to netcdf global attributes
  !               - all input params controllable in namelist input 
  !               - simple options added for fake climate change
  !               - simple migration trickle added
  ! v1.2 Jan 2012 - model structure rearranged to clean up and use directory structure
  !               - preprocessing upgraded orographic downscaling of temperature 
  ! v1.21         - mode 3 upgrades to use AR4.
  ! 
  ! v1.22-1.25    - various minor corrections to output files, naming changes, and script streamlines
  !               - ability to run with ECMWF seasonal forecasts.
  !
  ! v1.26***      - stable beta release - 
  !               - adjusted hydrology 
  !               - Version that is documented in Tompkins and Ermert (2013)***
  !
  ! v1.3 July 2012 - capability of exact restarts 
  !                - rain and temperature are read in one slice at a time (for AR5) and not necessarily in the output 
  !                - adapable output file
  !                - Long term support version
  !
  ! v1.3.1 July 2013 - Revised Larvae growth rates according to Craig et al. 1999 with 
  !                    new larvae mortality functions based on data from Bayoh and Lindsay 2004.
  !                 - indoor temperature included, with temperature experience a linear mix of indoor and outdoor
  !
  ! v1.3.2 Nov 2013 - Minor bug corrections
  !
  ! v1.3.3 April 2014 - Tidy up netcdf reading - all scale factors and latitude reversing now in pre-processing
  !                     All netcdf variables now controllable from namelist
  !
  ! v1.3.4 Oct 2014  - Bug correction on pond fraction initialization - 
  !
  ! v1.3.5 May 2015  - New transmission assumption that assumes EI categories receive more bites than susceptibles.
  !
  ! V1.3.6 March 2016 - Bayoh mortality added
  !
  ! V1.4   April 2016 - *** Major release ***
  !                     + SEIR Immunity added
  !                     + Clinical cases diagnostic added 
  !                     + New fit to larvae mortality data implemented using 
  !                        least squares fit to a logistic regression
  !                     + Revised stochastic term added to gonotrophic cycle to 
  !                        improve stability when ngono<100 
  !                     + Minor bug fixed to the transmission rates from host to vector 
  !                        that meant transmission was underestimated at temperatures<20C
  !                     + Minor bug fix in transmission to vector which neglect rbiteratio
  !                     + mortality of hosts included.
  !

  USE netcdf
  USE mo_constants
  USE mo_control
  USE mo_vectri
  USE mo_ncdf_tools

  IMPLICIT NONE

  ! local integers 
  INTEGER :: i,ix,iy,iday,idate,istep,iloop &
       &           ,iinfv,ivect,ihost,iinfh &
       &           ,icheck=0

! local real scalars
  REAL ::    zgonof, zsporof, zhostf &
       &           ,zsurvp &  ! adult vector survival probability
       &           ,zrain  &  ! dummy rain variables, daily
       &           ,zlarvmaturef, zdel, zmasslarv, zcapacity, zbiolimit &
       &           ,ztemp, ztempindoor, ztempwater & !dummy temperature, indoor temperature and water temperature
       &           ,zlimit, zcases &
       &           ,zvectinfectbite, zhost_SE, zhost_I, zhost_R &
       &           ,zprobhost2vect, zprobhost2vect_I, zprobhost2vect_R &
       &           ,zprobvect2hostS, zprobvect2hostEI & !transmission probabilities
       &           ,znnhbr, znndailyeir, znndailyeirS, znndailyeirEI &
       &           ,zpud1, zflushr & ! useful puddle variables
       &           ,zsurvpl

  ! fortran timer functions
  REAL :: time1=0.0,time2=0.0
  LOGICAL :: lspinup, &  ! spin up period - turn off output
       &            ltimer=.false.      ! turn the cpu timer on for the first timestep
  ! local real vectors
  REAL, ALLOCATABLE :: zarray1(:),zvector_density(:,:)
  REAL, ALLOCATABLE :: zpr(:,:) ! calculate of parasite ratio including ALL exposed (from day 1, i.e. not detectable)

! --------------------
! START OF VECTRI CODE
! --------------------
  CALL setup ! initial conditions for arrays

! allocate the 2d local arrays
  ALLOCATE(zvector_density(nlon,nlat))
  ALLOCATE(zpr(nlon,nlat))

  ! Spin up
  WRITE(iounit,*) 'starting run',nrun
  DO iloop=1,nrun
     WRITE(*,'(1a1,A5,I10,$)') char(13), 'step ',iloop
     IF (ltimer) CALL timer('go  ',icheck,time1,time2) !set up timer

     ! for spin up
     IF (iloop<=nyearspinup*INT(nlenspinup)/dt) THEN
        istep=MOD(iloop-1,INT(nlenspinup))+1
        IF (iloop==1) WRITE(iounit,*) 'spin up period switched on for ',nyearspinup*nlenspinup,' days' 
        lspinup=.TRUE. 
     ELSE 
        istep=iloop-nyearspinup*INT(nlenspinup/dt)
        IF (iloop==nyearspinup*INT(nlenspinup)/dt+1) write(iounit,*),'INTEGRATION STARTS' 
        lspinup=.FALSE. 
     ENDIF

     idate=ndate(istep)
     iday=istep*NINT(dt)

     !IF (MOD(iday,ndaydiag)==0) WRITE(iounit,*) 'date ',idate

     !---------------------------
     ! read in met data timeslice
     !---------------------------
     CALL getdata(istep)

     !-------------------------------------------------
     ! apply simple fixed climate or population changes
     !-------------------------------------------------
     CALL climate_pop_trends(istep)

     !------------------------------------------------
     ! 1. diagnostics for the present timestep
     !    Some of these are used in the analysis
     !    While some are for diagnostic purposes only 
     !------------------------------------------------

     !-----------------------
     ! 1.1 derive diagnostics
     !-----------------------

     !------------------------------------------------
     ! vector migration
     ! this will be a weighted random walk, weighted
     ! to allow clustering around populations
     ! for now we simply set a minimum  
     !------------------------------------------------
     rvect=MAX(rvect,0.0)
     rvect(:,0,:,:)=MAX(rvect(:,0,:,:),rvect_min)
     rlarv=MAX(rlarv,0.0)
     zvector_density(:,:)=SUM(SUM(rvect,DIM=1),DIM=1) ! total vector number = vector density    
     zpr=1.0-rhost(0,nadult_ni,:,:)/SUM(rhost(:,nadult_ni,:,:),DIM=1) ! parasite ratio

     ! store some diagnostics:
     IF (loutput_vector) rdiag2d(:,:,ncvar_vector(2))=zvector_density  ! store vector density
     IF (loutput_larvae) rdiag2d(:,:,ncvar_larvae(2))=SUM(rlarv,DIM=1) ! total larvae number
     IF (loutput_prd)    rdiag2d(:,:,ncvar_prd(2))= &                  ! parasite ratio
          &   SUM(rhost(INT(rhost_detectd):,1,:,:),DIM=1)/SUM(rhost(:,1,:,:),DIM=1)  !
     IF (loutput_pr)     rdiag2d(:,:,ncvar_pr(2))=zpr
     IF (loutput_cspr)   rdiag2d(:,:,ncvar_cspr(2))=SUM(rvect(:,ninfv,:,:),DIM=1)/MAX(zvector_density(:,:),reps)
     IF (loutput_vecthostratio) rdiag2d(:,:,ncvar_vecthostratio(2))=zvector_density(:,:)/MAX(rpopdensity(:,:),reps)

     ! zoophilic rates - well actually is anthropophilic rate.
     WHERE (rpopdensity(:,:)>=0.0) &
          & rzoophilic=1.0-(1.0-rzoophilic_min)*EXP(-rpopdensity(:,:)/rzoophilic_tau)
     rbitezoo(:,:)=1.0-(1.0-rbiteratio*rzoophilic(:,:))**dt ! product useful 

     !-----------------
     ! GRIDDED LOOP !!!
     !-----------------
     DO ix=1,nlon
        DO iy=1,nlat

           IF (ix==nxdg.and.iy==nydg) THEN
              zmasslarv=0.0 !dummy statement
           ENDIF

           !---------------------------------------------
           ! 1.2 point diagnostics for calculations below
           !---------------------------------------------
           ! total biomass of larvae per m2 of water surface
           zmasslarv=SUM(rlarv(:,ix,iy)*rmasslarv(:))
           if (loutput_lbiomass)rdiag2d(ix,iy,ncvar_lbiomass(2))=zmasslarv

           ! proportion of hosts in SEIR categories
           zhost_I=SUM(rhost(ninfh  ,:,ix,iy))/MAX(SUM(rhost(:,:,ix,iy)),reps) !+1
           zhost_R=SUM(rhost(ninfh+1,:,ix,iy))/MAX(SUM(rhost(:,:,ix,iy)),reps) !+1
           zhost_SE=1.0-zhost_I-zhost_R
           
           IF(lverbose)&
                &PRINT *,'host check ',zhost_SE,SUM(rhost(0:ninfh-1,:,ix,iy))/MAX(SUM(rhost(:,:,ix,iy)),reps)
           
           ! daily cspr: proportion of infected vectors that can tranmit and may bite in *this* timestep
           zvectinfectbite=rvect(0,ninfv,ix,iy)/MAX(SUM(rvect(0,:,ix,iy)),reps)

           !------------------------
           ! 2. meteorological data
           !------------------------
           zrain=MIN(MAX(rrain(ix,iy),0.0),200.0)  ! 1 day rainfall, with safety limits applied.

           ! indoor temperature - after Lunde et al. Malaria Journal 2013, 12:28
           ztempindoor=10.33+0.58*rtemp(ix,iy)

           ! water temperature 
           ztempwater=rtemp(ix,iy)+rwater_tempoffset ! water temperature

           IF(lverbose)&
           & PRINT *,'rtemp ',rtemp(ix,iy),rwater_tempoffset,ztemp
           ! temperature experience by vector is mix of indoor and outdoor 
           ! temperature on a daily timestep
           ztemp=MAX(rbeta_indoor*ztempindoor+(1.0-rbeta_indoor)*rtemp(ix,iy),-30.0)
           
           !---------------------------------------
           ! ONLY RUN MODEL IF NOT A LAKE/SEA POINT
           !---------------------------------------
           IF (rpopdensity(ix,iy)>=rpopdensity_min) THEN

              !---------------------
              ! 3. Biting treatment
              !---------------------

              ! --------------------------
              ! Transmission probabilities
              ! --------------------------
              !  
              ! - zbiteratio of vectors get a meal
              ! - rzoophilic proportion of these on humans 
              ! - zhostinfect of the people are infected
              ! - rpt* prob of transmission
              !   
              ! zeroth box movement depends on biting ratio
              ! - will possibly adjust the survival for this category 
              !   as feeding is dangerous.
              ! - this will depend on population of vector and host - 4now constant - 
              !  

              !----------------------------------------------
              ! 3.1 set up parameters needed in this timestep
              !----------------------------------------------

              ! ---------
              ! host2vect:
              ! ---------
              ! transmission from the host to vector - 
              ! Assume each biting vector gets one blood meal
              !
              ! v1.4 corrects a bug that was introduced with
              !      zhighrisk.  Categories I and R are considered to
              !      be at higher risk, S are more likely to be safer
              !      Previously this was accounted for for vect2host transmission
              !      But not host2vect.  Corrected in v1.4
              !   
              ! zprobhost2vect_I: probability of I getting a bite that leads host2vector trans
              !
              ! zprobhost2vect_R: probability of R getting a bite that leads host2vector trans
              ! 
              ! zprobhost2vect : overall probability of vector acquiring
              !  parasite accounting for highrisk, zoo factor and categories - sum of above two.
              !  categories as they already account for proportions and risk
              !
              zdel=zhost_SE+rbitehighrisk*(zhost_I+zhost_R) ! repeated RHS bot
              zprobhost2vect_I=rpthost2vect_I*rbitezoo(ix,iy)*zhost_I*rbitehighrisk/zdel
              zprobhost2vect_R=rpthost2vect_R*rbitezoo(ix,iy)*zhost_R*rbitehighrisk/zdel
              zprobhost2vect=zprobhost2vect_I+zprobhost2vect_R ! total prob

              IF(lverbose)PRINT *,'prob',rbitezoo(ix,iy),rpthost2vect_I,rpthost2vect_R,zprobhost2vect
              ! ---------
              ! vect2host:
              ! ---------
              !
              ! mean daily bites per host for non bednet - taking bednetuse into account
              ! recall nn=no net, as the calculation only for people not under net
              !
              znnhbr=rbiteratio*rzoophilic(ix,iy)* &
                & SUM(rvect(0,:,ix,iy))/(rpopdensity(ix,iy)*rnobednetuse)
              IF(loutput_hbr) rdiag2d(ix,iy,ncvar_hbr(2))=znnhbr

              ! daily EIR 
              znndailyeir=zvectinfectbite*znnhbr
              IF(loutput_eir) rdiag2d(ix,iy,ncvar_eir(2))=znndailyeir

              !
              ! Treat biting as a random event, therefore is distributed as a Poisson PDF.
              ! however, only need to call if vector infection rate greater than zero
              ! 
              zprobvect2hostS=0.0
              zprobvect2hostEI=0.0
              IF (znndailyeir>reps) THEN
                 ! we could speed this up if rbitehighrisk=1
                 ! highrisk refers to average number of bites on 
                 znndailyeirS=znndailyeir/(1.0+(rbitehighrisk-1.0)*zpr(ix,iy))

                 CALL transmission(znndailyeirS,zprobvect2hostS)

                 znndailyeirEI=znndailyeirS*rbitehighrisk
                 CALL transmission(znndailyeirEI,zprobvect2hostEI)
              ENDIF

              !--------------------------------
              ! Vector to human transmission
              !--------------------------------
              IF (ltimer) CALL timer('tran',icheck,time1,time2) ! cpu timer

              !--------------------------------------------
              ! SEIR model for Vector to human transmission
              !--------------------------------------------
              ! In all releases v1.X this is a simple SEI model 
              ! with multiple E boxes to model the delay correctly.
              ! From v1.4 the SEI model was extended to SEIR to add immunity 
              ! and cases diagnostic for E->I transition recorded.

              !-----------------------------------------
              ! 3.2 Evolution of disease in host S->E->I
              !-----------------------------------------
              ! record pre-transition I totals.
              zcases=SUM(rhost(ninfh,:,ix,iy)) ! this would need to be a weighted mean if nhost>1
              zhostf=dt/ninfh
              DO ihost=1,nhost
                 ! from v1.3.5 this prob is separate for S class...
                 ! immunity: Only advect first ninfh boxes
                 ! Transition S->E->I
                 CALL advection(zprobvect2hostS,zhostf,rhost(0:ninfh,ihost,ix,iy),ninfh) 
              ENDDO

              ! calculate new I arrivals.
              zcases=MAX(0.0,SUM(rhost(ninfh,:,ix,iy))-zcases) ! this would need to be a weighted mean if nhost>1
              IF (loutput_cases) rdiag2d(ix,iy,ncvar_cases(2))=zcases

              !-----------------------
              ! 3.5 Immunity Gain I->R 
              !-----------------------
              IF (rimmune_gain_tau>0.0) THEN ! if gain tau=0 immunity switched off
                 DO ihost=1,nhost
                    !    Pimmune=1.0-EXP(EIRd/EIRTAU)
                    zdel=MAX(0.0,1.0-EXP(-znndailyeir/rimmune_gain_tau))*rhost(ninfh,ihost,ix,iy)
                    ! Transition I->R
                    rhost(ninfh+1,ihost,ix,iy)=rhost(ninfh+1,ihost,ix,iy)+zdel
                    rhost(ninfh,ihost,ix,iy)=rhost(ninfh,ihost,ix,iy)-zdel
                 ENDDO

                 !-----------------------
                 ! 3.6 Immunity Loss R->S 
                 !-----------------------
                 ! As loss rate is >> timestep can use explicit numerics
                 DO ihost=1,nhost
                    zdel=rhost(ninfh+1,ihost,ix,iy)*dt/rimmune_loss_tau

                    ! Transition R->S
                    rhost(ninfh+1,ihost,ix,iy)=rhost(ninfh+1,ihost,ix,iy)-zdel
                    rhost(0,ihost,ix,iy)=rhost(0,ihost,ix,iy)+zdel 
                 ENDDO
              ENDIF ! immunity switch
              IF (loutput_immunity) rdiag2d(ix,iy,ncvar_immunity(2))=SUM(rhost(ninfh+1,:,ix,iy))/SUM(rhost(:,:,ix,iy))

              !-----------------------------------------------------
              ! 3.7 clearing rate: I->S
              !     from 1.4.0 only apply to ninfh
              !-----------------------------------------------------
              DO ihost=1,nhost
                 ! As clearing rate is >> timestep can use explicit numerics
                 zdel=dt/rhostclear*rhost(ninfh,ihost,ix,iy) 
                 ! Transition I->S
                 rhost(ninfh,ihost,ix,iy)=rhost(ninfh,ihost,ix,iy)-zdel
                 rhost(0,ihost,ix,iy)=rhost(0,ihost,ix,iy)+zdel
              ENDDO

              !-----------------------------------------------------
              ! 3.8 Death/Replacement of population: E,I,R->S
              !-----------------------------------------------------
              DO ihost=1,nhost
                 DO iinfh=1,ninfh+1
                    zdel=rpop_death_rate_daily*rhost(iinfh,ihost,ix,iy)
                    rhost(iinfh,ihost,ix,iy)=rhost(iinfh,ihost,ix,iy)-zdel
                    rhost(0,ihost,ix,iy)=rhost(0,ihost,ix,iy)+zdel
                 ENDDO
              ENDDO

              !---------------------
              ! 4. Sporogonic cycle
              !---------------------

              ! 4.1 transmission of parasite
              !     now done separately
              ! v1.4 bug correct - rbiteratio added to transmission probability:
              zdel=rbiteratio*zprobhost2vect*rvect(0,0,ix,iy)
              rvect(0,0,ix,iy)=rvect(0,0,ix,iy)-zdel
              rvect(0,1,ix,iy)=rvect(0,1,ix,iy)+zdel

              ! degree day concept of Detinova (1962):
              zsporof=dt*(ztemp-rtsporo)/dsporo
              zsporof=MIN(MAX(0.0,zsporof),1.0)

              IF (lverbose) THEN
                PRINT *,'prob ',zprobhost2vect,zdel,zsporof
                PRINT *,'vect tot',SUM(rvect)
                PRINT *,'cspr ',rdiag2d(:,:,ncvar_cspr(2))
                PRINT *,'PRd', rdiag2d(:,:,ncvar_prd(2))
                PRINT *,'EIRa ',365*znndailyeir
                PRINT *,'water frac ',rwaterfrac(ix,iy)
              ENDIF
             

              IF (zsporof>reps) THEN
                 ! need temprary array since using second index...
                 ALLOCATE(zarray1(0:ninfv)) ! before was 0:
                 ! NOTE: ivect=0 needs special treatment as is the 
                 !       only box where infection can start
                 !       IF wrap around is permitted this will also 
                 !       require special treatment
                 DO ivect=0,ngono
                    zarray1=rvect(ivect,0:ninfv,ix,iy)
!                    SELECT CASE(ivect)
!                    CASE(0)
                       !          CALL advection(zprobhost2vect,zsporof,zarray1,ninfv)
!                       CALL advection(0.0,zsporof,zarray1,ninfv)
!                    CASE DEFAULT
                       CALL advection(0.0,zsporof,zarray1,ninfv)
!                    END SELECT
                    rvect(ivect,0:ninfv,ix,iy)=zarray1
                 ENDDO
                 DEALLOCATE(zarray1)
              ENDIF

              IF (ltimer) CALL timer('spor',icheck,time1,time2) ! cpu timer

              !---------------------
              ! 5. Gonotrophic cycle
              !---------------------
              !
              ! note that we do not add a day for nighttime feeding (Ermert 10), NO LONGER TRUE 
              ! however, with a less than daily timestep, the gonotropic cycle should 
              ! always be rounded UP to the nearest integer number of days, i.e.
              ! egg-laying occurs in the nighttime timestep, since feeding and vector 
              ! oviposition only occurs at night... 
              !
              ! all temperature relationships will be modified to implement a 
              ! gaussian temperature spread to account for variations within the 
              ! grid cell due to topography and microclimates.  this effectively 
              ! adds a diffusive term to the parasite and eggs progressions.
              !
              ! Using degree day concept of Detinova (1962)
              zgonof=dt*(ztemp-rtgono)/dgono

              ! Here we want to add one resting day according to Hoshen and Morse
              ! But we also include a stochastic fluctuation of uniform distribution 
              ! and a width of one day.  This is due to the fact that the gonotrophic
              ! cycle length is poorly resolved by a daily timestep.  Without this term
              ! poor numerical behaviour can occur if ngono is set to less than 100
              ! With the stochastic term, the EIR-Temperature behaviour is smooth with ngono=25 
              zgonof=(2.0*rgono_stoch*RAND()+1.0-rgono_stoch)*zgonof/(zgonof+0.5+RAND())     
              zgonof=MIN(MAX(0.0,zgonof),1.0)
              IF(zgonof>reps) THEN
                 DO iinfv=0,ninfv
                    CALL advection(rbiteratio,zgonof,rvect(:,iinfv,ix,iy),ngono)
                 ENDDO
              ENDIF

              IF (ltimer) CALL timer('gono',icheck,time1,time2) ! cpu timer

              !-----------------------------------------
              ! 6. Vector temperature dependent survival
              !-----------------------------------------
              ! Martins I
              ! Martins II 
              ! Bayoh Scheme
              !
              ! See PhD Thesis of Ermert 2010
              !
              !-----------------------------------------

              IF (ztemp>rtvecsurvmin.AND.ztemp<rtvecsurvmax) THEN

                 SELECT CASE(nsurvival_scheme)
                 CASE(1) ! martins I
                    zsurvp=rmar1(0) + rmar1(1)*ztemp + rmar1(2)*ztemp**2
                 CASE(2) ! martins II
                    zsurvp=EXP(-1.0/(rmar2(0)+rmar2(1)*ztemp+rmar2(2)*ztemp**2))
                 CASE(3) ! Bayoh scheme 
                    zsurvp= -2.123e-7*ztemp**5 &
                         & +1.951e-5*ztemp**4 &
                         & -6.394e-4*ztemp**3 &
                         & +8.217e-3*ztemp**2 &
                         & -1.865e-2*ztemp + 7.238e-1
                 CASE DEFAULT
                    STOP 'invalid survival scheme'
                 END SELECT

                 ! from v1.3.5 base mortality rate added
                 zsurvp=zsurvp*rvecsurv
                 CALL RANDOM_NUMBER(zdel)
                 zsurvp=zsurvp+0.1*zdel-0.05 ! stochastic noise
                 zsurvp=MIN(MAX(0.0,zsurvp),1.0)
                 zsurvp=zsurvp**dt ! adjust by timestep:
              ELSE
                 zsurvp=0.0
              ENDIF

              rvect(:,:,ix,iy)=zsurvp*rvect(:,:,ix,iy)

              !-----------------------------------------
              ! 7. Vector oviposition
              !-----------------------------------------
              ! mean eggs per laying - eventually this will draw from a sample
              IF (ltimer) CALL timer('fuzz',icheck,time1,time2) ! cpu timer

              zlimit=1.0
              IF (nlayingmax>0) zlimit=MIN(1.0,nlayingmax/SUM(rvect(ngono,:,ix,iy)))

              rlarv(0,ix,iy)=0.0
              DO iinfv=0,ninfv
                 rlarv(0,ix,iy)=rlarv(0,ix,iy)+neggmn*zlimit*rvect(ngono,iinfv,ix,iy)

                 ! move the female to the meal searching box zero
                 rvect(0,iinfv,ix,iy)=rvect(0,iinfv,ix,iy)+rvect(ngono,iinfv,ix,iy)
                 rvect(ngono,iinfv,ix,iy)=0
              ENDDO

              IF (ltimer) CALL timer('eggs',icheck,time1,time2) ! cpu timer

              !-----------------------------------------
              ! 8. Puddle model
              !-----------------------------------------
              ! puddle model, increases surface area at rate related to rainfall

              ! This is not clean, but is a quick solution to a serious bug pre-1.3.3
              ! whereby waterfrac is the archived/initialized variable but rpuddle was the prognosed variable!
              rpuddle(ix,iy)=rwaterfrac(ix,iy)-rwaterperm(ix,iy)

              SELECT CASE(ipud_vers)
              CASE(125) !v1.25
                 STOP '125 removed'
              CASE(126)
                 zpud1=rwaterfrac_rate*dt
                 rpuddle(ix,iy)=(rpuddle(ix,iy)+zpud1*zrain*rwaterfrac_max)/(1.0+zpud1*(zrain+rwaterfrac_evap126))

              CASE(130) !v1.30
                 !-----------------------------------------          
                 ! dw/dt = K w^(-p/2) (fQ - w(E+fI))
                 !         K: pond geometry scale factor (related to S0 and h0 in geometry model)
                 !         p: pond geometry power factor (0.5-2 temporary ponds, 3-5 lakes)
                 !         f: proportion of maximum pond area factor 1-w/w_max
                 !         Q: run off - calculated from SCS formula Q=(P-0.2S)^2/(P+0.8S)
                 !         E: evaporation from ponds
                 !         I: maxmimum infiltration rate from ponds
                 !-----------------------------------------
                 STOP 'issues with 130 hydrology scheme'

              CASE DEFAULT
                 WRITE(iounit,*),&
              &'no default option for puddle - please set ipuddle correctly'
                 STOP
              END SELECT

              ! convert back to prognostic waterfrac:
              rwaterfrac(ix,iy)=rpuddle(ix,iy)+rwaterperm(ix,iy)
              rwaterfrac(ix,iy)=MIN(MAX(rwaterfrac(ix,iy),0.0),rwaterfrac_max) ! safety

              !-----------------------------------------------
              ! 9. Larvae maturation - limitation of resourses
              !-----------------------------------------------
              ! limitation of resources - negative feedback on numbers
              ! integrate biomass of larvae - bomblies instead slows grow rate 
              ! from v1.3.5 added rwateroccupancy, which gives fractional occupancy by LU type.
              zbiolimit=rbiocapacity*rwaterfrac(ix,iy) !*rwateroccupancy(ix,iy)
              zcapacity=MIN(MAX((zbiolimit-zmasslarv)/zbiolimit,0.01),1.0)

              !-------------------------------------------
              ! 10. Larvae maturation - larvae progression
              !-------------------------------------------
              ! larvae maturation progress rate  
              ! degree day concept of Detinova (1962) again:

              ! egg development and pupae development are all around one day, 
              ! therefore the immature phase length is mostly controlled by the 
              ! immature phase.
              IF (ztempwater>rlarv_tmin) THEN
!              IF (ztempwater>rlarv_tmin .AND. ztempwater <rlarv_tmax ) THEN
                 zlarvmaturef=rlarvmature(nlarv_scheme,1)*ztempwater+rlarvmature(nlarv_scheme,2)
                 !zlarvmaturef=zlarvmaturef/(zlarvmaturef*(rlarv_eggtime+rlarv_pupaetime))
                 zlarvmaturef=zlarvmaturef*dt ! timestep 
                 zlarvmaturef=MIN(MAX(0.0,zlarvmaturef),1.0)
              ELSE
                 zlarvmaturef=0.0
              ENDIF
              if (lverbose)PRINT *,'twater larvf',ztempwater,zlarvmaturef

              ! Bomblies slowed development due to overcrowding but 
              ! this causes unrealistically long development
              ! times, so instead we use the overcrowding factor 
              ! to alter the mortality rate for larvae
!!!! zlarvmaturef=zlarvmaturef*zcapacity
              IF (zlarvmaturef>reps) CALL advection(1.0, zlarvmaturef,rlarv(:,ix,iy),nlarv)

              IF (ltimer) CALL timer('larv',icheck,time1,time2) ! cpu timer

              !-----------------------------------------
              ! 11. Larvae mortality
              !-----------------------------------------
              ! larvae mortality is due to 
              ! 
              ! Dessication - if waterfrac reduces over a timestep, 
              !                  the larvae are reduced by the same frac
              !                  this of course assumes that the pool are independent
              !
              ! Water temperatures
              !
              ! Here for function of temperature we use Data from Bayoh and Lindsay 2003,2004
              !           Mean survival                         Proportion of terminal
              !            in days (95%         Range of larval  events occurring as    Equality of survival
              !Temperature
              !( C)       confidence interval) mortality (days) larval mortality (%)   distributions*
              !10           2.7 (2.6-2.8)        2-5             100.0                  a
              !12           3.7 (3.6-3.9)        1-6             100.0                  -
              !14          20.5 (19.3-21.8)      5-42            100.0                  -
              !16          25.5 (24.4-26.5)      9-39            100.0   *cycle survival*  b
              !18          24.9 (23.8-26.2)     10-38             58.0    30.9  0.972      b
              !20          24.9 (23.6-26.4)      3-31             24.7    23.0  0.987      b
              !22          18.1 (17.5-18.6)      5-20             24.0    18.3  0.985       -
              !24          16.4 (15.9-16.8)      6-18             20.7    15.3  0.985    -
              !26          13.5 (13.2-13.9)      5-15             27.3    13.0  0.976    -
              !28          11.0 (10.6-11.4)      3-14             33.3    11.4  0.965    c
              !30          11.2 (10.8-11.5)      4-16             72.7    10.1  0.879    c
              !32          10.2 (9.9-10.5)       5-13             70.0    9.1   0.876    -
              !34           8.9 (8.5-9.3)        4-14            100.0                  -
              !36           6.9 (6.8-7.2)        4-10            100.0                  -
              !38           4.8 (4.6-4.9)        3-7             100.0                  -
              !40           2.8 (2.6-2.9)        2-4             100.0                  a
              ! 
              ! The *cycle survival* is calculated from the survival and development length
              ! survival rate as a function of temperature is calculated.
              ! 
              ! We fit a scaled logistic curve  smin+(smax-smin)*(1.0-1.0/(1.0+exp((t0-tdata)/tau))) using the NLS package in R: 
              ! Smax, t0 and tau are fitted.  Smin is a fixed constant at 0.8 
              ! (you can achieve a good fit with smin at any value between 0 and 0.8, the data does not allow you to specify a first
              !  guess as it does not cover the high temperature range in a useful way).
              ! 
              ! The data is smoothed by a 3 point running mean, with the end points retained:
              ! 
              ! R code:
              ! sdata=c(0.972,0.987,0.985,0.985,0.976,0.965,0.879,0.876)
              ! sdata=c(0.972,rollmean(sdata,3),0.876)
              ! fit=nls(sdata ~ smin+(smax-smin)*(1.0-1.0/(1.0+exp((t0-tdata)/tau))) , start=list(t0=t0,tau=tau,smax=smax))
              ! 
              ! smin=0.0 gives:
              !   t0          tau        smax(=rlarvsurv)
              !   39.0608718  3.5141577  0.9871559 
              ! smin=0.8 gives:
              !     Estimate Std. Error t value Pr(>|t|)    
              ! t0   30.975279   0.310659  99.708 1.92e-09 ***
              ! tau   2.189333   0.385852   5.674  0.00237 ** 
              ! smax  0.983665   0.005228 188.164 8.04e-11 ***

              zlimit=0.7
              zsurvp=zlimit+(rlarvsurv-zlimit)*(1.0-1.0/(1.0+EXP((33.1-ztempwater)/2.75)))

              zsurvp=zsurvp*zcapacity ! Bomblies type capacity limitation
              !   cannibalism here if required
       
              !-----------------------------
              ! dessication by rain FLUSHING
              ! ----------------------------
              zflushr=(1.0-rlarv_flushmin)*EXP(-zrain/rlarv_flushtau)+rlarv_flushmin
              
              DO i=0,nlarv
                 ! greatest flushing to L1 larvae so apply as a linear function
                 zlimit=REAL(i)/REAL(nlarv)
                 zsurvpl=zsurvp*(zlimit*(1.0-zflushr)+zflushr)
                 zsurvpl=zsurvpl**dt ! timestep
                 zsurvpl=MIN(MAX(0.0,zsurvpl),1.0)
                 rlarv(i,ix,iy)=rlarv(i,ix,iy)*zsurvpl
              ENDDO

              !-----------------------------------------
              ! 12. Larvae hatching
              !-----------------------------------------
              rvect(0,0,ix,iy)=rvect(0,0,ix,iy)+rlarv(nlarv,ix,iy)
              rlarv(nlarv,ix,iy)=0.0    

              IF (ltimer) CALL timer('end1',icheck,time1,time2) ! cpu timer

              !------------------------------------------------------------------------------------
              ! 13. Migration of population and vectors
              !------------------------------------------------------------------------------------
              !     for the moment - consider a simple migration from "outside" as a trickle source
              !     a second option will be introduced to represent migration in the gridded world
              !     !------------------------------
              ! options:
              ! a. method 1: do not allow uninfected to exceed 99.5 % (remove in postprocessing)
              !
              ! b.  use the gravity model as in the cholera model 
              !     M=Pop(a).Pop(b)/(distance a-b) (set up matrix value before)
              !     but the matrix will limit maximum distance - small distances 
              !     will be set to zero since workers return to sleep in same 
              !     location each night when transmission occurs.
              !
              ! c. Full agent based model (for 2013)
              !
              !------------------------------------------------------------------------------------
              IF (rmigration>reps) THEN
                 zdel=MAX(rmigration-1.0+rhost(0,nadult_ni,ix,iy),0.0)
                 rhost(ninfh,nadult_ni,ix,iy)=rhost(ninfh,nadult_ni,ix,iy)+zdel
                 rhost(0,nadult_ni,ix,iy)    =rhost(0,nadult_ni,ix,iy)    -zdel
              ENDIF

              IF (lascii) THEN
                 WRITE(iounit,*)'ascii output no longer supported from v1.3.3'
                 lascii=.false.
              ENDIF
              IF (ltimer) CALL timer('diag',icheck,time1,time2) ! cpu timer
              !--------------------
              ! END OF SPATIAL LOOP
              !--------------------
           ELSE
              rdiag2d(ix,iy,:)=rfillvalue
           ENDIF !non-lake or sea point
        ENDDO !nlon
     ENDDO !nlat
     !------------
     ! 15. NCDF OUTPUT
     !------------
     IF (.NOT.lspinup) THEN
        CALL writedata(iday)
     ENDIF
  ENDDO ! date loop
  CALL setdown
  WRITE(iounit,*) 'integration finished'

  !---------------------------------------
CONTAINS
  SUBROUTINE timer(str,icheck,time1,time2)

    IMPLICIT NONE

    REAL :: time1,time2
    INTEGER :: icheck
    CHARACTER*4 :: str

    CALL cpu_time(time2)
    PRINT *,'Check point ',icheck,str,1000*(time2-time1)
    time1=time2
    icheck=icheck+1

    RETURN
  END SUBROUTINE timer

END PROGRAM VECTRI
