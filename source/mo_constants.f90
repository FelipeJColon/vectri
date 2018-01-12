MODULE mo_constants
  !--------------------------------------------------------- 
  ! VECTRI: VECtor borne disease infection model of TRIeste.  
  !
  ! Tompkins AM. 2011, ICTP
  ! tompkins@ictp.it
  !
  ! tunable constants for model
  ! 
  ! updates: v1.4.0 Added immunity
  !
  !---------------------------------------------------------

  IMPLICIT NONE

  ! larvae stage constants
  ! ----------------------
  REAL :: neggmn=80     ! E_p in Emert - mean FEMALE egg number per vector laying event
                           ! Bomblies assumes 150 - default range 60 to 80.

  INTEGER :: nlayingmax=-99 ! Ermert/LMM limit to maximum number of vector allowed to lay
                                       ! Set to negative number to switch this option off

  ! max and min temperature outside of which no larvae survive (c) from Bayoh and Lindsay (2004)
  ! need to combine with ***death rate*** - NOTE COMPATIBLE WITH SCHEME
  REAL :: rlarv_eggtime=1.0   ! days for egg hatching
  REAL :: rlarv_pupaetime=1.0 ! days for pupae stage

  ! vector survival constants
  ! maximum temperature Jepson et al., 1947; Depinay et al., 2004 as stated in Bomblies 2008
  !  REAL :: rtlarsurvmax=40.0 ! maximum temperature (C) above which larvae dies
  
  ! The rflushmin parameter is set to the SURVIVAL average of L1 to L4 larvae
  ! according to Krijn P. Paaijmans1,2*, Moses O. Wandago2, Andrew K. Githeko3, Willem Takken (2007)
  ! 0.89 = L1 larvae flushed out at 17.5% rate after heavy rainfall, L4 at 4.8%... 
  REAL :: rlarv_flushmin=0.4 ! minimum survivial with high rate 
  REAL :: rlarv_flushtau=20  ! (mm/day) describes how flushing increases with rainrate 


  ! Larvae growth rate scheme 
  !   1: constant as in Ermert 2010
  !   2: DD - Jepson 1947 as reported by HM04
  !   3: DD - BL2003 (approximated linearly)
  !   4: DD - Craig et al. 1999 
  !
  ! slope and offsets for larvae schemes - Set by nlarv_scheme in mo_control
  !
  REAL :: rlarv_tmin=12.16 ! 16 old constant - BL=18 slope of parameterization for fractional growth
  REAL :: rlarv_tmax=38.0  ! 37 BL=34 

  ! 1: This is a constant rate - Temperature *INdependent* - 12 days 
  REAL :: rlarv_ermert=0.08333  ! units day**-1 
 
  ! 2: Jepson 1948 cited (not used?) in HM04 - T-dependent rate
  REAL :: rlarv_jepson=0.011  ! units K**-1 day**-1 

  ! 3: from Bayoh and Lindsay 2003, (see p76 Ermert2010)
  ! the BL2003 paper uses an exponential function with very large factors
  ! we approximate this with a linear function
  REAL :: rlarv_bayoh=0.005 ! units K**-1 day**-1 

  ! 4: From Craig et al 1999 - this approximates Bayoh very closely
  REAL :: rlarv_craig=0.00554 !  units K**-1 day**-1  degree day constant for growth rate

  ! larval mass and ecocapacity constants - taken from Bomblies, we assume linear growth
  REAL :: rmasslarv_stage4=0.45    ! mass in mg of stage 4 adult 
  REAL :: rbiocapacity=300.0  ! max bio carry capacity of 1m**2 of pool. 

  ! base daily survival  due to preditors and real life (non-lab) factors
  ! before v1.2.9 was fixed to 0.825 
  !   HM04 and Bomblies use 0.88, Ermert 2010 use 0.825  
  ! After v1.3.0 this base rate is further reduced by Temperature-dependent Craig et al 1999/Bayoh and Lindsay 2004 
  ! All versions then reduced this further for earlier stage larvae (linear) for flushing effect.
  REAL :: rlarvsurv=0.987  ! The new parametrization fits with 0.987 to Craig/Bayoh 
                             
  ! -----------------
  ! Surface Hydrology 
  ! -----------------
  INTEGER :: ipud_vers=126    ! puddle model version from: 125=v1.2.5, 126=v1.2.6, 130=v1.30

  ! global constants:
  REAL :: rwaterperm_default=1.e-6    ! default water fraction of permanent water bodies used if
                                      ! Sentinel maps are not read in
  REAL :: rwaterfrac_rate=1./1000.      ! increase rate (m^-1)

  ! ipud_vers=125 constants
  REAL :: rwaterfrac_itau=1.0/4.0    ! 0.1 default decay timescale (days^-1)
  REAL :: rwaterfrac_infil125=4.e-6       ! fixed reduction rate due to infiltration and evaporation 

  ! ipud_vers=126 constants
  REAL :: rwaterfrac_evap126=250  !250 - evaporation+infiltration rate in mm/day (*Lv/86400 gives units of Wm^-2)
                                  ! 1.728=50w/m2
  REAL :: rwaterfrac_max=0.2      ! max default water fraction def=0.2 - also used in v130

  ! ipud_vers=130
  REAL :: rwaterfrac_CN=80         ! 25400/CN-254 curve number 85 
  REAL :: rwaterfrac_min=1.e-6     ! min default water fraction
  REAL :: rwaterfrac_evap=5        ! evaporation rate in mm/day (*Lv/86400 gives units of Wm^-2) 1.728=50w/m2
  REAL :: rwaterfrac_infil130=500  ! maximum infiltration rate over sandy soils. mm/day
  REAL :: rwaterfrac_shapep2=0.5   ! shape factor p DIVIDED BY 2 (default p=1)
  REAL :: rwaterfrac_S !=25400/rwaterfrac_CN-254 - defined in setup
  REAL :: rwaterfrac_ia=0.1        ! from v1.3.5 initial abstraction (0.05-0.2)

  ! offset of water temperature to air temperature
  REAL :: rwater_tempoffset=0.5 

  ! proportion of time spent indoor resting 
  REAL :: rbeta_indoor=0.0  ! Use with care - test carefully.

  !------------------------
  ! Intervention strategies - DO NOT USE YET...
  !------------------------
  !
  ! bednetuse: At the moment this is constant
  !    assumes bednet are distributed at constant rate to replace old ones
  !    a fraction of these are treated, i.e. deadly to vector, this fraction is simply the ratio
  !    lifetime of treatment effectiveness / lifetime of the bednet 
  !  
  REAL :: rnobednetuse=1.0  ! proportion of host NOT sleeping under bednet
  REAL :: rbednettreat=0.0  ! proportion of bednets that are treated and result in vector death

  !
  ! biting parameters
  ! -----------------
  ! daily bite rate of vectors to initiate gonotropic cycle
  ! this fraction manage to locate a host and get a meal each day! 
  REAL :: rbiteratio=0.5         ! Success rate of finding a blood meal per day 
                                 ! Note: values > 0.9 lead to "pulsing" as the gonotrophic cycle 
                                 ! is poorly resolved by the daily timestep

  REAL :: rzoophilic_tau=30.e-6  ! e-folding rate for zoophilicity (people/m**2) 
  REAL :: rzoophilic_min=0.1     ! minimum rate for low population densities 
                                                ! (rate lower since adjust to feed on cattle) 

  ! Distribution of bites in the post v1.3.5 version: rbitehighrisk is a factor between 1.0 and infinity...
  !
  ! rbitehighrisk: the ratio of risk between the higher vulnerable population (in the EI category) and the 
  !                lower vulnerable population (the susceptibles S).  Until we get an agent based approach, 
  !                we used prior bite history as a proxy for vilnerability in this way.
  ! Thus if b is the mean bites per person, k the factor below, bs bites per person in non-vulnerable S group
  ! b= pr k bs + (1-pr) bs, thus bs=b/(1+(k-1)pr)
  !   
  ! Within each categories the bite distribution is still distributed randomly, but the overall distribution of 
  ! bites per person across the population will have a wider distribution.
  !
  ! Setting rbitehighrisk=1.0: no difference in risk, this reproduces the default model prior to v.1.3.5
  ! 
  REAL :: rbitehighrisk=5.0 ! 1.0 no effect, >1.0 those already infected are assumed to be at high risk - 

! ---------------------------
! GARKI type model parameters
! ---------------------------
  ! This has to be defined here for dependencies
  ! v.1.4.0 DO NOT CHANGE FROM 1.  further age categories will be introduced with v2.X 
  INTEGER, PARAMETER :: nhost=1    ! number of host categories in garki type model
                               
  ! number of days after infection for clearance of disease in non-immune subjects
  ! If you are running with immunity switched off then I suggest to leave at high value 90 days
  ! 
  ! A much shorter value is preferable when immunity is switched on
  ! This is quite short as it is assumed that treatment is sought for most clinical cases
  ! With treatment clearance times are on the order of 2 days
  !  e.g. Brandts, Christian H., et al. "Effect of paracetamol on parasite clearance time in 
  !  Plasmodium falciparum malaria." The Lancet 350.9079 (1997): 704-709.)
  ! although it is assumed that it may take several days before treatment is sought.
  REAL :: rhostclear=20.0

  !  REAL, DIMENSION(nhost) :: rhostclear=(/rhostclear_ni,20.0/)      

  ! malaria transmission probabilities
  !   - these can be from 0.1 to 0.6 - depends on the immunity of the host
  !   - this will be generalized once the garki model is added and will be a function of host category
  ! these values are from LMM of ermert, which in turn are from    
  REAL :: rpthost2vect_I=0.2   ! I: Infectious (nonimmune)
  REAL :: rpthost2vect_R=0.04  ! R: Recovered  (immune)

  REAL :: rptvect2host=0.3

  ! Immunity module - Very simple single compartment immunity treatment to extend model from SEI to SEIR
  !
  ! So model for humans is a SE...EIR model, with n bins for E to get delay for clinical onset
  ! The assumption is that 95% of the population gets full immunity thus EXP(EIRa/TAU)=0.05
  ! The default assumes that an annual EIRa of 100 bites 
  !   is adequate for complete (=95%) immunity following (PAPER HERE)
  ! Thus TAU=-100/ln(0.05)=33.4  - if we use this on daily step the
  !   model is timestep independent as exp(-eir/(n*tau))^n = exp(-eir/tau)
  ! Set to ZERO to switch off immunity.
  REAL :: rimmune_gain_eira=100 ! default 100 - units EIR   [ -100.0/log(0.05) ]


  ! loss of immunity 
  ! tau set is set so that 95% of people lose immunity after X years.
  ! tau= -X*365/ln(0.05), so 365 gives 3 years loss...
  REAL :: rimmune_loss_tau=365  ! units days

  !
  ! vector parasite development constants - Detinova (1962)
  !
  ! Gonotropic cycle
  REAL :: rtgono=7.7 ! threshold temperature for egg development in vector
                    ! Detinova (1962) gives 4.5 9.9 and 7.7 at RH=30-40,70-80 and 90-100
                    ! this will be implemented in an interpolation scheme
  REAL :: dgono=37.1
                    ! Detinova (1962) gives 65.4,36.5,37.1 at RH=30-40,70-80 and 90-100
                    ! this will be implemented in an interpolation scheme
  REAL :: rgono_stoch=0.1 ! stochastic fractional perturbation on gonotrophic cycle

  ! Sporogonic cycle from Detinova (1962) - 18.0 used by HM04, 16.0 by Ermert
  REAL :: rtsporo=16.0 ! threshold temperature for parasite development
  REAL :: dsporo=111.0

  ! vector mortality functions
  REAL :: rvecsurv=0.96 ! underlying survival on which temperature based function is imposed.

  ! constants for martins schemes, indices are powers of T:
  REAL, PARAMETER :: rmar1(0:2)=(/0.45, 0.054, -0.0016/)   ! martins I
  REAL, PARAMETER :: rmar2(0:2)=(/-4.4, 1.31,  -0.03/)     ! martins II (used in Craig et al. 99)

  ! make sure these are consistent with Martin's functions
  REAL, PARAMETER :: rtvecsurvmin=5.0  ! minimum temperature (C) below which vector dies
  REAL, PARAMETER :: rtvecsurvmax=39.9 ! maximum temperature (C) above which vector dies

  ! ------------------------------------------
  ! garki model constants for host development
  ! ------------------------------------------
  ! adult/child indices, ni=non-imume, im=imune
  INTEGER, PARAMETER :: nadult_ni=1, nchild_ni=2, nadult_im=3, nchild_im=4

  ! number of days MINUS ONE! after biting before a human becomes infective 
  ! 20 days in the LMM2010 - reference Ermert2011a and references therein - range 10-26d
  ! *NOTE* this has to be a parameter as it is used in mo_control to set a dimension
  !        bizarrely this produces a segmentation fault in the compiler if you make a var
  REAL, PARAMETER :: rhost_infectd=20.0

  ! day number MINUS ONE at which malaria can be detected by blood tests
  REAL, PARAMETER :: rhost_detectd=9.0   

  REAL :: rhost_infect_init=0.1  ! initial value of host reservoir

  !---------------------
  ! population constants
  !---------------------

  ! Death rate - Units fraction per day.
  ! This simply causes the death of x% of the population, who are replaced by malaria naive 
  REAL :: rpop_death_rate=0.02  

  ! Population birth rate - at the moment 
  ! REAL :: rpop_birth_rate=rpop_death_rate ! at the moment, no population growth permitted.

  REAL :: rpopdensity_min=5.e-7 ! minimum population density (per m**2, i.e. 1e-6=1 per km2) 

  ! migration level as fraction of population arriving per year with malaria infection
  ! at the moment this is treated as a simple minimum threshold
  REAL :: rmigration=0.0

  ! latitude limits 
  REAL :: rlatmax=70. ! turn the model off outside this limit for efficiency (assume malaria free)

  ! quality control parameters
  REAL :: rqualtemp_min=10.0
  REAL :: rqualtemp_max=50.0

  ! math constants
  REAL, PARAMETER :: rpi=ACOS(-1.0)     ! pi
  REAL, PARAMETER :: reps=1.0e-10       ! small positive threshold number
  REAL, PARAMETER :: r0CtoK=273.15      ! freezing point in Kelvin

END MODULE mo_constants
