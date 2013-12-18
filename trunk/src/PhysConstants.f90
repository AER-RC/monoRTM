MODULE PhysConstants

! Physical constants
! Initialized with default reference values

  IMPLICIT NONE
  PRIVATE

  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module) 
  !---------------------------------------------------------------
  PUBLIC :: getPhysConst,setPhysConst

  !----------------------------------------------------------------------------------------
  !                       MATHEMATICAL CONSTANTS
  !----------------------------------------------------------------------------------------
  !TYPE    | Name        | Value               | Units
  !----------------------------------------------------------------------------------------
  REAL    :: PIref       = 3.1415926535898     ! [rad]

  !----------------------------------------------------------------------------------------
  !                   FUNDAMENTAL PHYSICAL CONSTANTS from NIST 01/11/2002
  !----------------------------------------------------------------------------------------
  !TYPE    | Name          | Value          |    Units         | Description             
  !----------------------------------------------------------------------------------------
  REAL    :: PLANCKref     = 6.62606876E-27 ! [g.cm^2/s]       | Planck constant          
  REAL    :: BOLTZref      = 1.3806503E-16  ! [g.cm^2/(s^2.K)] | Boltzman Constant        
  REAL    :: CLIGHTref     = 2.99792458E+10 ! [cm/s]           | Speed of Light           
  REAL    :: AVOGADref     = 6.02214199E+23 ! [molec./mole]    | Avogadro's Number        
  REAL    :: ALOSMTref     = 2.6867775E+19  !                  

  !----------------------------------------------------------------------------------------
  !                   DERIVED PHYSICAL CONSTANTS from NIST 01/11/2002
  !----------------------------------------------------------------------------------------
  !TYPE    | Name        | Value           |    Units               | Description
  !----------------------------------------------------------------------------------------
  REAL    :: GASCONref   = 8.314472E+07    !  [g.cm^2/(s^2.K.mole)] | Gas constant
  REAL    :: RADCN1ref   = 1.191042722E-12                          ! 1st radiation constant
  REAL    :: RADCN2ref   = 1.4387752                                ! 2nd radiation constant
  !       The first and second radiation constants are taken from NIST.
  !       They were previously obtained from the relations:
  !       RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07      
  !       RADCN2 = PLANCK*CLIGHT/BOLTZ 

CONTAINS

  subroutine getPhysConst(PI, PLANCK ,BOLTZ, CLIGHT ,AVOGAD ,ALOSMT, &
     GASCON, RADCN1, RADCN2)

  ! For access to constants

    !  Arguments:
    REAL, INTENT(OUT), OPTIONAL :: PI
    REAL, INTENT(OUT), OPTIONAL :: PLANCK 
    REAL, INTENT(OUT), OPTIONAL :: BOLTZ 
    REAL, INTENT(OUT), OPTIONAL :: CLIGHT 
    REAL, INTENT(OUT), OPTIONAL :: AVOGAD 
    REAL, INTENT(OUT), OPTIONAL :: ALOSMT 
    REAL, INTENT(OUT), OPTIONAL :: GASCON
    REAL, INTENT(OUT), OPTIONAL :: RADCN1
    REAL, INTENT(OUT), OPTIONAL :: RADCN2

    if (present(PI))         PI         = PIref
    if (present(PLANCK))     PLANCK     = PLANCKref
    if (present(BOLTZ))      BOLTZ      = BOLTZref
    if (present(CLIGHT))     CLIGHT     = CLIGHTref
    if (present(AVOGAD))     AVOGAD     = AVOGADref
    if (present(ALOSMT))     ALOSMT     = ALOSMTref
    if (present(GASCON))     GASCON     = GASCONref
    if (present(RADCN1))     RADCN1     = RADCN1ref
    if (present(RADCN2))     RADCN2     = RADCN2ref

  end subroutine getPhysConst

!-----------------------------------------------------------------------------

  subroutine setPhysConst(PI, PLANCK, BOLTZ ,CLIGHT ,AVOGAD ,ALOSMT, &
     GASCON, RADCN1, RADCN2)

  ! To override the default values

    !  Arguments:
    REAL, INTENT(IN), OPTIONAL :: PI
    REAL, INTENT(IN), OPTIONAL :: PLANCK 
    REAL, INTENT(IN), OPTIONAL :: BOLTZ 
    REAL, INTENT(IN), OPTIONAL :: CLIGHT 
    REAL, INTENT(IN), OPTIONAL :: AVOGAD 
    REAL, INTENT(IN), OPTIONAL :: ALOSMT 
    REAL, INTENT(IN), OPTIONAL :: GASCON
    REAL, INTENT(IN), OPTIONAL :: RADCN1
    REAL, INTENT(IN), OPTIONAL :: RADCN2

    if (present(PI))         PIref         = PI
    if (present(PLANCK))     PLANCKref     = PLANCK
    if (present(BOLTZ))      BOLTZref      = BOLTZ
    if (present(CLIGHT))     CLIGHTref     = CLIGHT
    if (present(AVOGAD))     AVOGADref     = AVOGAD
    if (present(ALOSMT))     ALOSMTref     = ALOSMT
    if (present(GASCON))     GASCONref     = GASCON
    if (present(RADCN1))     RADCN1ref     = RADCN1
    if (present(RADCN2))     RADCN2ref     = RADCN2

  end subroutine setPhysConst

END MODULE PhysConstants
