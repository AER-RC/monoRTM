MODULE RtmConstants

! Physical constants

  IMPLICIT NONE
  PRIVATE

  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module) 
  !---------------------------------------------------------------
  PUBLIC :: getRtmConst,setRtmConst

  !----------------------------------------------------------------------------------------
  !                   FUNDAMENTAL PHYSICAL CONSTANTS from NIST 01/11/2002
  !----------------------------------------------------------------------------------------
  !TYPE    | Name          | Value          | Description             |    Units
  !----------------------------------------------------------------------------------------
  REAL    :: PLANCKdfl     = 6.62606876E-27 !Planck constant          | [g.cm^2/s]
  REAL    :: CLIGHTdfl     = 2.99792458E+10 !Speed of Light           | [cm/s]
  REAL    :: BOLTZdfl      = 1.3806503E-16  !Boltzman Constant        | [g.cm^2/(s^2.K)]
  REAL    :: AVOGADdfl     = 6.02214199E+23 !Avogadro's Number        | [molec./mole]
  REAL    :: ALOSMTdfl     = 2.6867775E+19  !                         |
  REAL    :: WVMOLMASSdfl  = 18.016         !Molecular mass water     | [g/mole]
  REAL    :: DRYMOLMASSdfl = 28.964         !Molecular mass dry air   | [g/mole]

  !----------------------------------------------------------------------------------------
  !                   DERIVED PHYSICAL CONSTANTS from NIST 01/11/2002
  !----------------------------------------------------------------------------------------
  !TYPE    | Name        | Value          | Description             |    Units
  !----------------------------------------------------------------------------------------
  REAL    :: GASCONdfl   = 8.314472E+07   !Gas constant             | [g.cm^2/(s^2.K.mole)]
  REAL    :: RADCN1dfl   = 1.191042722E-12!1st radiation constant   |
  REAL    :: RADCN2dfl   = 1.4387752      !2nd radiation constant   |
  !       The first and second radiation constants are taken from NIST.
  !       They were previously obtained from the relations:
  !       RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07      
  !       RADCN2 = PLANCK*CLIGHT/BOLTZ 

  !----------------------------------------------------------------------------------------
  !                       MATHEMATICAL CONSTANTS
  !----------------------------------------------------------------------------------------
  !TYPE    | Name        | Value               | Description             |   Units
  !----------------------------------------------------------------------------------------
  REAL    :: PIdfl       = 3.1415926535897932  ! http:                   | [rad]
!                                              | //www.cecm.sfu.ca/pi9   |

CONTAINS

  subroutine getRtmConst(PLANCK, CLIGHT ,BOLTZ ,AVOGAD ,ALOSMT, &
     GASCON, WVMOLMASS, DRYMOLMASS, RADCN1, RADCN2, PI)

  ! For access to constants

    !  Arguments:
    REAL, INTENT(OUT), OPTIONAL :: PLANCK 
    REAL, INTENT(OUT), OPTIONAL :: CLIGHT 
    REAL, INTENT(OUT), OPTIONAL :: BOLTZ 
    REAL, INTENT(OUT), OPTIONAL :: AVOGAD 
    REAL, INTENT(OUT), OPTIONAL :: ALOSMT 
    REAL, INTENT(OUT), OPTIONAL :: WVMOLMASS
    REAL, INTENT(OUT), OPTIONAL :: DRYMOLMASS
    REAL, INTENT(OUT), OPTIONAL :: GASCON
    REAL, INTENT(OUT), OPTIONAL :: RADCN1
    REAL, INTENT(OUT), OPTIONAL :: RADCN2
    REAL, INTENT(OUT), OPTIONAL :: PI

    if (present(PLANCK))     PLANCK     = PLANCKdfl
    if (present(CLIGHT))     CLIGHT     = CLIGHTdfl
    if (present(BOLTZ))      BOLTZ      = BOLTZdfl
    if (present(AVOGAD))     AVOGAD     = AVOGADdfl
    if (present(ALOSMT))     ALOSMT     = ALOSMTdfl
    if (present(WVMOLMASS))  WVMOLMASS  = WVMOLMASSdfl
    if (present(DRYMOLMASS)) DRYMOLMASS = DRYMOLMASSdfl
    if (present(GASCON))     GASCON     = GASCONdfl
    if (present(RADCN1))     RADCN1     = RADCN1dfl
    if (present(RADCN2))     RADCN2     = RADCN2dfl
    if (present(PI))         PI         = PIdfl

  end subroutine getRtmConst

!-----------------------------------------------------------------------------

  subroutine setRtmConst(PLANCK, CLIGHT ,BOLTZ ,AVOGAD ,ALOSMT, &
     GASCON, WVMOLMASS, DRYMOLMASS, RADCN1, RADCN2, PI)

  ! To override the default values

    !  Arguments:
    REAL, INTENT(IN), OPTIONAL :: PLANCK 
    REAL, INTENT(IN), OPTIONAL :: CLIGHT 
    REAL, INTENT(IN), OPTIONAL :: BOLTZ 
    REAL, INTENT(IN), OPTIONAL :: AVOGAD 
    REAL, INTENT(IN), OPTIONAL :: ALOSMT 
    REAL, INTENT(IN), OPTIONAL :: WVMOLMASS
    REAL, INTENT(IN), OPTIONAL :: DRYMOLMASS
    REAL, INTENT(IN), OPTIONAL :: GASCON
    REAL, INTENT(IN), OPTIONAL :: RADCN1
    REAL, INTENT(IN), OPTIONAL :: RADCN2
    REAL, INTENT(IN), OPTIONAL :: PI

    if (present(PLANCK))     PLANCKdfl     = PLANCK
    if (present(CLIGHT))     CLIGHTdfl     = CLIGHT
    if (present(BOLTZ))      BOLTZdfl      = BOLTZ
    if (present(AVOGAD))     AVOGADdfl     = AVOGAD
    if (present(ALOSMT))     ALOSMTdfl     = ALOSMT
    if (present(WVMOLMASS))  WVMOLMASSdfl  = WVMOLMASS
    if (present(DRYMOLMASS)) DRYMOLMASSdfl = DRYMOLMASS
    if (present(GASCON))     GASCONdfl     = GASCON
    if (present(RADCN1))     RADCN1dfl     = RADCN1
    if (present(RADCN2))     RADCN2dfl     = RADCN2
    if (present(PI))         PIdfl         = PI

  end subroutine setRtmConst

END MODULE RtmConstants
