MODULE PlanetConstants

! Physical constants for Earth
! Initialized with default reference values

  IMPLICIT NONE
  PRIVATE

  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module) 
  !---------------------------------------------------------------
  PUBLIC :: getPlanetConst,setPlanetConst,gravConst

  !----------------------------------------------------------------------------------------
  !                   FUNDAMENTAL PHYSICAL CONSTANTS from NIST 01/11/2002
  !----------------------------------------------------------------------------------------
  !TYPE    | Name          | Value          | Description       |    Units
  !----------------------------------------------------------------------------------------
  REAL    :: WVMWTref  =     18.015         ![g/mole]           | Molecular mass water
  REAL    :: AIRMWTref =     28.964         ![g/mole]           | Molecular mass dry air

CONTAINS

  subroutine getPlanetConst(WVMWT, AIRMWT, XMASS_H2O, XMASS_DRY)

  ! For access to constants

    !  Arguments:
    REAL, INTENT(OUT), OPTIONAL :: WVMWT
    REAL, INTENT(OUT), OPTIONAL :: AIRMWT
    REAL, INTENT(OUT), OPTIONAL :: XMASS_H2O
    REAL, INTENT(OUT), OPTIONAL :: XMASS_DRY

    if (present(WVMWT))  WVMWT  = WVMWTref
    if (present(AIRMWT)) AIRMWT = AIRMWTref
    if (present(XMASS_H2O)) XMASS_H2O = WVMWTref*1.E-3
    if (present(XMASS_DRY)) XMASS_DRY = AIRMWTref*1.E-3

  end subroutine getPlanetConst

!-----------------------------------------------------------------------------

  subroutine setPlanetConst(WVMWT, AIRMWT)

  ! To override the default values

    !  Arguments:
    REAL, INTENT(IN), OPTIONAL :: WVMWT
    REAL, INTENT(IN), OPTIONAL :: AIRMWT

    if (present(WVMWT))  WVMWTref  = WVMWT
    if (present(AIRMWT)) AIRMWTref = AIRMWT

  end subroutine setPlanetConst

!-----------------------------------------------------------------------------

  function gravConst(LATITUDE)

!   Gravitational constant for Earth in meters/s^2

    USE PhysConstants, ONLY: getPhysConst

    !  Arguments:
    REAL, INTENT(IN), OPTIONAL  :: LATITUDE    ! Latitude (degrees) for which gravitational 
                                               ! constant is desired

    REAL                        :: gravConst   ! in meters/s^2
    REAL                        :: ref_lat
    REAL                        :: PI
    REAL, PARAMETER             :: DEFAULT_LAT= 45.   ! in degrees

    call getPhysConst(PI=PI)
!         Latitude for which gravitational constant desired
    if (present (LATITUDE) ) then
       ref_lat = LATITUDE
    else
       ref_lat = DEFAULT_LAT   
    endif

    gravConst = 9.80665 - 0.02586*COS(2.0*PI*REF_LAT/180.0)

  end function gravConst

END MODULE PlanetConstants
