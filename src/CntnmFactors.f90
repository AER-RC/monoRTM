MODULE CntnmFactors

! Processing continuum scale factors

  IMPLICIT NONE
  PRIVATE

  !---------------------------------------------------------------
  !  List of Public entities
  !---------------------------------------------------------------
  PUBLIC :: CntnmFactors_t,setCntnmFactors,pushCntnmFactors, &
            oneMolecCntnm,applyCntnmCombo

  !-----------------------------------------------------------------------------
  !  Derived types
  !-----------------------------------------------------------------------------
  TYPE CntnmFactors_t
    REAL :: xself, xfrgn, xco2c, xo3cn, xo2cn, xn2cn, xrayl
  END TYPE CntnmFactors_t

  !-----------------------------------------------------------------------------
  !  Global constants
  !-----------------------------------------------------------------------------
  TYPE(CntnmFactors_t), PARAMETER :: allOneCntnm &
                                      = CntnmFactors_t(1.,1.,1.,1.,1.,1.,1.)
  TYPE(CntnmFactors_t), PARAMETER :: allZeroCntnm &
                                      = CntnmFactors_t(0.,0.,0.,0.,0.,0.,0.)

CONTAINS

!-------------------------------------------------------------------------------

  subroutine setCntnmFactors(cFactors,XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL)

  ! For set all continuum factors equal to 1

    !  Arguments:
    TYPE(CntnmFactors_t), INTENT(INOUT) :: cFactors
    REAL, OPTIONAL,       INTENT(IN)    :: XSELF
    REAL, OPTIONAL,       INTENT(IN)    :: XFRGN
    REAL, OPTIONAL,       INTENT(IN)    :: XCO2C
    REAL, OPTIONAL,       INTENT(IN)    :: XO3CN
    REAL, OPTIONAL,       INTENT(IN)    :: XO2CN
    REAL, OPTIONAL,       INTENT(IN)    :: XN2CN
    REAL, OPTIONAL,       INTENT(IN)    :: XRAYL

    IF (PRESENT(XSELF)) THEN
      cFactors%xself = XSELF
    END IF
    IF (PRESENT(XFRGN)) THEN
      cFactors%xfrgn = XFRGN
    END IF
    IF (PRESENT(XCO2C)) THEN
      cFactors%xco2c = XCO2C
    END IF
    IF (PRESENT(XO3CN)) THEN
      cFactors%xo3cn = XO3CN
    END IF
    IF (PRESENT(XO2CN)) THEN
      cFactors%xo2cn = XO2CN
    END IF
    IF (PRESENT(XN2CN)) THEN
      cFactors%xn2cn = XN2CN
    END IF
    IF (PRESENT(XRAYL)) THEN
      cFactors%xrayl = XRAYL
    END IF

  end subroutine setCntnmFactors

!-------------------------------------------------------------------------------

  subroutine pushCntnmFactors(cFactors)

! Load continuum scale factors into the common block used by subroutine contnm 

    !  Arguments:
    TYPE(CntnmFactors_t), INTENT(IN) :: cFactors
     
    REAL ::         XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
    common /CNTSCL/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL

    XSELF = cFactors%xself
    XFRGN = cFactors%xfrgn
    XCO2C = cFactors%xco2c
    XO3CN = cFactors%xo3cn
    XO2CN = cFactors%xo2cn
    XN2CN = cFactors%xn2cn
    XRAYL = cFactors%xrayl

  end subroutine pushCntnmFactors

!-------------------------------------------------------------------------------

  subroutine oneMolecCntnm(molec,cFactors)

! Zeroes out all continua EXCEPT for molecule molec in the common block used by
! subroutine contnm 
! Continuua can subsequently be restored by calling pushCntnmFactors(cFactors), 
! where the actual argument corresponding to dummy argument cFactors is the same
! for both calls

    !  Arguments:
    INTEGER,              INTENT(IN) :: molec       ! index of molecule in HITRAN
    TYPE(CntnmFactors_t), INTENT(IN) :: cFactors

    !  Local constants:
    INTEGER, PARAMETER :: idH2O=1, idCO2=2, idO3=3, idO2=7, idN2=22, idRAYL=99  ! HITRAN indices

    !  Local variables:
    TYPE(CntnmFactors_t) :: myFactors

    myFactors = allZeroCntnm

    !! leave Rayleigh off to avoid double counting
    !myFactors%xrayl = cFactors%xrayl  ! Leave Rayleigh unchanged
    

    select case (molec)
       case (idH2O)
          myFactors%xself = cFactors%xself
          myFactors%xfrgn = cFactors%xfrgn
       case (idCO2)
          myFactors%xco2c = cFactors%xco2c
       case (idO3)
          myFactors%xo3cn = cFactors%xo3cn
       case (idO2)
          myFactors%xo2cn = cFactors%xo2cn
       case (idN2)
          myFactors%xn2cn = cFactors%xn2cn
       case (idRayl)
          myFactors%xrayl = cFactors%xrayl
        case default
          ! leave as all zero for molecules that have no modeled continuum
    end select

    call pushCntnmFactors(myFactors)

  end subroutine oneMolecCntnm

!-------------------------------------------------------------------------------

  subroutine applyCntnmCombo(ICNTNM,cFactors)

  ! Set continuum scale factors according to an index of factor combinations
  !
  !  Continuum index definitions
  !  ---------------------------
  !  ICNTNM Value      Self     Foreign    Rayleigh     Others
  !        0            no        no          no          no
  !        1            yes       yes         yes         yes
  !        2            no        yes         yes         yes
  !        3            yes       no          yes         yes
  !        4            no        no          yes         yes
  !        5            yes       yes         no          yes

    !  Arguments:
    INTEGER,              INTENT(IN)    :: ICNTNM
    TYPE(CntnmFactors_t), INTENT(INOUT) :: cFactors

    SELECT CASE (ICNTNM)
      CASE (0)
        cFactors = allZeroCntnm
      CASE (1)
        cFactors = allOneCntnm
      CASE (2)
        cFactors = allOneCntnm
        cFactors%xself = 0.
      CASE (3)
        cFactors = allOneCntnm
        cFactors%xfrgn = 0.
      CASE (4)
        cFactors = allOneCntnm
        cFactors%xself = 0.
        cFactors%xfrgn = 0.
      CASE (5)
        cFactors = allOneCntnm
        cFactors%xrayl = 0.
      CASE (6)
        ! null case
      CASE DEFAULT
        print *,'err:[CntnmFactors::applyCntnmCombo] Invalid ICNTNM:',ICNTNM
        STOP
    END SELECT

  end subroutine applyCntnmCombo

END MODULE CntnmFactors
