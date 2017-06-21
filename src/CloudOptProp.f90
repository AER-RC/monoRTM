MODULE CloudOptProp

! Subprograms for cloud optical properties

  USE PhysConstants, ONLY: getPhysConst

  IMPLICIT NONE
  PRIVATE

  !-----------------------------------------------------------------------------
  !  List of Public subprograms (accessible from outside module) 
  !-----------------------------------------------------------------------------
  PUBLIC :: ODCLW, ODCLW_LHM, ODCLW_TKC, Forward_TKC

  !-----------------------------------------------------------------------------
  !  Set the default version used when calling ODCLW
  !-----------------------------------------------------------------------------
  INTERFACE ODCLW
     MODULE PROCEDURE ODCLW_TKC
  END INTERFACE ODCLW


  ! Global constants
  REAL, PARAMETER :: Hz_per_GHz = 1.e9

CONTAINS

  !  OD function wrapper for "Turner-Kneifel-Cadeddu" liquid water absorption model
   REAL FUNCTION ODCLW_TKC(WN,TEMP,CLW)
      ! Output optical depth of cloud liquid water (unitless)
      REAL*8, INTENT(IN) :: WN    ! Wavenumber (cm-1)
      REAL,   INTENT(IN) :: TEMP  ! Temperature (K)
      REAL,   INTENT(IN) :: CLW   ! Cloud liquid water (kg/m2 = mm)
      ! Local constants
      REAL, PARAMETER :: K_at_0C=273.15
      ! Local variables
      REAL :: tempC
      REAL :: freq
      REAL :: speedOfLight ! (cm/s)
      REAL :: absCLW
      
      ! Convert units of inputs
      CALL getPhysConst(CLIGHT=speedOfLight)
      freq = WN*speedOfLight/Hz_per_GHz
      tempC = TEMP-K_at_0C
      
      ! Obtain absorption coefficient and convert to optical depth
      CALL Forward_TKC(freq,tempC,absCLW)
      ODCLW_TKC=absCLW*CLW

      RETURN

   END FUNCTION ODCLW_TKC


  !-----------------------------------------------------------------------------

  !  Abstract:
  !     The "Turner-Kneifel-Cadeddu" liquid water absorption model (submitted to JTECH 2015).
  !     It was built using both laboratory observations (primarily at warm temperatures) and 
  !     field data observed by MWRs at multiple frequencies at supercool temperatures. The field
  !     data were published in Kneifel et al. JAMC 2014.  The strength of the TKC model is the 
  !     use of an optimal estimation framework to determine the empirical coefficients of the 
  !     double-Debye model.  A full description of this model is given in 
  ! 
  !             Turner, D.D., S. Kneifel, and M.P. Cadeddu, 2015: An improved liquid
  !             water absorption model in the microwave for supercooled liquid clouds.
  !             J. Atmos. Oceanic Technol., submitted April 2015.
  ! 
  !      Note that the model is designed to operate over the frequency range from 0.5 to 500
  !      GHz, and temperatures from -40 degC to +50 degC.  
  ! 
  !  Authors:
  !     Dave Turner, National Severe Storms Laborotory / NOAA
  !     Stefan Kneifel, McGill University and the University of Cologne
  !     Maria Cadeddu, Argonne National Laboratory
  ! 
  !  Call:
   SUBROUTINE Forward_TKC(freq,temp,alpha,epsilon)
      ! Output mass absorption coefficient of cloud liquid water (m2 kg-1)
      REAL,    INTENT(IN)            :: freq    ! Frequency (GHz)
      REAL,    INTENT(IN)            :: temp    ! Temperature (degrees C)
      REAL,    INTENT(OUT)           :: alpha   ! Mass absorption coefficient of 
                                                ! cloud liquid water (m2 kg-1)
      COMPLEX, INTENT(OUT), OPTIONAL :: epsilon ! permitivity

      ! Local constants
      REAL, PARAMETER :: cm_per_m=100.
    
      !    Empirical coefficients for the TKC model. 
      REAL, PARAMETER :: a_1 = 8.110808E+01 
      REAL, PARAMETER :: b_1 = 4.433736E-03 
      REAL, PARAMETER :: c_1 = 1.301700E-13 
      REAL, PARAMETER :: d_1 = 6.627126E+02 
      REAL, PARAMETER :: a_2 = 2.025164E+00 
      REAL, PARAMETER :: b_2 = 1.072976E-02 
      REAL, PARAMETER :: c_2 = 1.011945E-14 
      REAL, PARAMETER :: d_2 = 6.089168E+02 
      REAL, PARAMETER :: t_c = 1.342433E+02 

      ! Local variables
      REAL :: frq
      REAL :: speedOfLight ! (cm/s)
      REAL :: pi
      REAL :: cl
      REAL :: eps_s
      REAL :: delta_1, delta_2
      REAL :: tau_1, tau_2
      REAL :: term1_p1, term2_p1
      REAL :: eps1, eps2
      COMPLEX :: epsilonLocal
      COMPLEX :: RE

  ! Convert the frequency from GHz to Hz
      frq = freq * Hz_per_GHz
    
      ! Some constants
      CALL getPhysConst(PI=pi,CLIGHT=speedOfLight) 
      cl =  speedOfLight/cm_per_m !speed of light in vacuum
    
      ! This helps to understand how things work below
    
      ! Compute the static dielectric permittivity (Eq 6)
      eps_s = 87.9144d0 - 0.404399d0 * temp + 9.58726D-4 * temp**2. - 1.32802D-6 * temp**3.
    
      ! Compute the components of the relaxation terms (Eqs 9 and 10)
              ! First Debye component
      delta_1 = a_1 * EXP(-b_1 * temp)
      tau_1   = c_1 * EXP(d_1 / (temp + t_c))
              ! Second Debye component
      delta_2 = a_2 * EXP(-b_2 * temp)
      tau_2   = c_2 * EXP(d_2 / (temp + t_c))
    
      ! Compute the relaxation terms (Eq 7) for the two Debye components
      term1_p1 = (tau_1**2.*delta_1) / (1.d0 + (2.d0*pi*frq*tau_1)**2.)
      term2_p1 = (tau_2**2.*delta_2) / (1.d0 + (2.d0*pi*frq*tau_2)**2.)
    
      ! Compute the real permittivitity coefficient (Eq 4)
      eps1 = eps_s - ((2.d0*pi*frq)**2.)*(term1_p1 + term2_p1) 
             
    
      ! Compute the relaxation terms (Eq 8) for the two Debye components
      term1_p1 = (tau_1 * delta_1) / (1.d0 + (2.d0*pi*frq*tau_1)**2.)
      term2_p1 = (tau_2 * delta_2) / (1.d0 + (2.d0*pi*frq*tau_2)**2.)
    
      ! Compute the imaginary permittivitity coefficient (Eq 5)
      eps2 = 2.d0*pi*frq * (term1_p1 + term2_p1)
                 
      epsilonLocal = CMPLX(eps1, eps2)
    
      ! Compute the mass absorption coefficient (Eq 1)
      RE = (epsilonLocal-1.)/(epsilonLocal+2.)
      alpha = 6.d0*pi*AIMAG(RE)*frq*1.d-3/cl

      RETURN

   END SUBROUTINE Forward_TKC


  !-----------------------------------------------------------------------------

   REAL FUNCTION ODCLW_LHM(WN,TEMP,CLW)
      !INPUTS: WN (WaveNUmber in cm-1)
      !        Temp (in K)
      !        CLW  (in mm or kg/m2)
      !OUTPUT: ODCLW: optical depth of the Cloud Liquid Water
      !FROM EQUATIONS OF LIEBE, HUFFORD AND MANABE
      !Ref:(INT. J. IR & MM WAVES V.12(17) JULY 1991
      COMPLEX EPS,RE
      REAL*8 WN
      REAL TEMP
      REAL CLW
      REAL PI,CLIGHT
      REAL :: EPS0,EPS1,EPS2
      REAL :: FREQ,FP,FS
      REAL :: THETA1
      CALL getPhysConst(PI=PI,CLIGHT=CLIGHT) 
      FREQ=WN*CLIGHT/1.E9
      IF ((FREQ.GT.3000.).AND.(CLW.GT.0.)) THEN
         WRITE(*,*) 'STOP: CLOUD IS PRESENT FOR SIMULATIONS'
         WRITE(*,*) 'IN A NON-MICROWAVE SPECTRAL REGION'
         STOP
      ENDIF
      THETA1 = 1.-300./TEMP
      EPS0 = 77.66 - 103.3*THETA1
      EPS1 = .0671*EPS0
      EPS2 = 3.52 + 7.52*THETA1
      FP = 20.1*EXP(7.88*THETA1)
      FS = 39.8*FP
      EPS = (EPS0-EPS1)/CMPLX(1.,FREQ/FP) + &
           (EPS1-EPS2)/CMPLX(1.,FREQ/FS) +EPS2
      RE = (EPS-1.)/(EPS+2.)
      ODCLW_LHM = -(6.*PI/299.792458)*CLW*AIMAG(RE)*FREQ
      RETURN
   END FUNCTION ODCLW_LHM

END MODULE CloudOptProp
