MODULE RTMmono

  PRIVATE :: RAD_UP_DN, bb_fn
  
  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module)
  !---------------------------------------------------------------
  PUBLIC :: RTM, calctmr 

  INTEGER, parameter :: NWNMX=10000

CONTAINS
  SUBROUTINE RTM(IOUT,IRT,NWN,WN,NLAY,T,TZ,O, &
       TMPSFC,  RUP,TRTOT,RDN,REFLC,EMISS,RAD,TB,IDU)

!
!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002 - 2009, Atmospheric & Environmental Research, Inc. (AER).|
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
!
!     PROGRAM:  RTM
!     -------
!
!     AUTHOR: Sid-Ahmed Boukabara
!     ------
!
!     AFFILIATION: ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!     -----------
!
!     DATE OF CREATION : May 1999
!     ----------------
!
!     AIM: This program is aimed at the simulation of the
!     ---- radiances using the Radiative Transfer equation.
!
!     INPUTS:
!     ------
!     - IOUT  : Flag to compute the radiances (IOUT=0) or both the
!               radiances and the brightness temperatures (IOUT=1)
!     - IRT   : Flag to compute the radiative transfer
!               IRT=1 from the surface to the space (satellite)
!               IRT=2 Limb measurements.
!               IRT=3 from the space to the surface (ground instrument)
!     - NWN   : Number of wavenumbers to be treated
!     - WN    : Vector of NWN wavenumbers [in cm-1], one should note that
!               this input could be a scalar (associated with NWN=1)
!     - NLAY  : Number of layers to be treated.
!     - P     : Vector of NLAY pressures (in mbar), one should note that
!               this input could be a scalar (associated with NLAY=1)
!     - T     : Vector of NLAY temperatures (layers averaged T) [in Kelvin]
!     - TZ    : Vector of NLAY+1 temperatures (levels T) [in Kelvin]
!     - TMPSFC: Surface temperatures [in Kelvin]
!     - O     : Total Optical depths in nepers, (NWNxNLAY)
!     - REFLC : Reflectivity vector of the Surface (NWN dimension)
!     - EMISS : Emissivity vector of the surface (NWN dimension)
!     - IDU   : Index for the Up/Down format of the profiles
!               IDU=0->the profiles are given from the top to the surface
!               IDU=1->the profiles are given from the surface to the top

!     Note:
!             If the surface was in thermodynamical equilibrium,
!             REFLC should be equal to (1-ESFC)
!
!     OUTPUTS:
!     -------
!     - RAD    : An array of NWN elts containing the radiances at the WN
!                frequencies.
!     - TB     : An array of NWN elts containing the brightness temperatures
!                at the WN frequencies.
!     - RUP    : Upwelling Contribution Radiance (NWN dimension)
!     - RDN    : DownWelling Contribution Radiance (NWN dimension)
!     - TRTOT  : Total Transmittance
!
!       Note:
!       -----
!       RTM takes into account the cosmic background contribution.
!       The cosmic radiation is hard coded (2.75 Kelvin).
!
!-------------------------------------------------------------------------------
  USE PhysConstants, ONLY: getPhysConst
  !include "declar.incl"
  USE lblparams, ONLY: MXLAY
  INTEGER NWN,NLAY,IRT,I,IOUT,IDU
  REAL RADCN1,RADCN2
  REAL*8 VV, WN(NWNMX)
  CHARACTER HVRSUB*15
  REAL TMPSFC,ESFC,RSFC,SURFRAD,ALPH,COSMOS,TSKY
  REAL O(:,:)
  REAL P(MXLAY), T(MXLAY), TZ(0:MXLAY)
  REAL, DIMENSION(:) :: RAD,EMISS,REFLC,RUP,TRTOT,TB,RDN
  REAL fbeta,beta,X
  COMMON /CVRSUB/ HVRSUB

  HVRSUB = '$Revision: 19812 $'

  call getPhysConst(RADCN1=RADCN1,RADCN2=RADCN2)

!---Up and Down radiances

  CALL RAD_UP_DN(T,NLAY,TZ,WN,rup,trtot,O,rdn,NWN,IDU,IRT)

!---RADIATIVE TRANSFER
  TSKY=2.75 !Cosmic background in Kelvin
  beta= RADCN2/TMPSFC
  alph= RADCN2/TSKY

  if (irt.eq.3) then
     print *, '     '
     print *, '***********************************'
     print *, &
          'NB: for Downwelling Radiance the Boundary is ', &
          'Internally Set to the Cosmic Value: 2.75K'

     print *, '***********************************'

  endif

  DO I=1,NWN
     vv = wn(i)
     SURFRAD  = bb_fn(vv,beta)
     COSMOS   = bb_fn(vv,alph)
     ESFC=EMISS(I)
     RSFC=REFLC(I)
!
!    Upwelling Case
     IF (IRT.EQ.1) RAD(I) = RUP(I) + &
         trtot(i) * (esfc*SURFRAD + rsfc*(rdn(i)+trtot(i)*COSMOS)) ! kcp 09/21/07

!
!    Limb Case      trtot is taken as the transmittance from the tangent point to h1 (SAC)
     IF (IRT.EQ.2) RAD(I) = RUP(I) + &
         trtot(i) * (rdn(i)+trtot(i)*COSMOS)
!
!    Downwelling Case
     IF (IRT.EQ.3) RAD(I)=RDN(I)+(trtot(i)*COSMOS)
!
     IF (IOUT.EQ.1) THEN
        X=RADCN1*(WN(I)**3)/RAD(I)+1.
        TB(I)=RADCN2*WN(I)/log(X)
     ENDIF
  ENDDO
  RETURN
  END

  SUBROUTINE RAD_UP_DN(T,nlayer,TZ,WN,rup,trtot,O,rdn,NWN,IDU,IRT)
  USE PhysConstants, ONLY: getPhysConst
  USE lblparams, ONLY: MXLAY
 
  IMPLICIT REAL*8 (V)
  !include "declar.incl"
  REAL T(MXLAY), TZ(0:MXLAY)
  REAL RADCN1,RADCN2
  INTEGER  layer,nlayer,NWN,IDU,lmin,lmax,nl
  REAL          beta,beta_a,bb,bba
!---local variables
  REAL          bbVEC(MXLAY),bbaVEC(0:MXLAY),ODTOT(NWN)
  REAL  O(:,:)
  REAL*8 WN(NWNMX)
  REAL, DIMENSION(:) :: RUP,RDN,TRTOT

  call getPhysConst(RADCN1=RADCN1,RADCN2=RADCN2)
  IF (IDU.NE.1) STOP 'ERROR IN IDU. OPTION NOT SUPPORTED YET'
  lmin=nlayer
  lmax=1
  nl=-1
  DO 60 I=1,NWN
     VV=WN(I)
     rup(I) = 0.
     rdn(I) = 0.
     trtot(I) = 1.
     ODTOT(I)=0.
     DO layer=1,nlayer
        ODTOT(I)=ODTOT(I)+O(I,layer)
        beta  = radcn2/t(layer)
        beta_a= radcn2/tz(layer)
        bbVEC(layer)  = bb_fn(VV,beta)
        bbaVEC(layer) = bb_fn(VV,beta_a)
        beta_a= radcn2/tz(layer-1)
        bbaVEC(layer-1) = bb_fn(VV,beta_a)
     ENDDO

     IF (IRT.NE.3) THEN !compute RUP only when IRT<>3
        ODT=ODTOT(I)
        DO 70 layer = 1,nlayer,1
           bb  = bbVEC(layer)
           bba = bbaVEC(layer)
           ODVI = O(I,layer)
           TRI = EXP(-ODVI)
           ODT=ODT-ODVI
           TRTOT(I)= EXP(-ODT)
           pade=0.193*ODVI+0.013*ODVI**2
           RUP(I)= RUP(I)+TRTOT(I)*(1.-TRI)*(bb+pade*bba)/(1.+pade)
 70     ENDDO
     ENDIF

     ODT=ODTOT(I)
     do 50 layer = nlayer,1,-1
        bb  = bbVEC(layer)
        bba = bbaVEC(layer-1)
        ODVI = O(I,layer)
        ODT=ODT-ODVI
        TRI = EXP(-ODVI)
        TRTOT(I)= EXP(-ODT)
        pade=0.193*ODVI+0.013*ODVI**2
        RDN(I)= RDN(I)+TRTOT(I)*(1.-TRI)*(bb+pade*bba)/(1.+pade)
 50  continue
     TRTOT(I)=EXP(-ODTOT(I))
 60 ENDDO
  RETURN
  END

        function bb_fn(v,fbeta)

          USE PhysConstants, ONLY: getPhysConst
          ! Arguments
          real*8, intent(in)  :: v
          real, intent(in)  :: fbeta
          real              :: bb_fn

          ! Variable
          real RADCN1

          call getPhysConst(RADCN1=RADCN1)
          bb_fn = radcn1*(v**3)/(exp(v*fbeta)-1.)

        end function bb_fn

      subroutine calctmr(nlayrs, nwn, wn, T, tz, O, tmr)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Author: Dave Turner, January 2008
!
!This routine computes the mean radiating temperature from the optical depth
! and temperature profiles.  The logic was provided by Vivienne Payne, AER,
! in an email on 10 Jan 2008.
!
! This routine is currently not connected to the body of the code, but can
! easily be called from the subroutine rtm.
!
! Results have been checked against results from a modified version of rtm.f
! that Tony !lough had supplied to Jim Liljegren and Maria Caddedu at ANL
! pre-2006 for the purposes of calculations for Jim''s statistical retrievals.
!
! Note that this routine is currently only applicable to downwelling
! calculations
!
! Inputs:
!    nlayrs:    The number of layers in the model atmosphere
!    nwn:       The number of spectral channels
!    wn:        The wavenumber array, in cm-1 (NWN)
!    T:         The layer-averaged temperature profile, in K (NLAYRS)
!    O:         The optical depth data, in nepers (NWN x NLAYRS)
! Output:
!    Tmr:       The mean radiating temperature spectrum, in K (NWN)
!
!  Vivienne Payne, AER Inc, 2008
!------------------------------------------------------------------------------
        USE PhysConstants, ONLY: getPhysConst
        !include "declar.incl"
        USE lblparams, ONLY: MXLAY

        real*8  wn(NWNMX), vv
        real    t(mxlay),tz(0:mxlay),o(:,:)
        integer nlayrs, nwn

        real    tmr(:)

        integer ifr, ilay
        real    sumtau, sumexp
        real    bbvec(MXLAY),bbavec(0:MXLAY),odtot(NWN),trtot(NWN)
        real    odt, odvi,  beta, beta_a
        real    radtmr, x
        REAL    RADCN1,RADCN2

        call getPhysConst(RADCN1=RADCN1,RADCN2=RADCN2)

        do ifr=1,nwn
            sumtau = 0.
            sumexp = 0.
            vv=wn(ifr)
            trtot(ifr) = 1.
            odtot(ifr)=0.
            do ilay=1,nlayrs
                odtot(ifr)=odtot(ifr)+O(ifr,ilay)
                beta  = radcn2/t(ilay)
                beta_a= radcn2/tz(ilay)
                bbvec(ilay)  = bb_fn(VV,beta)
                bbavec(ilay) = bb_fn(VV,beta_a)
                beta_a= radcn2/tz(ilay-1)
                bbavec(ilay-1) = bb_fn(VV,beta_a)
            enddo

            odt=odtot(ifr)
            do ilay = nlayrs,1,-1
                bb  = bbvec(ilay)
                bba = bbavec(ilay-1)
                odvi = O(ifr,ilay)
                odt=odt-odvi
                tri = exp(-odvi)
                trtot(ifr)= EXP(-odt)
! calculate the "effective emissivity" of the layer using "linear in tau"
! (see Clough et al 1992)
                pade=0.193*odvi+0.013*odvi**2
                beff = (bb + pade*bba)/(1.+pade)
                sumexp = sumexp + beff*trtot(ifr)*(1-tri)
            enddo

! this bit is based on Han & Westwater (2000) eq 14
            radtmr = sumexp / (1. - exp(-1*odtot(ifr)))
            x=radcn1*(wn(ifr)**3)/radtmr+1.
            tmr(ifr) = radcn2*wn(ifr) / log(x)

        enddo

        return
       end

END MODULE RTMmono

