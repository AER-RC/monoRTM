!     path:		$Source$
!     author:		$Author $
!     revision:	        $Revision: 11207 $
!     created:	        $Date: 2011-03-29 13:43:38 -0400 (Tue, 29 Mar 2011) $
MODULE ModmMod

  USE PhysConstants, ONLY: getPhysConst
  USE PlanetConstants, ONLY: getPlanetConst
  PRIVATE

  !---------------------------------------------------------------
  !  List of Public subroutines (accessible from outside module) 
  !---------------------------------------------------------------
  PUBLIC :: MODM

  INCLUDE 'isotope.incl'     !INCLUDE HITRSID DATA-STATEMENTS

CONTAINS

      SUBROUTINE MODM(IPR,ICP,NWN,WN,dvset,NLAY,P,T,CLW, &
                  O,O_BY_MOL, OC, O_CLW, ODXSEC, &
                       NMOL,WKL,WBRODL, &
                 SCLCPL,SCLHW,Y0RES,HFILE,cntnmScaleFac,ixsect)
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
!-------------------------------------------------------------------------------
!
!     PROGRAM:  MODM
!     -------
!
!     AUTHOR: Sid-Ahmed Boukabara 
!     ------
!
!     AFFILIATION: ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!     -----------
!
!     DATE OF CREATION : October 1998
!     ----------------
!
!     AIM: This program is aimed at the calculation of the
!     ---  atmospheric optical depths. The spectral validity depends
!          only on the region covered by the file:"spectral_lines.dat"
!          The components treated here are the water vapor, the
!          oxygen, the ozone, the nitrogen and nitrogen dioxide.
!
!     INPUTS:
!     ------
!     - ICP      : Flag to take(if =1) or not (if=0) the line coupling
!     - NWN      : Number of wavenumbers to be treated
!     - WN       : Vector of NWN wavenumbers [in cm-1], one should note that
!                  this input could be a scalar (associated with NWN=1)
!     - NLAY     : Number of layers to be treated.
!     - P        : Vector of NLAY pressures (in mbar), one should note that
!                  this input could be a scalar (associated with NLAY=1)
!     - T        : Vector of NLAY temperatures [in Kelvin]
!     - CLW      : Vector of NLAY Cloud Liquid Water amounts [in kg/m2 or mm]
!                  When Cloud is present, the frequency must be consistent
!                  with Rayleigh absorption (no scattering performed in 
!                  monortm). 
!     - SCLCPL   : Scaling factor of the Line Coupling parameters (usually SCLCPL=1)
!     - SCLHW    : Scaling factor of the pressure dependence of the halfwidth
!                  of the zero frequency band (usually SCLHW=1)
!     - Y0RES    : Y0 resonnance (usually Y0RES=0) to be added to the Yi
!     - HFILE    : Name of the spectral lines information file (HITRAN)
!     - cntnmScaleFac   : Structure containing scale factors to be applied to the continua
!
!
!     OUTPUTS:
!     -------
!     - O      : An array of NWNxNLAY elts containing the total optical depths
!                 due to all the active species [in nepers]
!
!     History of the modifications:
!     *****************************  
!     - written in 1999 by Sid Ahmed Boukabara, Ross Hoffman
!	and Tony Clough. 
!     - validated against ARM sondes in the
!	microwave spectrum (23.8 and 31.4 GHz). SAB, 2000.
!     - extended to more species by Sid Ahmed Boukabara in 03/2000.
!     - cleaned up and commented in 2001 for first public release.
!	Also put under CVS configuration management. SAB.
!     - Extended O2 lines to submillimeter. Extensive validation
!	by comparison to Rosenkranz model and MWR data.
!	Update of the LBLATM module (accepts inputs at pressure 
!	grid, along with altitude grid). 
!	Fixed the handling of N2 amount coming from LBLATM (which
!	depends on the number of molecules NMOL). 
!	Adopted accurate constants values. 
!	Sid Ahmed Boukabara. Dec 14th 2001.
!     - Updated on January 7th 2002. ARM option (INP=2) updated and
!       made more efficient after Jim's comments. (INP=3) option optimized.
!       WV line intensities modified in the microwave (see Tony's email).
!     - Updated on October 2nd 2002. SAB. Speed option implemented.
!       We added also the possibility to run Voigt or Lorentz
!       line shape (speed up process) depending on the current
!       condition (parameter zeta). The pressure induced
!       shifted frequency is also passed to the line shape 
!       computation (instead of the spectroscopic wavenumber).
!     - September 2003: Modified spectral lines file to improve agreement 
!       with SGP MWRP data (provided by Nico Cimini). Scaled O2 line coupling
!       parameters: Y * 0.87, G* 0.
!     - 2006: Implementation of Tretyakov et al (2005) O2 line coupling.
!       Validated against ARM MWRP data (see Cadeddu et al, 2007)
!     - 2007: Updated spectral line file to change the widths and 
!       temperature dependences of the widths for the 22 and 183 GHz lines
!     - 2008: Extensive update to enable the use of MonoRTM beyond 
!       the microwave region and to use the MT_CKD continuum.
!       Updates to self and foreign broadened microwave continuum based
!       on ARM SGP MWR data at 31.4 GHz and ARM FKB (COPS) data at 150 GHz.  
!
!     Comments should be forwarded to Karen Cady-Pereira (cadyp@aer.com)
!     or Vivienne Payne (vpayne@aer.com).
!
!-------------------------------------------------------------------------------
      USE CntnmFactors, ONLY: CntnmFactors_t,oneMolecCntnm,pushCntnmFactors
      USE lnfl_mod, ONLY : GET_LNFL
      USE RTMmono,  ONLY : NWNMX
      USE lblparams, ONLY: MXLAY,MXMOL
      !include "declar.incl"

      parameter (n_absrb=5050,ncont=6)
      INTEGER, INTENT(IN) :: IPR
      real *8  v1abs,v2abs
      real*8 v1, v2
      real   O(NWNMX,MXLAY),OC(NWNMX,MXMOL,MXLAY), &
             O_BY_MOL(NWNMX,MXMOL,MXLAY),O_CLW(NWNMX,MXLAY), &
      	     odxsec(nwnmx,mxlay),CLW(MXLAY),P(MXLAY),T(MXLAY)
      REAL*8 WN(NWNMX)
      REAL WKL(MXMOL,MXLAY),WBRODL(MXLAY)
      real oc_rayl(nwnmx,mxlay)
      real scor(42,9)
      integer index_cont(ncont), imol
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)                
                                                                         
      TYPE(CntnmFactors_t) :: cntnmScaleFac
      CHARACTER*8      XID,       HMOLID,      YID     
      REAL*8               SECANT,       XALTZ
                                                                         
      COMMON /CVRCNT/ HNAMCNT,HVRCNT
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),        &
                      WKC(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   &
                      EMISIV,FSCDID(17),NMOLC,LAYER ,YI1,YID(10),LSTWDF

      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN


      REAL RADCN2
      CHARACTER HFILE*80,HVRMODM*15
      COMMON /CVRMODM/ HVRMODM
      LOGICAL INIT
      DATA INIT/.TRUE./
      SAVE INIT

! h2o, co2, o3, o2,n2,rayleigh/
      data index_cont/1,2,3,7,22,99/
      HVRMODM = '$Revision: 11207 $' 


! Set up useful constants
      CALL getPhysConst(RADCN2=RADCN2)
      ONEPL = 1.001                                                       
      ONEMI = 0.999                                                      
      ARGMIN = 34.                                                      
      EXPMIN = EXP(-ARGMIN) 

! Set up constants for call to contnm
         jrad = 0
         nmolc = nmol
         v1 = wn(1)
         v2 = wn(nwn)
         dvabs = 1.0
         v1abs = int(v1)-3.*DVABS
         V2ABS = INT(V2+3.*DVABS+0.5)
         NPTABS = (V2ABS-V1ABS)/DVABS+1.5

      IF(INIT)THEN
         CALL GET_LNFL(IPR,ICP,HFILE,v1,v2) !reads the HITRAN data
         INIT=.FALSE.
      ENDIF
!  Initialize
        oc(1:nwn,1:mxmol,1:nlay) = 0.
        odxsec(1:nwn,1:nlay) = 0.
        oc_rayl(1:nwn,1:nlay) = 0.
        o(1:nwn,1:nlay) = 0.

         if (ixsect .eq. 1)  &
             call monortm_xsec_sub(wn,nwn,p,t,nlay,odxsec)

      DO K=1,NLAY            !Loop over the Temp/Press/Amount profile 
! Set up variables for call to contnm
         pave = p(k)
         tave = t(k)
         wbroad = wbrodl(k)
         xkt = tave/radcn2

! Call CONTNM for each molecule and for Rayleigh
        wkc(1:nmol) = wkl(1:nmol,k)
        if (nmol.lt.22) wkc(index_cont(ncont-1)) = wbroad    ! set n2 to wbroad if nmol < 22
        do icount=1,ncont
              im = index_cont(icount) 
	      absrb(:) = 0.
	      call oneMolecCntnm(im,cntnmScaleFac)
	      call contnm(jrad)
              !print *, ' '
	      call pushCntnmFactors(cntnmScaleFac)   ! Restore factors
              if (icount.LT.ncont) then 
! Interpolate for gridded spectral resolution in one step
		 if (dvset.ne.0) call xint(v1abs,v2abs,dvabs,absrb,1.0,v1, &
			      dvset,oc(1:nwn,im,k),1,nwn) 
   ! Interpolate for specific wavenumbers one at a time
		 if (dvset.eq.0) then 
		    do iw=1,nwn
		       call xint(v1abs,v2abs,dvabs,absrb,1.0,wn(iw),1.0, &
			    oc(iw,im,k),1,1)
		    end do
		 end if
! Multiply by radiation term
		 do iw=1,nwn
		    oc(iw,im,k) = oc(iw,im,k)*radfn(wn(iw),xkt)
		 end do
              else   !Rayleigh !!
! Interpolate for gridded spectral resolution in one step
		 if (dvset.ne.0) call xint(v1abs,v2abs,dvabs,absrb,1.0,v1, &
			      dvset,oc_rayl(1:nwn,k),1,nwn) 
   ! Interpolate for specific wavenumbers one at a time
		 if (dvset.eq.0) then 
		    do iw=1,nwn
		       call xint(v1abs,v2abs,dvabs,absrb,1.0,wn(iw),1.0, &
			    oc_rayl(iw,k),1,1)
		    end do
		 end if
! Multiply by radiation term
		 do iw=1,nwn
		    oc_rayl(iw,k) = oc_rayl(iw,k)*wn(iw)/1.0e4
		 end do
              endif
         end do                    ! end molecule loop

! calculate TIPS using Gamache routine rather that QOFT
         call tips_2003(nmol,t(k),scor)
 

	 DO M=1,NWN                !loop over the wavenumbers

            CALL INITI(P(K),T(K),RADCT,T0,P0,NMOL, &
              WKL(1:NMOL,K),WBRODL(K),XN0,Xn,Xn_WV)         !INITIALIZATION
            RFT=WN(M)*TANH((RADCT*WN(M))/(2*T(K)))	              !RAD_FIELD_TERM

            CALL LINES(Xn,WN(M),T(K),NMOL,WKL(1:nmol,k)      , &       !PROCESS_LINES
            wbrodl(K),RADCT,T0,o_by_mol(m,1:nmol,k),XN0,RFT, &
            P(K),P0,SCLCPL,SCLHW,Y0RES,scor)
   
            O_CLW(M,K)=ODCLW(WN(M),T(K),CLW(K))                       !OPTDEPTH CLW
            do imol = 1,nmol
                O(M,K) = o(m,k) + O_BY_MOL(M,imol,K) 
                !print *, imol,k,o_by_mol(m,imol,k)
                !print *, oc(m,1:index_cont(5),k)
            enddo
            o(m,k) = o(m,k) + odxsec(m,k) +  oc_rayl(m,k) +  &
                    sum(oc(m,1:index_cont(5),k))+O_CLW(M,K)

         ENDDO

      ENDDO                     ! end layer loop
      RETURN
      END SUBROUTINE MODM


      SUBROUTINE LINES(Xn,WN,T,NMOL,WK, &
           wbrod,RADCT, T0,o_by_mol, &
           XN0,RFT,P,P0,SCLCPL,SCLHW,Y0RES,scor)


      USE lnfl_mod, ONLY : NBLM,ISO,XNU0,DELTNU,E,ALPS,ALPF,X,XG,S0,Rmol, &
             brd_mol_flg,brd_mol_tmp,brd_mol_hw,brd_mol_shft,sdep

      PARAMETER (NNM=  39,IIM= 75000,MXBRDMOL=7)
      REAL WK(NMOL),o_by_mol(nmol)
      REAL*8 WN,XNU
      REAL A(4),B(4),TEMPLC(4)
      real scor(42,9)

      real, dimension(mxbrdmol) ::  rhoslf,tmpcor_arr,alfa_tmp

      DATA TEMPLC /200.0,250.0,296.0,340.0 /

      deltnuC=25.          !cm-1
      WTOT=sum(wk)+wbrod
      RP=P/P0                   !ratio of pressure
      RP2=RP*RP                 !square of the ratio of pressure
      DO IL=1,3                 !find correct temp. interval for interpolation
         ILC=IL
         IF (T.LT.TEMPLC(ILC+1)) GOTO 20
      ENDDO
 20   RECTLC=1.0/(TEMPLC(ILC+1)-TEMPLC(ILC))
      TMPDIF=T-TEMPLC(ILC)                                             
      RT=T/T0                   !ratio of temperature
      RHORAT=(Xn/XN0)               !ratio of number density
      rhoslf(:) = rhorat*wk(1:mxbrdmol)/wtot
      o_by_mol(:)  = 0.         !initialization

      DO I=1,NMOL
         W_SPECIES = WK(i)
         IF (W_SPECIES.EQ.0.) THEN
            OL = 0.
            GOTO 10
         ENDIF
         SF=0.
         J=0
         RAT=(Xn/XN0)*(W_SPECIES/WTOT)
         DO WHILE (J.LT.NBLM(I))
            J=J+1
            JJ=J
            IF ((XG(I,J).EQ.-1).OR.(XG(I,J).EQ.-3).OR.(XG(I,J).EQ.-5))THEN
               JJ=J+1 !the LCC are stored in XNU0(J+1),DELTNU(J+1),etc..
               A(1)=XNU0(I,JJ)
               B(1)=S0(I,JJ)
               A(2)=alpf(I,JJ)
               B(2)=E(I,JJ)
               A(3)=RMOL(I,JJ)   
               B(3)=ALPS(I,JJ)
               A(4)=X(I,JJ)
               B(4)=deltnu(I,JJ)
               AIP=A(ILC)+((A(ILC+1)-A(ILC))*RECTLC)*TMPDIF
               BIP=B(ILC)+((B(ILC+1)-B(ILC))*RECTLC)*TMPDIF
            ENDIF
            IF ((XG(I,J).EQ.-5).AND.(XG(I,J-1).EQ.-5)) THEN   !Self LC!!
               rho_for = (rhorat-rhoslf(i))/rhorat
               rho_sel = rhoslf(i)/rhorat
               A(1) = rho_for*A(1)+rho_sel*XNU0(I,JJ)
               B(1) = rho_for*B(1)+rho_sel*S0(I,JJ)
               A(2) = rho_for*A(2)+rho_sel*alpf(I,JJ)
               B(2) = rho_for*B(2)+rho_sel*E(I,JJ)
               A(3) = rho_for*A(3)+rho_sel*RMOL(I,JJ)
               B(3) = rho_for*B(3)+rho_sel*ALPS(I,JJ)
               A(4) = rho_for*A(4)+rho_sel*X(I,JJ)
               B(4) = rho_for*B(4)+rho_sel*deltnu(I,JJ)
            ENDIF

          
            !---application of the scaling factors
            IF ((XG(I,J).EQ.-1)) THEN !Scaling of the Line coupling parameters
               AIP=AIP*SCLCPL+Y0RES
               BIP=BIP*SCLCPL+Y0RES
            ENDIF
            IF ((XG(I,J).EQ.-3)) THEN !Press depend of the hwidth of the 0 band 
               AIP=AIP*SCLHW
               BIP=BIP*SCLHW
            ENDIF
            ! convert S0 back to HITRAN form
	    S0_adj = S0(I,J)*(xnu0(i,j)*(1.0-exp(-(RADCT*Xnu0(i,j)/T0))))

            ! get shift
             Xnu=Xnu0(I,J)+(deltnu(I,J)*(Xn/xn0))
             !print *,'xnu0,xnu ',xnu0(i,j),xnu
            

            ! modify shift due to specific broadening by other molecules if information is available
            xnu = xnu+sum(rhoslf(:)*brd_mol_flg(i,:,j)*(brd_mol_shft(i,:,j)-deltnu(i,j)))
               
            !check line within 25cm-1 from WN, (except for O2, cause line coupling)
            IF ((ABS(WN-Xnu).GT.deltnuC).and.(I.NE.7))  &
                 GOTO 30   

            XIPSF = scor(i,iso(i,j))

            !CALL INTENS(T,S0(I,J),E(I,J),RADCT,T0,Xnu,STILD,XIPSF)
            CALL INTENS(T,S0_adj,E(I,J),RADCT,T0,Xnu,STILD,XIPSF)

            ! since parameters are now coming from binary file tdep (here called X) comes in correctly
            XTILD=X(I,J)

            ! calculate  Lorentz halfwidth
            ! if specific broadening by other molecules available, recalculate halfwidth
            HWHM_C=HALFWHM_C(alpf(I,J),alps(I,J),RT,XTILD,RHORAT,I, rhoslf,  &
                 brd_mol_flg(i,:,j),brd_mol_hw(i,:,j),brd_mol_tmp(i,:,j))

            !stop
            ! calculate Doppler width
            HWHM_D=HALFWHM_D(I,ISO(I,J),Xnu,T)
         
            IF(XG(I,J).EQ.-3.) THEN
               HWHM_C=HWHM_C*(1-(AIP*(RP))-(BIP*(RP2)))
            ENDIF
            zeta=HWHM_C/(HWHM_C+HWHM_D)
            ilshp=1                            !=0->Lorentz, =1->Voigt

            !MJA 20130517 Assuming that even for speed dependence 
            !we should go to lorentz at high zeta and in wings
            !Seems consistent with Figure 1 of Boone et al., JQSRT, 105, 525-532, 2007.
            if ((ABS(WN-Xnu).GT.(100.*HWHM_D)).or.(zeta.GT.0.99)) ilshp=0
            IF (ilshp.eq.0) CALL LSF_LORTZ(XG(I,J),RP,RP2,AIP,BIP, &
                 HWHM_C,WN,Xnu,SLS,I)
            !MJA 20130517 New speed dependent voigt line shape
            !IF (ilshp.eq.1) CALL LSF_VOIGT(XG(I,J),RP,RP2,AIP,BIP, &
            !     HWHM_C,WN,Xnu,SLS,HWHM_D,I)
!            IF (sdep(I,J) .ne. 0.0) THEN
!               print *, "WVN",  Xnu0(I,J)
!               print *, "SDEP", sdep(I,J)
!               print *, "RP", RP
!            ENDIF
            IF (ilshp.eq.1) CALL LSF_SDVOIGT(XG(I,J),RP,RP2,AIP,BIP, &
                 HWHM_C,WN,Xnu,SLS,HWHM_D,I, sdep(I,J))

            SF=SF+(STILD*SLS)
 30         CONTINUE
            J=JJ
         END DO
         SPSD=W_species*SF
         OL =RFT*SPSD
 10      o_by_mol(i)  = OL
      ENDDO
      END SUBROUTINE LINES

      FUNCTION HALFWHM_D(MOL,ISO,XNU,T)
      PARAMETER (NMOL=39,Nspeci=85)
      REAL BOLTZ,CLIGHT,AVOGAD
      REAL C,K,T,M
      REAL*8 XNU  
      INTEGER ILOC,ISO
      COMMON /ISVECT/ ISO_MAX(NMOL),SMASS(nmol,9)
      common /iso_id/ iso_82(98)
      
      call getPhysConst(BOLTZ=BOLTZ,CLIGHT=CLIGHT,AVOGAD=AVOGAD)
      M=SMASS(mol,iso)
      HALFWHM_D=(XNU/CLIGHT)*SQRT(2.*LOG(2.)*((BOLTZ*T)/(M/AVOGAD)))
      END FUNCTION HALFWHM_D


      SUBROUTINE LSF_VOIGT(XF,RP,RP2,AIP,BIP,HWHM,WN,Xnu,SLS, &
           AD,MOL)
      REAL*8 WN,Xnu,deltXNU,deltnuC,CHI
      DATA MOL_WV/1/,MOL_CO2/2/,MOL_O3/3/,MOL_O2/7/,MOL_N2/22/,&
           MOL_N2O/4/

      deltnuC=25.          !cm-1
      DIFF=(WN+Xnu)-deltnuC
      SLS = 0.
      chi = 1.

! JULY 2008 VHP: 
!                Note that there is currently no pedestal subtraction here for O2.
!                This choice was made in order to avoid discontinuities due to O2 line coupling.
!                We could get around this by generating an O2 continuum in the same way
!                that we generate the CO2 continuum.
      
      IF ((MOL.NE.MOL_O2).AND.(MOL.NE.MOL_CO2)) THEN ! no possibility of line coupling
                                                     ! check for line within 25cm-1 has already
                                                     ! been performed in modm.f

          deltXNU=(WN-Xnu)
          XL1=VOIGT(deltXNU,HWHM,AD) !VOIGT for (+) osc.
          XL3=VOIGT(deltnuC,HWHM,AD) !VOIGT for 25cm-1 wn    
          IF (DIFF .LE. 0.) THEN
             deltXNU=(WN+Xnu)
              XL2=VOIGT(deltXNU,HWHM,AD) !VOIGT for (-) osc.
              SLS = (XL1 + XL2 - (2 * XL3)) 
          ELSE
              SLS = (XL1 - XL3) 
          ENDIF
      ELSE ! O2 or CO2 (check for line within 25 cm-1 has to be performed here for O2)
          IF ((ABS(WN-Xnu).LE.deltnuC).and.(XF.NE.-1).and. &
                 (XF.NE.-3).AND.(XF.NE.-5)) THEN    !no line coupling 
              deltXNU=(WN-Xnu)
              XL1=VOIGT(deltXNU,HWHM,AD) !VOIGT for (+) osc.
              IF (MOL.EQ.MOL_O2) THEN ! O2, no line coupling
                  IF (DIFF .LE. 0.) THEN
                      deltXNU=(WN+Xnu)
                      XL2=VOIGT(deltXNU,HWHM,AD) !VOIGT for (-) osc.
                      SLS = (XL1+XL2) ! no pedestal subtraction for O2
                  ELSE
                      SLS = (XL1)
                  ENDIF
              ELSE ! CO2, no line coupling (no CO2 lines within 25cm-1 of zero: don't need (-) osc)
                  deltXNU = (WN-Xnu)
                  CALL CHI_FN(deltXNU,CHI)
                  XL3 = VOIGT(deltnuC,HWHM,AD) !VOIGT for 25cm-1 wn
!           multiply the co2 pedestal contribution by the chi factor
                  XL3 = XL3*(2.-(deltXNU**2/(deltXNU**2+HWHM**2))) ! co2 pedestal
                  SLS = CHI*(XL1-XL3)
              ENDIF
          ELSE ! line has line coupling
              IF (MOL.EQ.MOL_O2) THEN  
                  IF ((XF.EQ.-1).or.(XF.EQ.-3)) THEN !O2 line coupling: don't implement 25cm-1 check
                      deltXNU=(WN-Xnu)
                      XL1=VOIGT(deltXNU,HWHM,AD) !VOIGT for (+) osc.
                      deltXNU=(WN+Xnu)
                      XL2=VOIGT(deltXNU,HWHM,AD) !VOIGT for (-) osc.
                      IF (XF.EQ.-1) THEN
                          Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
                          Y2=(1.-(AIP*(1/HWHM)*RP*(WN+Xnu))+(BIP*RP2))
                          SLS = (XL1*(Y1)+XL2*(Y2))
                      ELSE
                          Y1=1.
                          Y2=1.
                          SLS=(XL1+XL2)
                      ENDIF
                  ENDIF
              ELSE ! co2
                  IF ((XF.EQ.-1).or.(XF.EQ.-3).OR.(XF.EQ.-5)) THEN ! CO2 line coupling
                                ! For CO2 (unlike O2) contributions beyond 25 cm-1 are in the cntnm
                                ! The "within 25cm-1" check for CO2 has already been performed in modm.f
                                ! calculate pedestal contribution without line coupling (impact)
                                ! multiply this by the chi factor
                                ! calculate the pedestal contribution from line coupling
                                ! multiply this by the chi factor
                                ! add the pedestal contributions from impact and line coupling
                      deltXNU = (WN-Xnu)

                      CALL CHI_FN(deltXNU,CHI)

                      XL1=VOIGT(deltXNU,HWHM,AD) !VOIGT for (+) osc.
                      deltXNU=(WN+Xnu)
!                     no negative oscillation, since no CO2 lines within 25cm-1 of zero cm-1
                      XL3 = VOIGT(deltnuC,HWHM,AD) !VOIGT for 25cm-1 wn
                      IF ((XF.EQ.-1).OR.(XF.EQ.-5)) THEN
                          Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
                          XP4=XL3* &
                              (1./((deltnuC)**2+HWHM**2)) * &
                              (2.-((WN-Xnu)**2/((deltnuC)**2+HWHM**2))) ! co2 impact pedestal
                          YP1=(Y1-1.)* &
                              (2.-(WN-Xnu)**2/((deltnuC)**2+HWHM**2)) ! line coupling contributions to pedestal
                          SLS=CHI*(XL1* &
                              (Y1)-XP4-XL3*(YP1))
                      ELSE
                          deltXNU = (WN-Xnu)
                          CALL CHI_FN(deltXNU,CHI)
                          Y1 = 1.
                          Y2 = 1.
                          XP4=XL3* &
                              (1./((deltnuC)**2+HWHM**2)) * &
                              (2.-(WN-Xnu)**2/((deltnuC)**2+HWHM**2)) ! co2 impact pedestal
                          SLS = CHI*(XL1-XP4)
                      ENDIF
                  ENDIF ! end test for CO2 line coupling
              ENDIF ! end test to distinguish between CO2 and O2
          ENDIF ! end test for line coupled lines within CO2 OR O2
      ENDIF ! end test for any possibility of line coupling

      END SUBROUTINE LSF_VOIGT

      SUBROUTINE LSF_SDVOIGT(XF,RP,RP2,AIP,BIP,HWHM,WN,Xnu,SLS, &
           AD,MOL, SDEP)
      REAL*8 WN,Xnu,deltXNU,deltnuC,CHI
      DATA MOL_WV/1/,MOL_CO2/2/,MOL_O3/3/,MOL_O2/7/,MOL_N2/22/, &
           MOL_N2O/4/
      deltnuC=25.          !cm-1
      DIFF=(WN+Xnu)-deltnuC
      SLS = 0.
      chi = 1.

! JULY 2008 VHP: 
!                Note that there is currently no pedestal subtraction here for O2.
!                This choice was made in order to avoid discontinuities due to O2 line coupling.
!                We could get around this by generating an O2 continuum in the same way
!                that we generate the CO2 continuum.
      !print *, SDEP
      IF ((MOL.NE.MOL_O2).AND.(MOL.NE.MOL_CO2)) THEN ! check for line within 25cm-1 has already
                                                     ! been performed in modm.f
          IF ((XF.EQ.-1).or.(XF.EQ.-3).or.(XF.EQ.-5)) THEN !line coupling
               deltXNU=(WN-Xnu)
               XL1=SDVOIGT(deltXNU,HWHM,AD, SDEP) !VOIGT for (+) osc.
               XL3=SDVOIGT(deltnuC,HWHM,AD, SDEP) !VOIGT for 25cm-1 wn (pedestal) 
               Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))! line coupling for (+) osc
               Y1P=(1.+(AIP*(1/HWHM)*RP*(deltnuC-Xnu))+(BIP*RP2))! line coupling for (+) osc pedestal            
               IF (DIFF .LE. 0.) THEN !Within 25 cm-1 of 0 cm-1
                    deltXNU=(WN+Xnu)
                    XL2=SDVOIGT(deltXNU,HWHM,AD, SDEP) !VOIGT for (-) osc.                    
                    Y2=(1.-(AIP*(1/HWHM)*RP*(WN+Xnu))+(BIP*RP2)) ! line coupling for (-) osc
                    Y2P=(1.-(AIP*(1/HWHM)*RP*(WN+Xnu))+(BIP*RP2)) ! line coupling contributions to (-)pedestal                  
                    SLS = (Y1*(XL1)-Y1P*(XL3)+Y2*(XL2)-Y2P*(XL3))
                    !SLS = (XL1 + XL2 - (2 * XL3)) 
               ELSE
                    SLS = Y1*(XL1) - Y1P*(XL3)
                    !SLS = (XL1 - XL3) 
               ENDIF
          ELSE !No line coupling
               deltXNU=(WN-Xnu)
               XL1=SDVOIGT(deltXNU,HWHM,AD, SDEP) !VOIGT for (+) osc.
               XL3=SDVOIGT(deltnuC,HWHM,AD, SDEP) !VOIGT for 25cm-1 wn    
               IF (DIFF .LE. 0.) THEN !Within 25 cm-1 of 0 cm-1
                    deltXNU=(WN+Xnu)
                    XL2=SDVOIGT(deltXNU,HWHM,AD, SDEP) !VOIGT for (-) osc.
                    SLS = (XL1 + XL2 - (2 * XL3)) 
               ELSE
                    SLS = (XL1 - XL3) 
               ENDIF
          ENDIF


      ELSE ! O2 or CO2 (check for line within 25 cm-1 has to be performed here for O2)
          IF ((ABS(WN-Xnu).LE.deltnuC).and.(XF.NE.-1).and. &
                 (XF.NE.-3).and.(XF.NE.-5)) THEN    !no line coupling 
              deltXNU=(WN-Xnu)
!              print *, "WN", WN, "SDEP", SDEP
              XL1=SDVOIGT(deltXNU,HWHM,AD, SDEP) !VOIGT for (+) osc.
              IF (MOL.EQ.MOL_O2) THEN ! O2, no line coupling
                  IF (DIFF .LE. 0.) THEN
                      deltXNU=(WN+Xnu)
                      XL2=SDVOIGT(deltXNU,HWHM,AD, SDEP) !VOIGT for (-) osc.
                      SLS = (XL1+XL2) ! no pedestal subtraction for O2
                  ELSE
                      SLS = (XL1)
                  ENDIF
              ELSE ! CO2, no line coupling (no CO2 lines within 25cm-1 of zero: don't need (-) osc)
                  deltXNU = (WN-Xnu)
                  CALL CHI_FN(deltXNU,CHI)
                  XL3 = SDVOIGT(deltnuC,HWHM,AD, SDEP) !VOIGT for 25cm-1 wn
!           multiply the co2 pedestal contribution by the chi factor
                  XL3 = XL3*(2.-(deltXNU**2/(deltXNU**2+HWHM**2))) ! co2 pedestal
                  SLS = CHI*(XL1-XL3)
              ENDIF

          ELSE ! line has line coupling
              IF (MOL.EQ.MOL_O2) THEN  
                  IF ((XF.EQ.-1).or.(XF.EQ.-3).or.(XF.EQ.-5)) THEN !O2 line coupling: don't implement 25cm-1 check
                      deltXNU=(WN-Xnu)
                      XL1=SDVOIGT(deltXNU,HWHM,AD, SDEP) !VOIGT for (+) osc.
                      deltXNU=(WN+Xnu)
                      XL2=SDVOIGT(deltXNU,HWHM,AD, SDEP) !VOIGT for (-) osc.
                      IF (XF.EQ.-1) THEN
                          Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
                          Y2=(1.-(AIP*(1/HWHM)*RP*(WN+Xnu))+(BIP*RP2))
                          SLS = (XL1*(Y1)+XL2*(Y2))
                      ELSE
                          Y1=1.
                          Y2=1.
                          SLS=(XL1+XL2)
                      ENDIF
                  ENDIF
              ELSE ! co2
                  IF ((XF.EQ.-1).or.(XF.EQ.-3).or.(XF.NE.-5)) THEN ! CO2 line coupling
                                ! For CO2 (unlike O2) contributions beyond 25 cm-1 are in the cntnm
                                ! The "within 25cm-1" check for CO2 has already been performed in modm.f
                                ! calculate pedestal contribution without line coupling (impact)
                                ! multiply this by the chi factor
                                ! calculate the pedestal contribution from line coupling
                                ! multiply this by the chi factor
                                ! add the pedestal contributions from impact and line coupling
                      deltXNU = (WN-Xnu)

                      CALL CHI_FN(deltXNU,CHI)
!                      print *, "CO2 pos osc"
                      XL1=SDVOIGT(deltXNU,HWHM,AD, SDEP) !VOIGT for (+) osc.
                      deltXNU=(WN+Xnu)
!                     no negative oscillation, since no CO2 lines within 25cm-1 of zero cm-1
!                      print *, "CO2 pedestal"
                      XL3 = SDVOIGT(deltnuC,HWHM,AD, SDEP) !VOIGT for 25cm-1 wn

!                      IF (XL1 .LT. 0.0 .OR. XL3 .LT. 0.0  &
!                          .OR. (XL1-XL3) .LT. 0.0) THEN
!                         print *, "XL1-XL3", XL1, XL3, XL1-XL3
!                      ENDIF

                      IF (XF.EQ.-1.or.XF.EQ.-5) THEN
                          Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
                          XP4=XL3* &
                              (1./((deltnuC)**2+HWHM**2)) * &
                              (2.-((WN-Xnu)**2/((deltnuC)**2+HWHM**2))) ! co2 impact pedestal
                          YP1=(Y1-1.)* &
                              (2.-(WN-Xnu)**2/((deltnuC)**2+HWHM**2)) ! line coupling contributions to pedestal
                          SLS=CHI*(XL1* &
                              (Y1)-XP4-XL3*(YP1))
                      ELSE
                          deltXNU = (WN-Xnu)
                          CALL CHI_FN(deltXNU,CHI)
                          Y1 = 1.
                          Y2 = 1.
                          XP4=XL3* &
                              (1./((deltnuC)**2+HWHM**2)) * &
                              (2.-(WN-Xnu)**2/((deltnuC)**2+HWHM**2)) ! co2 impact pedestal
                          SLS = CHI*(XL1-XP4)
                      ENDIF
                  ENDIF ! end test for CO2 line coupling
              ENDIF ! end test to distinguish between CO2 and O2
          ENDIF ! end test for line coupled lines within CO2 OR O2
      ENDIF ! end test for any possibility of line coupling

      END SUBROUTINE LSF_SDVOIGT

      
      SUBROUTINE LSF_LORTZ(XF,RP,RP2,AIP,BIP,HWHM,WN,Xnu,SLS, &
          MOL)
      REAL*8 WN,Xnu,deltnuC,deltXNU,CHI
      DATA MOL_WV/1/,MOL_CO2/2/,MOL_O3/3/,MOL_O2/7/,MOL_N2/22/, &
           MOL_N2O/4/

      deltnuC=25.          !cm-1
      DIFF=(WN+Xnu)-deltnuC
       SLS = 0.
       chi= 1.

! JULY 2008 VHP: 
!                Note that there is currently no pedestal subtraction here for O2.
!                This choice was made in order to avoid discontinuities due to O2 line coupling.
!                In the future we could generate an O2 continuum in the same way
!                that we generate the CO2 continuum.

      IF ((MOL.NE.MOL_O2).AND.(MOL.NE.MOL_CO2)) THEN ! no possibility of line coupling
                                                     ! check for line within 25cm-1 has already
                                                     ! been performed in modm.f

          IF ((XF.EQ.-1).or.(XF.EQ.-3).or.(XF.EQ.-5)) THEN !line coupling
               deltXNU=(WN-Xnu)
               XL1=XLORENTZ((deltXNU)/HWHM) !LORENTZ for (+) osc.
               XL3=XLORENTZ((deltnuC)/HWHM) !LORENTZ for 25cm-1 wn               
               Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))! line coupling for (+) osc
               Y1P=(1.+(AIP*(1/HWHM)*RP*(deltnuC-Xnu))+(BIP*RP2))! line coupling for (+) osc pedestal            
               IF (DIFF .LE. 0.) THEN !Within 25 cm-1 of 0 cm-1
                    deltXNU=(WN+Xnu)
                    XL2=XLORENTZ((deltXNU)/HWHM) !LORENTZ for (-) osc.
                    Y2=(1.-(AIP*(1/HWHM)*RP*(WN+Xnu))+(BIP*RP2)) ! line coupling for (-) osc
                    Y2P=(1.-(AIP*(1/HWHM)*RP*(WN+Xnu))+(BIP*RP2)) ! line coupling contributions to (-)pedestal                  
                    SLS = (Y1*(XL1)-Y1P*(XL3)+Y2*(XL2)-Y2P*(XL3))
                    !SLS = (XL1 + XL2 - (2 * XL3)) 
               ELSE
                    SLS = Y1*(XL1) - Y1P*(XL3)
                    !SLS = (XL1 - XL3) 
               ENDIF
          ELSE !No line coupling
              deltXNU=(WN-Xnu)
              XL1=XLORENTZ((deltXNU)/HWHM) !LORENTZ for (+) osc.
              XL3=XLORENTZ((deltnuC)/HWHM) !LORENTZ for 25cm-1 wn    
              IF (DIFF .LE. 0.) THEN
                  deltXNU=(WN+Xnu)
                  XL2=XLORENTZ((deltXNU)/HWHM) !LORENTZ for (-) osc.
                  SLS = (XL1 + XL2 - (2 * XL3)) / HWHM
              ELSE
                  SLS = (XL1 - XL3) / HWHM
              ENDIF
          ENDIF
      ELSE                      ! O2 or CO2 (check for line within 25 cm-1 has to be performed for O2)
          IF ((ABS(WN-Xnu).LE.deltnuC).and.(XF.NE.-1).and. &
                 (XF.NE.-3).AND.(XF.NE.-5)) THEN    !no line coupling 
              deltXNU=(WN-Xnu)
              XL1=XLORENTZ((deltXNU)/HWHM) !LORENTZ for (+) osc.
              IF (MOL.EQ.MOL_O2) THEN ! O2, no line coupling
                  IF (DIFF .LE. 0.) THEN
                      deltXNU=(WN+Xnu)
                      XL2=XLORENTZ((deltXNU)/HWHM) !LORENTZ for (-) osc.
                      SLS = (XL1+XL2 ) / HWHM
                  ELSE
                      SLS = (XL1) / HWHM
                  ENDIF
              ELSE ! CO2, no line coupling (no CO2 lines within 25cm-1 of zero: don't need (-) osc)
                  deltXNU = (WN-Xnu)
                  CALL CHI_FN(deltXNU, CHI)
                  XL3=XLORENTZ((deltnuC)/HWHM) !LORENTZ for 25cm-1 wn
!           multiply the co2 pedestal contribution by the chi factor
                  XL3=XL3*(2.-(deltXNU**2/(deltXNU**2+HWHM**2))) ! co2 pedestal
                  SLS = CHI*(XL1-XL3) / HWHM
              ENDIF
          ELSE ! line has line coupling
              IF (MOL.EQ.MOL_O2) THEN ! for O2 line-coupled lines the 25 cm-1 check is not implemented
                  IF ((XF.EQ.-1).or.(XF.EQ.-3).or.(XF.EQ.-5)) THEN !O2 line coupling
                      deltXNU=(WN-Xnu)
                      XL1=XLORENTZ((deltXNU)/HWHM) !LORENTZ for (+) osc.
                      deltXNU=(WN+Xnu)
                      XL2=XLORENTZ((deltXNU)/HWHM) !LORENTZ for (-) osc.
                      IF (XF.EQ.-1) THEN
                          Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
                          Y2=(1.-(AIP*(1/HWHM)*RP*(WN+Xnu))+(BIP*RP2))
                          SLS = (XL1*(Y1)+XL2*(Y2)) / HWHM
                      ELSE
                          Y1=1.
                          Y2=1.
                          SLS=(XL1+XL2) / HWHM
                      ENDIF
                  ENDIF
              ELSE
                  IF ((XF.EQ.-1).or.(XF.EQ.-3).OR.(XF.EQ.-5)) THEN ! CO2 line coupling
                                ! For CO2 (unlike O2) contributions beyond 25 cm-1 are in the cntnm
                                ! The "within 25cm-1" check for CO2 has already been performed in modm.f
                                ! calculate pedestal contribution without line coupling (impact)
                                ! multiply this by the chi factor
                                ! calculate the pedestal contribution from line coupling
                                ! multiply this by the chi factor
                                ! add the pedestal contributions from impact and line coupling
                      deltXNU = (WN-Xnu)
                      CALL CHI_FN(deltXNU, CHI)
                      XL1=XLORENTZ((deltXNU)/HWHM) !LORENTZ for (+) osc.
                      deltXNU=(WN+Xnu)
!                     no negative oscillation, since no co2 lines within 25cm-1 of zero
                      XL3 = XLORENTZ((deltnuC)/HWHM) !LORENTZ for 25cm-1 wn
                      IF ((XF.EQ.-1).OR.(XF.EQ.-5)) THEN
                          Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
                          XP4=XL3* &
                              (1./((deltnuC)**2+HWHM**2)) * &
                              (2.-((WN-Xnu)**2/((deltnuC)**2+HWHM**2))) ! co2 impact pedestal
                          YP1=(Y1-1.)* &
                              (2.-(WN-Xnu)**2/((deltnuC)**2+HWHM**2))
                          SLS=CHI*(XL1* &
                              (Y1)-XP4-XL3*(YP1))/HWHM
                      ELSE
                          Y1 = 1.
                          Y2 = 1.
                          XP4=XL3* &
                              (1./((deltnuC)**2+HWHM**2)) * &
                              (2.-((WN-Xnu)**2/((deltnuC)**2+HWHM**2))) ! co2 impact pedestal
                          SLS = CHI*(XL1 - XP4) / HWHM
                      ENDIF
                  ENDIF ! end test for CO2 line coupling
              ENDIF ! end test to distinguish between CO2 and O2
          ENDIF ! end test for line coupled lines within CO2 OR O2
      ENDIF ! end test for any possibility of line coupling

      END SUBROUTINE LSF_LORTZ


      
      FUNCTION HALFWHM_C(AF,AS,RT,XTILD,RHORAT,MOL,rhoslf,brd_flg,brd_hw,brd_tmp)
      parameter (mxbrdmol=7)
      integer*4, dimension(mxbrdmol) :: brd_flg
      real, dimension(mxbrdmol)      :: brd_hw,brd_tmp,rhoslf

      real, dimension(mxbrdmol)  :: tmpcor,alfa_tmp

      IF ((MOL.EQ.1).AND.(AS.EQ.0.)) AS=5*AF
      alfa0i = AF*(RT**XTILD)
      hwhmsi = AS*(RT**XTILD)
      !HALFWHM_C=AF*(RT**XTILD)*(RHORAT-RAT)+AS*(RT**XTILD)*RAT
      HALFWHM_C= alfa0i*(RHORAT-rhoslf(mol)) + hwhmsi*rhoslf(mol)

      !recalculate halfwith if information on broadening by specific molecules is available
	 if(sum(brd_flg(:)).gt.0) then
	    tmpcor = RT**brd_tmp(:)
	    alfa_tmp = brd_hw(:)*tmpcor
	    alfsum = sum(rhoslf(:)*brd_flg(:)*alfa_tmp)
	    HALFWHM_C = (rhorat-sum(rhoslf(:)*brd_flg(:)))    &
	      *alfa0i + alfsum
  ! if no new self info, need to add standard self to total half width
	    if(brd_flg(mol).eq.0)                                   &
		 HALFWHM_C = HALFWHM_C + rhoslf(mol)*(hwhmsi-alfa0i)
	 end if
      END FUNCTION HALFWHM_C

      
      SUBROUTINE INTENS(T,S0s,Es,RADCT,T0,Xnus,STILD,XIPSF)
      REAL*8 Xnus
      S=S0s*(EXP(-RADCT*Es/T)/EXP(-RADCT*Es/T0))*XIPSF
      STILD=S*((1+EXP(-(RADCT*Xnus/T)))/ &
      (Xnus*(1-EXP(-(RADCT*Xnus/T0)))))
      END SUBROUTINE INTENS

      SUBROUTINE INITI(P,T,RADCT,T0,P0,NMOL,WK,WBROD,XN0, &
                       Xn,Xn_WV)
      REAL PLANCK,BOLTZ,CLIGHT,WVMOLMASS,DRYMOLMASS
      REAL wk(39),wbrod
      CALL getPhysConst(PLANCK=PLANCK,BOLTZ=BOLTZ,CLIGHT=CLIGHT)
      CALL getPlanetConst(WVMWT=WVMOLMASS,AIRMWT=DRYMOLMASS)
      RADCT=PLANCK*CLIGHT/BOLTZ !in K/cm-1
      T0=296.                   !in K
      P0=1013.25                !in HPa
      XN0=(P0/(BOLTZ*T0))*1.E+3
      XN =(P /(BOLTZ*T ))*1.E+3
      WDRY=sum(wk(2:nmol))+wbrod
      RATIOMIX=(Wk(1)*WVMOLMASS)/(WDRY*DRYMOLMASS)
      WVPRESS=(RATIOMIX/(RATIOMIX+(WVMOLMASS/DRYMOLMASS)))*P
      Xn_WV=(WVPRESS/(BOLTZ*T))*1.E+3
      END SUBROUTINE INITI
      
        

      FUNCTION ODCLW(WN,TEMP,CLW)
      !INPUTS: WN (WaveNUmber in cm-1)
      !        Temp (in K)
      !        CLW  (in mm or kg/m2)
      !OUTPUT: ODCLW: optical depth of the Cloud Liquid Water
      !FROM EQUATIONS OF LIEBE, HUFFORD AND MANABE
      !Ref:(INT. J. IR & MM WAVES V.12(17) JULY 1991
      COMPLEX EPS,RE
      REAL*8 WN
      REAL PI,CLIGHT
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
      ODCLW = -(6.*PI/299.792458)*CLW*AIMAG(RE)*FREQ
      RETURN
      END FUNCTION ODCLW


      FUNCTION XLORENTZ(Z)
      REAL*8 Z
      REAL PI
      CALL getPhysConst(PI=PI)
      XLORENTZ=1./(PI*(1.+(Z**2))) 
      RETURN
      END FUNCTION XLORENTZ



      FUNCTION VOIGT(DELTNU,ALPHAL,ALPHAD)
      COMPLEX v
      REAL*8 DELTNU
      REAL avc(0:101)
      REAL ALPHAL,ALPHAD,ZETA,ALPHAV,AVCINTERP,DNU
      REAL AL,AD,X,Y,ANORM1,VOIGT,DZETA
      INTEGER IZETA2,IZETA1
      REAL PI
      DATA AVC/                                                    &
        .10000E+01,.99535E+00,.99073E+00,.98613E+00,.98155E+00,    &
        .97700E+00,.97247E+00,.96797E+00,.96350E+00,.95905E+00,    &
        .95464E+00,.95025E+00,.94589E+00,.94156E+00,.93727E+00,    &
        .93301E+00,.92879E+00,.92460E+00,.92045E+00,.91634E+00,    &
        .91227E+00,.90824E+00,.90425E+00,.90031E+00,.89641E+00,    &
        .89256E+00,.88876E+00,.88501E+00,.88132E+00,.87768E+00,    &
        .87410E+00,.87058E+00,.86712E+00,.86372E+00,.86039E+00,    &
        .85713E+00,.85395E+00,.85083E+00,.84780E+00,.84484E+00,    &
        .84197E+00,.83919E+00,.83650E+00,.83390E+00,.83141E+00,    &
        .82901E+00,.82672E+00,.82454E+00,.82248E+00,.82053E+00,    &
        .81871E+00,.81702E+00,.81547E+00,.81405E+00,.81278E+00,    &
        .81166E+00,.81069E+00,.80989E+00,.80925E+00,.80879E+00,    &
        .80851E+00,.80842E+00,.80852E+00,.80882E+00,.80932E+00,    &
        .81004E+00,.81098E+00,.81214E+00,.81353E+00,.81516E+00,    &
        .81704E+00,.81916E+00,.82154E+00,.82418E+00,.82708E+00,    &
        .83025E+00,.83370E+00,.83742E+00,.84143E+00,.84572E+00,    &
        .85029E+00,.85515E+00,.86030E+00,.86573E+00,.87146E+00,    &
        .87747E+00,.88376E+00,.89035E+00,.89721E+00,.90435E+00,    &
        .91176E+00,.91945E+00,.92741E+00,.93562E+00,.94409E+00,    &
        .95282E+00,.96179E+00,.97100E+00,.98044E+00,.99011E+00,    &
        .10000E+01,.10000E+01/                             
      call getPhysConst(PI=PI)
      !---computes zeta
      zeta=alphal/(alphal+alphad)
      !---interpolation of the AVC
      IZETA1=ZETA*100
      IZETA2=IZETA1+1
      DZETA=((ZETA*100.)-INT(ZETA*100.))
      AVCINTERP=AVC(IZETA1)+DZETA*(AVC(IZETA2)-AVC(IZETA1))
      !---voigt width 
      ALPHAV = AVCINTERP*(ALPHAD + ALPHAL)
      !-----GENERATE VOIGT PROFILE SUCH THAT THE VOIGT HALFWIDTH = 1.
      if (zeta .lt. 1.00) then 
         AL=ALPHAL/ALPHAD
         AD= 1.
         dnu=deltnu/alphad
      end if
      !---case of pure lorentz
      if (zeta .eq. 1.00) then 
         VOIGT=(ALPHAL/(PI*(ALPHAL**2+(DELTNU)**2)))
         RETURN
      end if
      !---SETUP PARAMETERS FOR CALL TO VOIGT GENERATOR (HUMLICEK)
      x = sqrt(log(2.))*(dnu)
      y = 1000.
      if (zeta .lt. 1.000) then
         y = sqrt(log(2.))*AL
      end if
      !---CALL the Humlicek subroutine
      v=W4(x,y)
      anorm1 = sqrt(log(2.)/PI)/alphad
      VOIGT=REAL(v)*ANORM1
      RETURN
      end FUNCTION VOIGT


      FUNCTION SDVOIGT(DELTNU,ALPHAL,ALPHAD,SDEP)
      COMPLEX v, lm_factor
      REAL*8 DELTNU
      REAL avc(0:101)
      REAL ALPHAL,ALPHAD,ZETA,ALPHAV,AVCINTERP,DNU
      REAL AL,AD,ANORM1,SDVOIGT,DZETA
      !MJA 20130517 Speed Dependence
      REAL SDEP, alfa, beta, delta, temp, x1, x2, y1, y2, sign, TINY
      REAL alfadelta
      INTEGER LM_flag
      INTEGER IZETA2,IZETA1
      REAL PI
      DATA AVC/                                                    &
        .10000E+01,.99535E+00,.99073E+00,.98613E+00,.98155E+00,    &
        .97700E+00,.97247E+00,.96797E+00,.96350E+00,.95905E+00,    &
        .95464E+00,.95025E+00,.94589E+00,.94156E+00,.93727E+00,    &
        .93301E+00,.92879E+00,.92460E+00,.92045E+00,.91634E+00,    &
        .91227E+00,.90824E+00,.90425E+00,.90031E+00,.89641E+00,    &
        .89256E+00,.88876E+00,.88501E+00,.88132E+00,.87768E+00,    &
        .87410E+00,.87058E+00,.86712E+00,.86372E+00,.86039E+00,    &
        .85713E+00,.85395E+00,.85083E+00,.84780E+00,.84484E+00,    &
        .84197E+00,.83919E+00,.83650E+00,.83390E+00,.83141E+00,    &
        .82901E+00,.82672E+00,.82454E+00,.82248E+00,.82053E+00,    &
        .81871E+00,.81702E+00,.81547E+00,.81405E+00,.81278E+00,    &
        .81166E+00,.81069E+00,.80989E+00,.80925E+00,.80879E+00,    &
        .80851E+00,.80842E+00,.80852E+00,.80882E+00,.80932E+00,    &
        .81004E+00,.81098E+00,.81214E+00,.81353E+00,.81516E+00,    &
        .81704E+00,.81916E+00,.82154E+00,.82418E+00,.82708E+00,    &
        .83025E+00,.83370E+00,.83742E+00,.84143E+00,.84572E+00,    &
        .85029E+00,.85515E+00,.86030E+00,.86573E+00,.87146E+00,    &
        .87747E+00,.88376E+00,.89035E+00,.89721E+00,.90435E+00,    &
        .91176E+00,.91945E+00,.92741E+00,.93562E+00,.94409E+00,    &
        .95282E+00,.96179E+00,.97100E+00,.98044E+00,.99011E+00,    &
        .10000E+01,.10000E+01/                             
      call getPhysConst(PI=PI)
      TINY = 1.0e-4
      !---computes zeta
      zeta=alphal/(alphal+alphad)
      !---interpolation of the AVC
      IZETA1=ZETA*100
      IZETA2=IZETA1+1
      DZETA=((ZETA*100.)-INT(ZETA*100.))
      AVCINTERP=AVC(IZETA1)+DZETA*(AVC(IZETA2)-AVC(IZETA1))
      !---voigt width 
      ALPHAV = AVCINTERP*(ALPHAD + ALPHAL)
      !-----GENERATE VOIGT PROFILE SUCH THAT THE VOIGT HALFWIDTH = 1.
      if (zeta .lt. 1.00) then 
         AL=ALPHAL/ALPHAD
         AD= 1.
         dnu=deltnu/alphad
      end if
      !---case of pure lorentz
      if (zeta .eq. 1.00 .and. ABS(SDEP) .LT. TINY) then 
         SDVOIGT=(ALPHAL/(PI*(ALPHAL**2+(DELTNU)**2)))
         RETURN
      end if
      !---SETUP PARAMETERS FOR CALL TO VOIGT GENERATOR (HUMLICEK)
      IF (ABS(SDEP) .GT. TINY) THEN !SD Voigt 
      !---Speed Dependence follows Boone et al., 2011, An efficient analytical 
      !---approach for calculating line mixing in atmospheric remote sensing 
      !---applications, JQSRT, 112, 980-989.

          !SDEP is Benner's speed dependence definition
          !Therefore Boone's gamma2 = alphal*SDEP and his alfa = (1/SDEP) - 1.5
          gamma2 = alphal*SDEP
          alfa = (alphal/gamma2) - 1.5 !Boone et al., 2011 Eq 13
          beta = (deltnu/gamma2) !Boone et al., 2011 Eq 13
          delta = (1.0/4.0/log(2.))*(alphad*alphad/gamma2/gamma2) !Boone et al., 2011 Eq 13   
          alfadelta = alfa+delta

          !Boone et al., 2011, Eq. 12
          temp = sqrt(alfadelta*alfadelta + beta*beta) !"temp" means temporary, used below
          !print *, "temp", temp
          !Real part of z1, z2, Boone et al., 2011, Eq. 12
          x1 = (1.0/sqrt(2.0))*sqrt(temp+alfadelta)-sqrt(delta)
          x2 = x1+2.0*sqrt(delta)
!          print *, "x1", x1, "x2", x2, "beta", beta, "deltnu",deltnu 
          !Imag part of z1, z2, Boone et al., 2011, Eq. 12
          if (beta .gt. 0.0) then
               sign = 1 
          else if (beta .eq. 0.0) then 
               sign = 0 
          else 
               sign = -1
          endif
!          print *, "sign", sign
          y1 = sign*sqrt((temp-delta-alfa)/2.0)
!          y1 = (1.0/sqrt(2.0))*sqrt(temp-alfadelta)
          y2 = y1 
!          print *, "y1", y1, "y2", y2


          !---CALL the modified Humlicek subroutine to calc speed dependent Voigt
          v=SD_Humlicek(y1,x1,y2,x2)
          vtemp1 = W4(y1,x1)
          vtemp2 = W4(y2,x2)
!          print *, "SD_Humlicek", REAL(v), REAL(vtemp1-vtemp2), REAL(vtemp1), REAL(vtemp2)
          IF (REAL(v).LT.0.0) STOP

          anorm1 = sqrt(log(2.)/PI)/alphad !Boone et al., 2011 Eq 10 
          !print *, 'anorm', anorm1, 'alphad', alphad, 'PI', PI
           SDVOIGT=REAL(v)*ANORM1

      ELSE !Normal Voigt
          !---SETUP PARAMETERS FOR CALL TO VOIGT GENERATOR (HUMLICEK)
          x = sqrt(log(2.))*(dnu)
          y = 1000.
          if (zeta .lt. 1.000) then
             y = sqrt(log(2.))*AL
          end if
          !---CALL the Humlicek subroutine
          v=W4(x,y)
          !print *, "Humlicek", REAL(v)
          anorm1 = sqrt(log(2.)/PI)/alphad
          !print *, 'anorm', anorm1, 'alphad', alphad, 'PI', PI
          SDVOIGT=REAL(v)*ANORM1
      ENDIF
      anorm1 = sqrt(log(2.)/PI)/alphad
      !print *, 'anorm', anorm1, 'alphad', alphad, 'PI', PI
      SDVOIGT=REAL(v)*ANORM1
      RETURN
    
      end FUNCTION SDVOIGT


!*************************************************************************
!     FORTRAN function for the complex probability function w(z).
!     COMPUTES THE COMPLEX PROBABILITY FUNCTION W(Z)=EXP(-Z**2)*ERFC(-I*Z)
!     IN THE UPPER HALF-PLANE Z=X+I*Y (I.E. FOR Y>=0)
!     MAXIMUM RELATIVE ERROR OF BOTH REAL AND IMAGINARY PARTS IS <1*10**(-4)
!     Reference:
!     Humlicek, J., 1982; Optimized Computation of the Voigt and Complex
!     Probability Functions, J. Quant. Spectrosc. Radiat. Transfer, 27, 437-444.
!     S-A Boukabara AER Inc, 2000
!*************************************************************************
      FUNCTION W4(x,y)
      COMPLEX W4,T,U
      
      T=CMPLX(Y,-X)
      S=ABS(X)+Y
      IF(S.LT.15.)GOTO 1
!     ***   REGION I
      W4=T*.5641896/(.5+T*T)
      RETURN
 1    IF(S.LT.5.5)GOTO 2
!     ***   REGION II
      U=T*T
      W4=T*(1.410474+U*.5641896)/(.75+U*(3.+U))
      RETURN
 2    IF(Y.LT. 0.195*ABS(X)-0.176)GOTO 3
!     ***   REGION III
      W4=(16.4955+T*(20.20933+T*(11.96482+ &
           T*(3.778987+T*.5642236))))/ &
           (16.4955+T*(38.82363+T*(39.27121+ &
           T*(21.69274+T*(6.699398+T)))))
      RETURN
!     ***   REGION IV
 3    U=T*T
      W4=CEXP(U)-T*(36183.31-U*(3321.9905- &
           U*(1540.787-U*(219.0313-U* &
           (35.76683-U*(1.320522-U*.56419))))))/ &
           (32066.6-U*(24322.84-U* &
           (9022.228-U*(2186.181-U*(364.2191- &
           U*(61.57037-U*(1.841439-U)))))))
      RETURN
      END FUNCTION W4

!*************************************************************************
!     FORTRAN function for the difference of two complex probability 
!     functions w(z1)-w(z2)..
!     COMPUTES THE COMPLEX PROBABILITY FUNCTION W(Z)=EXP(-Z**2)*ERFC(-I*Z)
!     IN THE UPPER HALF-PLANE Z=X+I*Y (I.E. FOR Y>=0)
!     MAXIMUM RELATIVE ERROR OF BOTH REAL AND IMAGINARY PARTS IS <1*10**(-4)
!     Modified for speed dependece as recommended in Boone et al., 2011.
!     References:
!     Humlicek, J., 1982; Optimized Computation of the Voigt and Complex
!     Probability Functions, J. Quant. Spectrosc. Radiat. Transfer, 27, 437-444.
!     Boone et al., 2011, An efficient analytical 
!     approach for calculating line mixing in atmospheric remote sensing 
!     applications, JQSRT, 112, 980-989.
!     Original Humlicek Implementation: 
!     S-A Boukabara AER Inc, 2000
!     Modifications for speed dependence:
!     M. J. Alvarado, AER, 2013
!*************************************************************************
      FUNCTION SD_Humlicek(x1,y1,x2,y2)
      COMPLEX SD_Humlicek,T1,T2,U1,U2, W1, W2
      REAL S1, S2, X1, Y1, X2, Y2
      INTEGER Region1, Region2, Region
      
      T1=CMPLX(Y1,-X1)
      T2=CMPLX(Y2,-X2)
      S1=ABS(X1)+Y1
      S2=ABS(X2)+Y2

!     Find correct Humlicek region for each line
      Region1 = 1
      IF(S1.GE.15.0) THEN
          Region1 = 1
      ELSE IF (S1.GE.6.0 .AND. S1.LT.15.0) THEN
          Region1 = 2
      ELSE !S1 .LT. 6.0
          Region1 = 3
          IF(Y1.LT. 0.195*ABS(X1)-0.176) Region1 = 4
      ENDIF

      Region2 = 1
      IF(S2.GE.15.0) THEN
          Region2 = 1
      ELSE IF (S2.GE.6.0 .AND. S2.LT.15.0) THEN
          Region2 = 2
      ELSE !S2 .LT. 6.0
          Region2 = 3
          IF(Y2.LT. 0.195*ABS(X2)-0.176) Region2 = 4
      ENDIF 

      !Use Largest of two regions
      Region = MAX(Region1, Region2)
!      print *, "Region1" , Region1, "Region2", Region2, "Region", Region
      IF(Region .GT. 1)GOTO 1
!     ***   REGION I
!      print *, "Region I"
      W1=T1*.5641896/(.5+T1*T1)
      W2=T2*.5641896/(.5+T2*T2)
      SD_Humlicek = W1-W2
      RETURN
 1    IF(Region .GT. 2)GOTO 2 
!     Change in Region 2 and 3 boundary recommended on page 985 of Boone et al. 2011
!     (see 4th paragraph of second column)
! 1    IF(S.LT.5.5)GOTO 2
!     ***   REGION II
!      print *, "Region II"
      U1=T1*T1
      U2=T2*T2 
      W1=T1*(1.410474+U1*.5641896)/(.75+U1*(3.+U1))
      W2=T2*(1.410474+U2*.5641896)/(.75+U2*(3.+U2))
      SD_Humlicek = W1-W2
      RETURN
 2    IF(Region .GT. 3)GOTO 3
!     ***   REGION III
!      print *, "Region III"
      W1=(16.4955+T1*(20.20933+T1*(11.96482+ &
           T1*(3.778987+T1*.5642236))))/ &
           (16.4955+T1*(38.82363+T1*(39.27121+ &
           T1*(21.69274+T1*(6.699398+T1)))))
      W2=(16.4955+T2*(20.20933+T2*(11.96482+ &
           T2*(3.778987+T2*.5642236))))/ &
           (16.4955+T2*(38.82363+T2*(39.27121+ &
           T2*(21.69274+T2*(6.699398+T2)))))
      SD_Humlicek = W1-W2
      RETURN
!     ***   REGION IV
 3    U1=T1*T1
!      print *, "Region IV"
!     You get errors if you use Region IV approximations for lines outside
!     Region IV, so use Region III instead (MJA, 08062013)
      U2=T2*T2
      IF(Region1 .EQ. 4) THEN
        W1=CEXP(U1)-T1*(36183.31-U1*(3321.9905- &
           U1*(1540.787-U1*(219.0313-U1* &
           (35.76683-U1*(1.320522-U1*.56419))))))/ &
           (32066.6-U1*(24322.84-U1* &
           (9022.228-U1*(2186.181-U1*(364.2191- &
           U1*(61.57037-U1*(1.841439-U1)))))))
      ELSE
        W1=(16.4955+T1*(20.20933+T1*(11.96482+ &
           T1*(3.778987+T1*.5642236))))/ &
           (16.4955+T1*(38.82363+T1*(39.27121+ &
           T1*(21.69274+T1*(6.699398+T1)))))
      ENDIF
      IF(Region2 .EQ. 4) THEN
        W2=CEXP(U2)-T2*(36183.31-U2*(3321.9905- &
           U2*(1540.787-U2*(219.0313-U2* &
           (35.76683-U2*(1.320522-U2*.56419))))))/ &
           (32066.6-U2*(24322.84-U2* &
           (9022.228-U2*(2186.181-U2*(364.2191- &
           U2*(61.57037-U2*(1.841439-U2)))))))
      ELSE
        W2=(16.4955+T2*(20.20933+T2*(11.96482+ &
           T2*(3.778987+T2*.5642236))))/ &
           (16.4955+T2*(38.82363+T2*(39.27121+ &
           T2*(21.69274+T2*(6.699398+T2)))))
      ENDIF
!      print *, W1, W2
      SD_Humlicek = W1-W2
      RETURN
      END FUNCTION SD_Humlicek
      
      subroutine chi_fn (deltXNU,chi)

      real*8 chi, deltXNU

      DATA ASUBL / 0.800 /,BSUBL / 10.0 /
                                                                         
!     SET UP CHI SUB-LORENTZIAN FORM FACTOR FOR CARBON DIOXIDE            
!     POLYNOMIAL MATCHED TO AN EXPONENTIAL AT X0 = 10 CM-1                
                                                                         

! Commented lines here are lines contained, but NOT currently used in LBLRTM
! LBLRTM_v11.3 has chi=1 for all circumstances
!
!      X0 = 10.                                                            
!      Y0 = asubl*EXP(-x0/bsubl) 
!      F  = 1./bsubl
!      Y1 = -F*Y0                                                          
!      Y2 = Y1*((BSUBL-1)/X0-F)                                            
!      Z0 = (Y0-1)/X0**2                                                   
!      Z1 = Y1/(2*X0)                                                      
!      Z2 = Y2/2.                                                          
!      C6 = (Z0-Z1+(Z2-Z1)/4.)/X0**4                                       
!      C4 = (Z1-Z0)/X0**2-2.*X0**2*C6                                      
!      C2 = Z0-X0**2*C4-X0**4*C6                                           
!
      FI = abs(deltaXnu)
      IF (FI.LT.X0) THEN                                               
          CHI = 1.+C2*FI**2+C4*FI**4+C6*FI**6                   
      ELSE                                                            
          CHI = asubl*EXP(-FI/bsubl) 
      ENDIF                                                           

!**%%$$
      chi = 1.

      return
!
      end subroutine chi_fn

END MODULE ModmMod
