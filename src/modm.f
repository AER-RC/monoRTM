C     path:		$Source$
C     author:		$Author $
C     revision:	        $Revision$
C     created:	        $Date$

      SUBROUTINE MODM(ICP,NWN,WN,dvset,NLAY,P,T,
     &                 NMOL,WKL,WBRODL,
     &           SCLCPL,SCLHW,Y0RES,HFILE,ICNTNM,ixsect,ISPD)
C-------------------------------------------------------------------------------
C
C     PROGRAM:  MODM
C     -------
C
C     AUTHOR: Sid-Ahmed Boukabara 
C     ------
C
C     AFFILIATION: ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
C     -----------
C
C     DATE OF CREATION : October 1998
C     ----------------
C
C     AIM: This program is aimed at the calculation of the
C     ---  atmospheric optical depths. The spectral validity depends
C          only on the region covered by the file:"spectral_lines.dat"
C          The components treated here are the water vapor, the
C          oxygen, the ozone, the nitrogen and nitrogen dioxide.
C
C     INPUTS:
C     ------
C     - ICP      : Flag to take(if =1) or not (if=0) the line coupling
C     - NWN      : Number of wavenumbers to be treated
C     - WN       : Vector of NWN wavenumbers [in cm-1], one should note that
C                  this input could be a scalar (associated with NWN=1)
C     - NLAY     : Number of layers to be treated.
C     - P        : Vector of NLAY pressures (in mbar), one should note that
C                  this input could be a scalar (associated with NLAY=1)
C     - T        : Vector of NLAY temperatures [in Kelvin]
C     - CLW      : Vector of NLAY Cloud Liquid Water amounts [in kg/m2 or mm]
C                  When Cloud is present, the frequency must be consistent
C                  with Rayleigh absorption (no scattering performed in 
C                  monortm). 
C     - SCLCPL   : Scaling factor of the Line Coupling parameters (usually SCLCPL=1)
C     - SCLHW    : Scaling factor of the pressure dependence of the halfwidth
C                  of the zero frequency band (usually SCLHW=1)
C     - Y0RES    : Y0 resonnance (usually Y0RES=0) to be added to the Yi
C     - HFILE    : Name of the spectral lines information file (HITRAN)
C     - ICNTNM   : Flag to account for the continuum (if =1) or not (if =0)
C     - ISPD     : Flag to use the slow version (if 0) (all the lines present) 
C                  which is most accurate, or use the fast version (if 1)
C                  (only flagged lines used). valid only for the microwave.
C
C
C     OUTPUTS:
C     -------
C     - O      : An array of NWNxNLAY elts containing the total optical depths
C                 due to all the active species [in nepers]

C     History of the modifications:
C     *****************************  
C     - written in 1999 by Sid Ahmed Boukabara, Ross Hoffman
C	and Tony Clough. 
C     - validated against ARM sondes in the
C	microwave spectrum (23.8 and 31.4 GHz). SAB, 2000.
C     - extended to more species by Sid Ahmed Boukabara in 03/2000.
C     - cleaned up and commented in 2001 for first public release.
C	Also put under CVS configuration management. SAB.
C     - Extended O2 lines to submillimeter. Extensive validation
C	by comparison to Rosenkranz model and MWR data.
C	Update of the LBLATM module (accepts inputs at pressure 
C	grid, along with altitude grid). 
C	Fixed the handling of N2 amount coming from LBLATM (which
C	depends on the number of molecules NMOL). 
C	Adopted accurate constants values. 
C	Sid Ahmed Boukabara. Dec 14th 2001.
C     - Updated on January 7th 2002. ARM option (INP=2) updated and
C       made more efficient after Jim's comments. (INP=3) option optimized.
C       WV line intensities modified in the microwave (see Tony's email).
C     - Updated on October 2nd 2002. SAB. Speed option implemented.
C       We added also the possibility to run Voigt or Lorentz
C       line shape (speed up process) depending on the current
C       condition (parameter zeta). The pressure induced
C       shifted frequency is also passed to the line shape 
C       computation (instead of the spectroscopic wavenumber).
C     - September 2003: Modified spectral lines file to improve agreement 
C       with SGP MWRP data (provided by Nico Cimini). Scaled O2 line coupling
C       parameters: Y * 0.87, G* 0.
C     - 2006: Implementation of Tretyakov et al (2005) O2 line coupling.
C       Validated against ARM MWRP data (see Cadeddu et al, 2007)
C     - 2007: Updated spectral line file to change the widths and 
C       temperature dependences of the widths for the 22 and 183 GHz lines
C     - 2008: Extensive update to enable the use of MonoRTM beyond 
C       the microwave region and to use the MT_CKD continuum.
C       Updates to self and foreign broadened microwave continuum based
C       on ARM SGP MWR data at 31.4 GHz and ARM FKB (COPS) data at 150 GHz.  
C
C     Comments should be forwarded to Karen Cady-Pereira (cadyp@aer.com)
C     or Vivienne Payne (vpayne@aer.com).
C
C-------------------------------------------------------------------------------
      include "declar.incl"

      parameter (n_absrb=5050,ncont=5)
      real *8  v1abs,v2abs
      real*8 v1, v2
      real scor(42,9)
      integer index_cont(ncont), imol
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)                
C                                                                         
      CHARACTER*8      XID,       HMOLID,      YID     
      REAL*8               SECANT,       XALTZ
C                                                                         
      COMMON /CVRCNT/ HNAMCNT,HVRCNT
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       
     *                WKC(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   
     *                EMISIV,FSCDID(17),NMOLC,LAYER ,YI1,YID(10),LSTWDF

      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN

      COMMON /CNTSCL/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL

      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     $     RADCN1,RADCN2 

      CHARACTER HFILE*80,HVRMODM*15
      COMMON /CVRMODM/ HVRMODM
      LOGICAL INIT
      DATA INIT/.TRUE./
      SAVE INIT

      data index_cont/1,2,3,7,22/
      HVRMODM = '$Revision$' 

      IF(INIT)THEN
         CALL READ_HITR(ICP,HFILE,ISPD,MINWN,MAXWN) !reads the HITRAN data
         if(wn(1).lt.MINWN.OR.wn(nwn).gt.MAXWN) then
            print *,' '
            print *,'!!! WARNING !!!'
            print *,
     &         '   Spectral range of lines file does not cover requested 
     & spectral range.'
            print *,' '
            stop
         endif
         INIT=.FALSE.
      ENDIF

! Set up useful constants
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

!  Initialize
        oc(1:nwn,1:ncont,1:nlay) = 0.
        odxsec(1:nwn,1:nlay) = 0.

         if (ixsect .eq. 1) 
     &       call monortm_xsec_sub(wn,nwn,p,t,nlay)

      DO K=1,NLAY            !Loop over the Temp/Press/Amount profile 
! Set up variables for call to contnm
         pave = p(k)
         tave = t(k)
         wbroad = wbrodl(k)
         xkt = tave/radcn2

! Call CONTNM for each molecule
        wkc(1:nmol) = wkl(1:nmol,k)
        if (nmol.lt.22) wkc(index_cont(ncont)) = wbroad    ! set n2 to wbroad if nmol < 22
        do icount=1,ncont
              im = index_cont(icount) 
	      absrb(:) = 0.
	      call zero_cntnm(im,xself,xfrgn,xco2c,xo3cn,xo2cn,xn2cn)
	      call contnm(jrad)
	      call zero_cntnm(99,xself,xfrgn,xco2c,xo3cn,xo2cn,xn2cn)
c Interpolate for gridded spectral resolution in one step
	      if (dvset.ne.0) call xint(v1abs,v2abs,dvabs,absrb,1.0,v1,
     &                     dvset,oc(1:nwn,im,k),1,nwn) 
c Interpolate for specific wavenumbers one at a time
	      if (dvset.eq.0) then 
		 do iw=1,nwn
		    call xint(v1abs,v2abs,dvabs,absrb,1.0,wn(iw),1.0,
     &                   oc(iw,im,k),1,1)
		 end do
	      end if
c Multiply by radiation term
	       do iw=1,nwn
		  oc(iw,im,k) = oc(iw,im,k)*radfn(wn(iw),xkt)
	       end do
         end do                    ! end molecule loop

c calculate TIPS using Gamache routine rather that QOFT
         call tips_2003(nmol,t(k),scor)

	 DO M=1,NWN                !loop over the wavenumbers

            CALL INITI(P(K),T(K),RADCT,T0,P0,NMOL,
     &        WKL(1:NMOL,K),WBRODL(K),XN0,Xn,Xn_WV)         !INITIALIZATION
            RFT=WN(M)*TANH((RADCT*WN(M))/(2*T(K)))	              !RAD_FIELD_TERM

            CALL LINES(Xn,WN(M),T(K),NMOL,WKL(1:nmol,k)      ,         !PROCESS_LINES
     &      wbrodl(K),RADCT,T0,o_by_mol(m,1:nmol,k),XN0,RFT,
     &      P(K),P0,SCLCPL,SCLHW,Y0RES,scor)
   
            O_CLW(M,K)=ODCLW(WN(M),T(K),CLW(K))                       !OPTDEPTH CLW
            do imol = 1,nmol
                O(M,K) = o(m,k) + O_BY_MOL(M,imol,K) 
            enddo
            o(m,k) = o(m,k) + odxsec(m,k) + 
     &              +sum(oc(m,1:index_cont(5),k))+O_CLW(M,K)

         ENDDO

      ENDDO                     ! end layer loop
      RETURN
      END 


      SUBROUTINE LINES(Xn,WN,T,NMOL,WK,
     &     wbrod,RADCT, T0,o_by_mol,
     &     XN0,RFT,P,P0,SCLCPL,SCLHW,Y0RES,scor)
      PARAMETER (NNM=  39,IIM= 75000)
      REAL WK(NMOL),o_by_mol(nmol)
      REAL*8 WN,XNU,XNU0(NNM,IIM)
      REAL DELTNU(NNM,IIM),E(NNM,IIM),ALPS(NNM,IIM),XG(NNM,IIM)
      REAL ALPF(NNM,IIM),X(NNM,IIM),A(4),S0(NNM,IIM),B(4),TEMPLC(4)
      INTEGER NBLM(NNM),ISO(NNM,IIM)
      real scor(42,9)

      COMMON/HITR/NBLM,ISO,XNU0,DELTNU,S0,E,ALPS,ALPF,X,XG,NMOLEC
      DATA TEMPLC /200.0,250.0,296.0,340.0 /
      DATA MOL_WV/1/,MOL_O3/3/,MOL_O2/7/,MOL_N2/22/,MOL_N2O/4/
      DATA MOL_CO/5/,MOL_SO2/9/,MOL_NO2/10/,MOL_OH/13/

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
      RN=(Xn/XN0)               !ratio of number density
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
            IF ((XG(I,J).EQ.-1).OR.(XG(I,J).EQ.-3)) THEN
               JJ=J+1 !the LCC are stored in XNU0(J+1),DELTNU(J+1),etc..
               A(1)=XNU0(I,JJ)
               B(1)=DELTNU(I,JJ)
               A(2)=S0(I,JJ)
               B(2)=E(I,JJ)
               A(3)=ALPS(I,JJ)
               B(3)=ALPF(I,JJ)
               A(4)=X(I,JJ)
               B(4)=XG(I,JJ)
               AIP=A(ILC)+((A(ILC+1)-A(ILC))*RECTLC)*TMPDIF
               BIP=B(ILC)+((B(ILC+1)-B(ILC))*RECTLC)*TMPDIF
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
            Xnu=Xnu0(I,J)+(deltnu(I,J)*(Xn/xn0))
            !check line within 25cm-1 from WN, (except for O2, cause line coupling)
            IF ((ABS(WN-Xnu).GT.deltnuC).and.(I.NE.7)) 
     &           GOTO 30   

            XIPSF = scor(i,iso(i,j))

            CALL INTENS(T,S0(I,J),E(I,J),RADCT,T0,Xnu,STILD,XIPSF)
            XTILD=1-X(I,J)
            HWHM_C=HALFWHM_C(alpf(I,J),alps(I,J),RT,XTILD,RN,I,
     &           RAT)
            HWHM_D=HALFWHM_D(I,ISO(I,J),Xnu,T)
            IF(XG(I,J).EQ.-3.) THEN
               HWHM_C=HWHM_C*(1-(AIP*(RP))-(BIP*(RP2)))
            ENDIF
            zeta=HWHM_C/(HWHM_C+HWHM_D)
            ilshp=1                            !=0->Lorentz, =1->Voigt
            if ((ABS(WN-Xnu).GT.(10.*HWHM_D)).or.(zeta.GT.0.99)) ilshp=0
            IF (ilshp.eq.0) CALL LSF_LORTZ(XG(I,J),RP,RP2,AIP,BIP,
     &           HWHM_C,WN,Xnu,SLS,I)
            IF (ilshp.eq.1) CALL LSF_VOIGT(XG(I,J),RP,RP2,AIP,BIP,
     &           HWHM_C,WN,Xnu,SLS,HWHM_D,I)
            SF=SF+(STILD*SLS)
 30         CONTINUE
            J=JJ
         END DO
         SPSD=W_species*SF
         OL =RFT*SPSD
 10      o_by_mol(i)  = OL
      ENDDO
      END 

      FUNCTION HALFWHM_D(MOL,ISO,XNU,T)
      PARAMETER (NMOL=39,Nspeci=85)
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     $     RADCN1,RADCN2 
      REAL C,K,T,M,AVOG
      REAL*8 XNU  
      INTEGER ILOC,ISO
      COMMON /ISVECT/ ISO_MAX(NMOL),SMASS(nmol,9)
      common /iso_id/ iso_82(98)
      
      M=SMASS(mol,iso)
      HALFWHM_D=(XNU/CLIGHT)*SQRT(2.*LOG(2.)*((BOLTZ*T)/(M/AVOGAD)))
      END


      SUBROUTINE LSF_VOIGT(XF,RP,RP2,AIP,BIP,HWHM,WN,Xnu,SLS,
     &     AD,MOL)
      REAL*8 WN,Xnu,deltXNU,deltnuC,CHI
      DATA MOL_WV/1/,MOL_CO2/2/,MOL_O3/3/,MOL_O2/7/,MOL_N2/22/,
     &     MOL_N2O/4/
      deltnuC=25.          !cm-1
      DIFF=(WN+Xnu)-deltnuC
      SLS = 0.
      chi = 1.

C JULY 2008 VHP: 
C                Note that there is currently no pedestal subtraction here for O2.
C                This choice was made in order to avoid discontinuities due to O2 line coupling.
C                We could get around this by generating an O2 continuum in the same way
C                that we generate the CO2 continuum.
      
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
          IF ((ABS(WN-Xnu).LE.deltnuC).and.(XF.NE.-1).and.
     &           (XF.NE.-3)) THEN    !no line coupling 
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
C           multiply the co2 pedestal contribution by the chi factor
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
                  IF ((XF.EQ.-1).or.(XF.EQ.-3)) THEN ! CO2 line coupling
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
C                     no negative oscillation, since no CO2 lines within 25cm-1 of zero cm-1
                      XL3 = VOIGT(deltnuC,HWHM,AD) !VOIGT for 25cm-1 wn
                      IF (XF.EQ.-1) THEN
                          Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
                          XP4=XL3*
     &                        (1./((deltnuC)**2+HWHM**2)) *
     &                        (2.-((WN-Xnu)**2/((deltnuC)**2+HWHM**2))) ! co2 impact pedestal
                          YP1=(Y1-1.)*
     &                        (2.-(WN-Xnu)**2/((deltnuC)**2+HWHM**2)) ! line coupling contributions to pedestal
                          SLS=CHI*(XL1*
     &                        (Y1)-XP4-XL3*(YP1))
                      ELSE
                          deltXNU = (WN-Xnu)
                          CALL CHI_FN(deltXNU,CHI)
                          Y1 = 1.
                          Y2 = 1.
                          XP4=XL3*
     &                        (1./((deltnuC)**2+HWHM**2)) *
     &                        (2.-(WN-Xnu)**2/((deltnuC)**2+HWHM**2)) ! co2 impact pedestal
                          SLS = CHI*(XL1-XP4)
                      ENDIF
                  ENDIF ! end test for CO2 line coupling
              ENDIF ! end test to distinguish between CO2 and O2
          ENDIF ! end test for line coupled lines within CO2 OR O2
      ENDIF ! end test for any possibility of line coupling

      END

      
      SUBROUTINE LSF_LORTZ(XF,RP,RP2,AIP,BIP,HWHM,WN,Xnu,SLS,
     &    MOL)
      REAL*8 WN,Xnu,deltnuC,deltXNU,CHI
      DATA MOL_WV/1/,MOL_CO2/2/,MOL_O3/3/,MOL_O2/7/,MOL_N2/22/,
     &     MOL_N2O/4/
      deltnuC=25.          !cm-1
      DIFF=(WN+Xnu)-deltnuC
       SLS = 0.
       chi= 1.

C JULY 2008 VHP: 
C                Note that there is currently no pedestal subtraction here for O2.
C                This choice was made in order to avoid discontinuities due to O2 line coupling.
C                In the future we could generate an O2 continuum in the same way
C                that we generate the CO2 continuum.

      IF ((MOL.NE.MOL_O2).AND.(MOL.NE.MOL_CO2)) THEN ! no possibility of line coupling
                                                     ! check for line within 25cm-1 has already
                                                     ! been performed in modm.f

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
      ELSE                      ! O2 or CO2 (check for line within 25 cm-1 has to be performed for O2)
          IF ((ABS(WN-Xnu).LE.deltnuC).and.(XF.NE.-1).and.
     &           (XF.NE.-3)) THEN    !no line coupling 
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
C           multiply the co2 pedestal contribution by the chi factor
                  XL3=XL3*(2.-(deltXNU**2/(deltXNU**2+HWHM**2))) ! co2 pedestal
                  SLS = CHI*(XL1-XL3) / HWHM
              ENDIF
          ELSE ! line has line coupling
              IF (MOL.EQ.MOL_O2) THEN ! for O2 line-coupled lines the 25 cm-1 check is not implemented
                  IF ((XF.EQ.-1).or.(XF.EQ.-3)) THEN !O2 line coupling
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
                  IF ((XF.EQ.-1).or.(XF.EQ.-3)) THEN ! CO2 line coupling
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
c                     no negative oscillation, since no co2 lines within 25cm-1 of zero
                      XL3 = XLORENTZ((deltnuC)/HWHM) !LORENTZ for 25cm-1 wn
                      IF (XF.EQ.-1) THEN
                          Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
                          XP4=XL3*
     &                        (1./((deltnuC)**2+HWHM**2)) *
     &                        (2.-((WN-Xnu)**2/((deltnuC)**2+HWHM**2))) ! co2 impact pedestal
                          YP1=(Y1-1.)*
     &                        (2.-(WN-Xnu)**2/((deltnuC)**2+HWHM**2))
                          SLS=CHI*(XL1*
     &                        (Y1)-XP4-XL3*(YP1))/HWHM
                      ELSE
                          Y1 = 1.
                          Y2 = 1.
                          XP4=XL3*
     &                        (1./((deltnuC)**2+HWHM**2)) *
     &                        (2.-((WN-Xnu)**2/((deltnuC)**2+HWHM**2))) ! co2 impact pedestal
                          SLS = CHI*(XL1 - XP4) / HWHM
                      ENDIF
                  ENDIF ! end test for CO2 line coupling
              ENDIF ! end test to distinguish between CO2 and O2
          ENDIF ! end test for line coupled lines within CO2 OR O2
      ENDIF ! end test for any possibility of line coupling

      END


      
      FUNCTION HALFWHM_C(AF,AS,RT,XTILD,RN,MOL,RAT)
      IF ((MOL.EQ.1).AND.(AS.EQ.0.)) AS=5*AF
      HALFWHM_C=AF*(RT**XTILD)*(RN-RAT)+AS*(RT**XTILD)*RAT
      END

      
      SUBROUTINE INTENS(T,S0s,Es,RADCT,T0,Xnus,STILD,XIPSF)
      REAL*8 Xnus
      S=S0s*(EXP(-RADCT*Es/T)/EXP(-RADCT*Es/T0))*XIPSF
      STILD=S*((1+EXP(-(RADCT*Xnus/T)))/
     &(Xnus*(1-EXP(-(RADCT*Xnus/T0)))))
      END 

      SUBROUTINE INITI(P,T,RADCT,T0,P0,NMOL,WK,WBROD,XN0,
     $                 Xn,Xn_WV)
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     $     RADCN1,RADCN2 
      DATA WVMOLMASS /18.016 /, DRYMOLMASS/28.97/   
      REAL wk(39),wbrod
      RADCT=PLANCK*CLIGHT/BOLTZ !in K/cm-1
      T0=296.                   !in K
      P0=1013.25                !in HPa
      XN0=(P0/(BOLTZ*T0))*1.E+3
      XN =(P /(BOLTZ*T ))*1.E+3
      WDRY=sum(wk(2:nmol))+wbrod
      RATIOMIX=(Wk(1)*WVMOLMASS)/(WDRY*DRYMOLMASS)
      WVPRESS=(RATIOMIX/(RATIOMIX+(WVMOLMASS/DRYMOLMASS)))*P
      Xn_WV=(WVPRESS/(BOLTZ*T))*1.E+3
      END 
      
        
      SUBROUTINE READ_HITR(ICP,HFILE,ISPD,MINWN,MAXWN)
      PARAMETER (NNM=  39,IIM=  75000)
      REAL*8 XNU0(NNM,IIM),nu0
      REAL DELTNU(NNM,IIM),E(NNM,IIM),ALPS(NNM,IIM),ALPF(NNM,IIM)
      REAL X(NNM,IIM),XG(NNM,IIM),S0(NNM,IIM)
      INTEGER NBLM(NNM),ISO(NNM,IIM)
      Character Q1*9,Q2*9,HFILE*80,CXID*1,HVRSPEC*15
      character(len=100) cxidline
      character*25 cdate
      COMMON /CVRSPEC/ HVRSPEC
      COMMON/HITR/NBLM,ISO,XNU0,DELTNU,S0,E,ALPS,ALPF,X,XG,NMOLEC
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,  
     $     NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,       
     $     NLTEFL,LNFIL4,LNGTH4                                 
      DATA MOL_WV/1/,MOL_O3/3/,MOL_O2/7/,MOL_N2/22/,MOL_N2O/4/
      DATA MOL_CO/5/,MOL_SO2/9/,MOL_NO2/10/,MOL_OH/13/
      EQUIVALENCE (CXID,CXIDLINE)                    

      NMOLEC=NNM

      OPEN(9,FILE=HFILE,form='formatted',status='old',ERR=20)
      !---We read comments /put them in LOG file
      IKOUNT = 0
      WRITE(IPR,'(a80)') ' '
      WRITE(IPR,'(a80)') ' '
      WRITE(IPR,'(a80)') '********************************'
      WRITE(IPR,'(a80)') ' SPECTROSCOPIC FILE INFORMATION '
      WRITE(IPR,'(a80)') '********************************'
      WRITE(*,'(a33)') '  '
      WRITE(*,'(a33)') '********************************'
      WRITE(*,'(a33)') ' SPECTROSCOPIC FILE INFORMATION '
      WRITE(*,'(a33)') '********************************'
 23   READ (9,'(A100)',END=30,ERR=30) CXIDLINE
 
      WRITE(IPR,'(a100)') CXIDLINE
      if (ikount.eq.2) then
         read (CXIDLINE,'(12X,A15)') HVRSPEC
         write(*,*) HVRSPEC
      end if
      if (ikount.eq.3) then
         read (CXIDLINE,'(12X,A25)') CDATE
         write(*,*) CDATE
         write(*,*) ' '
      endif
c      if (ikount.ge.12 .AND. ikount.le.14) then
c         write(*,'(a100)') cxidline(2:100)
c      endif

      ikount = ikount + 1

      IF (CXID.EQ.'$'.OR.CXID.EQ.'>'.OR.CXID.EQ.'%'.OR.CXID.EQ.'C'.OR.
     &    CXID.EQ.'c') then
         GO TO 23                            
      else
         backspace (unit=9)
      end if

      ILINE=0
      ILINE_O2=0

 22   READ(9,100,END=3,ERR=30) mo,iso_scal,nu0,S0_scal,R,
     &     alpf_scal,alps_scal,E_scal,X_scal,deltnu_scal,
     &     iv1,iv2,Q1,Q2,ier1,ier2,ier3,iref1,iref2,iref3,iflg
      !---test to limit the number of lines used-----
      !if ((nu0.ge.60.)) goto 3
      !---test to limit the strength of the WV lines used
      !if(((S0_scal.lt.1.E-25)).and.(mo.eq.1)) goto 22
      !---test to limit the isotopes of WV 
      !if((iso_scal.ne.1).and.(mo.eq.1)) goto 22
      !---test to remove the wavenumber shift effect
      !deltnu_scal=0.
      !---test to use only those lines that are flagged
      IF (ISPD .eq. 1) then
         IF (iflg.ne.mo) then 
            IF ((iref3.EQ.-1).OR.(iref3.EQ.-3)) THEN 
               READ(9,2,END=3,ERR=30)JF,XF,DF,SF,EF,AS,AF,XS,XF
            ENDIF
            goto 22
         ENDIF
      ENDIF
      !-----------------------------------------------
      IF (mo.EQ.mol_o2) ILINE_O2=ILINE_O2+1
      ILINE=ILINE+1
      IF (ILINE.EQ.1) MINWN=nu0
      jj = mo
      nblm(jj) = nblm(jj)+1
      ii = nblm(jj)

      ISO(JJ,II)=iso_scal
      XNU0(JJ,II)=nu0
      DELTNU(JJ,II)=deltnu_scal
      S0(JJ,II)=S0_scal!provoques underflow message
      E(JJ,II)=E_scal
      ALPS(JJ,II)=alps_scal
      ALPF(JJ,II)=alpf_scal
      X(JJ,II)=X_scal
      XG(JJ,II)=iref3
      IF (((iref3.EQ.-1).OR.(iref3.EQ.-3)).AND.(ICP.EQ.1)) THEN !Lines Coupling
         nblm(jj) = nblm(jj)+1
         ii = nblm(jj)
         READ(9,2,END=3,ERR=30)J,XNU0(JJ,II),DELTNU(JJ,II),
     &        S0(JJ,II),E(JJ,II),ALPS(JJ,II),ALPF(JJ,II),
     &        X(JJ,II),XG(JJ,II)
      ENDIF
      IF (((iref3.EQ.-1).OR.(iref3.EQ.-3)).AND.(ICP.EQ.0)) THEN !Lines Cplng ignored
         READ(9,2,END=3,ERR=30)JF,XF,DF,SF,EF,AS,AF,XS,XF
         XG(JJ,II)=0.
      ENDIF
      GOTO 22
 3    CONTINUE
      MAXWN=NU0
      WRITE(*,'(a33)') '***********************************'
      WRITE(*,*) ' '
      IF (ISPD.EQ.1) THEN
         WRITE(*,*) '****************************************'
         WRITE(*,*) '*            W A R N I N G             *'
         WRITE(*,*) '****************************************'
         WRITE(*,*) 'FAST VERSION IS RUNNING.'
         WRITE(*,*) 'CURRENTLY VALID IN THE MICROWAVE ONLY.'
      ENDIF
      WRITE(*,*) ' '
      WRITE(*,*) '****************************************'
      WRITE(*,*) '* SPECTRAL LINES INFORMATION AVAILABLE *'
      WRITE(*,*) '****************************************'
      WRITE(*,*) 'Minimum Wavenumber:',MINWN,' cm-1'
      WRITE(*,*) 'Maximum Wavenumber:',MAXWN,' cm-1'
      WRITE(*,*)
      WRITE(*,'(2x,a8,8x,a8,3x,a8)') 'Molecule','# lines' 
      DO J=1,NMOLEC
          WRITE(IPR,'(2x,i8,5x,i8,5x,i8)') J, NBLM(J)
      ENDDO
      WRITE(*,*)
      WRITE(*,*) '****************************************'
      CLOSE(9)
 100  FORMAT (I2,I1,F12.6,1P,2E10.3,0P,2F5.4,F10.4,F4.2,F8.6,
     &     2I3,2A9,3I1,3I2,I2)
1000  FORMAT (I2,I1,F12.6,1P,2E10.3)
 2    FORMAT (I2,1P,4(E13.6,E11.4),0P,I2)                           
      RETURN
 20   PRINT *, 'ERROR OPENING HITRAN FILE:',HFILE
      STOP
 30   PRINT *, 'ERROR READING HITRAN FILE:',HFILE
      STOP
      END   

      INCLUDE 'isotope.dat'     !INCLUDE HITRSID DATA-STATEMENTS

      FUNCTION ODCLW(WN,TEMP,CLW)
      !INPUTS: WN (WaveNUmber in cm-1)
      !        Temp (in K)
      !        CLW  (in mm or kg/m2)
      !OUTPUT: ODCLW: optical depth of the Cloud Liquid Water
      !FROM EQUATIONS OF LIEBE, HUFFORD AND MANABE
      !Ref:(INT. J. IR & MM WAVES V.12(17) JULY 1991
      COMPLEX EPS,RE
      REAL*8 WN
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     $     RADCN1,RADCN2 
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
      EPS = (EPS0-EPS1)/CMPLX(1.,FREQ/FP) +
     $     (EPS1-EPS2)/CMPLX(1.,FREQ/FS) +EPS2
      RE = (EPS-1.)/(EPS+2.)
      !ODCLW = -.06286057*CLW*AIMAG(RE)*FREQ
      ODCLW = -(6.*PI/299.792458)*CLW*AIMAG(RE)*FREQ
      RETURN
      END


      FUNCTION XLORENTZ(Z)
      REAL*8 Z
      REAL XLORENTZ
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     $     RADCN1,RADCN2 
      XLORENTZ=1./(PI*(1.+(Z**2))) 
      RETURN
      END



      FUNCTION VOIGT(DELTNU,ALPHAL,ALPHAD)
      COMPLEX v,w4
      REAL*8 DELTNU
      REAL avc(0:101)
      REAL ALPHAL,ALPHAD,ZETA,ALPHAV,AVCINTERP,DNU
      REAL AL,AD,X,Y,ANORM1,VOIGT,DZETA
      INTEGER IZETA2,IZETA1
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     $     RADCN1,RADCN2 
      DATA AVC/                                       
     *  .10000E+01,.99535E+00,.99073E+00,.98613E+00,.98155E+00,   
     *  .97700E+00,.97247E+00,.96797E+00,.96350E+00,.95905E+00,   
     *  .95464E+00,.95025E+00,.94589E+00,.94156E+00,.93727E+00,   
     *  .93301E+00,.92879E+00,.92460E+00,.92045E+00,.91634E+00,   
     *  .91227E+00,.90824E+00,.90425E+00,.90031E+00,.89641E+00,   
     *  .89256E+00,.88876E+00,.88501E+00,.88132E+00,.87768E+00,   
     *  .87410E+00,.87058E+00,.86712E+00,.86372E+00,.86039E+00,   
     *  .85713E+00,.85395E+00,.85083E+00,.84780E+00,.84484E+00,   
     *  .84197E+00,.83919E+00,.83650E+00,.83390E+00,.83141E+00,   
     *  .82901E+00,.82672E+00,.82454E+00,.82248E+00,.82053E+00,   
     *  .81871E+00,.81702E+00,.81547E+00,.81405E+00,.81278E+00,   
     *  .81166E+00,.81069E+00,.80989E+00,.80925E+00,.80879E+00,   
     *  .80851E+00,.80842E+00,.80852E+00,.80882E+00,.80932E+00,   
     *  .81004E+00,.81098E+00,.81214E+00,.81353E+00,.81516E+00,   
     *  .81704E+00,.81916E+00,.82154E+00,.82418E+00,.82708E+00,   
     *  .83025E+00,.83370E+00,.83742E+00,.84143E+00,.84572E+00,   
     *  .85029E+00,.85515E+00,.86030E+00,.86573E+00,.87146E+00,   
     *  .87747E+00,.88376E+00,.89035E+00,.89721E+00,.90435E+00,   
     *  .91176E+00,.91945E+00,.92741E+00,.93562E+00,.94409E+00,   
     *  .95282E+00,.96179E+00,.97100E+00,.98044E+00,.99011E+00,   
     *  .10000E+01,.10000E+01/                             
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
      end



C*************************************************************************
C     FORTRAN function for the complex probability function w(z).
C     COMPUTES THE COMPLEX PROBABILITY FUNCTION W(Z)=EXP(-Z**2)*ERFC(-I*Z)
C     IN THE UPPER HALF-PLANE Z=X+I*Y (I.E. FOR Y>=0)
C     MAXIMUM RELATIVE ERROR OF BOTH REAL AND IMAGINARY PARTS IS <1*10**(-4)
C     Reference:
C     Humlicek, J., 1982; Optimized Computation of the Voigt and Complex
C     Probability Functions, J. Quant. Spectrosc. Radiat. Transfer, 27, 437-444.
C     S-A Boukabara AER Inc, 2000
C*************************************************************************
      FUNCTION W4(x,y)
      COMPLEX W4,T,U
      
      T=CMPLX(Y,-X)
      S=ABS(X)+Y
      IF(S.LT.15.)GOTO 1
C     ***   REGION I
      W4=T*.5641896/(.5+T*T)
      RETURN
 1    IF(S.LT.5.5)GOTO 2
C     ***   REGION II
      U=T*T
      W4=T*(1.410474+U*.5641896)/(.75+U*(3.+U))
      RETURN
 2    IF(Y.LT. 0.195*ABS(X)-0.176)GOTO 3
C     ***   REGION III
      W4=(16.4955+T*(20.20933+T*(11.96482+
     &     T*(3.778987+T*.5642236))))/
     &     (16.4955+T*(38.82363+T*(39.27121+
     &     T*(21.69274+T*(6.699398+T)))))
      RETURN
C     ***   REGION IV
 3    U=T*T
      W4=CEXP(U)-T*(36183.31-U*(3321.9905-
     &     U*(1540.787-U*(219.0313-U*
     &     (35.76683-U*(1.320522-U*.56419))))))/
     &     (32066.6-U*(24322.84-U*
     &     (9022.228-U*(2186.181-U*(364.2191-
     &     U*(61.57037-U*(1.841439-U)))))))
      RETURN
      END

      subroutine zero_cntnm(molec,xself,xfrgn,xco2c,xo3cn,xo2cn,xn2cn)

c zeroes out all continua EXCEPT for molecule molec, but saves initial scaling factors
c if molec is 99 then restores initial factors
      
      save xself_safe, xfrgn_safe, xco2c_safe,xo3cn_safe, xo2cn_safe, 
     &     xn2cn_safe

      if(molec.ne.99) then
	 xself_safe = xself
	 xfrgn_safe = xfrgn
	 xco2c_safe = xco2c
	 xo2cn_safe = xo2cn
	 xo3cn_safe = xo3cn
	 xn2cn_safe = xn2cn
      end if
      select case (molec)
         case (1)
            xco2c = 0.
            xo2cn = 0.
            xo3cn = 0.
            xn2cn = 0.
         case (2)
            xself = 0.
            xfrgn = 0.
            xo2cn = 0.
            xo3cn = 0.
            xn2cn = 0.
         case (3)
            xself = 0.
            xfrgn = 0.
            xco2c = 0.
            xo2cn = 0.
            xn2cn = 0.
         case (7)
            xself = 0.
            xfrgn = 0.
            xco2c = 0.
            xo3cn = 0.
            xn2cn = 0.
         case (22)
            xself = 0.
            xfrgn = 0.
            xco2c = 0.
            xo2cn = 0.
            xo3cn = 0.
         case (99)
            xself= xself_safe
            xfrgn= xfrgn_safe
            xco2c= xco2c_safe
            xo2cn= xo2cn_safe
            xo3cn= xo3cn_safe
            xn2cn= xn2cn_safe
          case default
            xself = 0.
            xfrgn = 0.
            xco2c = 0.
            xo2cn = 0.
            xo3cn = 0.
            xn2cn = 0.
       end select

       return
       end
      
      subroutine chi_fn (deltXNU,chi)
c
      real chi, deltXNU
c
      DATA ASUBL / 0.800 /,BSUBL / 10.0 /
C                                                                         
C     SET UP CHI SUB-LORENTZIAN FORM FACTOR FOR CARBON DIOXIDE            
C     POLYNOMIAL MATCHED TO AN EXPONENTIAL AT X0 = 10 CM-1                
C                                                                         

C Commented lines here are lines contained, but NOT currently used in LBLRTM
C LBLRTM_v11.3 has chi=1 for all circumstances

C      X0 = 10.                                                            
C      Y0 = asubl*EXP(-x0/bsubl) 
C      F  = 1./bsubl
C      Y1 = -F*Y0                                                          
C      Y2 = Y1*((BSUBL-1)/X0-F)                                            
C      Z0 = (Y0-1)/X0**2                                                   
C      Z1 = Y1/(2*X0)                                                      
C      Z2 = Y2/2.                                                          
C      C6 = (Z0-Z1+(Z2-Z1)/4.)/X0**4                                       
C      C4 = (Z1-Z0)/X0**2-2.*X0**2*C6                                      
C      C2 = Z0-X0**2*C4-X0**4*C6                                           
C
      FI = abs(deltaXnu)
      IF (FI.LT.X0) THEN                                               
          CHI = 1.+C2*FI**2+C4*FI**4+C6*FI**6                   
      ELSE                                                            
          CHI = asubl*EXP(-FI/bsubl) 
      ENDIF                                                           

c**%%$$
      chi = 1.

      return
c
      end
