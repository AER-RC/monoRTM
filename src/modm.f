      SUBROUTINE MODM(IVC,ICP,NWN,WN,NLAY,P,T,W_wv,
     &     W_o2,W_o3,W_n2,W_n2O,W_co,W_so2,W_no2,
     &     W_oh,CLW,O,OL_WV,OS_WV,OF_WV,OL_O2,OL_O3,
     &     OL_N2,OC_N2,OL_N2O,OL_CO,OL_SO2,OL_NO2,
     &     OL_OH,O_CLW,XSLF,XFRG,XCN2,SCLCPL,SCLHW,
     &     Y0RES,HFILE,ICNTNM)

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
C     - IVC   : Flag of contin. vers.: CKD2.4(if=2)  MPMf87/s93 (if=3)
C     - ICP   : Flag to take(if =1) or not (if=0) the line coupling
C     - NWN   : Number of wavenumbers to be treated
C     - WN    : Vector of NWN wavenumbers [in cm-1], one should note that
C               this input could be a scalar (associated with NWN=1)
C     - NLAY  : Number of layers to be treated.
C     - P     : Vector of NLAY pressures (in mbar), one should note that
C               this input could be a scalar (associated with NLAY=1)
C     - T     : Vector of NLAY temperatures [in Kelvin]
C     - W_WV  : Vector of NLAY water vapor column amounts [in molecules/cm2]
C     - W_O2  : Vector of NLAY oxygen column amounts [in molecules/cm2]
C     - W_O3  : Vector of NLAY ozone column amounts [in molecules/cm2]
C     - W_N2  : Vector of NLAY nitrogen column amounts [in molecules/cm2]
C     - W_N2O : Vector of NLAY N2O column amounts [in molecules/cm2]
C     - W_CO  : Vector of NLAY CO column amounts [in molecules/cm2]
C     - W_SO2 : Vector of NLAY SO2 column amounts [in molecules/cm2]
C     - W_NO2 : Vector of NLAY NO2 column amounts [in molecules/cm2]
C     - W_OH  : Vector of NLAY OH column amounts [in molecules/cm2]
C     - CLW   : Vector of NLAY Cloud Liquid Water amounts [in kg/m2 or mm]
C               When Cloud is present, the frequency must be consistent
C               with Rayleigh absorption (no scattering performed in 
C               monortm). 
C     - XSLF  : Scaling factor of the self WV continuum (usually XSLF=1)
C     - XFRG  : Scaling factor of the foreign WV continuum (usually XFRG=1)
C     - XCN2  : Scaling factor of the N2 continuum (usually XCN2=1)
C     - SCLCPL: Scaling factor of the Line Coupling parameters (usually SCLCPL=1)
C     - SCLHW : Scaling factor of the pressure dependence of the halfwidth
C               of the zero frequency band (usually SCLHW=1)
C     - Y0RES : Y0 resonnance (usually Y0RES=0) to be added to the Yi
C     - HFILE : Name of the spectral lines information file (HITRAN)
C     - ICNTNM: Flag to account for the continuum (if =1) or not (if =0)
C
C
C     OUTPUTS:
C     -------
C     - O      : An array of NWNxNLAY elts containing the total optical depths
C                 due to all the active species [in nepers]
C     - OL_WV  : An array of NWNxNLAY elts containing the water vapor optical
C                 depth (due to lines only), [in Nepers]
C     - OS_WV  : An array of NWNxNLAY elts containing the water vapor optical
C                 depth (due to self continuum), [in Nepers]
C     - OF_WV  : An array of NWNxNLAY elts containing the water vapor optical
C                 depth (due to foreign continuum), [in Nepers]
C     - OL_O2  : An array of NWNxNLAY elts containing the oxygen optical
C                 depth (due to lines only), [in Nepers]
C     - OL_O3  : An array of NWNxNLAY elts containing the ozone optical
C                 depth (due to lines only), [in Nepers]
C     - OL_N2  : An array of NWNxNLAY elts containing the nitrogen optical
C                 depth (due to lines only), [in Nepers]
C     - OC_N2  : An array of NWNxNLAY elts containing the nitrogen optical
C                 depth (due to continuum), [in Nepers]
C     - OL_N2O : An array of NWNxNLAY elts containing the N2O optical
C                 depth (due to lines only), [in Nepers]
C     - OL_CO  : An array of NWNxNLAY elts containing the CO optical
C                 depth (due to lines only), [in Nepers]
C     - OL_SO2 : An array of NWNxNLAY elts containing the SO2 optical
C                 depth (due to lines only), [in Nepers]
C     - OL_NO2 : An array of NWNxNLAY elts containing the NO2 optical
C                 depth (due to lines only), [in Nepers]
C     - OL_OH  : An array of NWNxNLAY elts containing the OH optical
C                 depth (due to lines only), [in Nepers]
C     - O_CLW  : An array of NWNxNLAY elts containing the CLW optical
C                 depth , [in Nepers]
C
C-------------------------------------------------------------------------------
      include "declar.incl"

      CHARACTER HFILE*80
      LOGICAL INIT
      DATA INIT/.TRUE./
      SAVE INIT
      IF(INIT)THEN
         CALL VECISO               !set-up the isotopes
         CALL READ_HITR(ICP,HFILE) !reads the HITRAN data
         INIT=.FALSE.
      ENDIF
      DO M=1,NWN                !loop over the wavenumbers
         DO K=1,NLAY            !Loop over the Temp/Press/Amount profile 
            CALL INITI(P(K),T(K),RADCT,T0,P0,W_wv(K),W_o2(K),         !INITIALIZATION
     &      W_o3(K),W_n2(K),W_n2O(K),W_co(K),W_so2(K),
     &      W_no2(K),W_oh(K),XN0,Xn,Xn_WV,RHOFAC)
            RFT=WN(M)*TANH((RADCT*WN(M))/(2*T(K)))	              !RAD_FIELD_TERM
            CALL LINES(Xn,WN(M),T(K),W_wv(K),W_o2(K),W_o3(K),         !PROCESS_LINES
     &      W_n2(K),W_n2O(K),W_co(K),W_so2(K),W_no2(K),
     &      W_oh(K),RADCT,T0,OL_WV(M,K),OL_O2(M,K),
     &      OL_O3(M,K),OL_N2(M,K),OL_N2O(M,K),OL_CO(M,K),
     &      OL_SO2(M,K),OL_NO2(M,K),OL_OH(M,K),XN0,RFT,
     &      P(K),P0,SCLCPL,SCLHW,Y0RES)
            O_CLW(M,K)=ODCLW(WN(M),T(K),CLW(K))                           !OPTDEPTH CLW
            IF (ICNTNM.EQ.1) THEN 
               OS_WV(M,K)=SWV(IVC,WN(M),T(K),T0,W_wv(K),RFT,Xn,           !CONT_SELF_WV
     &              Xn_WV,XN0,XSLF) 
               OF_WV(M,K)=FWV(IVC,WN(M),W_wv(K),RFT,Xn,Xn_WV,XN0,XFRG)    !CONT_FRGN_WV 
               OC_N2(M,K)=CONTI_N2(WN(M),T(K),T0,W_n2(K),RFT,RHOFAC,
     &              XCN2)                                                 !CONT_N2
               O(M,K)=OL_WV(M,K)+OL_O2(M,K)+OL_O3(M,K)+OL_N2(M,K)+        !TOTAL_OPTDEPTH
     &              OL_N2O(M,K)+OL_CO(M,K)+OL_SO2(M,K)+OL_NO2(M,K)+
     &              OL_OH(M,K)+OS_WV(M,K)+OF_WV(M,K)+OC_N2(M,K)+
     &              O_CLW(M,K)
            ENDIF
            IF (ICNTNM.NE.1) THEN 
               O(M,K)=OL_WV(M,K)+OL_O2(M,K)+OL_O3(M,K)+OL_N2(M,K)+        !TOTAL_OPTDEPTH
     &              OL_N2O(M,K)+OL_CO(M,K)+OL_SO2(M,K)+OL_NO2(M,K)+
     &              OL_OH(M,K)+O_CLW(M,K)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END 

      FUNCTION FWV(IVC,WN,W_wv,RFT,Xn,Xn_WV,XN0,XFRG)
      REAL*8 WN
      IF (W_wv.EQ.0.) THEN
         FWV=0.
         RETURN
      ENDIF
      IF (IVC.EQ.2) THEN                               !CKD2.4 CONTINUUM
         FWV=FWV24(WN,W_wv,RFT,Xn,Xn_WV,XN0,XFRG)      !CNT_FRG_WV CKD2.4
      ENDIF
      IF (IVC.EQ.3) THEN                                 !MPMf87s93 CONTINUUM
         FWV=FWV_MPMf87s93(WN,W_wv,RFT,Xn,Xn_WV,XN0,XFRG)!CNT_FRG_WV 
      ENDIF
      RETURN
      END 

      FUNCTION SWV(IVC,WN,T,T0,W_wv,RFT,Xn,Xn_WV,XN0,XSLF)
      REAL*8 WN
      IF (W_wv.EQ.0.) THEN
         SWV=0.
         RETURN
      ENDIF
      IF (IVC.EQ.2) THEN                               !CKD2.4 CONTINUUM
         SWV=SWV24(WN,T,T0,W_wv,RFT,Xn,Xn_WV,XN0,XSLF) !CNT_SLF_WV CKD2.4 
      ENDIF
      IF (IVC.EQ.3) THEN                                !MPMf87s93 CKD2.4 CONT.
         SWV=SWV_MPMf87s93(WN,T,T0,W_wv,RFT,Xn,Xn_WV,XN0,XSLF)!CNT_SLF_WV
      ENDIF
      RETURN
      END 



      FUNCTION SWV_MPMf87s93(WN,T,T0,W_wv,RFT,Xn,Xn_WV,XN0,XSLF)
      REAL*8 WN,X(4)
      COMMON /SH2O/ V1,V2,DV,NPTSLFWV,SWV296(2003)!---UNITS(CM**3/MOL)*1.E-20   
      COMMON /S260/ V1_,V2_,DV_,NPTSLFWV_,SWV260(2003)
      J=INT((WN-V1)/DV)+1
      X(1)=((SWV296(J-1))*(SWV260(J-1)/SWV296(J-1))**((T-T0)/(260.-T0)))
      X(2)=((SWV296(J))*(SWV260(J)/SWV296(J))**((T-T0)/(260.-T0)))
      X(3)=((SWV296(J+1))*(SWV260(J+1)/SWV296(J+1))**((T-T0)/(260.-T0)))
      X(4)=((SWV296(J+2))*(SWV260(J+2)/SWV296(J+2))**((T-T0)/(260.-T0)))
      XF     = ( WN- (V1 + DV * FLOAT(J-1)) ) / DV
      SFAC = 3. 
      SWV_MPMf87s93=(W_wv)*RFT*
     &     (Xn_WV/XN0)*XLGR(XF,X)*1.E-20*SFAC*XSLF
      RETURN
      END 



      FUNCTION FWV_MPMf87s93(WN,W_wv,RFT,Xn,Xn_WV,XN0,XFRG)
      REAL*8 WN,X(4)
      COMMON /FH2O/ V1,V2,DV,NPTFH2O,FH2O(2003)                              
      J=INT((WN-V1)/DV)+1
      DO I=1,4
         X(I)=FH2O(J+I-2)
      ENDDO
      XF     = ( WN- (V1 + DV * FLOAT(J-1)) ) / DV
      FSCAL = 0.8
      FWV_MPMf87s93=XLGR(XF,X)*1.E-20*
     &     ((W_wv)*RFT*((Xn-Xn_WV)/XN0))*FSCAL*XFRG  
      RETURN
      END 



      FUNCTION FWV24(WN,W_wv,RFT,Xn,Xn_WV,XN0,XFRG)
      REAL*8 V0F1,V0F1a,V0F2,V0F3,VF2,VF4,VF6
      REAL*8 WN,X(4)
      COMMON /FH2O/ V1,V2,DV,NPTFH2O,FH2O(2003)     
      DATA V0F1,HWSQF1,BETAF1 /350.,40000.,5.e-09 /
      DATA FACTRF1,V0F1a,HWSQF1a /-0.7,630.,4225/
      DATA BETAF1a,FACTRF1a,V0F2 /2.e-08,+0.75,1130./
      DATA HWSQF2,BETAF2,FACTRF2 /108900.,8.E-11,-0.97/
      DATA V0F3,HWSQF3,BETAF3,FACTRF3 /1975.,62500.,5.E-06,-0.65/
      J=INT((WN-V1)/DV)+1
      DO I=1,4
         X(I)=FH2O(J+I-2)
      ENDDO
      XF     = ( WN- (V1 + DV * FLOAT(J-1)) ) / DV
      !---added correction to the forgn continuum
      VF2 = (WN-V0F1)**2
      VF6 = VF2 * VF2 * VF2
      FSCAL = (1.+FACTRF1*(HWSQF1/(VF2+(BETAF1*VF6)+HWSQF1)))
      VF2 = (WN-V0F1a)**2
      VF6 = VF2 * VF2 * VF2
      FSCAL = FSCAL* 
     *     (1.+FACTRF1a*(HWSQF1a/(VF2+(BETAF1a*VF6)+HWSQF1a)))
      VF2 = (WN-V0F2)**2
      VF6 = VF2 * VF2 * VF2
      FSCAL = FSCAL* 
     *     (1.+FACTRF2*(HWSQF2/(VF2+(BETAF2*VF6)+HWSQF2)))
      VF2 = (WN-V0F3)**2
      VF4 = VF2*VF2
      FSCAL = FSCAL* 
     *     (1.+FACTRF3*(HWSQF3/(VF2+BETAF3*VF4+HWSQF3)))
      FWV24=XLGR(XF,X)*1.E-20*
     &     ((W_wv)*RFT*((Xn-Xn_WV)/XN0))*FSCAL*XFRG  
      RETURN
      END 





      FUNCTION SWV24(WN,T,T0,W_wv,RFT,Xn,Xn_WV,XN0,XSLF)
      REAL*8 WN,X(4)
      COMMON /SH2O/ V1,V2,DV,NPTSLFWV,SWV296(2003)!---UNITS(CM**3/MOL)*1.E-20   
      COMMON /S260/ V1_,V2_,DV_,NPTSLFWV_,SWV260(2003)
      DATA V0S1,HWSQ1,BETAS1 /0.,10000.,1.E-04/
      DATA FACTRS1,V0S2,HWSQ2 /0.688,1050.,40000/
      DATA FACTRS2,V0S3,HWSQ3 /-0.2333,1310.,14400./
      DATA BETAS3,FACTRS3 /5.E-06,-0.15/
      J=INT((WN-V1)/DV)+1
      X(1)=((SWV296(J-1))*(SWV260(J-1)/SWV296(J-1))**((T-T0)/(260.-T0)))
      X(2)=((SWV296(J))*(SWV260(J)/SWV296(J))**((T-T0)/(260.-T0)))
      X(3)=((SWV296(J+1))*(SWV260(J+1)/SWV296(J+1))**((T-T0)/(260.-T0)))
      X(4)=((SWV296(J+2))*(SWV260(J+2)/SWV296(J+2))**((T-T0)/(260.-T0)))
      XF     = ( WN- (V1 + DV * FLOAT(J-1)) ) / DV
      SFAC = 1.
      VS2 = (WN-V0S1)**2
      VS4 = VS2*VS2
      SFAC = SFAC * 
     *     (1.+FACTRS1*(HWSQ1/(WN**2+(BETAS1*VS4)+HWSQ1)))  
      VS2 = (WN-V0S2)**2
      SFAC = SFAC *
     *     (1.+FACTRS2*(HWSQ2/(VS2+HWSQ2)))
      VS2 = (WN-V0S3)**2
      VS4 = VS2*VS2
      SFAC = SFAC *
     *     (1.+FACTRS3*(HWSQ3/(VS2+(BETAS3*VS4)+HWSQ3))) 
      SWV24=(W_wv)*RFT*(Xn_WV/XN0)*XLGR(XF,X)*1.E-20*SFAC*XSLF
      RETURN
      END 




      FUNCTION CONTI_N2(WN,T,T0,W_n2,RFT,RHOFAC,XCN2)
      REAL*8 WN,X(4)
      COMMON /N2RT0/ V1,V2,DV,NPTCONTN2,CT296(73)
      COMMON /N2RT1/ V1_,V2_,DV_,NPTCONTN2_,CT220(73)
      IF (W_n2.EQ.0.) THEN
         CONTI_N2=0.
         RETURN
      ENDIF
      J=INT((WN-V1)/DV)+1
      X(1)=(CT296(J-1))*(CT296(J-1)/CT220(J-1))**((T-T0)/(220.-T0))
      X(2)=(CT296(J))*(CT296(J)/CT220(J))**((T-T0)/(220.-T0))
      X(3)=(CT296(J+1))*(CT296(J+1)/CT220(J+1))**((T-T0)/(220.-T0))
      X(4)=(CT296(J+2))*(CT296(J+2)/CT220(J+2))**((T-T0)/(220.-T0))
      XF     = ( WN- (V1 + DV * FLOAT(J-1)) ) / DV
      CONTI_N2=XLGR(XF,X)*1.E-20*((W_n2/0.268675)*RFT*(RHOFAC))*XCN2     
      RETURN
      END 

      FUNCTION XLGR(XF,X) !4 points Lagrange interpolation
      REAL*8 A(4),X(4)    !with continous derivatives
      B=0.5*XF*(1.-XF)
      A(1)   = -B*(1.-XF)
      A(2)   = 1-((3.-2.*XF)*XF*XF)+B*XF
      A(3)   = ((3.-2.*XF)*XF*XF)+B*(1.-XF)
      A(4)   = -(B*XF)
      XLGR  = A(1)*X(1)+A(2)*X(2)+A(3)*X(3)+A(4)*X(4)
      END

      SUBROUTINE LINES(Xn,WN,T,W_wv,W_o2,W_o3,W_n2,
     &     W_n2O,W_co,W_so2,W_no2,W_oh,RADCT,T0,
     &     OL_WV,OL_O2,OL_O3,OL_N2,OL_N2O,
     &     OL_CO,OL_SO2,OL_NO2,OL_OH,
     &     XN0,RFT,P,P0,SCLCPL,SCLHW,Y0RES)
      PARAMETER (NNM=   9,IIM=   5000)
      REAL*8 WN,XNU,XNU0(NNM,IIM)
      REAL DELTNU(NNM,IIM),E(NNM,IIM),ALPS(NNM,IIM),XG(NNM,IIM)
      REAL ALPF(NNM,IIM),X(NNM,IIM),A(4),S0(NNM,IIM),B(4),TEMPLC(4)
      INTEGER MOL(NNM),NBLM(NNM),ISO(NNM,IIM),ICF(22)
      COMMON/HITR/MOL,NBLM,ISO,XNU0,DELTNU,S0,E,ALPS,ALPF,X,XG,NMOLEC
      DATA TEMPLC /200.0,250.0,296.0,340.0 /
      DATA ICF /1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/
      DATA MOL_WV/1/,MOL_O3/3/,MOL_O2/7/,MOL_N2/22/,MOL_N2O/4/
      DATA MOL_CO/5/,MOL_SO2/9/,MOL_NO2/10/,MOL_OH/13/
      deltnuC=25.          !cm-1
      WTOT=W_wv+W_o2+W_o3+W_n2+W_n2O+W_co+W_so2+W_no2+W_oh
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
      OL_WV  = 0.               !initialization
      OL_O2  = 0.               !initialization
      OL_O3  = 0.               !initialization
      OL_N2  = 0.               !initialization
      OL_N2O = 0.               !initialization
      OL_CO  = 0.               !initialization
      OL_SO2 = 0.               !initialization
      OL_NO2 = 0.               !initialization
      OL_OH  = 0.               !initialization
      DO I=1,NMOLEC
         IF (MOL(I).EQ.MOL_WV)  W_SPECIES=W_WV
         IF (MOL(I).EQ.MOL_O2)  W_SPECIES=W_O2
         IF (MOL(I).EQ.MOL_O3)  W_SPECIES=W_O3
         IF (MOL(I).EQ.MOL_N2)  W_SPECIES=W_N2
         IF (MOL(I).EQ.MOL_N2O) W_SPECIES=W_N2O
         IF (MOL(I).EQ.MOL_CO)  W_SPECIES=W_CO
         IF (MOL(I).EQ.MOL_SO2) W_SPECIES=W_SO2
         IF (MOL(I).EQ.MOL_NO2) W_SPECIES=W_NO2
         IF (MOL(I).EQ.MOL_OH)  W_SPECIES=W_OH
         IF (W_SPECIES.EQ.0.) THEN
            OL=0.
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
            IF ((ABS(WN-Xnu).GT.deltnuC).and.(MOL(I).NE.MOL_O2)) GOTO 30   
            CALL QOFT (MOL(I),ISO(I,J),296.,QT_296) 
            CALL QOFT (MOL(I),ISO(I,J),T,QT) 
            XIPSF = QT_296/QT
            CALL INTENS(T,S0(I,J),E(I,J),RADCT,T0,Xnu,STILD,XIPSF)
            XTILD=1-X(I,J)
            HWHM_C=HALFWHM_C(alpf(I,J),alps(I,J),RT,XTILD,RN,MOL(I),RAT)
            HWHM_D=HALFWHM_D(MOL(I),ISO(I,J),Xnu,T)
            IF(XG(I,J).EQ.-3.) THEN
               HWHM_C=HWHM_C*(1-(AIP*(RP))-(BIP*(RP2)))
            ENDIF
            CALL LSF_VOIGT(XG(I,J),RP,RP2,AIP,BIP,HWHM_C,WN,
     &           Xnu,ICF(MOL(I)),SLS,HWHM_D,MOL(I))
            SF=SF+(STILD*SLS)
 30         CONTINUE
            J=JJ
         END DO
         SPSD=W_species*SF
         OL=RFT*SPSD
 10      IF (MOL(I).EQ.MOL_WV)  OL_WV  = OL
         IF (MOL(I).EQ.MOL_O2)  OL_O2  = OL
         IF (MOL(I).EQ.MOL_O3)  OL_O3  = OL
         IF (MOL(I).EQ.MOL_N2)  OL_N2  = OL
         IF (MOL(I).EQ.MOL_N2O) OL_N2O = OL
         IF (MOL(I).EQ.MOL_CO)  OL_CO  = OL
         IF (MOL(I).EQ.MOL_SO2) OL_SO2 = OL
         IF (MOL(I).EQ.MOL_NO2) OL_NO2 = OL
         IF (MOL(I).EQ.MOL_OH)  OL_OH  = OL
      ENDDO
      END 

      FUNCTION HALFWHM_D(MOL,ISO,XNU,T)
      PARAMETER (NMOL=36,Nspeci=85)
      REAL C,K,T,M,AVOG,M
      REAL*8 XNU  
      INTEGER ILOC,ISO
      COMMON /ISVECT/ ISOVEC(NMOL),ISO82(Nspeci),ISONM(NMOL),
     *     smassi(Nspeci)
      C=2.997925E10      !light speed
      K=1.380662E-16     !Boltzman constant
      AVOG= 6.022045E23  !Avogadro number
      ILOC = ISOVEC(MOL)+ISO                                       
      M=SMASSI(ILOC)
      HALFWHM_D=(XNU/C)*SQRT(2.*ALOG(2.)*((k*T)/(M/AVOG)))
      END


      SUBROUTINE LSF_VOIGT(XF,RP,RP2,AIP,BIP,HWHM,WN,Xnu,ICF,SLS,AD,MOL)
      REAL*8 WN,Xnu,deltXNU,deltnuC
      DATA MOL_WV/1/,MOL_O3/3/,MOL_O2/7/,MOL_N2/22/,MOL_N2O/4/
      deltnuC=25.          !cm-1
      DIFF=(WN+Xnu)-deltnuC
      IF (ICF .EQ. 1) THEN        !continuum modeled:substract the 25cm-1
         deltXNU=(WN-Xnu)
         XL1=VOIGT(deltXNU,HWHM,AD)  !VOIGT for (+) osc.
         XL3=VOIGT(deltnuC,HWHM,AD)  !VOIGT for 25cm-1 wn
         IF (DIFF .LE. 0.) THEN
            deltXNU=(WN+Xnu)
            XL2=VOIGT(deltXNU,HWHM,AD) !VOIGT for (-) osc.
            SLS = (XL1 + XL2 - (2 * XL3)) 
         ELSE
            SLS = (XL1 - XL3) 
         ENDIF
      ENDIF
      IF (ICF .EQ. 0) THEN      !continuum not modeled
         IF (MOL.NE.MOL_O2) THEN !all cases except  (O2) 
            deltXNU=(WN-Xnu)
            XL1=VOIGT(deltXNU,HWHM,AD) !VOIGT for (+) osc.
            IF (DIFF .LE. 0.) THEN
               deltXNU=(WN+Xnu)
               XL2=VOIGT(deltXNU,HWHM,AD) !VOIGT for (-) osc.
               IF (XF.EQ.-1) THEN
                  Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
                  Y2=(1.-(AIP*(1/HWHM)*RP*(WN+Xnu))+(BIP*RP2))
                  SLS = (XL1*(Y1)+XL2*(Y2))
               ELSE
                  SLS = (XL1+XL2)
               ENDIF
            ELSE
               IF (XF.EQ.-1) THEN
                  Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
                  SLS = (XL1*(Y1)) 
               ELSE
                  SLS = (XL1)
               ENDIF
            ENDIF
         ENDIF
         IF (MOL.EQ.MOL_O2) THEN !case of the O2 
            SLS=0.
            IF ((ABS(WN-Xnu).LE.deltnuC).and.(XF.NE.-1).and.
     &           (XF.NE.-3)) THEN    !no O2 line coupling 
               deltXNU=(WN-Xnu)
               XL1=VOIGT(deltXNU,HWHM,AD) !VOIGT for (+) osc.
               IF (DIFF .LE. 0.) THEN
                  deltXNU=(WN+Xnu)
                  XL2=VOIGT(deltXNU,HWHM,AD) !VOIGT for (-) osc.
                  SLS = (XL1+XL2)
               ELSE
                  SLS = (XL1)
               ENDIF
            ENDIF
            IF ((XF.EQ.-1).or.(XF.EQ.-3)) THEN !O2 line coupling
               deltXNU=(WN-Xnu)
               XL1=VOIGT(deltXNU,HWHM,AD) !VOIGT for (+) osc.
               deltXNU=(WN+Xnu)
               XL2=VOIGT(deltXNU,HWHM,AD) !VOIGT for (-) osc.
               Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
               Y2=(1.-(AIP*(1/HWHM)*RP*(WN+Xnu))+(BIP*RP2))
               SLS = (XL1*(Y1)+XL2*(Y2))
            ENDIF
         ENDIF
      ENDIF
      END 

      
      SUBROUTINE LSF_LORTZ(XF,RP,RP2,AIP,BIP,HWHM,WN,Xnu,ICF,SLS)
      REAL*8 WN,Xnu,deltnuC
      deltnuC=25.          !cm-1
      DIFF=(WN+Xnu)-deltnuC
      IF (ICF .EQ. 1) THEN      !continuum modeled:substract the 25cm-1
         XL1=XLORENTZ((WN-Xnu)/HWHM) !Lorentz for (+) osc.
         XL3=XLORENTZ(deltnuC/HWHM)  !Lorentz for 25cm-1 wn
         IF (DIFF .LE. 0.) THEN
            XL2=XLORENTZ((WN+Xnu)/HWHM) !Lorentz for (-) osc.
            SLS = (XL1 + XL2 - (2 * XL3)) / HWHM
         ELSE
            SLS = (XL1 - XL3) / HWHM
         ENDIF
      ENDIF
      IF (ICF .EQ. 0) THEN      !continuum not modeled
         XL1=XLORENTZ((WN-Xnu)/HWHM) !Lorentz for (+) osc.
         IF (DIFF .LE. 0.) THEN
            XL2=XLORENTZ((WN+Xnu)/HWHM) !Lorentz for (-) osc.
            IF (XF.EQ.-1) THEN
               Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
               Y2=(1.-(AIP*(1/HWHM)*RP*(WN+Xnu))+(BIP*RP2))
               SLS = (XL1*(Y1)+XL2*(Y2))/HWHM
            ELSE
               SLS = (XL1+XL2)/HWHM
            ENDIF
         ELSE
            IF (XF.EQ.-1) THEN
               Y1=(1.+(AIP*(1/HWHM)*RP*(WN-Xnu))+(BIP*RP2))
               SLS = (XL1*(Y1)) / HWHM
            ELSE
               SLS = (XL1) / HWHM
            ENDIF
         ENDIF
      ENDIF
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

      SUBROUTINE INITI(P,T,RADCT,T0,P0,W_wv,W_o2,
     &     W_o3,W_n2,W_n2O,W_co,W_so2,W_no2,W_oh,
     &     XN0,Xn,Xn_WV,RHOFAC)
      DATA OXYGEN /2.090E+05/
      DATA AN2    /7.81E+05/
      DATA WVMOLMASS /18.016 /, DRYMOLMASS/28.97/   
      RADCT=1.4387865           !in K/cm-1
      T0=296.                   !in K
      P0=1013.                  !in HPa
      XN0=(P0/(1.380662E-16*T0))*1.E+3
      XN=(P/(1.380662E-16*T))*1.E+3
      WDRY=W_o2+W_o3+W_n2+W_n2O+W_co+W_so2+W_no2+W_oh
      RATIOMIX=(W_wv*WVMOLMASS)/(WDRY*DRYMOLMASS)
      WVPRESS=(RATIOMIX/(RATIOMIX+0.622))*P
      Xn_WV=(WVPRESS/(1.380662E-16*T))*1.E+3
      RHOFAC=(W_n2/(WDRY+W_wv))*(P/P0)*(273./T)
      END 
      
      SUBROUTINE QofT(Mol, Iso, Tout, QT) !version of October 28th 1998
      PARAMETER (NMOL=36,Nspeci=85)       !OPTIMIZED ON 10/28/1998 S.A.BOUKABARA
      COMMON /ISVECT/ ISOVEC(NMOL),ISO82(Nspeci),ISONM(NMOL),
     *     sdum(Nspeci)
      common/Qtot/ Qcoef(Nspeci,3,5), Q296(Nspeci), aQ(Nspeci), 
     + bQ(Nspeci), gj(Nspeci)
      ivec = isovec(Mol) + iso
      irange = 1
      QT = (((Qcoef(ivec,irange,5)*Tout+Qcoef(ivec,irange,4))*Tout+
     &     Qcoef(ivec,irange,3))*Tout+Qcoef(ivec,irange,2))*Tout+
     &     Qcoef(ivec,irange,1)
      RETURN
      END

      SUBROUTINE vecIso          !version got on october 28th 1998
      PARAMETER (NMOL=36,Nspeci=85)
      COMMON /ISVECT/ ISOVEC(NMOL),ISO82(Nspeci),ISONM(NMOL),
     *     sdum(Nspeci)
         ISOVEC(1) = 0
         DO 20 I = 2,NMOL
          ISOVEC(I) = 0
          DO 10 J = 1,I-1
           ISOVEC(I) = ISOVEC(I)+ISONM(J)
   10     CONTINUE
   20    CONTINUE
      RETURN
      END
        
      SUBROUTINE READ_HITR(ICP,HFILE)
      PARAMETER (NNM=   9,IIM=    5000)
      REAL*8 XNU0(NNM,IIM),nu0
      REAL DELTNU(NNM,IIM),E(NNM,IIM),ALPS(NNM,IIM),ALPF(NNM,IIM)
      REAL X(NNM,IIM),XG(NNM,IIM),S0(NNM,IIM)
      INTEGER MOL(NNM),NBLM(NNM),ISO(NNM,IIM)
      Character Q1*9,Q2*9,HFILE*80
      COMMON/HITR/MOL,NBLM,ISO,XNU0,DELTNU,S0,E,ALPS,ALPF,X,XG,NMOLEC
      DATA MOL_WV/1/,MOL_O3/3/,MOL_O2/7/,MOL_N2/22/,MOL_N2O/4/
      DATA MOL_CO/5/,MOL_SO2/9/,MOL_NO2/10/,MOL_OH/13/
      I_WV=0
      I_O3=0
      I_O2=0
      I_N2=0
      I_N2O=0
      I_CO=0
      I_SO2=0
      I_NO2=0
      I_OH=0
      NMOLEC=0
      OPEN(9,FILE=HFILE,form='formatted',status='old',ERR=20)
      ILINE=0
 22   READ(9,100,END=3,ERR=30) mo,iso_scal,nu0,S0_scal,R,
     &     alpf_scal,alps_scal,E_scal,X_scal,deltnu_scal,
     &     iv1,iv2,Q1,Q2,ier1,ier2,ier3,iref1,iref2,iref3
      ILINE=ILINE+1
      IF (ILINE.EQ.1) MINWN=nu0
      CALL S_INDX(mo,mol_wv,I_WV,NMOLEC,NMOL_WV,MOL,II,JJ)
      CALL S_INDX(mo,mol_o2,I_O2,NMOLEC,NMOL_O2,MOL,II,JJ)
      CALL S_INDX(mo,mol_o3,I_O3,NMOLEC,NMOL_O3,MOL,II,JJ)
      CALL S_INDX(mo,mol_n2,I_N2,NMOLEC,NMOL_N2,MOL,II,JJ)
      CALL S_INDX(mo,mol_n2o,I_N2O,NMOLEC,NMOL_N2O,MOL,II,JJ)
      CALL S_INDX(mo,mol_co,I_CO,NMOLEC,NMOL_CO,MOL,II,JJ)
      CALL S_INDX(mo,mol_so2,I_SO2,NMOLEC,NMOL_SO2,MOL,II,JJ)
      CALL S_INDX(mo,mol_no2,I_NO2,NMOLEC,NMOL_NO2,MOL,II,JJ)
      CALL S_INDX(mo,mol_oh,I_OH,NMOLEC,NMOL_OH,MOL,II,JJ)
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
         CALL COUPL_INDX(mo,mol_wv,I_WV,II)
         CALL COUPL_INDX(mo,mol_o2,I_O2,II)
         CALL COUPL_INDX(mo,mol_o3,I_O3,II)
         CALL COUPL_INDX(mo,mol_n2,I_N2,II)
         CALL COUPL_INDX(mo,mol_n2o,I_N2O,II)
         CALL COUPL_INDX(mo,mol_co,I_CO,II)
         CALL COUPL_INDX(mo,mol_so2,I_SO2,II)
         CALL COUPL_INDX(mo,mol_no2,I_NO2,II)
         CALL COUPL_INDX(mo,mol_oh,I_OH,II)
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
      WRITE(*,*) '****************************************'
      WRITE(*,*) '* SPECTRAL LINES INFORMATION AVAILABLE *'
      WRITE(*,*) '****************************************'
      WRITE(*,*) 'Minimum Wavenumber:',MINWN,' cm-1'
      WRITE(*,*) 'Maximum Wavenumber:',MAXWN,' cm-1'
      WRITE(*,*)
      WRITE(*,'(2x,a8,8x,a8,3x,a8)') '','Molecule','# lines' 
      DO J=1,NMOLEC
         IF (MOL(J).EQ.MOL_WV)  NBLM(J)=I_WV
         IF (MOL(J).EQ.MOL_O2)  NBLM(J)=I_O2
         IF (MOL(J).EQ.MOL_O3)  NBLM(J)=I_O3
         IF (MOL(J).EQ.MOL_N2)  NBLM(J)=I_N2
         IF (MOL(J).EQ.MOL_N2O) NBLM(J)=I_N2O
         IF (MOL(J).EQ.MOL_CO)  NBLM(J)=I_CO
         IF (MOL(J).EQ.MOL_SO2) NBLM(J)=I_SO2
         IF (MOL(J).EQ.MOL_NO2) NBLM(J)=I_NO2
         IF (MOL(J).EQ.MOL_OH)  NBLM(J)=I_OH
         WRITE(*,'(2x,i8,5x,i8,5x,i8)') J,MOL(J),NBLM(J)
      ENDDO
      WRITE(*,*)
      WRITE(*,*) '****************************************'
      CLOSE(9)
 100  FORMAT (I2,I1,F12.6,1P,2E10.3,0P,2F5.4,F10.4,F4.2,F8.6,
     &     2I3,2A9,3I1,3I2)
 2    FORMAT (I2,1P,4(E13.6,E11.4),0P,I2)                           
      RETURN
 20   PRINT *, 'ERROR OPENING HITRAN FILE:',HFILE
      STOP
 30   PRINT *, 'ERROR READING HITRAN FILE:',HFILE
      STOP
      END   

      SUBROUTINE S_INDX(mo,molSP,I_SP,NMOLEC,NMOL_SP,MOL,II,JJ)
      DIMENSION MOL(*)
      IF(mo.eq.molsp) THEN
         I_SP=I_SP+1
         IF (I_SP.EQ.1) THEN
            NMOLEC=NMOLEC+1
            NMOL_SP=NMOLEC
            MOL(NMOLEC)=molSP
         ENDIF
         II=I_SP
         JJ=NMOL_SP
      ENDIF
      END

      SUBROUTINE COUPL_INDX(mo,mol_sp,I_SP,II)
      IF(mo.eq.mol_sp) THEN
         I_SP=I_SP+1
         II=I_SP
      ENDIF
      END
      INCLUDE 'isotope.dat'     !INCLUDE HITRSID DATA-STATEMENTS
      INCLUDE 'continuumDATASTAT.dat'!INCLUDE CONTINUUM DATA-STATEMENTS

      FUNCTION ODCLW(WN,TEMP,CLW)
      !INPUTS: WN (WaveNUmber in cm-1)
      !        Temp (in K)
      !        CLW  (in mm or kg/m2)
      !OUTPUT: ODCLW: optical depth of the Cloud Liquid Water
      !FROM EQUATIONS OF LIEBE, HUFFORD AND MANABE
      !Ref:(INT. J. IR & MM WAVES V.12(17) JULY 1991
      COMPLEX EPS,RE
      REAL*8 WN
      FREQ=WN*29.98
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
      ODCLW = -.06286057*CLW*AIMAG(RE)*FREQ
      RETURN
      END


      FUNCTION XLORENTZ(Z)
      REAL*8 Z
      REAL XLORENTZ
         XLORENTZ=1./(3.14159*(1.+(Z**2))) 
      RETURN
      END



      FUNCTION VOIGT(DELTNU,ALPHAL,ALPHAD)
      IMPLICIT NONE
      COMPLEX v,w4
      REAL*8 DELTNU
      REAL avc(0:101)
      REAL ALPHAL,ALPHAD,ZETA,ALPHAV,AVCINTERP,DNU
      REAL AL,AD,X,Y,ANORM1,VOIGT,DZETA
      INTEGER IZETA2,IZETA1
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
         VOIGT=(ALPHAL/(3.14159*(ALPHAL**2+(DELTNU)**2)))
         RETURN
      end if
      !---SETUP PARAMETERS FOR CALL TO VOIGT GENERATOR (HUMLICEK)
      x = sqrt(alog(2.))*(dnu)
      y = 1000.
      if (zeta .lt. 1.000) then
         y = sqrt(alog(2.))*AL
      end if
      !---CALL the Humlicek subroutine
      v=W4(x,y)
      anorm1 = sqrt(alog(2.)/3.14159)/alphad
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
 2    IF(Y.LT..195*ABS(X)-.176)GOTO 3
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
      
