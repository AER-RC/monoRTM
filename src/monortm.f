	PROGRAM MONORTM
	!IMPLICIT NONE
	INTEGER NWN,I,ICPL,IS,IOUT,IRT,J,MXLAY
	INTEGER NWNMX,IATM,INP,IVC
	PARAMETER (NWNMX=60000,MXLAY=203,NPROFMX=5000)
	REAL*8 WN(NWNMX),WNMIN,WNMAX,V1,V2,wnlblrtm
	REAL W_wv(MXLAY),W_o2(MXLAY),W_o3(MXLAY),W_n2(MXLAY),
     &     W_n2O(MXLAY),W_co(MXLAY),W_so2(MXLAY),
     &     W_no2(MXLAY),W_oh(MXLAY),CLW(MXLAY)
	REAL O(NWNMX,MXLAY),
     &     OL_WV(NWNMX,MXLAY),OS_WV(NWNMX,MXLAY),
     &     OF_WV(NWNMX,MXLAY),OL_O2(NWNMX,MXLAY),
     &     OL_O3(NWNMX,MXLAY),OL_N2(NWNMX,MXLAY),
     &     OC_N2(NWNMX,MXLAY),OL_N2O(NWNMX,MXLAY),
     &     OL_CO(NWNMX,MXLAY),OL_SO2(NWNMX,MXLAY),
     &     OL_NO2(NWNMX,MXLAY),OL_OH(NWNMX,MXLAY),
     &     O_CLW(NWNMX,MXLAY),W_wv0(MXLAY)
	REAL TMPSFC,RAD(NWNMX),EMISS(NWNMX),REFLC(NWNMX)
	REAL RUP(NWNMX),TRTOT(NWNMX),RDN(NWNMX),TB(NWNMX)
	character*4 ht1,ht2
	CHARACTER fic0*64,fileARMlist*64,filearm*80,fileprof*80
	CHARACTER*80 filearmTAB(NPROFMX)
	character*8 XID,HMOLID,YID,HDATE,HTIME         
	real*8      SECANT,XALTZ 
	COMMON /PATHD/ PP(MXLAY),TT(MXLAY),WKL(35,MXLAY),
	1    WBRODL(MXLAY),DVL(MXLAY),WTOTL(MXLAY),ALBL(MXLAY),
	2    ADBL(MXLAY),AVBL(MXLAY),H2OSL(MXLAY),IPTH(MXLAY),
	3    ITYL(MXLAY),SECNTA(MXLAY),HT1,HT2,ALTZ(0:MXLAY),
	4    PZ(0:MXLAY),TZ(0:MXLAY) 
	COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,  
	1    AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      
	2    DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,      
	3    ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,     
	4    EXTID(10)    
	COMMON /BNDPRP/ TMPSFC,BNDEMI(3),BNDRFL(3),IBPROP          
	COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       
	1    WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   
	2    EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    


	!---INPUTS & GENERAL CHARACTERISTICS
	IVC=2    !if=2->CKD2.4  MPMf87/s93 (if=3)
	ICPL=1   !=1->cpl =0->nocpl
	IDU=1    !IDU=0->the prof are from top to surface, IDU=1->the opposite
	INP=3    !INP=1->the inputs read by LBLATM input (LBLRTM.IN)
		 !INP=2->ARM inputs transformed to LBLATM input
		 !INP=3->ARM inputs transformed to layer inputs (these are in
	         !       MONORTM_PROF.IN from createMONOinputfromARMfiles.f
	XSLF=1.  !scaling factor (SLF cont)
	XFRG=1.  !scaling factor (FRG cont)
	XCN2=1.  !scaling factor (N2 cont)
	SCLCPL=1.!scaling factor (Line Coupling Parameters)
	SCLHW=1. !scaling factor (Pressure Dependence of the halfwidth of the 0 band)
	Y0RES=0. !Y0RES of the line coupling coeffs
	SCALWV=1.!scaling of the WV profile
	NPROF=1  !number of profiles (by default, LBLRTM.IN contains 1)
	IOUT=1   !=0->radiances / =1->TB+RAD
	CLOUD=0. !cloud liquid water

	IF (INP.EQ.2) THEN
	   fileARMlist='ARM.IN'
	   open(23,file=fileARMlist,status='old',form='formatted')
	   NPROF=0
	   DO WHILE (.true.)
	      read(23,'(a)',end=70) filearm
	      NPROF=NPROF+1
	      filearmTAB(NPROF)=filearm
	   ENDDO
 70	   close(23)
	   IF (NPROF.EQ.0) THEN
	      PRINT *, ' NO FILES FOUND IN THE ARM FILES LIST'
	      STOP
	   ENDIF
	ENDIF
	IF (INP.EQ.3) THEN
	   fileprof='in/MONORTM_PROF.IN'
	   !fileprof='MONORTM_PROF.IN_TOT_RHCORR_31LEV'
	   !fileprof='MONORTM_PROF.IN_TOT_RHORIG_31LEV'
	   open(53,file=fileprof,status='old',form='formatted')
	   read(53,12) NPROF
	ENDIF
c---------------------
c
c	NWN=4
c	WN(1)=0.789344
c	WN(2)=0.79828
c	WN(3)=1.043027
c	WN(4)=1.051763
c---------------------
	NWN=1000
	!WNMIN=3.33566
	!WNMAX=6.67111
	WNMIN=0.1/30.
	WNMAX=300./30.
	DELTAWN=(WNMAX-WNMIN)/(NWN-1)
	DO I=1,NWN
	   EMISS(I)=0.6
	   REFLC(I)=0.4
	   WN(I)=WNMIN+(I-1)*DELTAWN
	ENDDO

	NREC=0

	DO 111 NPR=1,NPROF
	   IF (INP.EQ.2) THEN
	      CALL ARM2LBLATM(filearmTAB(NPR),IFLAG,ilaunchdate,
	1	   ilaunchtime,ibasetime,iserialnumber,isondeage)
	      IF (IFLAG.EQ.2) GOTO 111
	   ENDIF
	   IF (INP.EQ.1.OR.INP.EQ.2) CALL RDLBLINP(IATM,IOUT,IRT)
	   IF (INP.EQ.3) THEN
	      READ(53,13,END=110) NPR1,NLAYRS,TMPSFC,IRT
	      READ(53,15,END=110) ilaunchdate,ilaunchtime,
	1	   ibasetime,iserialnumber,isondeage
	      DO IL=1,NLAYRS
		 IF (IL.EQ.1) THEN                                                 
		    READ(53,910,END=110) PP(IL),TT(IL),SECNTA(IL),  
	1		 PZ(IL-1),TZ(IL-1),PZ(IL),TZ(IL) 
		 ELSE                                                             
		    READ(53,915,END=110) PP(IL),TT(IL),SECNTA(IL),
	1		 PZ(IL),TZ(IL)  
		 ENDIF                                                            
		 READ(53,14,END=110) WKL(1,IL),WKL(3,IL),WKL(4,IL),
	1	      WKL(7,IL),WBRODL(IL),WKL(5,IL),WKL(9,IL),
	2	      WKL(10,IL),WKL(13,IL),WKL(22,IL)
	      ENDDO   
	   ENDIF
	   DO IL=1,NLAYRS
	      W_wv0(IL)=WKL(1,IL)
	      W_o2(IL)=WKL(7,IL)
	      !W_n2(IL)=WBRODL(IL)
	      W_n2(IL)=WKL(22,IL)

	      !W_o3(IL)=WKL(3,IL)
	      !W_n2O(IL)=WKL(4,IL)
	      !W_co(IL)=WKL(5,IL)
	      !W_so2(IL)=WKL(9,IL)
	      !W_no2(IL)=WKL(10,IL)
	      !W_oh(IL)=WKL(13,IL)

	      W_o3(IL)=0.
	      W_n2O(IL)=0.
	      W_co(IL)=0.
	      W_so2(IL)=0.
	      W_no2(IL)=0.
	      W_oh(IL)=0.
	   ENDDO

	   !---GET THE INSTRUMENTAL CONFIGURATION
	   !IF (INP.EQ.1.OR.INP.EQ.2) CALL INSTRCONF(V1,V2,DVSET,NWN,WN)

	   !---SET-UP THE EMISS/REFLEC VECTORS
	   !CALL EMISS_REFLEC(NWN,EMISS,REFLC,WN)

	   NREC=NREC+1
				!---SCALE THE WATER VAPOR PROFILE
	   DO L=1,NLAYRS
	      W_wv(L)=W_wv0(L)*SCALWV
	   ENDDO
				!---COLUMN WATER VAPOR/LIQUID WATER
	   CLW(1)=cloud
	   CALL INTEGR(W_wv,CLW,NLAYRS,WVCOLMN,CLWCOLMN)
				!---OPTICAL DEPTHS
           CALL MODM(IVC,ICPL,NWN,WN,NLAYRS,PP,TT,W_wv,
	1	W_o2,W_o3,W_n2,W_n2O,W_co,W_so2,W_no2,
	2	W_oh,CLW,O,OL_WV,OS_WV,OF_WV,OL_O2,
	3	OL_O3,OL_N2,OC_N2,OL_N2O,OL_CO,OL_SO2,
	4	OL_NO2,OL_OH,O_CLW,XSLF,XFRG,XCN2,
	5	SCLCPL,SCLHW,Y0RES)
				!---CORRECT THE OPTDEPTHS FOR THE SLANT PATH
	   CALL CORR_OPTDEPTH(INP,NLAYRS,SECNTA,NWN,ANGLE,O)
				!---RADIATIVE TRANSFER
	   CALL RTM(IOUT,IRT,NWN,WN,NLAYRS,TT,TZ,
	1	TMPSFC,O,RUP,TRTOT,RDN,REFLC,EMISS,RAD,TB,IDU)
				!---WRITE OUT THE RESULTS
	   CALL STOREOUT(NWN,WN,RAD,TB,TRTOT,SCLCPL,SCLHW,NREC,
	1	WVCOLMN,0,XSLF,0,Y0RES,0,INP,NPR,ilaunchdate,
	1	ilaunchtime,ibasetime,iserialnumber,isondeage,
	2	CLWCOLMN,TMPSFC,REFLC,EMISS,O,OL_WV,OS_WV,OF_WV,
	3	OL_O2,OL_O3,OL_N2,OC_N2,OL_N2O,OL_CO,OL_SO2,
	4	OL_NO2,OL_OH,O_CLW,NLAYRS,PP)


 111	ENDDO
 110	CLOSE(1)
 910	FORMAT (E15.7,F10.4,F10.4,1X,2(F8.3,F7.2))
 915	FORMAT (E15.7,F10.4,F10.4,F8.3,F7.2)
 14	FORMAT (10E15.7)
 12	FORMAT (19x,I7)
 13	format(9x,I7,20x,I7,7x,f9.2,6x,I4)
 15	FORMAT (8x,i12,9x,i12,7x,i12,8x,i12,10x,i12)
	CLOSE(53)
	STOP
	END

