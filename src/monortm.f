	PROGRAM MONORTM
	!***************************************************************
	!        Monochromatic Radiative Transfer Model 
	!***************************************************************
	!  
	!   This code is a forward model. It generates the brightness
	!   temperatures by:
	!
	!   - calling MODM (Monochromatic Optical Depth Model) for 
	!     the optical depth computation and
	!   - running the radiative transfer code RTM to simulate  
	!     the radiances/brightness temperatures.
	!   - it optionally uses LBLATM to generate the internal inputs 
	!
	!                           MONORTM
	!                       _______|______
	!                      |       |      |
	!                      |       |      |
	!                     MODM   LBLATM  RTM
	!
	!   MONORTM:
	!   ********
	!   MonoRTM is designed to be very flexible. We can either use 
	!   it as a black box and control everything from the LBLRTM.IN
	!   input file, or one can modify the code itself and
	!   recompile it. In the latter, it is structured in such a way 
	!   that the changes should always be done in "monortm.f" (the  
	!   driver program). The other auxillary files (monortm_sub.f, 
	!   modm.f, lblatm.f and declar.incl, continuumDATASTAT.dat) 
	!   should normally not be touched, except in rare situations.
	!   In case there is an update in the continuum calculations
	!   a new file "continuumDATASTAT.dat" will be generated
	!   and sent to the users (or made available on the WEB/ftp
	!   site). In the same way, if the spectroscopic database is
	!   to be updated, only the file "spectral_lines.dat" file
	!   should be replaced. The code would not need to be 
	!   modified.
	!   
	!   MODM:
	!   *****
	!   The core of MonoRTM is the computation of the optical depths.
	!   It is designed as a subroutine for flexibility and could be
	!   easily plugged in a different radiative transfer model.
	!
	!
	!   Several options are possible in this code. We can:
	!   **************************************************  
	!
	!   - modify the version of the continuum model (see IVC param)
	!     IVC=1  -> not used 
	!     IVC=2  -> CKD version 2.4 (by default)
	!     IVC=3  -> Rosenkranz 1998's suggestions (MPMf87/s93)
	!
	!   - turn on/off the line coupling of the O2 lines (see ICPL)
	!     ICPL=1 -> Line coupling turned on  
	!     ICPL=0 -> Line coupling turned off (even if we have line 
        !               couplng coefficients in the spectroscopic 
	!               selected_lines)
	!   - assume that the atmospheric profiles are organized
	!     from the top to the surface (or the opposite), see IDU.
	!     IDU=0  -> profiles from top (1st elt) to surface (last elt)
	!     IDU=1  -> profiles from surface to top (default)
	!
	!   - read the inputs from different sources/formats. see INP
	!     INP=1  ->the inputs read by LBLATM input (LBLRTM.IN). This 
	!              input file is consistent with LBLRTM's TAPE5.
	!     INP=2  ->ARM files transformed to LBLATM input (needs ARM.IN 
	!              that will contain the list of ARM sondes files with 
	!              their absolute path). The sondes are quality-flagged
	!              See subroutine ARM2LBLATM for the flags.
	!              Only sondes with FLAG=0 are simulated.
	!     INP=3  ->inputs are read from a TAPE7-like file called
	!              MONORTM_PROF.IN. It contains the layers information
	!              WARNING: In this case, the surface temperature is 
	!              taken to be the surface level temperature.
	!     INP comes from the LBLRTM.IN file (see instructions)
	!   - scale/tune several parameters for the line coupling and
	!     the continuum computation, see :
	!     XSLF   : scaling factor for self continuum (default=1)
	!     XFRG   : scaling factor for foreign continuum (default=1)
	!     XCN2   : scaling factor for N2 continuum (default =1)
	!     SCLCPL : scaling factor for Line Coupling Parameters
	!     SCLHW  : scaling factor for the ressure Dependence of the 
	!              halfwidth of the O2 0-zero band
	!     Y0RES  : Y0RES of the line coupling coeffs.
	!
	!   - scale the water vapor profile. Useful when we try to scale
	!     the sondes profiles or when we perform the retrievals.
	!     SCALWV : factor by which the WV profile will be multiplied.
	!
	!   - output the simulations in radiances or in brightness
	!     temperatures. Note that the internal calculations are first
	!     done in radiances. IOUT comes from LBLRTM.IN
	!     IOUT   : =0 -> radiances 
	!     IOUT   : =1 -> TB+RAD
	!  
	!   History of the modifications:
	!   *****************************  
	!   - written in 1999 by Sid Ahmed Boukabara, Ross Hoffman
	!     and Tony Clough. 
	!   - validated against ARM sondes in the
	!     microwave spectrum (23.8 and 31.4 GHz). SAB, 2000.
	!   - extended to more species by Sid Ahmed Boukabara in 03/2000.
	!   - cleaned up and commented in 2001 for first public release.
	!     Also put under CVS configuration management. SAB.
	! 
	!   Comments or Questions:
	!   **********************
	!   Comments and/or questions are appreciated. 
	!   They should be forwarded to:
	!   AER Inc. (Sid Boukabara or Tony Clough)
        !   131 Hartwell Avenue, Lexington, MA02421, USA
	!   By phone  : 1 781 761 2213
	!   By Fax    : 1 781 761 2299
	!   By E-mail : sboukaba@aer.com, clough@aer.com
	!
	!
	!   Sid Ahmed Boukabara, AER Inc, April, 2001
	!
	!***************************************************************
	include "declar.incl"
	INTEGER NWN,I,ICPL,IS,IOUT,IRT,J,ICNTNM,IATM,INP,IVC
	REAL*8 V1,V2,SECANT,XALTZ 
	REAL TMPSFC
	CHARACTER fileARMlist*64,hmod*60,CTYPE*3
	CHARACTER fileprof*80,HFILE*80,FILEOUT*60,ht1*4,ht2*4
	CHARACTER*60 FILEIN,FILESONDE
	character*8 XID,HMOLID,YID,HDATE,HTIME
	CHARACTER*48 CFORM1,CFORM2                                         
	CHARACTER*4 PZFORM(5)
	CHARACTER*7 PAFORM(2)                                              
	DATA PAFORM / '1PE15.7','  G15.7'/                                 
	DATA PZFORM / 'F8.6','F8.5','F8.4','F8.3','F8.2'/                  
	DATA CFORM1 / '(1PE15.7,0PF10.2,10X,A3,I2,1X,2(F7.3,F8.3,F7.2))'/  
	DATA CFORM2 / '(  G15.7,0PF10.2,10X,A3,I2,23X,(F7.3,F8.3,F7.2))'/  
	COMMON /PATHD/ P,T,WKL,WBRODL,DVL,WTOTL,ALBL,ADBL,AVBL,
	2    H2OSL,IPTH,ITYL,SECNTA,HT1,HT2,ALTZ,PZ,TZ
	COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,  
	1    AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      
	2    DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,      
	3    ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,     
	4    EXTID(10)    
	COMMON /BNDPRP/ TMPSFC,BNDEMI(3),BNDRFL(3),IBPROP          
	COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       
	1    WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   
	2    EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    


	!---INPUTS & GENERAL CONTROL PARAMETERS
	IVC=2    !if=2->CKD2.4  MPMf87/s93 (if=3)
	ICPL=1   !=1->cpl =0->nocpl
	IDU=1    !IDU=0->the prof are from top to surface, IDU=1->the opposite

	XSLF=1.  !scaling factor (SLF cont)
	XFRG=1.  !scaling factor (FRG cont)
	XCN2=1.  !scaling factor (N2 cont)
	SCLCPL=1.!scaling factor (Line Coupling Parameters)
	SCLHW=1. !scaling factor (Pressure Dependence of the halfwidth of the 0 band)
	Y0RES=0. !Y0RES of the line coupling coeffs
	SCALWV=1.!scaling of the WV profile

	!---FILES NAMES ALL PUT HERE FOR CONVENIENCE
	FILEIN      ='in/LBLRTM.IN'
	FILESONDE   ='in/SONDE.IN'
	fileARMlist ='in/ARM.IN'
	fileprof    ='in/MONORTM_PROF.IN'
	HFILE       ='in/spectral_lines.dat'
	FILEOUT     ='out/MONORTM.OUT'


	!---GET THE WAVENUMBERS AND INP,IRT INFORMATION
	CALL RDLBLINP(0,IATM,IOUT,IRT,NWN,WN,FILEIN,ICNTNM,CLW,INP)                

	!---GET THE PROFILES NUMBER
	CALL GETPROFNUMBER(INP,FILEIN,fileARMlist,fileprof,
	1    NPROF,filearmTAB)

	!---PRINT OUT MONORTM VERSION AND PROFILES NUMBER
	CALL start(NPROF,INP)

	!---CHECK INPUTS AND THEIR CONSISTENCY WITH MONORTM
	CALL CHECKINPUTS(NWN,NPROF,NWNMX,NPROFMX)

	!---Loop over the number of profiles to be processed
	NREC=0
	DO 111 NPR=1,NPROF
	   !*********************************************
	   !* First Step: Read in the inputs
	   !*********************************************	   
	   !---Inputs: LBLRTM.IN (TAPE5-type of file)
	   IF (INP.EQ.1) THEN
	      CALL RDLBLINP(1,IATM,IOUT,IRT,NWN,WN,FILEIN,ICNTNM,CLW,INP)
	   ENDIF
	   !---Inputs: wavenumbers from LBLRTM.IN, profiles from ARM sondes
	   IF (INP.EQ.2) THEN
	      CALL ARM2LBLATM(filearmTAB(NPR),IFLAG,ilaunchdate, 
	1	   ilaunchtime,ibasetime,iserialnumber,isondeage,
	2	   NWN,WN,V1,V2,DVSET,FILESONDE)
	      IF (IFLAG.EQ.2) THEN
		 WRITE(*,'(a30,i5,a8)') 'PROCESSING PROFILE NUMBER:',NPR,' FLAG=2'
		 GOTO 111
	      ENDIF
	      CALL RDLBLINP(0,IATM,IOUT,IRT,NWN,WN,FILESONDE,ICNTNM,CLW,INP)
	   ENDIF
	   !---Inputs: MONORTM_PROF.IN (TAPE7 consistent)
	   IF (INP.EQ.3) THEN
	      READ (53,924,END=110,ERR=50) IFORM,NLAYRS,NMOL,SECNT0,HMOD
	      DO IL=1,NLAYRS
		 SECNTA(IL)=SECNT0
		 LTST = IL                                                     
		 IF (IL.EQ.1) LTST = 0                                         
		 PTST = ALOG10(PZ(LTST))                                      
		 NPTST = PTST+2                                               
		 IF (PTST.LT.0.0) NPTST = 1                                   
		 CFORM1(38:41) = PZFORM(NPTST)                                
		 CFORM2(38:41) = PZFORM(NPTST)                                
		 NPTST = 1                                                    
		 IF (P(IL).GE.0.1) NPTST = 2                                
		 CFORM1(2:8) = PAFORM(NPTST)                                  
		 CFORM2(2:8) = PAFORM(NPTST)                                  
		 IF (IL.EQ.1) THEN                                                 
		    READ (53,CFORM1,END=110,ERR=50) P(IL),T(IL),  
	1		 CTYPE,IPATH,ALTZ(IL-1),PZ(IL-1),        
	2		 TZ(IL-1),ALTZ(IL),  PZ(IL),  TZ(IL)  
		    TMPSFC=TZ(IL-1)
		    IF (IPATH.EQ.1) IRT=1 !space-based observer
		    IF (IPATH.EQ.3) IRT=3 !ground-based observer
		    IF (IPATH.EQ.2) IRT=2 !limb measurement
		 ELSE                                                             
		    READ (53,CFORM2,END=110,ERR=50) P(IL),T(IL),  
	1		 CTYPE,IPATH,ALTZ(IL),PZ(IL),TZ(IL)             
		 ENDIF                                                            
		 READ(53,978,END=110,ERR=50) (WKL(K,IL),K=1,7),WBRODL(IL)          
		 IF (NMOL.GT.7) READ(53,978) (WKL(K,IL),K=8,NMOL)
	      ENDDO 
	   ENDIF

	   !---PREPARE THE INPUTS FOR MODM
	   DO IL=1,NLAYRS
	      W_wv0(IL)=WKL(1,IL)
	      W_o2(IL)=WKL(7,IL)
	      W_n2(IL)=WBRODL(IL)
	      !W_n2(IL)=WKL(22,IL)
	      W_o3(IL)=WKL(3,IL)
	      W_n2O(IL)=WKL(4,IL)
	      W_co(IL)=WKL(5,IL)
	      W_so2(IL)=WKL(9,IL)
	      W_no2(IL)=WKL(10,IL)
	      W_oh(IL)=WKL(13,IL)
	   ENDDO

	   !***********************************************
	   !* Second Step: SET-UP THE EMISS/REFLEC VECTORS
	   !***********************************************	   
	   CALL EMISS_REFLEC(NWN,EMISS,REFLC,WN) 

	   NREC=NREC+1

	   !---SCALE THE WATER VAPOR PROFILE
	   DO L=1,NLAYRS
	      W_wv(L)=W_wv0(L)*SCALWV
	   ENDDO

	   !---COLUMN WATER VAPOR/LIQUID WATER
	   CALL INTEGR(W_wv,CLW,NLAYRS,WVCOLMN,CLWCOLMN)

	   !***********************************************
	   !* Third Step: OPTICAL DEPTHS COMPUTATION
	   !***********************************************	   
           CALL MODM(IVC,ICPL,NWN,WN,NLAYRS,P,T,W_wv,
	1	W_o2,W_o3,W_n2,W_n2O,W_co,W_so2,W_no2,
	2	W_oh,CLW,O,OL_WV,OS_WV,OF_WV,OL_O2,
	3	OL_O3,OL_N2,OC_N2,OL_N2O,OL_CO,OL_SO2,
	4	OL_NO2,OL_OH,O_CLW,XSLF,XFRG,XCN2,
	5	SCLCPL,SCLHW,Y0RES,HFILE,ICNTNM)

	   !***********************************************
	   !* Fourth Step: CORRECT FOR THE SLANT PAT
	   !***********************************************	   
	   !---CORRECT THE OPTDEPTHS FOR THE SLANT PATH
	   CALL CORR_OPTDEPTH(INP,NLAYRS,SECNTA,NWN,ANGLE,O,IRT)

	   !***********************************************
	   !* Fifth Step: RADIATIVE TRANSFER
	   !***********************************************	   
	   CALL RTM(IOUT,IRT,NWN,WN,NLAYRS,T,TZ,
	1	TMPSFC,O,RUP,TRTOT,RDN,REFLC,EMISS,RAD,TB,IDU)
				
	   !***********************************************
	   !* Sixth Step: WRITE OUT THE RESULTS
	   !***********************************************	   
	   CALL STOREOUT(NWN,WN,RAD,TB,TRTOT,SCLCPL,SCLHW,NREC,
	1	WVCOLMN,0,XSLF,0,Y0RES,0,INP,NPR,ilaunchdate,
	1	ilaunchtime,ibasetime,iserialnumber,isondeage,
	2	CLWCOLMN,TMPSFC,REFLC,EMISS,O,OL_WV,OS_WV,OF_WV,
	3	OL_O2,OL_O3,OL_N2,OC_N2,OL_N2O,OL_CO,OL_SO2,
	4	OL_NO2,OL_OH,O_CLW,NLAYRS,P,FILEOUT,ANGLE)

	   WRITE(*,'(a30,i5)') 'PROCESSING PROFILE NUMBER:',NPR

 111	ENDDO			!Loop over the profiles
 110	CLOSE(1)		!closes the output file
 924	FORMAT (1X,I1,I3,I5,F10.6,3A8) 
 926	FORMAT (E15.7,F10.4,10X,I5,1X,F7.3,15X,F7.3,/,(1P8E15.7))          
 978	FORMAT (1P8E15.7)                                                  
	CLOSE(53)		!closes the MONORTM_PROF.IN file
	CLOSE(55)		!closes the LBLRTM.IN file
	CLOSE(66)		!closes the LBLATM.LOG file
	STOP
 50	WRITE(*,*) 'ERROR OPENING/READING FILE:',fileprof
	STOP
	END

