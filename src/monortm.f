C     path:		$Source$
C     author:		$Author $
C     revision:	        $Revision$
C     created:	        $Date$

	PROGRAM MONORTM
                                                                         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C             ALGORITHM AUTHORS:    
C                                     S. BOUKABARA
C                                     S.A. CLOUGH                      
C                                     R. HOFFMAN                     
C                                                                      
C                                                                      
C                      Atmospheric and Environmental Research Inc. (AER)
C                      131 Hartwell Avenue
C                      Lexington, MA, 02421      
C                                                                      
C----------------------------------------------------------------------
C                                                                      
C               WORK SUPPORTED BY:    THE ARM PROGRAM                  
C                                     OFFICE OF ENERGY RESEARCH        
C                                     DEPARTMENT OF ENERGY                          
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C  Copyright 2002, 2003, Atmospheric & Environmental Research, Inc. (AER).
C  This software may be used, copied, or redistributed as long as it is not sold and this copyright notice is
C  reproduced on each copy made.  This model is provided as is without any express or implied warranties.
C  (http://www.rtweb.aer.com/)


C                                                                      
C**********************************************************************
C   Comments and/or questions are appreciated. 
C   They should be forwarded to:
C   AER Inc. (Sid Boukabara or Tony Clough)
C   131 Hartwell Avenue, Lexington, MA02421, USA
C   By phone  : 1 781 761 2213
C   By Fax    : 1 781 761 2299
C   By E-mail : sboukaba@aer.com, clough@aer.com
C
C
C   Sid Ahmed Boukabara, AER Inc, April, 2001
C                                                                     
C**********************************************************************

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
	!   it as a black box and control everything from the MONORTM.IN
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
	!     INP=1  ->the inputs read by LBLATM input (MONORTM.IN). This 
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
	!     INP comes from the MONORTM.IN file (see instructions)
	!   - scale/tune several parameters for the line coupling and
	!     the continuum computation, see :
	!     XSLF   : scaling factor for self continuum (default=1)
	!     XFRG   : scaling factor for foreign continuum (default=1)
	!     XCN2   : scaling factor for N2 continuum (default =1)
	!     SCLCPL : scaling factor for Line Coupling Parameters
	!     SCLHW  : scaling factor for the pressure Dependence of the 
	!              halfwidth of the O2 0-zero band
	!     Y0RES  : Y0RES of the line coupling coeffs.
	!
	!   - output the simulations in radiances or in brightness
	!     temperatures. Note that the internal calculations are first
	!     done in radiances. IOUT comes from MONORTM.IN
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
	!   - Extended O2 lines to submillimeter. Extensive validation
	!     by comparison to Rosenkranz model and MWR data.
	!     Update of the LBLATM module (accepts inputs at pressure 
	!     grid, along with altitude grid). 
	!     Fixed the handling of N2 amount coming from LBLATM (which
	!     depends on the number of molecules NMOL). Added version 
	!     numbers and comments in the output file.
	!     Adopted accurate constants values. SAB.
	!     Sid Ahmed Boukabara. Dec 14th 2001.
        !   - Updated on January 7th 2002. ARM option (INP=2) updated and
	!     made more efficient after Jim's comments. (INP=3) option optimized.
	!     WV line intensities modified in the microwave (see Tony's email).
	!     Sid Ahmed Boukabara AER Inc. 2002. 
	!   - Updated on October 2nd 2002. SAB. Major changes: fix in the 
	!     way we read radiosondes, no correction of slant path 
	!     (lblatm does that). Speed option implemented.
	!     We added also the possibility to run Voigt or Lorentz
	!     line shape (speed up process) depending on the current
	!     condition (parameter zeta). The pressure induced
	!     shifted frequency is also passed to the line shape 
	!     computation (instead of the spectroscopic wavenumber).
	!     Also, in the uplooking configuration we compute only the 
	!     downwelling radiance (again, for speed purposes).

        !     September 2003: Modified spectral lines file to improve agreement 
        !     with SGP MWRP data (provided by Nico Cimini). Scaled O2 line coupling
        !     parameters: Y * 0.87, G* 0.
        !     
	!
	!***************************************************************
	include "declar.incl"
	INTEGER NWN,I,ICPL,IS,IOUT,IRT,J,ICNTNM,IATM,INP,IVC
	REAL*8 V1,V2,SECANT,XALTZ 
	REAL TMPSFC,TPROF(mxlay),qprof(mxlay),press(mxlay)
        REAL zvec(mxlay),dzvec(mxlay),zbnd(mxfsc),zbnd2(mxfsc)
        CHARACTER HVRATM*15,HVRMODM*15,HVRSUB*15,HVRMON*15
        CHARACTER HVRREL*15, HVRSPEC*15
	CHARACTER fileARMlist*64,hmod*60,CTYPE*3
	CHARACTER fileprof*80,HFILE*80,FILEOUT*60,ht1*4,ht2*4
	CHARACTER*60 FILEIN,FILESONDE,FILELOG
	character*8 XID,HMOLID,YID,HDATE,HTIME
	character*1 hmol_scal
	character*10 holn2
	real tarray(2)
	COMMON /CVRMON  / HVRMON
        COMMON /CVRATM  / HVRATM
        COMMON /CVRMODM / HVRMODM
	COMMON /CVRSUB  / HVRSUB
	COMMON /CVRREL  / HVRREL
	COMMON /CVRSPEC  / HVRSPEC
	COMMON /PATHD/ P,T,WKL,WBRODL,DVL,WTOTL,ALBL,ADBL,AVBL,
     &     H2OSL,IPTH,ITYL,SECNTA,HT1,HT2,ALTZ,PZ,TZ
	common /profil_scal/ nmol_scal,hmol_scal(64),xmol_scal(64)

	COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,
     1    ADBAR,AVBAR,  
     1    AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      
     2    DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,      
     3    ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,     
     4    EXTID(10)    
	COMMON /BNDPRP/ TMPSFC,BNDEMI(3),BNDRFL(3),IBPROP          
	COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4), 
     1    WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   
     2    EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    
	COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,  
     1    NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,       
     2    NLTEFL,LNFIL4,LNGTH4                                 

	DIMENSION WMT(64)

	!---INPUTS & GENERAL CONTROL PARAMETERS
	IVC=2    !if=2->CKD2.4  MPMf87/s93 (if=3)
	ICPL=1   !=1->cpl =0->nocpl
	IDU=1    !IDU=0->from top 2 surface (not supported), IDU=1->the opposite
	XSLF=1.  !scaling factor (SLF cont)
	XFRG=1.  !scaling factor (FRG cont)
	XCN2=1.  !scaling factor (N2 cont)
	SCLCPL= 1. !scaling factor (Line Coupling Parameters)
	SCLHW=1. !scaling factor (Pressure Dependence of the halfwidth of the 0 band)
	Y0RES=0. !Y0RES of the line coupling coeffs
	IPUNCH=1 !flag to create (1) or not (0) TAPE7 in case INP=2

	!---FILES NAMES ALL PUT HERE FOR CONVENIENCE
	FILEIN      ='MONORTM.IN'
	FILESONDE   ='SONDE.IN'
	fileARMlist ='ARM.IN'
	fileprof    ='MONORTM_PROF.IN'
	HFILE       ='spectral_lines.dat'
	FILEOUT     ='MONORTM.OUT'
	FILELOG     ='MONORTM.LOG'

	!---Initializations of the version numbers
	HVRMON      ='NOT USED       ' 
	HVRATM      ='NOT USED       ' 
	HVRMODM     ='NOT USED       ' 
	HVRSUB      ='NOT USED       ' 
	HVRREL      ='NOT USED       ' 

	!---Version number of MonoRTM
	HVRMON = '$Revision$' 

	!---Release number of MonoRTM
	HVRREL = 'Release  3.2'

	!---GET THE PROFILES NUMBER

	CALL GETPROFNUMBER(INP,FILEIN,fileARMlist,fileprof,
     1    NPROF,filearmTAB)
	!---File Unit numbers/Open files
	IPU  =7 
	IPR  =66
	IOT  =1
	IRD  =55
	IPF  =53
	OPEN (IPR,FILE=FILELOG,STATUS='UNKNOWN',ERR=2000)        
	OPEN (IRD,FILE=FILEIN,STATUS='UNKNOWN',ERR=3000)        
	OPEN (IOT,file=FILEOUT,status='unknown',
     1    form='formatted',ERR=4000)
	IF (INP.EQ.3) OPEN(IPF,file=fileprof,status='old',
     1    form='formatted',err=6000)

	!---Get info about IBMAX/ZBND/H1/H2...
	CALL RDLBLINP(IATM,IOUT,IRT,NWN,WN,FILEIN,
     1    ICNTNM,    INP,IBMAX,ZBND,H1f,H2f,ISPD)

	!---Write header in output file
	WRITE(IOT,'(a)') 'MONORTM RESULTS:'
	WRITE(IOT,'(a)') '----------------' 
	WRITE(IOT,'(a5,I8)') 'NWN :',NWN 

	!---PRINT OUT MONORTM VERSION AND PROFILES NUMBER
	CALL start(NPROF,INP)

	!---CHECK INPUTS AND THEIR CONSISTENCY WITH MONORTM
	CALL CHECKINPUTS(NWN,NPROF,NWNMX,NPROFMX,INP)

	!---Loop over the number of profiles to be processed
	NREC=0
	DO 111 NPR=1,NPROF
	   !*********************************************
	   !* First Step: Read in the inputs
	   !*********************************************	   


	   !---Inputs: MONORTM.IN (TAPE5-type of file)
	   IF (INP.EQ.1) THEN
	      IF (NPR.EQ.1) THEN
		 REWIND(IPU)
		 REWIND(IRD)
		 IPASS=0
	      ENDIF
	      CALL RDLBLINP(IATM,IOUT,IRT,NWN,WN,FILEIN,
     1	         ICNTNM,    INP,IBMAX,ZBND,H1f,H2f,ISPD)
	   ENDIF


	   !---Inputs: wave#/path/angle from MONORTM.IN, profiles from ARM sondes
	   IF (INP.EQ.2) THEN
	      IF (NPR.EQ.1) REWIND(IPU)
	      CALL ARM2LBLATM(filearmTAB(NPR),IFLAG,ilaunchdate, 
     1  	   ilaunchtime,ibasetime,iserialnumber,isondeage,
     2	           NWN,WN,V1,V2,DVSET,FILESONDE,NLAYRS,IBMAX,ZBND,
     3  	   ANGLE,H1F,H2F,NMOL,IPUNCH)

	      IF (IFLAG.EQ.2) THEN
		 WRITE(*,'(a30,i5,a8)') 'PROCESSING PROFILE NUMBER:',
     1	         NPR,' FLAG=2'
		 GOTO 111
	      ENDIF
	      OPEN (IRD,FILE=FILESONDE,STATUS='UNKNOWN',ERR=5000)        
	      CALL RDLBLINP(IATM,IOUT,IRT,NWN,WN,FILESONDE,
     1             ICNTNM,    INP,IBMAX2,ZBND2,H1,H2,ISPD)
	      CLOSE(IRD)
	   ENDIF


	   !---Inputs: MONORTM_PROF.IN (TAPE7 consistent)
	   IF (INP.EQ.3) THEN

c	      READ (IPF,'(1x,i5,10a8)') ipass, xid

	      READ (IPF,925,END=110,ERR=50) IFORM,NLAYRS,NMOL,SECNT0,
     1             HMOD,HMOD,H1,H2,ANGLE,LEN 

	      IF (ANGLE.GT.90.) IRT = 1 !space-based observer (looking down) 
	      IF (ANGLE.LT.90.) IRT = 3 !ground-based observer (looking up)
	      IF (ANGLE.EQ.90.) IRT = 2 !limb measurements
	      DO IL=1,NLAYRS
		 SECNTA(IL)=SECNT0
		 IF (IL.EQ.1) THEN  
		    READ (IPF,*,END=110,ERR=50) P(IL),T(IL),  
     1                   IPATH,ALTZ(IL-1),PZ(IL-1),        
     2	                 TZ(IL-1),ALTZ(IL),  PZ(IL),  TZ(IL)  
		    TMPSFC=TZ(IL-1)
		 ELSE                                             
		    READ (IPF,*,END=110,ERR=50) P(IL),T(IL),  
     1	                 IPATH,ALTZ(IL),PZ(IL),TZ(IL)   
		 ENDIF                                              
		 READ(IPF,978,END=110,ERR=50) (WKL(K,IL),K=1,7),
     1	              WBRODL(IL)          
		 IF (NMOL.GT.7) READ(IPF,978) (WKL(K,IL),K=8,NMOL)
	      ENDDO 
	   ENDIF
 
c_______________________________________________________________________
c
c     at this point scale profile if option selected
c
	   if (nmol_scal .gt. 0 ) call profil_scal_sub(nlayrs)

c_______________________________________________________________________
C
	   !---PREPARE THE INPUTS FOR MODM
	   DO IL=1,NLAYRS
	      W_wv(IL)=WKL(1,IL)
	      W_o2(IL)=WKL(7,IL)
	      IF (NMOL.GE.22) THEN
		 W_n2(IL)=WKL(22,IL)
		 W_other(IL)=WBRODL(IL)
	      ENDIF
	      IF (NMOL.LT.22) THEN
		 W_n2(IL)=WBRODL(IL)
		 W_other(IL)=0.
	      ENDIF
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

	   !---COLUMN WATER VAPOR/LIQUID WATER
	   CALL INTEGR(W_wv,CLW,NLAYRS,WVCOLMN,CLWCOLMN)

	   !***********************************************
	   !* Third Step: OPTICAL DEPTHS COMPUTATION
	   !***********************************************	

           CALL MODM(IVC,ICPL,NWN,WN,NLAYRS,P,T,
     4	                             XSLF,XFRG,XCN2,
     5	        SCLCPL,SCLHW,Y0RES,HFILE,ICNTNM,ISPD)
	   
	   !***********************************************
	   !* Fourth Step: CORRECT FOR THE SLANT PATH
	   !***********************************************	   
	   !---CORRECT THE OPTDEPTHS FOR THE SLANT PATH
	   !CALL CORR_OPTDEPTH(INP,NLAYRS,SECNTA,NWN,ANGLE,O,IRT)

	   !***********************************************
	   !* Fifth Step: RADIATIVE TRANSFER
	   !***********************************************	   

	   CALL RTM(IOUT,IRT,NWN,WN,NLAYRS,T,TZ,
     1          TMPSFC,  RUP,TRTOT,RDN,REFLC,EMISS,RAD,TB,IDU)				

	   !***********************************************
	   !* Sixth Step: WRITE OUT THE RESULTS
	   !***********************************************	   
	   CALL STOREOUT(NWN,WN,RAD,TB,TRTOT,SCLCPL,SCLHW,NREC,
     1          WVCOLMN,0,XSLF,0,Y0RES,0,INP,NPR,ilaunchdate,
     1          ilaunchtime,ibasetime,iserialnumber,isondeage,
     2          CLWCOLMN,TMPSFC,REFLC,EMISS,
     4	                              NLAYRS,P,FILEOUT,ANGLE)

	   WRITE(*,'(a30,i5)') 'PROCESSING PROFILE NUMBER:',NPR

 111	ENDDO			!Loop over the profiles

	!---Write out tail of the output file
	WRITE(IPR,'(a)') 
	WRITE(IPR,'(a)') '--------------------------------------'
	WRITE(IPR,1000) HVRREL,HVRSPEC,HVRMON,HVRMODM,HVRSUB,HVRATM  

	!---Different formats
 924	FORMAT (1X,I1,I3,I5,F10.6,3A8) 
 925	FORMAT(1X,I1,I3,I5,F10.6,2A8,4X,F8.2,4X,F8.2,5X,F8.3,5X,I2) 
 926	FORMAT (E15.7,F10.4,10X,I5,1X,F7.3,15X,F7.3,/,(1P8E15.7))    
 972	FORMAT (64a1)
 973	FORMAT (7e15.7,/,(8e15.7,/))
 978	FORMAT (1P8E15.7)                                             
 1000	FORMAT ('Modules and versions used in this calculation:',/,/,
     &    A15,/,/,
     &    5x,  'spectral file :',          5X,A15,/, 
     &    5x,  'monortm.f     : ',         4X,A15,10X,
     &             'modm.f           :  ', 4X,A15,/,
     &    5x,'monortm_sub.f : ',           4X,A15,10X,
     &             'lblatm_monortm.f :  ', 4X,A15)
 110	CONTINUE

	!---Close all files
        CLOSE(IPF)		!closes the MONORTM_PROF.IN file
	CLOSE(IOT)		!closes the OUTPUT file
	CLOSE(IRD)		!closes the SONDE.IN file
	CLOSE(IPR)		!closes the LBLATM.LOG file
	CLOSE(IPU)		!closes the MONORTM.IN file
	STOP

	!---Error messages
 50	WRITE(*,*) ' EXIT; ERROR READING :',fileprof
	STOP
 2000	WRITE(*,*) ' EXIT; ERROR OPENING :',FILELOG          
	STOP
 3000	WRITE(*,*) ' EXIT; ERROR OPENING :',FILEIN           
	STOP
 4000	WRITE(*,*) ' EXIT; ERROR OPENING :',FILEOUT
	STOP
 5000	WRITE(*,*) ' EXIT; ERROR OPENING :',FILESONDE           
	STOP
 6000	WRITE(*,*) ' EXIT; ERROR OPENING :',FILEPROF           
	STOP
	END

	!---Block data to be consistent with LBLRTM/LBLATM
	Block Data phys_consts
	COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     1       RADCN1,RADCN2 
        DATA PI / 3.1415926535897932 /   ! from http://www.cecm.sfu.ca/pi9

c---------------------------------------------                
c       Constants from NIST 01/11/2002
c---------------------------------------------                
	DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
     1       CLIGHT / 2.99792458E+10 /, 
     2       AVOGAD / 6.02214199E+23 /, ALOSMT / 2.6867775E+19 /,
     3       GASCON / 8.314472  E+07 /,
     4       RADCN1 / 1.191042722E-12 /, RADCN2 / 1.4387752    /
c---------------------------------------------                
c       units are generally cgs
c       The first and second radiation constants are taken from NIST.
c       They were previously obtained from the relations:
c       RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07      
c       RADCN2 = PLANCK*CLIGHT/BOLTZ 
c---------------------------------------------                
	end
c___________________________________________________________________
c___________________________________________________________________
c___________________________________________________________________


