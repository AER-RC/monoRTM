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
C                      Atmospheric and Environmental Research Inc. (AER)
C                      131 Hartwell Avenue
C                      Lexington, MA, 02421      
C                                                                      
C----------------------------------------------------------------------
C                                                                      
C        WORK SUPPORTED BY:    THE ARM PROGRAM                  
C                              OFFICE OF ENERGY RESEARCH        
C                              DEPARTMENT OF ENERGY   
C                              THE JOINT CENTER FOR SATELLITE DATA ASSIMILATION
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).
C  This software may be used, copied, or redistributed as long as it is not sold and this copyright notice is
C  reproduced on each copy made.  This model is provided as is without any express or implied warranties.
C  (http://www.rtweb.aer.com/)

C                                                                      
C**********************************************************************
C   Comments and/or questions are appreciated. 
C   They should be forwarded to:
C   AER Inc. (Karen Cady-Pereira or Vivienne Payne)
C   131 Hartwell Avenue, Lexington, MA02421, USA
C   By phone  : 1 781 761 2288
C   By Fax    : 1 781 761 2299
C   By E-mail : cadyp@aer.com, vpayne@aer.com
C
C
C   Vivienne H. Payne, AER Inc, August, 2008
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
	!   - turn on/off the line coupling (see ICPL)
	!     ICPL=1 -> Line coupling turned on  
	!     ICPL=0 -> Line coupling turned off (even if we have line 
        !               couplng coefficients in the spectroscopic 
	!               selected_lines)
	!   - assume that the atmospheric profiles are organized
	!     from the top to the surface (or the opposite), see IDU.
	!     IDU=1  -> profiles from surface to top (default)
	!
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
        !   - output the optical depths for each layer
	!     IOD comes from MONORTM.IN
	!     IOD    : =0 -> do not write out layer optical depths
	!     IOD    : =1 -> write out layer optical depths
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

        !   - September 2003: Modified spectral lines file to improve agreement 
        !     with SGP MWRP data (provided by Nico Cimini). Scaled O2 line coupling
        !     parameters: Y * 0.87, G* 0.
	!   - 2006: Implementation of Tretyakov et al (2005) O2 line coupling.
        !     Validated against ARM MWRP data (see Cadeddu et al, 2007)
	!   - 2007: Updated spectral line file to change the widths and 
        !     temperature dependences of the widths for the 22 and 183 GHz lines
        !     (see Payne et al., 2008)
	!   - 2008: Extensive update to enable the use of MonoRTM beyond 
        !     the microwave region and to use the MT_CKD continuum.
	!   - June 2009: Updates to self and foreign broadened H2O microwave 
        !     continuum based on ARM SGP MWR and GVRP data at 31.4 and 170 GHz 
        !     and ARM FKB (COPS) data at 150 GHz.   
	!
	!
	!***************************************************************
        USE ModmMod, ONLY: MODM
        USE CntnmFactors, ONLY: CntnmFactors_t
	IMPLICIT REAL*8           (V) ! for consistency with LBLRTM routines
	include "declar.incl"
	INTEGER NWN,I,ICPL,IS,IOUT,IOD,IRT,J,IATM
	REAL*8 V1,V2,SECANT,XALTZ 
	REAL TMPSFC,TPROF(mxlay),qprof(mxlay),press(mxlay)
        REAL zvec(mxlay),dzvec(mxlay),zbnd(mxfsc),zbnd2(mxfsc)
	REAL O(NWNMX,MXLAY),OC(NWNMX,MXMOL,MXLAY),
     &     O_BY_MOL(NWNMX,MXMOL,MXLAY),O_CLW(NWNMX,MXLAY),
     &	   odxsec(nwnmx,mxlay)
        CHARACTER HVRATM*15,HVRMODM*15,HVRSUB*15,HVRMON*15
        CHARACTER HVRREL*15, HVRSPEC*15
	CHARACTER fileARMlist*64,hmod*60,CTYPE*3
	CHARACTER fileprof*80,HFILE*80,FILEOUT*60,ht1*4,ht2*4
	CHARACTER*60 FILEIN,FILELOG
	character*8 XID,HMOLID,YID,HDATE,HTIME
	character*1 hmol_scal
	character*10 holn2
	CHARACTER*10 XSFILE,XSNAME,ALIAS,XNAME,XFILS(6),BLANK  
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
        
        COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(mx_xs),XAMNT(mx_xs,MXLAY)

        COMMON /XSECTF/ XSFILE(6,5,mx_xs),XSNAME(mx_xs),ALIAS(4,mx_xs) 
        COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS),     
     *                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),      
     *                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS), 
     *                NUMXS,IXSBIN

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
	DIMENSION Wdrair(mxlay)
        
        TYPE(CntnmFactors_t) :: cntnmScaleFac

c------------------------------------
c Variables for analytic derivative calculation
c These are not used in MonoRTM, but are initialized here in order
c to maintain consistency between contnm.f for MonoRTM and LBLRTM
c note: ipts  = same dimension as ABSRB
c       ipts2 = same dimension as C
	parameter (ipts=5050,ipts2=6000)
	common /CDERIV/ icflg,iuf,v1absc,v2absc,dvabsc,nptabsc,delT_pert,
     &    dqh2oC(ipts),dTh2oC(ipts),dUh2o

	real cself(ipts),cfrgn_aj(ipts)

	icflg = -999
	iuf = 0
	v1absc = 0
	v2absc = 0
	dvabsc = 0
	nptabsc = 0
	delT_pert = 0
c------------------------------------

	!---INPUTS & GENERAL CONTROL PARAMETERS
	ICPL=1   !=1->cpl =0->nocpl
	IDU=1    ! ->profiles input fromn surface to TOA
	SCLCPL= 1. !scaling factor (Line Coupling Parameters)
	SCLHW=1. !scaling factor (Pressure Dependence of the halfwidth of the 0 band)
	Y0RES=0. !Y0RES of the line coupling coeffs

	!---FILES NAMES ALL PUT HERE FOR CONVENIENCE
	FILEIN      ='MONORTM.IN'
	fileARMlist ='ARM.IN'
	fileprof    ='MONORTM_PROF.IN'
	!HFILE       ='spectral_lines.dat'
	HFILE       ='TAPE3'
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
	HVRREL = 'Release  4.6'

	!---GET THE PROFILES NUMBER

	CALL GETPROFNUMBER(IATM,FILEIN,fileARMlist,fileprof,
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
	IF (IATM.EQ.0) OPEN(IPF,file=fileprof,status='old',
     1    form='formatted',err=6000)

	IPASSATM = 0 ! flag to determine whether RDLBLINP has been called before

	!---Get info about IBMAX/ZBND/H1/H2...
	CALL RDLBLINP(IATM,IOUT,IOD,IRT,NWN,WN,FILEIN,
     1    cntnmScaleFac,IXSECT,IBMAX,ZBND,H1f,H2f,ISPD,IPASSATM)


	!---PRINT OUT MONORTM VERSION AND PROFILES NUMBER
	CALL start(NPROF,IATM)

	!---CHECK INPUTS AND THEIR CONSISTENCY WITH MONORTM
	CALL CHECKINPUTS(NWN,NPROF,NWNMX,NPROFMX)

	!---Loop over the number of profiles to be processed
	NREC=0
	DO 111 NPR=1,NPROF
	   !*********************************************
	   !* First Step: Read in the inputs
	   !*********************************************	   


	   !---Inputs: MONORTM.IN (TAPE5-type of file)
	   IF (IATM.EQ.1) THEN
	      IF (NPR.EQ.1) THEN
		 REWIND(IPU)
		 REWIND(IRD)
		 IPASS=0
	      ENDIF
	      CALL RDLBLINP(IATM,IOUT,IOD,IRT,NWN,WN,FILEIN,
     1	         cntnmScaleFac,IXSECT,IBMAX,ZBND,H1f,H2f,ISPD,IPASSATM)
	   ENDIF


	   !---Inputs: MONORTM_PROF.IN (TAPE7 consistent)
	   IF (IATM.EQ.0) THEN

c	      READ (IPF,'(1x,i5,10a8)') ipass, xid

	      READ (IPF,925,END=110,ERR=50) IFORM,NLAYRS,NMOL,SECNT0,
     1             HMOD,HMOD,H1,H2,ANGLE,LEN 

	      IF (ANGLE.GT.90.) IRT = 1 !space-based observer (looking down) 
	      IF (ANGLE.LT.90.) IRT = 3 !ground-based observer (looking up)
	      IF (ANGLE.EQ.90.) IRT = 2 !limb measurements
	      DO IL=1,NLAYRS
		 SECNTA(IL)=SECNT0
		 IF (IFORM.EQ.0) THEN  
                    if (il.eq.1) then 
		       READ (IPF,974,END=110,ERR=50) P(IL),T(IL),secnt,ipath,  
     1                   ALTZ(IL-1),PZ(IL-1), TZ(IL-1),       
     2	                 ALTZ(IL),  PZ(IL),  TZ(IL), CLW(IL)  
                    else
		       READ (IPF,9742,END=110,ERR=50) P(IL),T(IL),secnt,ipath,  
     2	                 ALTZ(IL),  PZ(IL),  TZ(IL), CLW(IL)  
                    endif
		 ELSE                                             
                    if (il.eq.1) then 
		       READ (IPF,975,END=110,ERR=50) P(IL),T(IL),secnt,ipath,  
     1                   ALTZ(IL-1),PZ(IL-1),TZ(IL-1),        
     2	                 ALTZ(IL),PZ(IL),TZ(IL) , CLW(IL)  
                     else
		       READ (IPF,9752,END=110,ERR=50) P(IL),T(IL),secnt,ipath,  
     2	                 ALTZ(IL),PZ(IL),TZ(IL) , CLW(IL)  
                     endif
		 ENDIF                                              
		 READ(IPF,978,END=110,ERR=50) (WKL(K,IL),K=1,7),
     1	              WBRODL(IL)          
		 IF (NMOL.GT.7) READ(IPF,978) (WKL(K,IL),K=8,NMOL)
C     --------------------------------------------------------------
C
C                     MIXING RATIO INPUT
C
C
C     First calculate the column amount of dry air ("WDRAIR")
C     Initialize WDNSTY to WBRODL(L) (always in column density)
C     Determine if each molecule is in column density.
C        - if so, just add to WDNSTY
C        - if not, add to WMXRAT
C
C     NOTE that if WKL is greater than one, then column density 
C               if WKL is less than one, then mixing ratio
C
         WDNSTY = WBRODL(IL)
         WMXRAT = 0.0
         WDRAIR(IL) = 0.0

         DO 22 M = 2,NMOL
            IF (WKL(M,IL).GT.1) THEN
               WDNSTY = WDNSTY + WKL(M,IL)
            ELSE
               WMXRAT = WMXRAT + WKL(M,IL)
            ENDIF
 22      CONTINUE

C
C        EXECUTE TESTS TO ENSURE ALL COMBINATION OF COLUMN DENSITIES
C        AND MIXING RATIOS FOR EACH LAYER HAVE BEEN PROPERLY SPECIFIED.

C        IF THE LAYER SUM OF MIXING RATIOS IS LESS THAN ONE (WHICH
C        IT SHOULD BE, GIVEN THAT WBROAD CONTRIBUTES TO THE DRY AIR 
C        MIXING RATIO), THEN COMPUTE DRY AIR BY DIVIDING THE TOTAL
C        MOLECULAR AMOUNTS GIVEN IN DENSITY BY THE FRACTION OF DRY
C        AIR (MIXING RATIO) THOSE MOLECULES COMPRISE.
C
C        IF THE LAYER SUM OF MIXING RATIOS IS GREATER THAN OR EQUAL
C        TO ONE, THAN AN ERROR HAS OCCURRED, SO STOP THE PROGRAM.
C        WBROAD IS ALWAYS LISTED IN COLUMN DENSITY, SO THE SUM OF
C        THE GIVEN MIXING RATIOS MUST ALWAYS BE LESS THAN ONE.
C

         IF (WBRODL(IL).LT.1.0 .AND. WBRODL(IL).NE.0.0) THEN
c            WRITE(IPR,918) IL
c            WRITE(*,918) IL
            STOP
         ENDIF

         IF (WDNSTY.EQ.0.0 .AND. WMXRAT.NE.0.0) THEN
c            WRITE(IPR,921) IL,WDNSTY,WMXRAT
c            WRITE(*,921) IL,WDNSTY,WMXRAT
            STOP 'WMXRAT AND/OR WDNSTY NOT PROPERLY SPECIFIED IN PATH'
         ENDIF
         IF (WMXRAT.LT.1.0) THEN
            WDRAIR(IL) = WDNSTY/(1.0-WMXRAT)
         ELSE
c           WRITE(IPR,921) IL,WMXRAT, WDNSTY
c           WRITE(*,921) IL,WMXRAT, WDNSTY
            STOP 'WMXRAT EXCEEDS 1.0'
         ENDIF

         IF (WKL(1,IL).LE.1.0 .AND. WKL(1,IL) .NE. 0.0
     *        .AND. WDRAIR(IL).EQ.0.0) THEN
c           WRITE(IPR,921) IL,WKL(1,IL),WDRAIR(IL)
c           WRITE(*,921) IL,WKL(1,IL),WDRAIR(IL)
            STOP 'WMXRAT NOT PROPERLY SPECIFIED IN PATH'
         ENDIF

C
C     NOW CONVERT ALL OTHER MOLECULES WHICH MAY BE IN MIXING RATIO
C     TO MOLECULAR DENSITY USING WDRAIR(L)
C
         DO 25 M = 1,NMOL
            IF (WKL(M,IL).LT.1.) WKL(M,IL) = WKL(M,IL)*WDRAIR(IL)
 25      CONTINUE
C
C     --------------------------------------------------------------
C

      END DO


	      
	      IF (IXSECT.GE.1) THEN
	         READ(IPF,930) IXMOLS,IXSBIN
		 XV1 = MINIMUM(WN,nwn)
		 XV2 = MAXIMUM(WN,nwn)

		 CALL XSREAD (ipf,XV1,XV2)                                           
		 WRITE (IPR,932) (I,XSNAME(I),I=1,IXMOLS)
		 READ (IPF,900) IFRMX,NLAYXS,IXMOL,SECNTX,HEDXS                   
		 IF (IXMOL.EQ.0) THEN                                             
		    WRITE (IPR,935) IXMOL                                         
		    STOP ' PATH - IXMOL 0 '                                       
	         ENDIF                                                            
	         IF (IXMOL.NE.IXMOLS) THEN                                        
	            WRITE (IPR,937) IXMOL,IXMOLS                                  
	            STOP ' PATH - IXMOL .NE. IXMOLS '                             
	         ENDIF                                                            
	         IF (NLAYRS.NE.NLAYXS) THEN                                       
	            WRITE (IPR,940) NLAYRS,NLAYXS                                 
	            STOP ' PATH - NLAYRS .NE. NLAYXS '                            
	         ENDIF  
              
		 SECNTX = ABS(SECNTX)                                             
		 WRITE (IPR,942) SECNTX,NLAYXS,IXMOLS,HEDXS                       
C       
		 DO 40 L = 1, NLAYXS                                              
C                                                                         
		     IF (L.EQ.1) THEN                                              
			 READ (IPF,910) PAVX,TAVX,SECKXS,CINPX,IPTHKX,ALTXB,
     &			     PZXB,TZXB,ALTXT,PZXT,TZXT                          
		     ELSE                                                          
			 READ (IPF,915) PAVX,TAVX,SECKXS,CINPX,IPTHKX,ALTXT, 
     &		              PZXT,TZXT                                              
		     ENDIF                                                         
			 READ (IPF,978) 
     *                    (XAMNT(M,L),M=1,7),WBRODX                   
			 IF (IXMOL.GT.7) 
     *                    READ (IPF,978) (XAMNT(M,L),M=8,IXMOL)      

  40		 CONTINUE

	     ENDIF ! test on IXSECT=1

	 ENDIF ! test on IATM=0
 
c_______________________________________________________________________
c
c     at this point scale profile if option selected
c
	   if (nmol_scal .gt. 0 ) call profil_scal_sub(nlayrs)

c_______________________________________________________________________
C
	   !***********************************************
	   !* Second Step: SET-UP THE EMISS/REFLEC VECTORS
	   !***********************************************	   
	   CALL EMISS_REFLEC(NWN,EMISS,REFLC,WN) 

	   NREC=NREC+1

	   !---COLUMN WATER VAPOR/LIQUID WATER
	   CALL INTEGR(WKL(1,:),CLW,NLAYRS,WVCOLMN,CLWCOLMN)

	   !***********************************************
	   !* Third Step: OPTICAL DEPTHS COMPUTATION
	   !***********************************************	

           CALL MODM(IPR,ICPL,NWN,WN,dvset,NLAYRS,P,T,
     1                 O,O_BY_MOL, OC, O_CLW, ODXSEC,
     4	               NMOL,WKL,WBRODL,
     5	        SCLCPL,SCLHW,Y0RES,HFILE,cntnmScaleFac,ixsect,ISPD)
	   
	   !***********************************************
	   !* Fifth Step: RADIATIVE TRANSFER
	   !***********************************************	   

	   CALL RTM(IOUT,IRT,NWN,WN,NLAYRS,T,TZ,O,
     1          TMPSFC,  RUP,TRTOT,RDN,REFLC,EMISS,RAD,TB,IDU)				

	   !***********************************************
	   !* Sixth Step: WRITE OUT THE RESULTS
	   !***********************************************	   
	   CALL STOREOUT(NWN,WN,WKL,WBRODL,RAD,TB,TRTOT,NPR,
     1          O,O_BY_MOL, OC, O_CLW, ODXSEC,
     2          WVCOLMN,CLWCOLMN,TMPSFC,REFLC,EMISS,
     4	        NLAYRS,NMOL,ANGLE,IOT,IOD,FILEOUT)

	   WRITE(*,'(a30,i5)') 'PROCESSING PROFILE NUMBER:',NPR

 111	ENDDO			!Loop over the profiles

	!---Write out tail of the output file
	WRITE(IPR,'(a)') 
	WRITE(IPR,'(a)') '--------------------------------------'
	WRITE(IPR,1000) HVRREL,HVRSPEC,HVRMON,HVRMODM,HVRSUB,HVRATM  

	!---Different formats
 900	FORMAT (1X,I1,I3,I5,F10.2,15A4) 
 910	FORMAT (E15.7,F10.4,F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))
 911	FORMAT (3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))                          
 915	FORMAT (E15.7,F10.4,F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))
 916	FORMAT (3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))              
 924	FORMAT (1X,I1,I3,I5,F10.6,3A8) 
 925	FORMAT(1X,I1,I3,I5,F10.6,2A8,4X,F8.2,4X,F8.2,5X,F8.3,5X,I2) 
 926	FORMAT (E15.7,F10.4,10X,I5,1X,F7.3,15X,F7.3,/,(1P8E15.7))    
 927	FORMAT (8E10.3) 
 930	FORMAT (I5,5X,I5)
 932	FORMAT (/,'  THE CROSS-SECTION MOLECULES SELECTED ARE: ',/,/,(5X,   
     *        I5,3X,A))                                                   
 935	FORMAT (/,'***** IXMOL = ',I5,' *****',/)                           
 937	FORMAT (/,'***** IXMOL = ',I5,' .NE. IXMOLS = ',I5,' *****',/)     
 940	FORMAT (/,'***** NLAYRS = ',I5,' .NE. NLAYXS = ',I5,' *****',/)     
 942	FORMAT (/,'0 SECANTX  =',F13.4,/'0 NLAYXS=',I4,/'0 ISMOLS=',I4,/,   
     &   '0',15A4)                                                   
 972	FORMAT (64a1)
 973	FORMAT (7e15.7,/,(8e15.7,/))
 974    FORMAT (3f10.4,3x,i2,1x,2(f7.2,f8.3,f7.2),f7.2)
 9742   FORMAT (3f10.4,3x,i2,1x,22x,1(f7.2,f8.3,f7.2),f7.2)
 975    FORMAT (e15.7,2f10.4,3x,i2,1x,2(f7.2,f8.3,f7.2),f7.2)
 9752   FORMAT (e15.7,2f10.4,3x,i2,1x,22x,1(f7.2,f8.3,f7.2),f7.2)
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
 6000	WRITE(*,*) ' EXIT; ERROR OPENING :',FILEPROF           
	STOP
	END

c___________________________________________________________________
c___________________________________________________________________
c___________________________________________________________________


