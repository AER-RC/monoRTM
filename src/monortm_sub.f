C     path:		$Source$
C     author:		$Author $
C     revision:	        $Revision$
C     created:	        $Date$

	SUBROUTINE RTM(IOUT,IRT,NWN,WN,NLAY,T,TZ,
     &    TMPSFC,O,RUP,TRTOT,RDN,REFLC,EMISS,RAD,TB,IDU)
C-------------------------------------------------------------------------------
C
C     PROGRAM:  RTM
C     -------
C
C     AUTHOR: Sid-Ahmed Boukabara 
C     ------
C
C     AFFILIATION: ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
C     -----------
C
C     DATE OF CREATION : May 1999
C     ----------------
C
C     AIM: This program is aimed at the simulation of the 
C     ---- radiances using the Radiative Transfer equation. 
C
C     INPUTS:
C     ------
C     - IOUT  : Flag to compute the radiances (IOUT=0) or both the 
C               radiances and the brightness temperatures (IOUT=1)
C     - IRT   : Flag to compute the radiative transfer
C               IRT=1 from the surface to the space (satellite)
C               IRT=2 Limb measurements.
C               IRT=3 from the space to the surface (ground instrument)
C     - NWN   : Number of wavenumbers to be treated
C     - WN    : Vector of NWN wavenumbers [in cm-1], one should note that
C               this input could be a scalar (associated with NWN=1)
C     - NLAY  : Number of layers to be treated.
C     - P     : Vector of NLAY pressures (in mbar), one should note that
C               this input could be a scalar (associated with NLAY=1)
C     - T     : Vector of NLAY temperatures (layers averaged T) [in Kelvin]
C     - TZ    : Vector of NLAY+1 temperatures (levels T) [in Kelvin]
C     - TMPSFC: Surface temperatures [in Kelvin]
C     - O     : Total Optical depths in nepers, (NWNxNLAY)
C     - REFLC : Reflectivity vector of the Surface (NWN dimension)
C     - EMISS : Emissivity vector of the surface (NWN dimension)
C     - IDU   : Index for the Up/Down format of the profiles
C               IDU=0->the profiles are given from the top to the surface
C               IDU=1->the profiles are given from the surface to the top
C
C     Note:
C             If the surface was in thermodynamical equilibrium, 
C             REFLC should be equal to (1-ESFC)
C
C     OUTPUTS:
C     -------
C     - RAD    : An array of NWN elts containing the radiances at the WN 
C                frequencies.
C     - TB     : An array of NWN elts containing the brightness temperatures 
C                at the WN frequencies.
C     - RUP    : Upwelling Contribution Radiance (NWN dimension)
C     - RDN    : DownWelling Contribution Radiance (NWN dimension)
C     - TRTOT  : Total Transmittance
C
C       Note:
C       -----
C       RTM takes into account the cosmic background contribution.
C       The cosmic radiation is hard coded (2.75 Kelvin). But this should
C       be read from the input file (to be updated in next version).
C
C-------------------------------------------------------------------------------
	include "declar.incl"
	INTEGER NWN,NLAY,IRT,I,IOUT,IDU
	REAL RADCN1,RADCN2
	REAL*8 V
	CHARACTER HVRSUB*15
	REAL TMPSFC,ESFC,RSFC,SURFRAD,ALPH,COSMOS,TSKY
	REAL fbeta,beta,bb_fn,X
	COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
	1    RADCN1,RADCN2 
	COMMON /CVRSUB/ HVRSUB
	BB_fn(V,fbeta)  = RADCN1*(V**3)/(EXP(V*fbeta)-1.)   
	HVRSUB = '$Revision$' 
	!---Up and Down radiances
	CALL RAD_UP_DN(T,NLAY,TZ,WN,rup,trtot,rdn,O,NWN,IDU,IRT)
	!---RADIATIVE TRANSFER
	TSKY=2.75 !Cosmic background in Kelvin
	beta= RADCN2/TMPSFC
	alph= RADCN2/TSKY
	DO I=1,NWN
	   SURFRAD  = bb_fn(WN(I),beta)
	   COSMOS   = bb_fn(WN(I),alph)
	   ESFC=EMISS(I)
	   RSFC=REFLC(I)
	   IF (IRT.EQ.1) RAD(I)=RUP(I)+
	1   trtot(i)*( (rsfc*(trtot(i)*COSMOS)+rdn(i))+esfc*SURFRAD ) ! sac 06/04/02
c
	   IF (IRT.EQ.3) RAD(I)=RDN(I)+(trtot(i)*COSMOS)
c
	   IF (IRT.EQ.2) RAD(I)=RUP(I)+ 
	1   trtot(i)*( (trtot(i)*COSMOS)+rdn(i) )  !sac 06/04/02
c
	   IF (IOUT.EQ.1) THEN
	      X=RADCN1*(WN(I)**3)/RAD(I)+1.
	      TB(I)=RADCN2*WN(I)/log(X)
	   ENDIF
	ENDDO
 	RETURN
	END 


	SUBROUTINE RAD_UP_DN(T,nlayer,TZ,WN,rup,trtot,rdn,O,NWN,
	1    IDU,IRT)
	IMPLICIT REAL*8 (V)      
	include "declar.incl"
	COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
	1    RADCN1,RADCN2 
	INTEGER  layer,nlayer,NWN,IDU,lmin,lmax,nl
	REAL          beta,beta_a,bb,bba
	!---local variables
	REAL          bbVEC(MXLAY),bbaVEC(0:MXLAY),ODTOT(NWNMX)
	BB_fn(V,fbeta)  = RADCN1*(V**3)/(EXP(V*fbeta)-1.)
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

	   IF (IRT.NE.3) THEN	!compute RUP only when IRT<>3
	      ODT=ODTOT(I)
	      DO 70 layer = 1,nlayer,1
		 bb  = bbVEC(layer)
		 bba = bbaVEC(layer)
		 ODVI = O(I,layer)
		 TRI = EXP(-ODVI)
		 ODT=ODT-ODVI
		 TRTOT(I)= EXP(-ODT)
		 pade=0.193*ODVI+0.013*ODVI**2
		 RUP(I)= RUP(I)+
	1	      TRTOT(I)*(1.-TRI)*(bb+pade*bba)/(1.+pade)
 70	      ENDDO
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
 50	   continue
	   TRTOT(I)=EXP(-ODTOT(I))
 60	ENDDO
	RETURN                                                   
	END                                                       



	SUBROUTINE READEM(ICOEF)
	!Reads in emission function values directly from file "EMISSION"
	IMPLICIT REAL*8           (V)
	PARAMETER (NMAXCO=4040)
	COMMON /EMSFIN/ V1EMIS,V2EMIS,DVEMIS,NLIMEM,ZEMIS(NMAXCO)
	READ (ICOEF,900,END=20,ERR=20) V1EMIS,V2EMIS,DVEMIS,NLIMEM
	DO 100 NGNU = 1,NLIMEM
	   READ (ICOEF,910,END=20,ERR=20) ZEMIS(NGNU)
 100	CONTINUE
	RETURN
 900	FORMAT (3E10.3,5X,I5)
 910	FORMAT (E15.7)
 20	WRITE(*,*) 'INCONSISTENT DATA OR ERROR OPENING IN READEM'
	STOP
	END

	SUBROUTINE READRF(ICOEF)
	IMPLICIT REAL*8           (V)
	PARAMETER (NMAXCO=4040)
	COMMON /RFLTIN/ V1RFLT,V2RFLT,DVRFLT,NLIMRF,ZRFLT(NMAXCO)
	READ (ICOEF,900,END=20,ERR=20) V1RFLT,V2RFLT,DVRFLT,NLIMRF
	DO 100 NGNU = 1,NLIMRF
	   READ (ICOEF,910,END=20,ERR=20) ZRFLT(NGNU)
 100	CONTINUE
	RETURN
 900	FORMAT (3E10.3,5X,I5)
 910	FORMAT (E15.7)
 20	WRITE(*,*) 'INCONSISTENT DATA OR ERROR OPENING IN READRF'
	END



	SUBROUTINE RDLBLINP(IATM,IPLOT,IRT,NWN,WN,
	1    FILEIN,ICNTNM,CLW,INP,IBMAXOUT,ZBNDOUT,
	2    H1fout,H2fout,ISPD)
C-------------------------------------------------------------------------------
C
C     SUBROUTINE:  RDLBLINP
C     -----------
C
C     AUTHOR: Sid-Ahmed Boukabara 
C     ------
C
C     AFFILIATION: ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
C     -----------
C
C     DATE OF CREATION : April 2001 
C     Modified on Jan 7th 2002.
C
C
C     AIM: This subroutine reads the control parameters
C     ---- from the input file MONORTM.IN (former TAPE5) 
C          and makes available some of them via common 
C          blocks. This subroutine has been added to make 
C          MONORTM compatible with LBLRTM inputs.
C          S.A. Boukabara AER INC. 1999
C	  IRT=1  
C	 			!=1->Space
C	 			!=2->limb
C	 			!=3->ground
C     UPDATES:
C     --------
C	Extension of the MONORTM.IN option for MonoRTM
C	If V1 or V2 is negative, then we expect a 
C	finite number of wavenumbers to be included
C	in the input file. These wavenumbers do not have 
C	to be equally spaced. This is a special option
C	specific for MonoRTM. This is particularly
C	useful for those interested in simulating 
C	ARM monochromatic frequencies for instance (23.8 
C	and 31.4 GHz)
C	
C	 Sid Ahmed Boukabara, April 2001
C	
C-------------------------------------------------------------------------------
	include "declar.incl"
	REAL*8           V1,V2,SECANT,XALTZ 
	character*4 ht1,ht2
	CHARACTER*1 CMRG(2),CONE,CTWO,CTHREE,CFOUR,CXIDLINE*80
	CHARACTER CDOL*1,CPRCNT*1,CXID*1,CA*1,CB*1,CC*1
	INTEGER IRD,IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,L,           
	1    ISCAN,IFILTR,IPLOT,ITEST,IATM,ILAS,ILNFLG,      
	2    IOD,IXSECT,IRAD,MPTS,NPTS,INFLAG,IOTFLG,JULDAT
	REAL SAMPLE,DVSET,ALFAL0,AVMASS,DPTMIN,DPTFAC,DVOUT 
	REAL TMPBND,XVMID,EMITST,REFTST
	INTEGER IBPROP,IBND,ICOEF,IMRG,LAYTOT,IFORM,NLAYRS,NMOL
	REAL PATH1,PTHODT,SECNT0,ZH1,ZH2,ZANGLE,PAVE,TAVE,SECNTK
	CHARACTER HEAD20*20,HEAD7*7,HEAD5*5,HEAD4*4,CINP*3
	CHARACTER*60 FILEOUT,FILEIN
	REAL SECL(64),WDNSTY,WMXRAT,WDRAIR(MXLAY)
	REAL ZBNDOUT(MXFSC)
	INTEGER IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,        
	1    NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,       
	2    NLTEFL,LNFIL4,LNGTH4,IPTHRK,IPATHL,M
	character*8 XID,HMOLID,YID      
	COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,NOZERO,NOP,H1F,H2F, 
	1    ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,        
	2    XVBAR, HMINF,PHIF,IERRF,HSPACE                     
 	COMMON /BNDRY/ ZBND(MXFSC),PBND(MXFSC),TBND(MXFSC),
	1    ALORNZ(MXFSC),ADOPP(MXFSC),AVOIGT(MXFSC)       
	COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,
	1    IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,             
	2    IPHMID,IPDIM,KDIM,KMXNOM,KMAX                      
	COMMON /PATHD/ P,T,WKL,WBRODL,DVL,WTOTL,ALBL,
	1    ADBL,AVBL,H2OSL,IPTH,ITYL,SECNTA,HT1,HT2,ALTZ,         
	2    PZ,TZ                          
	COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,  
	1    NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,       
	2    NLTEFL,LNFIL4,LNGTH4                                 
	COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,
	1    AVBAR,  
	1    AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      
	2    DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,      
	3    ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,     
	4    EXTID(10)    
	COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP            
	COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),  
	1    WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   
	2    EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    
	DATA CDOL / '$'/,CPRCNT / '%'/
	DATA CONE / '1'/,CTWO / '2'/,CTHREE / '3'/,CFOUR / '4'/,       
	1    CA / 'A'/,CB / 'B'/,CC / 'C'/                   
	CHARACTER CEX*2,CEXST*2,CPRGID*60    
	DATA CEXST/'EX'/
	EQUIVALENCE (CXID,CXIDLINE)                    


				!---record 1.1
 20	READ (IRD,905,END=80,ERR=6000) CXIDLINE 
	IF (CXID.EQ.CPRCNT) THEN
	   WRITE(*,*) '-END OF FILE:',FILEIN
	   RETURN 
	ENDIF                       
	IF (CXID.NE.CDOL) GO TO 20                            
	READ (CXIDLINE,'(1x,10A8)') (XID(I),I=1,10)     
				!---record 1.2
	READ(IRD,925,END=80,ERR=6000) IHIRAC,ILBLF4,             
	1    ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,    
	2    ITEST,IATM, CMRG,ILAS, IOD,IXSECT,IRAD,
	3    MPTS,NPTS,INP,ISPD  
	
	IF ((INP.LE.1).OR.(INP.GE.4)) INP=1
	!---IF INP=1 !MONORTM.IN INPUT
	!---IF INP=2 !ARM SONDES INPUTS
	!---IF INP=3 !MONORTM_PROF.IN INPUT FILE
	
	IF (INP.EQ.2 .AND. IATM.NE.1) STOP 'INP=2 => IATM=1'

	!----CHECKING THE INPUTS FROM RECORD 1.2
	IF (IAERSL.GT.0) THEN
	   WRITE(*,*) 'CURRENTLY MONORTM DOES NOT HANDLE AEROSOLS'
	   WRITE(*,*) 'CHECK YOUR INPUT FILE: ',FILEIN
	   STOP
	ENDIF
	IF (ILBLF4.GT.0) THEN
	   WRITE(*,*) '----------------------------------------'
	   WRITE(*,*) 'WARNING: ILBLF4 IS IGNORED IN MONORTM'
	   WRITE(*,*) 'IN ',FILEIN,' ILBLF4=',ILBLF4
	ENDIF
	IF (IEMIT.NE.1) THEN
	   WRITE(*,*) '----------------------------------------'
	   WRITE(*,*) 'WARNING: IEMIT IS IGNORED IN MONORTM'
	   WRITE(*,*) 'IT IS SET INTERNALLY TO ONE'
	   WRITE(*,*) 'IN ',FILEIN,' IEMIT=',IEMIT
	ENDIF
	IF (ISCAN.NE.0) THEN
	   PRINT *, 'MONORTM DOES NOT SCANNING/INTERPOL/FFT'
	   PRINT *, 'PLEASE CHECK YOUR FILE:',FILEIN
	   STOP
	ENDIF
	IF (IFILTR.NE.0) THEN
	   PRINT *, 'MONORTM DOES NOT ANY FILTERING'
	   PRINT *, 'PLEASE CHECK YOUR FILE:',FILEIN
	   STOP
	ENDIF
	IF (ILAS.GT.0) THEN
	   WRITE(*,*) '----------------------------------------'
	   WRITE(*,*) 'WARNING: ILAS IS IGNORED IN MONORTM'
	   WRITE(*,*) 'IN ',FILEIN,' ILAS=',ILAS
	ENDIF
	IF (IOD.GT.0) THEN
	   WRITE(*,*) '----------------------------------------'
	   WRITE(*,*) 'WARNING: IOD IS IGNORED IN MONORTM'
	   WRITE(*,*) 'IN ',FILEIN,' IOD=',IOD
	ENDIF
	IF (IXSECT.NE.0) THEN
	   PRINT *, 'MONORTM DOES NOT ACCEPT CROSS SECTIONS INPUTS'
	   PRINT *, 'PLEASE CHECK YOUR FILE:',FILEIN
	   STOP
	ENDIF
	IF (MPTS.GT.0) THEN
	   WRITE(*,*) '----------------------------------------'
	   WRITE(*,*) 'WARNING: MPTS IS IGNORED IN MONORTM'
	   WRITE(*,*) 'IN ',FILEIN,' MPTS=',MPTS
	ENDIF
	IF (NPTS.GT.0) THEN
	   WRITE(*,*) '----------------------------------------'
	   WRITE(*,*) 'WARNING: NPTS IS IGNORED IN MONORTM'
	   WRITE(*,*) 'IN ',FILEIN,' NPTS=',NPTS
	ENDIF
	IF (CMRG(2).EQ.CA) THEN                               
	   IMRG = 12                                          
	ELSEIF (CMRG(2).EQ.CB) THEN                              
	   IMRG = 22                                             
	ELSEIF (CMRG(2).EQ.CC) THEN                           
	   IMRG = 32                                            
	ELSE                                                    
	   READ (CMRG(2),930) IMRG                              
	   IF (CMRG(1).EQ.CONE) IMRG = IMRG+10                    
	   IF (CMRG(1).EQ.CTWO) IMRG = IMRG+20                      
	   IF (CMRG(1).EQ.CTHREE) IMRG = IMRG+30               
	   IF (CMRG(1).EQ.CFOUR) IMRG = IMRG+40
	ENDIF                                                  
	!------END CHECKING RECORD 1.2
				!---record 1.2.1
	IF (IEMIT.EQ.2) THEN
	   READ(IRD,1010,ERR=6000) INFLAG,IOTFLG,JULDAT
	ENDIF
	IF (IEMIT.EQ.3) THEN
	   WRITE(*,*) 'CURRENTLY MONORTM DOES NOT HANDLE DERIVATIVES'
	   WRITE(*,*) 'CHECK YOUR INPUT FILE: ',FILEIN
	   STOP
	ENDIF	
				!---record 1.3
	IF ((IHIRAC+IAERSL+IEMIT+IATM+ILAS).GT.0) THEN   
	   READ (IRD,970,END=80,ERR=6000) V1,V2,SAMPLE,DVSET, 
	1	ALFAL0,AVMASS,DPTMIN,DPTFAC,ILNFLG,DVOUT 
	   !---CHECKING RECORD 1.3
	   IF (SAMPLE.GT.0) THEN
	      WRITE(*,*) '----------------------------------------'
	      WRITE(*,*) 'WARNING: SAMPLE IS IGNORED IN MONORTM'
	      WRITE(*,*) 'IN ',FILEIN,' SAMPLE=',SAMPLE
	   ENDIF
	   IF (ALFAL0.GT.0) THEN
	      WRITE(*,*) '----------------------------------------'
	      WRITE(*,*) 'WARNING: ALFAL0 IS IGNORED IN MONORTM'
	      WRITE(*,*) 'IN ',FILEIN,' ALFAL0=',ALFAL0
	   ENDIF
	   IF (AVMASS.GT.0) THEN
	      WRITE(*,*) '----------------------------------------'
	      WRITE(*,*) 'WARNING: AVMASS IS IGNORED IN MONORTM'
	      WRITE(*,*) 'IN ',FILEIN,' AVMASS=',AVMASS
	   ENDIF
	   IF (DPTMIN.GT.0) THEN
	      WRITE(*,*) '----------------------------------------'
	      WRITE(*,*) 'WARNING: DPTMIN IS IGNORED IN MONORTM'
	      WRITE(*,*) 'IN ',FILEIN,' DPTMIN=',DPTMIN
	   ENDIF
	   IF (DPTFAC.GT.0) THEN
	      WRITE(*,*) '----------------------------------------'
	      WRITE(*,*) 'WARNING: DPTFAC IS IGNORED IN MONORTM'
	      WRITE(*,*) 'IN ',FILEIN,' DPTFAC=',DPTFAC
	   ENDIF
	   IF (ILNFLG.GT.0) THEN
	      WRITE(*,*) 'STOP: ILNFLG MUST BE 0 FOR MONORTM'
	      WRITE(*,*) 'IN ',FILEIN,' ILNFLG=',ILNFLG
	      STOP
	   ENDIF
	   IF (DVOUT.GT.0) THEN
	      WRITE(*,*) '----------------------------------------'
	      WRITE(*,*) 'WARNING: DVOUT IS IGNORED IN MONORTM'
	      WRITE(*,*) 'IN ',FILEIN,' DVOUT=',DVOUT
	   ENDIF
	   IF ((DVSET.LE.0.).AND.(V1.NE.V2).AND.(V1.GT.0.).AND.
	1	(V2.GT.0.)) THEN
	      WRITE(*,*) 'MONORTM REQUIRES POSITIVE DVSET,'
	      WRITE(*,*) 'OR (V1=V2 AND DVSET=0) IF YOU WANT ONLY '
	      WRITE(*,*) 'ONE FREQ PROCESSED'
	      WRITE(*,*) 'CHECK YOUR INPUT FILE: ',FILEIN,V1,V2,DVSET
	      STOP
	   ENDIF
	   !---record 1.3.1
	   IF ((V1.LT.0.).OR.(V2.LT.0.)) THEN
	      READ(IRD,'(I8)',ERR=6000)NWN
	      IF (NWN.GT.NWNMX) THEN
		 WRITE(*,*) 'STOP: NUMBER OF WAVENUMBERS ',
	1	      'EXCEEDS LIMIT. ',NWN,NWNMX
		 WRITE(*,*) 'FIX: EXTEND NWNMX IN DECLAR.INCL'
		 STOP
	      ENDIF
	      !---record 1.3.2
	      DO IWN=1,NWN
		 READ(IRD,'(E19.7)',ERR=6000) WN(IWN)
	      ENDDO
	   ELSE
	      IF (DVSET.NE.0.) THEN 
		 NWN=NINT(((V2-V1)/DVSET)+1.)
		 IF (NWN.GT.NWNMX) THEN
		    WRITE(*,*) 'STOP: NUMBER OF WAVENUMBERS ',
	1		 'EXCEEDS LIMIT. ',NWN,NWNMX
		    WRITE(*,*) 'FIX: EXTEND NWNMX IN DECLAR.INCL'
		    STOP
		 ENDIF
		 DO J=1,NWN
		    WN(J)=V1+(J-1)*DVSET
		 ENDDO
	      ENDIF
	      IF (DVSET.EQ.0.) THEN 
		 IF (V1.NE.V2) THEN
		    PRINT *, 'AMBIGUITY IN THE WAVENUMBER',V1,V2,DVSET
		    STOP
		 ENDIF
		 IF (V1.EQ.V2) THEN
		    NWN=1
		    WN(1)=V1
		 ENDIF
	      ENDIF
	   ENDIF
	ENDIF                              
	!---END CHECKING RECORD 1.3
				!---record 1.4
	IEMIT=1
	IF (IEMIT.GT.0) THEN                                   
	   READ (IRD,970,END=80,ERR=6000) TMPBND,
	1	(BNDEMI(IBND),IBND=1,3),            
	1	(BNDRFL(IBND),IBND=1,3)                    
	ENDIF   
				!---record 1.4 (continued: manual emissivities)
	ICOEF = 13
	IF (BNDEMI(1).LT.0) THEN
	   OPEN (UNIT=ICOEF,FILE='in/EMISSION',STATUS='OLD',ERR=4000)
	   CALL READEM(ICOEF)
	   CLOSE (ICOEF)
	ELSE
	   XVMID = (V1+V2)/2.                                
	   EMITST = BNDEMI(1)+BNDEMI(2)*XVMID+BNDEMI(3)*XVMID*XVMID   
	   IF (EMITST.LT.0..OR.EMITST.GT.1.) THEN                     
	      STOP 'BNDEMI OUTSIDE PHYSICAL RANGE'                  
	   ENDIF                                                     
	ENDIF
				!---record 1.4 (continued: manual reflectivities)
	IF (BNDRFL(1).LT.0) THEN
	   OPEN (UNIT=ICOEF,FILE='in/REFLECTION',STATUS='OLD',ERR=5000)
	   CALL READRF(ICOEF)
	   CLOSE (ICOEF)
	ELSE
	   REFTST = BNDRFL(1)+BNDRFL(2)*XVMID+BNDRFL(3)*XVMID*XVMID   
	   IF (REFTST.LT.0..OR.REFTST.GT.1.) THEN                     
	      STOP 'BNDRFL OUTSIDE PHYSICAL RANGE'                    
	   ENDIF                                     
	ENDIF     
				!---record 1.6a
	IF (IMRG.GE.35) THEN
	   READ (IRD,945,ERR=6000) PATH1,LAYTOT
	   IF ((IMRG.GE.40).AND.(IEMIT.EQ.3)) THEN
	      READ (IRD,946) PTHODT
	   ENDIF
	ENDIF
				!---record 2.1
	IF (IATM.EQ.0) THEN 
	   READ (IRD,901,ERR=6000) IFORM,NLAYRS,NMOL,SECNT0,HEAD20,ZH1,
	1	HEAD4,ZH2,HEAD5,ANGLE,HEAD7    
				!---record 2.1.1
	   DO 30 L = 1, NLAYRS                                      
	      IF (L.EQ.1) THEN                                       
		 IF (IFORM.EQ.1) THEN
		    READ (IRD,910,ERR=6000) PAVE,TAVE,SECNTK,
	1		 CINP,IPTHRK,ALTZ(L-1),     
	1		 PZ(L-1),TZ(L-1),ALTZ(L),PZ(L),TZ(L),CLW(L)  
		 ELSE
		    READ (IRD,911,ERR=6000) PAVE,TAVE,SECNTK,
	1		 CINP,IPTHRK,ALTZ(L-1),
	1		 PZ(L-1),TZ(L-1),ALTZ(L),PZ(L),TZ(L),CLW(L) 
		 ENDIF
	      ELSE                                                   
		 IF (IFORM.EQ.1) THEN
		    READ (IRD,915,ERR=6000) PAVE,TAVE,SECNTK,
	1		 CINP,IPTHRK,ALTZ(L),PZ(L),TZ(L),CLW(L)     
		 ELSE
		    READ (IRD,916,ERR=6000) PAVE,TAVE,SECNTK,
	1		 CINP,IPTHRK,ALTZ(L),PZ(L),TZ(L),CLW(L) 
		 ENDIF
	      ENDIF                                                  
	      IF (TZ(L).EQ.0.) TZ(L) = TAVE
	      P(L) = PAVE                                         
	      T(L) = TAVE                                        
	      SECANT = SECNT0                                        
	      SECL(L) = SECANT                                        
	      SECNTA(L) = SECANT                                     
				!---record 2.1.2
	      IF (IFORM.EQ.1) THEN
		 READ (IRD,9255,ERR=6000) 
	1	      (WKL(M,L),M=1,7),WBRODL(L)
				!---record 2.1.3
		 IF (NMOL.GT.7) 
	1	      READ (IRD,9255,ERR=6000) (WKL(M,L),M=8,NMOL)   
	      ELSE
		 READ (IRD,927,ERR=6000) (WKL(M,L),M=1,7),WBRODL(L)
				!---record 2.1.3
		 IF (NMOL.GT.7) 
	1	      READ (IRD,927,ERR=6000) (WKL(M,L),M=8,NMOL)
	      ENDIF
				!---MIXING RATIO INPUT
	      WDNSTY = WBRODL(L)
	      WMXRAT = 0.0
	      DO 22 M = 2,NMOL
		 IF (WKL(M,L).GT.1) THEN
		    WDNSTY = WDNSTY + WKL(M,L)
		 ELSE
		    WMXRAT = WMXRAT + WKL(M,L)
		 ENDIF
 22	      CONTINUE
	      IF (WDNSTY.EQ.0.0) THEN
		 WRITE(*,920) L
		 STOP 'WDNSTY ERROR IN PATH'
	      ENDIF
	      IF (WMXRAT.LT.1.0) THEN
		 WDRAIR(L) = WDNSTY/(1.0-WMXRAT)
	      ELSE
		 WRITE(*,1000) L,WMXRAT, WDNSTY
		 STOP 'WMXRAT ERROR IN PATH'
	      ENDIF
	      DO 25 M = 1,NMOL
		 IF (WKL(M,L).LT.1) WKL(M,L) = WKL(M,L)*WDRAIR(L)
 25	      CONTINUE
 30	   CONTINUE        
	ENDIF 
	IF (INP.EQ.3) RETURN !In case inp=3 we need only the wave# info.
	IF (IATM.EQ.1) CALL LBLATM

	!---assignment of output variables
	IBMAXOUT=IBMAX
	H1Fout=H1F
	H2Fout=H2F
	do i=1,ibmax
	   ZBNDOUT(i)=ZBND(i)
	enddo
	IF ((INP.EQ.2).AND.(H1F.GE.H2F)) STOP 'INP=2 -> H1<H2'
	!---CHECKING RECORD 2.1
	IF (NLAYRS.GT.200) THEN
	   WRITE(*,*) '----------------------------------------'
	   WRITE(*,*) 'WARNING: NLAYRS MUST BE LESS THAN 200'
	   WRITE(*,*) 'IN ',FILEIN,' NLAYRS=',NLAYRS
	   WRITE(*,*) 'NLAYRS IS SET INTERNALLY TO 200'
	   NLAYRS=200
	ENDIF
	IF (NMOL.GT.35) THEN
	   WRITE(*,*) '----------------------------------------'
	   WRITE(*,*) 'WARNING: NMOL MUST BE LESS THAN 35'
	   WRITE(*,*) 'IN ',FILEIN,' NMOL=',NMOL
	   WRITE(*,*) 'NMOL IS SET INTERNALLY TO 35'
	   NMOL=35
	ENDIF
	IF (NMOL.EQ.0) NMOL = 7                             
	IF ((SECNT0.NE.1.).AND.(SECNT0.NE.-1.)) THEN
	   WRITE(*,*) '----------------------------------------'
	   WRITE(*,*) 'WARNING: SECNT0 MUST BE EITHER 1 OR -1'
	   WRITE(*,*) 'IN ',FILEIN,' SECNT0=',SECNT0
	ENDIF
	IF (ANGLE.GT.90.) IRT = 1 !space-based observer (looking down) 
	IF (ANGLE.LT.90.) IRT = 3 !ground-based observer (looking up)
	IF (ANGLE.EQ.90.) IRT = 2 !limb measurements
	!---END OF RECORD 2.1 CHECKING

 10	FORMAT (i4,5f19.13)
 901	FORMAT (1X,I1,I3,I5,F10.2,A20,F8.2,A4,F8.2,A5,F8.3,A7)
 905	FORMAT (A80)                                                   
 910	FORMAT (E15.7,F10.4,F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2),E15.7)
 911	FORMAT (3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2),E15.7)           
 915	FORMAT (E15.7,F10.4,F10.4,A3,I2,23X,(F7.2,F8.3,F7.2),E15.7)
 916	FORMAT (3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2),E15.7)         
 920	FORMAT (I3)                                                 
 925	FORMAT (10(4X,I1),3X,2A1,3(4X,I1),I1,I4,1X,I4,1X,I4,1X,I4)    
 9255	FORMAT (8E15.7)
 927	FORMAT (8E10.3)                                          
 930	FORMAT (I1)                                               
 945	FORMAT (A55,1X,I4)
 946	FORMAT (A55)
 970	FORMAT (8E10.3,4X,I1,5x,e10.3)                           
 990	FORMAT (F20.8)                                   
 1000	FORMAT ('Layer',I2,': Changing molecule ',I2,' from ',E10.3,
	1    ' to 1.000E+20.')
 1010	FORMAT (2I5,2X,I3)
 1015	FORMAT (I5)
	RETURN
 80	WRITE(*,*) ' EXIT; EOF ON :',FILEIN             
	STOP
 4000	WRITE(*,*) ' EXIT; ERROR OPENING EMISSION FILE'          
	STOP
 5000	WRITE(*,*) ' EXIT; ERROR OPENING REFLECTION FILE'          
	STOP
 6000	WRITE(*,*) ' EXIT; ERROR READING :',FILEIN          
	STOP
	END
	
	
	FUNCTION EMISFN (VI)                         
	IMPLICIT REAL*8           (V)                        
        !  FUNCTION EMISFN CALCULATES BOUNDARY EMISSIVITY FOR WAVE NUMBER      
        !  VALUE CORRESPONDING TO       
	PARAMETER (NMAXCO=4040)
	COMMON /EMSFIN/ V1EMIS,V2EMIS,DVEMIS,NLIMEM,ZEMIS(NMAXCO)
	COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP       
	EQUIVALENCE (BNDEMI(1),A) , (BNDEMI(2),B) , (BNDEMI(3),C)     
	!---Check for A < 0->use inputs in file "EMISSION"
	IF (A.LT.0.) THEN
	   NELMNT = INT((VI-V1EMIS)/DVEMIS)
	   IF ((NELMNT.LE.0).OR.(NELMNT.GE.NLIMEM)) THEN
	      WRITE(*,*) 'Frequency range of calculation exceeded',
	1	   ' emissivity input.'
	      WRITE(*,*) ' VI = ',VI,' V1EMIS = ',V1EMIS,' V2EMIS = ',
	1	   V2EMIS
	      STOP 'ERROR IN EMISFN'
	   ENDIF
	   V1A = V1EMIS+DVEMIS*NELMNT
	   V1B = V1EMIS+DVEMIS*(NELMNT+1)
	   CALL LINTCO(V1A,ZEMIS(NELMNT),V1B,ZEMIS(NELMNT+1),VI,ZINT,
	1	ZDEL)
	   EMISFN = ZINT
	   RETURN
	ENDIF
	IF (B.EQ.0..AND.C.EQ.0.) THEN                               
	   EMISFN = A                                               
	   RETURN                                                
	ENDIF                                                   
	XVI = VI                                                 
	EMISFN = A+B*XVI+C*XVI*XVI    
	RETURN                                                   
	END  
                                                               
	FUNCTION REFLFN (VI)                         
	IMPLICIT REAL*8           (V)                           
	!---FUNCTION REFLFN CALCULATES BOUNDARY REFLECTIVITY FOR WAVE NUMBER    
	!   VALUE CORRESPONDING TO VI       
	PARAMETER (NMAXCO=4040)
	COMMON /RFLTIN/ V1RFLT,V2RFLT,DVRFLT,NLIMRF,ZRFLT(NMAXCO)
	COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP         
	EQUIVALENCE (BNDRFL(1),A) , (BNDRFL(2),B) , (BNDRFL(3),C)     
	!---Check for A < 0->use values in from file"REFLECTION"
	IF (A.LT.0.) THEN
	   NELMNT = INT((VI-V1RFLT)/DVRFLT)
	   IF ((NELMNT.LE.0).OR.(NELMNT.GE.NLIMRF)) THEN
	      WRITE(*,*) 'Frequency range of calculation exceeded',
	1	   ' reflectivity input.'
	      WRITE(*,*) ' VI = ',VI,' V1RFLT = ',V1RFLT,' V2RFLT = ',
	1	   V2RFLT
	      STOP 'ERROR IN REFLFN'
	   ENDIF
	   V1A = V1RFLT+DVRFLT*NELMNT
	   V1B = V1RFLT+DVRFLT*(NELMNT+1)
	   CALL LINTCO(V1A,ZRFLT(NELMNT),V1B,ZRFLT(NELMNT+1),VI,ZINT,
	1	ZDEL)
	   REFLFN = ZINT
	   RETURN
	ENDIF
	IF (B.EQ.0..AND.C.EQ.0.) THEN                         
	   REFLFN = A                                             
	   RETURN                                                 
	ENDIF                                                     
	XVI = VI                                                  
	REFLFN = A+B*XVI+C*XVI*XVI                                 
	RETURN                                                     
	END                                                        
	
	SUBROUTINE LINTCO(V1,Z1,V2,Z2,VINT,ZINT,ZDEL)
	!Linearly interpolates emission and reflection values which
	!are directly read in from ASCII files
	IMPLICIT REAL*8           (V)
	ZDEL = (Z2-Z1)/(V2-V1)	!--ZDEL is the slope of the line
	ZCEPT = Z1 - ZDEL*V1	!--ZCEPT is the intercept for V = 0.0
	ZINT = ZDEL*VINT + ZCEPT !--Calculate ZINT value at VINT
	RETURN
        END


	SUBROUTINE EMISS_REFLEC(NWN,EMISS,REFLC,WN)
	IMPLICIT NONE
	include "declar.incl"
	INTEGER NWN,J
	REAL REFLFN,EMISFN
	DO J=1,NWN
	   REFLC(J)=REFLFN(WN(J))
	   EMISS(J)=EMISFN(WN(J))
	ENDDO
	RETURN 
	END


	SUBROUTINE STOREOUT(NWN,WN,RAD,TB,TRTOT,SCLCPL,SCLHW,NR,
	1    WVCOLMN,M,XSLF,K,Y0RES,N,INP,NPR,ilaunchdate,
	1    ilaunchtime,ibasetime,iserialnumber,isondeage,
	2    CLWCOLMN,TMPSFC,REFLC,EMISS,O,OL_WV,OS_WV,OF_WV,
	3    OL_O2,OL_O3,OL_N2,OC_N2,OL_N2O,OL_CO,OL_SO2,
	4    OL_NO2,OL_OH,O_CLW,NLAY,P,FILEOUT,ANGLE)
	include "declar.incl"
	INTEGER I,NWN,NR,K,M,N,NLAY,J
	REAL OTOT,OTOT_WV,OTOT_O2,OTOT_N2,OTOT_O3,OTOT_N2O
	REAL FREQ,OTOT_CO,OTOT_SO2,OTOT_NO2,OTOT_OH,ANGLE
	INTEGER INP,NPR,ilaunchdate,ilaunchtime,ibasetime,
	1    iserialnumber,isondeage
	REAL SCLCPL,SCLHW,Y0RES,CLWCOLMN,TMPSFC,WVCOLMN,XSLF
	CHARACTER FILEOUT*60
	COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
	1    RADCN1,RADCN2 
	DO I=1,NWN
	   FREQ=WN(I)*CLIGHT/1.E9
	   !----Computation of the integrated optical depths
	   OTOT=0.
	   OTOT_WV=0.
	   OTOT_O2=0.
	   OTOT_N2=0.
	   OTOT_O3=0.
	   OTOT_N2O=0.
	   OTOT_CO=0.
	   OTOT_SO2=0.
	   OTOT_NO2=0.
	   OTOT_OH=0.
	   DO J=1,NLAY
	      OTOT=OTOT+O(I,J)
	      OTOT_WV=OTOT_WV+OL_WV(I,J)+OS_WV(I,J)+OF_WV(I,J)
	      OTOT_O2=OTOT_O2+OL_O2(I,J)
	      OTOT_N2=OTOT_N2+OL_N2(I,J)+OC_N2(I,J)
	      OTOT_O3=OTOT_O3+OL_O3(I,J)
	      OTOT_N2O=OTOT_N2O+OL_N2O(I,J)
	      OTOT_CO=OTOT_CO+OL_CO(I,J)
	      OTOT_SO2=OTOT_SO2+OL_SO2(I,J)
	      OTOT_NO2=OTOT_NO2+OL_NO2(I,J)
	      OTOT_OH=OTOT_OH+OL_OH(I,J)
	   ENDDO
	   WRITE(1,21) NPR,NR,i,FREQ,TB(I),RAD(I),TRTOT(I),
	2	WVCOLMN,CLWCOLMN,TMPSFC,REFLC(I),EMISS(I),
	3	ANGLE,
	2	OTOT,OTOT_WV,OTOT_O2,OTOT_N2,OTOT_O3,
	3	OTOT_N2O,OTOT_CO,OTOT_SO2,OTOT_NO2,OTOT_OH
	ENDDO
 21	format(3i5,2f9.3,E19.9,f9.5,2f8.4,3f8.2,f9.3,10E12.4)
	RETURN
 1000	WRITE(*,*) 'ERROR OPENING FILE:',FILEOUT
	STOP
	END


	SUBROUTINE CORR_OPTDEPTH(INP,NLAY,SECNTA,NWN,ANGLE,O,IRT)
	include "declar.incl"
	INTEGER INP,J,NLAY,NWN,I,IRT
	REAL PI,SECNT,ALPHA,ANGLE
	COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
	1    RADCN1,RADCN2 
	!----SANITY CHECK
	IF (IRT.EQ.3) alpha=(angle*PI)/180.
	IF ((IRT.EQ.1).and.(ANGLE.GT.90.)) alpha=((180.-angle)*PI)/180.
	IF ((IRT.EQ.1).and.(ANGLE.LE.90.)) then
	   WRITE(*,*) 'IRT AND ANGLE INCOMPATIBLE '
	   STOP
	ENDIF
	SECNT=1./(cos(alpha))
	IF (SECNT.NE.1.) THEN
	   DO J=1,NLAY
	      DO I=1,NWN
		 O(I,J)=O(I,J)*SECNT
	      ENDDO
	   ENDDO
	ENDIF
	RETURN
	END


	SUBROUTINE INTEGR(W,CLW,NLAYRS,WVCOLMN,CLWCOLMN)
	IMPLICIT NONE
	include "declar.incl"
	REAL W(MXLAY)
	REAL CLWCOLMN,WVCOLMN
	INTEGER I,NLAYRS
	WVCOLMN=0.
	CLWCOLMN=0.
	DO I=1,NLAYRS
	   WVCOLMN=WVCOLMN+W(I)
	   CLWCOLMN=CLWCOLMN+CLW(I)
	ENDDO 
	WVCOLMN=WVCOLMN*(18./6.022e+23)
	RETURN
	END




	SUBROUTINE ARM2LBLATM(filearm,IFLAG,ilaunchdate,
	1    ilaunchtime,ibasetime,iserialnumber,isondeage,
	2    NWN,WN,V1,V2,DVSET,FILEIN,NLAYRS,IBMAX,ZBND,
	3    angle,H1F,H2F,NMOL,IPUNCH)
	!This subroutine simply reads the file 'filearm'
	!containing the ARM sonde data (extracted as follow:
	!> > > ftp to vernlaw.er.anl.gov
	!> > > log in as user anonymous
	!> > > enter your e-mail address as a password
	!> > > cd to directory /pub/arm/sonde/YYYY/corr/asc or
	!> > /pub/arm/sonde/YYYY/corr/cdf
	!> > > get the files)
	!and puts them into the file MONORTM.IN to be read by
	!LBLATM and then put in common for MONORTM.
	!---Output
	!IF IFLAG=0 profile is valid
	!IF IFLAG=1 profile valid but levls with same press removed
	!IF IFLAG=2 profile is rejected
	!Sid Ahmed Boukabara
	include "declar.incl"	
	INTEGER ibmax
	PARAMETER (MXLEV=3400)
	REAL PRESS(MXLEV),TEMP(MXLEV),HUMID(MXLEV),ALTIT(MXLEV)
	real ZBND(IBMAX)
	REAL*8 V1,V2
	REAL DVSET
	character filearm*90,filein*60,ligne*80,hmod*60
	CHARACTER*1 JCHARP,JCHART,JCHAR(35)
	INTEGER ilaunchdate,ilaunchtime,iserialnumber,ibasetime
	INTEGER isondeage
	data 	JCHAR /'H','6','6','6','6','6','6','6','6','6','6','6',
	1    '6','6','6','6','6','6','6','6','6','6','6','6','6','6',
	2    '6','6','6','6','6','6','6','6','6'/

	open(33,file=filein,status='unknown',form='formatted',err=1000)
	open(44,file=filearm,status='old',form='formatted',err=2000)
	read(44,10) ilaunchdate
	read(44,10) ilaunchtime
	read(44,10) ibasetime
	read(44,10) iserialnumber
	read(44,10) isondeage
	read(44,11) PWVorig
	read(44,11) PWVcorr
	read(44,10) nlevels
	read(44,'(a)') ligne
	nlev=0
	Pr_old=2000.
	alt_old=-2000.
	IFLAG=0
	IF ((nlevels.le.5).or.(nlevels.gt.MXLEV)) THEN  !minimum levels number is 5
	   close(44)
	   close(33)
	   IFLAG=2
	   RETURN
	ENDIF

	DO 230 i=1,nlevels
	   read(44,12) iTOff,Pr,Tp,RHorig,RHcor,xLat,xLon,iAlt
	   !---To avoid negative values 
	   IF (Pr.lt.0..or.RHcor.lt.0.) THEN
	      GOTO 230
	   ENDIF
	   !---to avoid two levels with the same pressure/altitude
	   IF (Pr.lt.Pr_old .and. iAlt/1000. .gt. alt_old) THEN
	      nlev=nlev+1
	      PRESS(nlev)=Pr
	      TEMP(nlev)=Tp
	      HUMID(nlev)=RHcor
	      ALTIT(nlev)=iAlt/1000.
	      alt_old=ALTIT(nlev)
	   ENDIF
	   IF (Pr.ge.Pr_old) IFLAG=1
	   Pr_old=Pr
 230	ENDDO
	IF ((nlev.le.5).or.(nlev.gt.MXLEV)) THEN  !minimum levels number is 5
	   close(44)
	   close(33)
	   IFLAG=2
	   RETURN
	ENDIF
 100	IHIRAC=1
	ILBLF4=0
	ICNTNM=1
	IAERSL=0
	IEMIT=1
	ISCAN=0
	IFILTR=0
	IPLOT=1
	ITEST=0
	IATM=1
	ILAS= 0   
	IOD=0
	IXSECT=0
	MPTS=0
	NPTS=0 
	INP=2
	TMPBND=2.75
	model=0
	itype=2
	nozero=1
	noprnt=0
	immax=nlev
	H1=max(H1F,ALTIT(1))
	H2=min(H2F,ALTIT(IMMAX))
	IBMAXSELECT=ibmax
	IBMINSELECT=1
	DO i=1,ibmax
	   IF (ZBND(i).GT.ALTIT(IMMAX)) GOTO 123
	   IBMAXSELECT=I
	   IF (ZBND(i).LT.H1) IBMINSELECT=I
	ENDDO
 123	CONTINUE
	JCHARP='A'
	JCHART='B'
	write(hmod,'(i12.12)') iserialnumber
	write(33,'(a)') '$ Rundeck'
	write(33,13) IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,             
	1    ISCAN,IFILTR,IPLOT,ITEST,IATM,ILAS,      
	2    IOD,IXSECT,MPTS,NPTS,INP 
	write(33,14) V1,V2,0.,DVSET,0.,0.,0.,0.,0,0. 
	IF ((V1.LT.0.).OR.(V2.LT.0.)) THEN
	   write(33,'(I8)') NWN
	   DO I=1,NWN
	      WRITE(33,'(E19.7)') WN(I)
	   ENDDO
	ENDIF
	write(33,14) TMPBND,1.,0.,0.,0.,0.,0.  
	write(33,15) model,itype,ibmaxselect-IBMINSELECT+1,nozero,
	1    noprnt,nmol,ipunch
	write(33,16) H1,H2,angle
	write(33,17) H1,(ZBND(i),i=IBMINSELECT+1,ibmaxselect)
	write(33,18) immax,hmod
	Do i=1,immax
	   write(33,19) ALTIT(i),PRESS(i),TEMP(i),JCHARP,JCHART,
	1	(JCHAR(j),j=1,nmol)
	   write(33,20) HUMID(i),(0.,j=1,nmol-1)
	enddo
	write(33,'(a)') '-1'
	write(33,'(a)') '%%%%%%%%%%%%%%%%'
	close(44)
	close(33)
 22	format(2f10.4,2f10.3,20x,I5)
 21	format(4f10.4)
 10	format(19x,i12)
 11	format(19x,f12.3)
 12	format(i8,f13.1,f10.1,f9.1,f11.1,f13.5,f12.5,i8)
 13	FORMAT (10(4X,I1),3X,2X ,3(4X,I1),1X,I4,1X,I4,1X,I4)    
 14	FORMAT (8E10.3,4X,I1,5x,e10.3)                         
 15	FORMAT (7I5)                                      
 16	FORMAT (3f10.3)                                      
 17	FORMAT (8f10.3)                                      
 18	FORMAT (I5,' Serial#',3A8)                        
 19	FORMAT (3f10.3,5x,2A1,3x,28A1)                          
 20	FORMAT (8f10.3)                                      
	RETURN
 1000	PRINT *, 'ERROR OPENING in ARM2LBLATM:',filein
	STOP
 2000	PRINT *, 'ERROR OPENING in ARM2LBLATM:',filearm
	STOP
	END

	SUBROUTINE START(nprof,INP)
	CHARACTER*15 HVRMON
	COMMON /CVRMON  / HVRMON
	
	WRITE(*,'(a)') '**********************************'
	WRITE(*,'(a)') '*        M O N O R T M           *'
	WRITE(*,'(a)') '*        '//HVRMON//'         *'
	WRITE(*,'(a)') '**********************************'
	WRITE(*,*) 'NUMBER OF PROFILES:',nprof
	IF (INP.EQ.1) WRITE(*,*) 'INPUTS FROM MONORTM.IN'
	IF (INP.EQ.2) WRITE(*,*) 'INPUTS FROM ARM.IN'
	IF (INP.EQ.3) WRITE(*,*) 'INPUTS FROM MONORTM_PROF.IN'
	WRITE(*,*)
	WRITE(*,*)
	WRITE(*,*)
	RETURN
	END


	SUBROUTINE GETPROFNUMBER(INP,FILEIN,fileARMlist,fileprof,
	1    NPROF,filearmTAB)
	include "declar.incl"
	CHARACTER FILEIN*60,filearm*90,CXID*1
	CHARACTER CDOL*1,CPRCNT*1,CXID*1,fileprof*80
	CHARACTER fileARMlist*64
	CHARACTER*8      HMOD                           
	DATA CDOL / '$'/,CPRCNT / '%'/
	NPROF=0

	!---Get first the INP info
	OPEN (530,FILE=FILEIN,STATUS='OLD',ERR=1000) 
 40	READ (530,'(a1)',END=80) CXID
	IF (CXID.NE.CDOL) THEN
	   GO TO 40     
	ENDIF
	READ (530,'(81X,I4)',END=80,ERR=6000) INP  
	close(530)

	IF (INP.EQ.1) THEN
	   OPEN (530,FILE=FILEIN,STATUS='OLD',ERR=1000) 
 20	   READ (530,'(a1)',END=80) CXID  
	   IF (CXID.EQ.CDOL) NPROF=NPROF+1
	   GO TO 20
	ENDIF
	IF (INP.EQ.2) THEN
	   open(530,file=fileARMlist,status='old',
	1	form='formatted',err=1000)
	   DO WHILE (.true.)
	      read(530,'(a)',end=80,err=1000) filearm
	      NPROF=NPROF+1
	      IF (NPROF.GT.NPROFMX) STOP' ERROR: NPROF>NPROFMX'
	      filearmTAB(NPROF)=filearm
	   ENDDO
	ENDIF
	IF (INP.EQ.3) THEN
	   open(530,file=fileprof,status='old',
	1	form='formatted',err=1000)
	   DO WHILE (.true.) 
	      READ (530,972,end=80,ERR=22) IFORM,LMAX,NMOL,SECNT0,HMOD,
	1	   HMOD,H1,H2,ANGLE,LEN  
	      NPROF=NPROF+1
 22	      CONTINUE
	   ENDDO
	ENDIF
 80	IF (NPROF.EQ.0) THEN 
	   WRITE(*,*) 'NO PROFILE FOUND IN GETPROFNUMBER'
	   STOP
	ELSE
	   CLOSE(530)
	   RETURN
	ENDIF
 1000	WRITE(*,*) 'ERROR OPENING OR READING FILE in GETPROFNUMBER'
 924	FORMAT (1X,I1,I3,I5,F10.6,A24) 
 972	FORMAT(1X,I1,I3,I5,F10.6,2A8,4X,F8.2,4X,F8.2,5X,F8.3,5X,I2) 
 6000	WRITE(*,*) ' EXIT; ERROR READING :',FILEIN          
	STOP
	END

	SUBROUTINE CHECKINPUTS(NWN,NPROF,NWNMX,NPROFMX,INP)
	INTEGER NWN,NPROF,NWNMX,NPROFMX
	IF (NWN.GT.NWNMX) THEN
	   WRITE(*,*) 'Number of wavenumbers too big:',NWN
	   WRITE(*,*) 'Maximum Allowed:',NWNMX
	   WRITE(*,*) 'Must extend NWNMX in declar.incl AND RECOMPILE'
	   STOP
	ENDIF
	IF (NPROF.GT.NPROFMX .AND. INP.EQ.2) THEN
	   WRITE(*,*) 'Number of profiles too big:',NPROF
	   WRITE(*,*) 'Maximum Allowed:',NPROFMX
	   WRITE(*,*) 'Must extend NPROFMX  AND RECOMPILE'
	   STOP
	ENDIF
	RETURN
	END

	SUBROUTINE EXPINT (X,X1,X2,A)                                   
C********************************************************************
C       THIS SUBROUTINE EXPONENTIALLY INTERPOLATES X1 AND X2 TO X BY  
C       THE FACTOR A             . NEEDED BY LBLATM                   
C********************************************************************
	IF (X1.EQ.0.0.OR.X2.EQ.0.0) GO TO 10                          
	X = X1*(X2/X1)**A                                        
	RETURN                                                   
 10	X = X1+(X2-X1)*A                                         
	RETURN                                                   
	END                                                      

	!----TEST -----------------------------
	!  TO MAINTAIN COMPATIBILTY WITH LBLATM
	!  WE ADD THE FOLLOWING DUMMY ARGUMENTS
	!--------------------------------------
	SUBROUTINE LBLDAT(HDATE)  
        CHARACTER*8      HDATE
	RETURN
	END

        SUBROUTINE FTIME (HTIME)
	CHARACTER*8      HTIME   
	RETURN
	END

	SUBROUTINE XSREAD (XV1,XV2)   	
	RETURN
	END
