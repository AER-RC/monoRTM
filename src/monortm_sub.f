

	SUBROUTINE RTM(IOUT,IRT,NWN,WN,NLAY,TT,TZ,
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
C     AIM: This program is aimed at the calculation of the 
C          Radiative Transfer components. 
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
C     - PP    : Vector of NLAY pressures (in mbar), one should note that
C               this input could be a scalar (associated with NLAY=1)
C     - TT    : Vector of NLAY temperatures (layers averaged T) [in Kelvin]
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
C     - TRTOT  : Total Transmission
C
C       Note:
C       -----
C       RTM takes into account the cosmic background contribution.
C       The cosmic radiation is hard coded (2.7 Kelvin). 
C
C-------------------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER NWN,NLAY,IRT,I,IOUT,IDU
	REAL RADCN1,RADCN2
	PARAMETER (RADCN2=1.438786314,RADCN1=1.191061999E-12)
	REAL*8 WN(NWN),V
	REAL TT(NLAY),TZ(0:NLAY),O(NWN,NLAY),bb_fn,X,EMISS(NWN)
	REAL REFLC(NWN)
	REAL RUP(NWN),TRTOT(NWN),RDN(NWN),fbeta,beta,RAD(NWN)
	REAL TMPSFC,ESFC,RSFC,SURFRAD,ALPH,COSMOS,TSKY,TB(NWN)
	BB_fn(V,fbeta)  = RADCN1*(V**3)/(EXP(V*fbeta)-1.)               
	!---Up and Down radiances
	CALL RAD_UP_DN(TT,NLAY,TZ,WN,rup,trtot,rdn,O,NWN,IDU)
	!---RADIATIVE TRANSFER
	TSKY=2.75 !Cosmic background in Kelvin
	beta= RADCN2/TMPSFC
	alph= RADCN2/TSKY
	DO I=1,NWN
	   SURFRAD  = bb_fn(WN(I),beta)
	   COSMOS   = bb_fn(WN(I),alph)
	   ESFC=EMISS(I)
	   RSFC=REFLC(I)
	   IF (IRT.EQ.1) THEN
	      RAD(I)=RUP(I)+((trtot(i)**2)*COSMOS)+
	1	   trtot(i)*(rsfc*rdn(i)+esfc*SURFRAD)
	   ENDIF
	   IF (IRT.EQ.3) RAD(I)=RDN(I)+(trtot(i)*COSMOS)
	   IF (IRT.EQ.2) RAD(I)=RUP(I)+RDN(I)*trtot(i)+
	1	((trtot(i)**2)*COSMOS)
	   IF (IOUT.EQ.1) THEN
	      X=RADCN1*(WN(I)**3)/RAD(I)+1.
	      TB(I)=RADCN2*WN(I)/alog(X)
	   ENDIF
	   !print *, '****',I,TB(I),WN(I),RUP(I),RDN(I),trtot(i),esfc,TMPSFC
	ENDDO
 	RETURN
	END 

	SUBROUTINE RAD_UP_DN(TAVEL,nlayer,TZ,VI,rup,trtot,rdn,OD,NWN,IDU)
	IMPLICIT REAL*8 (V)      
	PARAMETER (RADCN2=1.438786314,RADCN1=1.191061999E-12,
	1    aa_inv=3.597122302)
	INTEGER  layer,nlayer,NWN,IDU,lmin,lmax,nl
	REAL          beta,beta_a,bb,bba,bbVEC(nlayer),bbaVEC(nlayer)
	DIMENSION     TAVEL(nlayer),TZ(0:nlayer),OD(NWN,nlayer),VI(NWN)
	DIMENSION     RUP(NWN),TRTOT(NWN),rdn(NWN),ODTOT(NWN)
	BB_fn(V,fbeta)  = RADCN1*(V**3)/(EXP(V*fbeta)-1.)
	lmin=1
	lmax=nlayer
	nl=1
	IF (IDU.EQ.1) then 
	   lmin=nlayer
	   lmax=1
	   nl=-1
	ENDIF
	DO 60 I=1,NWN
	   VV=VI(I)
	   rup(I) = 0.
	   rdn(I) = 0.
	   trtot(I) = 1.
	   ODTOT(I)=0.
	   DO layer=lmax,lmin,-nl
	      ODTOT(I)=ODTOT(I)+OD(I,layer)
	      beta  = radcn2/tavel(layer)
	      beta_a= radcn2/tz(layer)
	      bbVEC(layer)  = bb_fn(VV,beta)
	      bbaVEC(layer) = bb_fn(VV,beta_a)
	   ENDDO
	   ODT=ODTOT(I)
	   DO 70 layer = lmax,lmin,-nl
	      bb  = bbVEC(layer)
	      bba = bbaVEC(layer)
	      ODVI = OD(I,layer)
	      TRI = EXP(-ODVI)
	      ODT=ODT-ODVI
	      TRTOT(I)= EXP(-ODT)
	      pade=0.193*ODVI+0.013*ODVI**2
	      RUP(I)= RUP(I)+TRTOT(I)*(1.-TRI)*(bb+pade*bba)/(1.+pade)
 70	   ENDDO
	   ODT=ODTOT(I)
	   do 50 layer = lmin,lmax,nl
	      bb  = bbVEC(layer)
	      bba = bbaVEC(layer-1)
	      ODVI = OD(I,layer)
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
	READ (ICOEF,900) V1EMIS,V2EMIS,DVEMIS,NLIMEM
	DO 100 NGNU = 1,NLIMEM
	   READ (ICOEF,910) ZEMIS(NGNU)
 100	CONTINUE
	RETURN
 900	FORMAT (3E10.3,5X,I5)
 910	FORMAT (E15.7)
	END

	SUBROUTINE READRF(ICOEF)
	IMPLICIT REAL*8           (V)
	PARAMETER (NMAXCO=4040)
	COMMON /RFLTIN/ V1RFLT,V2RFLT,DVRFLT,NLIMRF,ZRFLT(NMAXCO)
	READ (ICOEF,900) V1RFLT,V2RFLT,DVRFLT,NLIMRF
	DO 100 NGNU = 1,NLIMRF
	   READ (ICOEF,910) ZRFLT(NGNU)
 100	CONTINUE
	RETURN
 900	FORMAT (3E10.3,5X,I5)
 910	FORMAT (E15.7)
	END



	SUBROUTINE RDLBLINP(IATM,JPLOT,IRT)
        !This subroutine reads some control parameters
        !from the input file (former TAPE5) and makes
        !available some of them via common blocks.
        !This subroutine has been added to make 
        !MONORTM compatible with LBLRTM inputs.
        !S.A. Boukabara AER INC. 1999
	!IRT=1  !=1->Space/=3->ground/=2->limb
	!Important: RDLBLINP supposes that the profile
	!of IPTH contains the same value (from which it
	!determines the IRT value.
	INTEGER NLAY
	PARAMETER (NLAY=203)
	REAL*8           V1,V2 
	character*4 ht1,ht2
	CHARACTER*1 CMRG(2),CONE,CTWO,CTHREE,CFOUR
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
	INTEGER IPTHRK,IPTH(NLAY),ITYL(NLAY),IPATHL,M
	REAL PZ(0:NLAY),TAVEL(NLAY),TZ(0:NLAY)
	REAL WKL(35,NLAY),PAVEL(NLAY),WBRODL(NLAY),
	1    WTOTL(NLAY),ALBL(NLAY),ADBL(NLAY),AVBL(NLAY),
	2    H2OSL(NLAY),SECNTA(NLAY),ALTZ(0:NLAY)
	REAL SECL(64),WDNSTY,WMXRAT,WDRAIR(NLAY),DVL(NLAY)
	INTEGER IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,        
	1    NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,       
	2    NLTEFL,LNFIL4,LNGTH4
	character*8 XID,HMOLID,YID      
	real*8      SECANT,XALTZ 
	COMMON /PATHD/ PAVEL,TAVEL,WKL,WBRODL,DVL,WTOTL,ALBL,
	1    ADBL,AVBL,H2OSL,IPTH,ITYL,SECNTA,HT1,HT2,ALTZ,         
	2    PZ,TZ                          
	COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,        
	1    NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,       
	2    NLTEFL,LNFIL4,LNGTH4                                 
	COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,  
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
				!--OUTPUT FILE 
	FILEOUT='out/LBLATM.LOG'
	IPR=66
	OPEN (IPR,FILE=FILEOUT,STATUS='UNKNOWN')        
				!---INPUT FILE 
	FILEIN='in/LBLRTM.IN'
	IRD=55
	OPEN (IRD,FILE=FILEIN,STATUS='UNKNOWN')        
				!---record 1.1
 20	READ (IRD,905,END=80) CXID  
	IF (CXID.EQ.CPRCNT) THEN
	   WRITE(*,*) 'ERROR STOP -END OF FILE:',FILEIN
	   STOP 
	ENDIF                       
	IF (CXID.NE.CDOL) GO TO 20                            
				!---record 1.2
				!print *, 'rec 1.2'
	READ(IRD,925,END=80) IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,             
	1    ISCAN,IFILTR,IPLOT,ITEST,IATM,CMRG,ILAS,      
	2    IOD,IXSECT,IRAD,MPTS,NPTS  
	IF (IAERSL.GT.0) THEN
	   WRITE(*,*) 'CURRENTLY MONORTM DOES NOT HANDLE AEROSOLS'
	   WRITE(*,*) 'CHECK YOUR INPUT FILE: ',FILEIN
	   STOP
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
				!---record 1.2.1
	IF (IEMIT.EQ.2) THEN
				!print *, 'rec 1.2.1'
	   READ(IRD,1010) INFLAG,IOTFLG,JULDAT
	ENDIF
	IF (IEMIT.EQ.3) THEN
	   WRITE(*,*) 'CURRENTLY MONORTM DOES NOT HANDLE DERIVATIVES'
	   WRITE(*,*) 'CHECK YOUR INPUT FILE: ',FILEIN
	   STOP
	ENDIF
	
				!---record 1.3
	IF ((IHIRAC+IAERSL+IEMIT+IATM+ILAS).GT.0) THEN   
				!print *, 'rec 1.3'
	   READ (IRD,970,END=80) V1,V2,SAMPLE,DVSET,ALFAL0,AVMASS,DPTMIN, 
	1	DPTFAC,ILNFLG,DVOUT 
	   IF ((DVSET.LE.0.).AND.(V1.NE.V2)) THEN
	      WRITE(*,*) 'MONORTM REQUIRES A POSITIVE DVSET VALUE'
	      WRITE(*,*) 'CHECK YOUR INPUT FILE: ',FILEIN
	      STOP
	   ENDIF
	ENDIF                              
				!---record 1.4
	IF (IEMIT.GT.0) THEN                                                
				!print *, 'rec 1.4'
	   READ (IRD,970,END=80) TMPBND,(BNDEMI(IBND),IBND=1,3),            
	1	(BNDRFL(IBND),IBND=1,3)                    
	ENDIF   
				!---record 1.4 (continued: manual emissivities)
	ICOEF = 13
	IF (BNDEMI(1).LT.0) THEN
	   OPEN (UNIT=ICOEF,FILE='in/EMISSION',STATUS='OLD')
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
	   OPEN (UNIT=ICOEF,FILE='in/REFLECTION',STATUS='OLD')
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
				!print *, 'rec 1.6.a'
	   READ (IRD,945) PATH1,LAYTOT
	   IF ((IMRG.GE.40).AND.(IEMIT.EQ.3)) THEN
	      READ (IRD,946) PTHODT
	   ENDIF
	ENDIF
				!---record 2.1
	IF (IATM.EQ.0) THEN 
	   READ (IRD,901) IFORM,NLAYRS,NMOL,SECNT0,HEAD20,ZH1,
	1	HEAD4,ZH2,HEAD5,ZANGLE,HEAD7    
	   IF (NMOL.EQ.0) NMOL = 7                                          
	   IF (SECNT0.LT.0.) THEN                                           
	      IPATHL = 1                                                    
	   ELSE                                                             
	      IPATHL = 3                                                    
	   ENDIF                                                            
				!---record 2.1.1
	   DO 30 L = 1, NLAYRS                                                 
	      IF (L.EQ.1) THEN                                                 
		 IF (IFORM.EQ.1) THEN
		    READ (IRD,910) PAVE,TAVE,SECNTK,CINP,IPTHRK,ALTZ(L-1),     
	1		 PZ(L-1),TZ(L-1),ALTZ(L),PZ(L),TZ(L)                   
		 ELSE
		    READ (IRD,911) PAVE,TAVE,SECNTK,CINP,IPTHRK,ALTZ(L-1),
	1		 PZ(L-1),TZ(L-1),ALTZ(L),PZ(L),TZ(L)
		 ENDIF
	      ELSE                                                             
		 IF (IFORM.EQ.1) THEN
		    READ (IRD,915) PAVE,TAVE,SECNTK,CINP,IPTHRK,
	1		 ALTZ(L),PZ(L),TZ(L)                                                 
		 ELSE
		    READ (IRD,916) PAVE,TAVE,SECNTK,CINP,IPTHRK,
	1		 ALTZ(L),PZ(L),TZ(L)
		 ENDIF
	      ENDIF                                                            
	      IF (TZ(L).EQ.0.) TZ(L) = TAVE
	      PAVEL(L) = PAVE                                                  
	      TAVEL(L) = TAVE                                                  
	      SECANT = SECNT0                                                  
	      IF (SECNTK.GT.0.) SECANT = SECNTK                                
	      SECL(L) = SECANT                                                 
	      SECNTA(L) = SECANT                                               
	      IF (IPTHRK.NE.0) IPATHL = IPTHRK                                 
	      IPTH(L) = IPATHL                                                 
	      IF (SECANT.EQ.0.) STOP 'RDLBL SECANT = 0'                        
				!---record 2.1.2
	      IF (IFORM.EQ.1) THEN
		 READ (IRD,9255) (WKL(M,L),M=1,7),WBRODL(L)
				!---record 2.1.3
		 IF (NMOL.GT.7) READ (IRD,9255) (WKL(M,L),M=8,NMOL)             
	      ELSE
		 READ (IRD,927) (WKL(M,L),M=1,7),WBRODL(L)
				!---record 2.1.3
		 IF (NMOL.GT.7) READ (IRD,927) (WKL(M,L),M=8,NMOL)
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
	   IF (IPTH(1).EQ.1) IRT=1 !space-based observer
	   IF (IPTH(1).EQ.3) IRT=3 !ground-based observer
	   IF (IPTH(1).EQ.2) IRT=2 !limb measurement
	   IF (IXSECT.NE.0) THEN
	      PRINT *, 'MONORTM DOES NOT ACCEPT CROSS SECTIONS INPUTS'
	      PRINT *, 'PLEASE CHECK YOUR FILE:',FILEIN
	      STOP
	   ENDIF
	ENDIF 
	IF (IATM.EQ.1) CALL LBLATM
	IF (IPTH(1).EQ.1) IRT=1 !space-based observer
	IF (IPTH(1).EQ.3) IRT=3 !ground-based observer
	IF (IPTH(1).EQ.2) IRT=2 !limb measurement
	IF (ISCAN.NE.0) THEN
	   PRINT *, 'MONORTM DOES NOT SCANNING/INTERPOL/FFT'
	   PRINT *, 'PLEASE CHECK YOUR FILE:',FILEIN
	   STOP
	ENDIF
	IF (IFILTR.NE.0) THEN
	   PRINT *, 'MONORTM DOES NOT FILTERING'
	   PRINT *, 'PLEASE CHECK YOUR FILE:',FILEIN
	   STOP
	ENDIF
	JPLOT=0
	IF (IPLOT.NE.0) THEN
	   READ (IRD,900) CPRGID,CEX            
	   IEXTRC=0                                                            
	   IF (CEX.EQ.CEXST) IEXTRC = 1                                        
	   READ (IRD,9055) XV1,XV2,XSIZ,DELV,NUMSBX,NOENDX,LFILE,         
	1	LSKIPF, SCALE,IOPT,I4P,IXDEC                                 
	   READ (IRD,9100) YMINR,YMAX,YSIZ,DELY,NUMSBY,NOENDY,IDEC,     
	1	JEMIT,JPLOT,LOGPLT,JHDR,JDUMMY,JOUT,JPLTFL    
	ENDIF
 10	FORMAT (i4,5f19.13)
 900	FORMAT (A60,18X,A2)                                                 
 901	FORMAT (1X,I1,I3,I5,F10.2,A20,F8.2,A4,F8.2,A5,F8.3,A7)
 905	FORMAT (A1)                                                        
 9055	FORMAT (4F10.4,4I5,F10.3,I2,I3,I5)                                  
 910	FORMAT (E15.7,F10.4,F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))
 9100	FORMAT (2G10.4,2G10.3,6I5,2(I2,I3))                                 
 911	FORMAT (3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))                          
 915	FORMAT (E15.7,F10.4,F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))
 916	FORMAT (3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))                          
 920	FORMAT (I3)                                                         
 925	FORMAT (10(4X,I1),3X,2A1,3(4X,I1),I1,I4,1X,I4)    
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
	CLOSE(IRD)
	CLOSE(IPR)
	RETURN
 80	WRITE(*,*) ' MONORTM EXIT; EOF ON :',FILEIN,' (fmr TAPE5)'             
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

	SUBROUTINE INSTRCONF(V1,V2,DVSET,NWN,WN)
	IMPLICIT NONE
	REAL*8 V1,V2,WN(NWN),WNMIN,WNMAX
	REAL DVSET,DWN
	INTEGER NWN,J
	NWN=INT((V2-V1)/DVSET)
	WNMIN=V1
	WNMAX=V2
	DWN=DVSET
	DO J=1,NWN
	   WN(J)=WNMIN+(J-1)*DWN
	ENDDO
	RETURN
	END

	SUBROUTINE EMISS_REFLEC(NWN,EMISS,REFLC,WN)
	IMPLICIT NONE
	INTEGER NWN,J
	REAL REFLC(NWN),EMISS(NWN),REFLFN,EMISFN
	REAL*8 WN(NWN)
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
	4    OL_NO2,OL_OH,O_CLW,NLAY,PRESS)
	IMPLICIT NONE
	INTEGER I,NWN,NR,K,M,N,NLAY,J
	REAL OTOT,OTOT_WV,OTOT_O2,OTOT_N2,OTOT_O3,OTOT_N2O
	REAL OTOT_CO,OTOT_SO2,OTOT_NO2,OTOT_OH
	INTEGER INP,NPR,ilaunchdate,ilaunchtime,ibasetime,
	1    iserialnumber,isondeage
	REAL*8 WN(NWN)
	REAL O(NWN,NLAY),OL_WV(NWN,NLAY),
	1    OS_WV(NWN,NLAY),OF_WV(NWN,NLAY),OL_O2(NWN,NLAY),
	2    OL_O3(NWN,NLAY),OL_N2(NWN,NLAY),OC_N2(NWN,NLAY),
	3    OL_N2O(NWN,NLAY),O_CLW(NWN,NLAY),OL_CO(NWN,NLAY),
	4    OL_SO2(NWN,NLAY),OL_NO2(NWN,NLAY),OL_OH(NWN,NLAY)
	REAL PRESS(NLAY)
	REAL TB(NWN),TRTOT(NWN),FREQ,RAD(NWN),WVCOLMN,XSLF
	REAL SCLCPL,SCLHW,Y0RES,CLWCOLMN
	REAL EMISS(NWN),REFLC(NWN),TMPSFC
	CHARACTER FILEOUT*60
	LOGICAL INIT
	DATA INIT/.TRUE./
	IF(INIT)THEN
	   FILEOUT='out/MONORTM.OUT'
	   OPEN(1,file=FILEOUT,status='unknown',form='formatted')
	   WRITE(1,'(a)') 'MONORTM RESULTS:'
	   WRITE(1,'(a)') ' PROF#,CHANNEL#,FREQ,TB,TRANSMITT,RAD' 
	   INIT=.FALSE.
	ENDIF
	DO I=1,NWN
	   FREQ=WN(I)*29.98
	   IF (INP.EQ.1) THEN  !LBLRTM.IN input
	      WRITE(1,10) NPR,NR,i,FREQ,TB(I),TRTOT(I),RAD(I),SCLCPL,
	1	   SCLHW,WVCOLMN,M,XSLF,K,Y0RES,N,
	2	   CLWCOLMN,TMPSFC,REFLC(I),EMISS(I)
	   ENDIF
	   IF (INP.EQ.2) THEN  !ARM.IN input
	      WRITE(1,20) NPR,NR,i,FREQ,TB(I),TRTOT(I),RAD(I),SCLCPL,
	1	   SCLHW,WVCOLMN,M,XSLF,K,Y0RES,N,ilaunchdate,
	1	   ilaunchtime,ibasetime,iserialnumber,isondeage,
	2	   CLWCOLMN
	   ENDIF
	   IF (INP.EQ.3) THEN  !pre-stored ARM data input(MONRTM_PROF.IN)
c	      WRITE(1,20) NPR,NR,i,FREQ,TB(I),TRTOT(I),RAD(I),SCLCPL,
c	1	   SCLHW,WVCOLMN,M,XSLF,K,Y0RES,N,ilaunchdate,
c	1	   ilaunchtime,ibasetime,iserialnumber,isondeage,
c	2	   CLWCOLMN
	      !----specific output
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
	      !---------------------
	      WRITE(1,21) NPR,NR,i,FREQ,TB(I),TRTOT(I),WVCOLMN,
	2	   CLWCOLMN,TMPSFC,REFLC(I),EMISS(I),OTOT,
	2	   OTOT_WV,OTOT_O2,OTOT_N2,OTOT_O3,OTOT_N2O,
	3	   OTOT_CO,OTOT_SO2,OTOT_NO2,OTOT_OH
	   ENDIF
	ENDDO

 10	format(3i5,2f9.3,f9.5,E15.7,2f6.3,3(f9.4,i4),f7.3,3f8.2)
 20	format(3i5,2f9.3,f9.5,E15.7,2f6.3,f9.4,i4,f9.4,i4,f9.4,i4,
	1    2I7,2I10,I4,f7.3)
 21	format(3i5,2f9.3,f9.5,2f6.3,3f8.2,10E12.4)
	RETURN
	END


	SUBROUTINE CORR_OPTDEPTH(INP,NLAY,SECNTA,NWN,ANGLE,O)
	IMPLICIT NONE
	INTEGER INP,J,NLAY,NWN,I
	REAL PI,SECNTA(NLAY),O(NWN,NLAY),SECNT,ALPHA,ANGLE
	PI=3.14159265
	IF ((INP.EQ.1).OR.(INP.EQ.2).OR.(INP.EQ.3)) THEN
	   DO J=1,NLAY
	      DO I=1,NWN
		 O(I,J)=O(I,J)*SECNTA(J)
	      ENDDO
	   ENDDO
	ENDIF
	IF ((INP.EQ.0).AND.(ANGLE.NE.0.)) THEN
	   alpha=(angle*PI)/180.
	   SECNT=1./(cos(alpha))
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
	REAL W(NLAYRS),WVCOLMN
	REAL CLW(NLAYRS),CLWCOLMN
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
	1    ilaunchtime,ibasetime,iserialnumber,isondeage)
	!This subroutine simply reads the file 'filearm'
	!containing the ARM sonde data (extracted as follow:
	!> > > ftp to vernlaw.er.anl.gov
	!> > > log in as user anonymous
	!> > > enter your e-mail address as a password
	!> > > cd to directory /pub/arm/sonde/YYYY/corr/asc or
	!> > /pub/arm/sonde/YYYY/corr/cdf
	!> > > get the files)
	!and puts them into the file LBLRTM.IN to be read by
	!LBLATM and then put in common for MONORTM.
	!---Output
	!IF IFLAG=0 profile is valid
	!IF IFLAG=1 profile valid but levls with same press removed
	!IF IFLAG=2 profile is rejected
	!Sid Ahmed Boukabara
	INTEGER NLEVMAX,ibmax
	PARAMETER (NLEVMAX=10000,ibmax=32)
	REAL PRESS(NLEVMAX),TEMP(NLEVMAX),HUMID(NLEVMAX),
	1    ALTIT(NLEVMAX)
	real alti(ibmax)
	character filearm*80,filein*60,ligne*80,hmod*60
	CHARACTER*1 JCHARP,JCHART,JCHAR(7)
	INTEGER ilaunchdate,ilaunchtime,iserialnumber,ibasetime
	INTEGER isondeage
	data alti /0.320,0.374,0.500,0.750,1.000,1.500,2.000,2.500,
     &       3.000,3.500,4.000,4.500,5.000,5.500,6.000,6.500,
     &       7.000,7.500,8.000,8.500,9.000,9.500,10.000,11.000,
     &       12.000,13.000,14.000,15.000,16.500,18.000,19.5,21./
	data 	JCHAR /'H','6','6','6','6','6','6'/

	FILEIN='in/LBLRTM.IN'
	open(33,file=filein,status='unknown',form='formatted')
	open(44,file=filearm,status='old',form='formatted')
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
	ialtmax=21000
	df=-9999.
	IFLAG=0
	IF ((nlevels.le.5).or.(nlevels.gt.3400)) THEN  !minimum levels number is 5
	   IFLAG=2
	   RETURN
	ENDIF

	DO i=1,nlevels
	   read(44,12) iTOff,Pr,T,RHorig,RHcor,xLat,xLon,iAlt
	   !---To avoid -9999. below 21 kms (profile rejected)
	   IF ((Pr.eq.df.or.t.eq.df.or.RHcor.eq.df).and.(ialt.le.ialtmax)) THEN
	      IFLAG=2
	      RETURN
	   ENDIF
	   !---To avoid -9999. above 21 kms (profile kept)
	   IF ((Pr.eq.df.or.t.eq.df.or.RHcor.eq.df).and.(ialt.gt.ialtmax)) THEN
	      IFLAG=0
	      GOTO 100
	   ENDIF
	   !---to avoid two levels with the same pressure
	   IF (Pr.lt.Pr_old) THEN
	      nlev=nlev+1
	      PRESS(nlev)=Pr
	      TEMP(nlev)=T
	      HUMID(nlev)=RHcor
	      ALTIT(nlev)=iAlt/1000.
	   ENDIF
	   IF (Pr.ge.Pr_old) IFLAG=1
	   Pr_old=Pr
	ENDDO
 100	IHIRAC=1
	ILBLF4=1
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
	IRAD=0
	MPTS=0
	NPTS=0 
	TMPBND=2.75
	model=0
	itype=2
	nozero=1
	noprnt=1
	nmol=7
	ipunch=1
	H1=0.320
	H2=21.000
	angle=0.000
	immax=nlev
	JCHARP='A'
	JCHART='B'
	write(hmod,'(i12.12)') iserialnumber
	write(33,'(a)') '$ Rundeck'
	write(33,13) IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,             
	1    ISCAN,IFILTR,IPLOT,ITEST,IATM,ILAS,      
	2    IOD,IXSECT,IRAD,MPTS,NPTS  
	write(33,14) 0.,0.,0.,0.,0.,0.,0.,0.,0,0. 
	write(33,14) TMPBND,1.,0.,0.,0.,0.,0.           
	write(33,15) model,itype,ibmax,nozero,noprnt,nmol,ipunch
	write(33,16) H1,H2,angle
	write(33,17) (alti(i),i=1,ibmax)
	write(33,18) immax,hmod
	Do i=1,immax
	   write(33,19) ALTIT(i),PRESS(i),TEMP(i),JCHARP,JCHART,JCHAR
	   write(33,20) HUMID(i),0.,0.,0.,0.,0.,0.
	enddo
	write(33,'(a)') 'TEST ARM RADIOSONDE: INPUT FOR MONORTM'
	write(33,21) 0.1,10.,10.2000,100.00
	write(33,22) 200.,300.,7.02000,20.,1
	write(33,'(a)') '-1'
	write(33,'(a)') '%%%%%%%%%%%%%%%%'
	close(44)
	close(33)
 22	format(2f10.4,2f10.3,20x,I5)
 21	format(4f10.4)
 10	format(19x,i12)
 11	format(19x,f12.3)
 12	format(i8,f13.1,f10.1,f9.1,f11.1,f13.5,f12.5,i8)
 13	FORMAT (10(4X,I1),3X,2X,3(4X,I1),I1,1X,I4,1X,I4)    
 14	FORMAT (8E10.3,4X,I1,5x,e10.3)                                      
 15	FORMAT (7I5)                                      
 16	FORMAT (3f10.3)                                      
 17	FORMAT (8f10.3)                                      
 18	FORMAT (I5,' Serial#',3A8)                                      
 19	FORMAT (3f10.3,5x,2A1,3x,28A1)                                      
 20	FORMAT (8f10.3)                                      
	RETURN
	END

