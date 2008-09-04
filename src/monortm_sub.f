C     path:		$Source$
C     author:		$Author $
C     revision:	        $Revision$
C     created:	        $Date$

	SUBROUTINE RTM(IOUT,IRT,NWN,WN,NLAY,T,TZ,
     &    TMPSFC,  RUP,TRTOT,RDN,REFLC,EMISS,RAD,TB,IDU)
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
C       The cosmic radiation is hard coded (2.75 Kelvin). 
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
     1       RADCN1,RADCN2 
	COMMON /CVRSUB/ HVRSUB
	BB_fn(V,fbeta)  = RADCN1*(V**3)/(EXP(V*fbeta)-1.)   
	HVRSUB = '$Revision$' 

	!---Up and Down radiances

	CALL RAD_UP_DN(T,NLAY,TZ,WN,rup,trtot,rdn,  NWN,IDU,IRT)

	!---RADIATIVE TRANSFER
	TSKY=2.75 !Cosmic background in Kelvin
	beta= RADCN2/TMPSFC
	alph= RADCN2/TSKY

	if (irt.eq.3) then
	   print *, '     '
           print *, '***********************************'
	   print *,
     *        'NB: for Downwelling Radiance the Boundary is ',
     *        'Internally Set to the Cosmic Value: 2.75K'

           print *, '***********************************'

	endif

	DO I=1,NWN
	   SURFRAD  = bb_fn(WN(I),beta)
	   COSMOS   = bb_fn(WN(I),alph)
	   ESFC=EMISS(I)
	   RSFC=REFLC(I)
c
c     Upwelling Case
	   IF (IRT.EQ.1) RAD(I) = RUP(I) +
     1       trtot(i) * (esfc*SURFRAD + rsfc*(rdn(i)+trtot(i)*COSMOS)) ! kcp 09/21/07
c
c     Limb Case      trtot is taken as the transmittance from the tangent point to h1 (SAC)
	   IF (IRT.EQ.2) RAD(I) = RUP(I) +
     1       trtot(i) * (rdn(i)+trtot(i)*COSMOS)
c
c     Downwelling Case
	   IF (IRT.EQ.3) RAD(I)=RDN(I)+(trtot(i)*COSMOS)
c
	   IF (IOUT.EQ.1) THEN
	      X=RADCN1*(WN(I)**3)/RAD(I)+1.
	      TB(I)=RADCN2*WN(I)/log(X)
	   ENDIF
	ENDDO
 	RETURN
	END 


	SUBROUTINE RAD_UP_DN(T,nlayer,TZ,WN,rup,trtot,rdn,  NWN,
     1       IDU,IRT)
	IMPLICIT REAL*8 (V)      
	include "declar.incl"
	COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     1       RADCN1,RADCN2 
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
     1                TRTOT(I)*(1.-TRI)*(bb+pade*bba)/(1.+pade)
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
     1       FILEIN,ICNTNM,IXSECT,IBMAXOUT,ZBNDOUT,
     2       H1fout,H2fout,ISPD,IPASSATM)
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
     1       ISCAN,IFILTR,IPLOT,ITEST,IATM,ILAS,ILNFLG,      
     2       IOD,IXSECT,MPTS,NPTS,INFLAG,IOTFLG,JULDAT
	REAL SAMPLE,DVSET,ALFAL0,AVMASS,DPTMIN,DPTFAC,DVOUT 
	REAL TMPBND,XVMID,EMITST,REFTST
	INTEGER IBPROP,IBND,ICOEF,IMRG,LAYTOT,IFORM,NLAYRS,NMOL
	REAL PATH1,PTHODT,SECNT0,ZH1,ZH2,ZANGLE,PAVE,TAVE,SECNTK
	CHARACTER HEAD20*20,HEAD7*7,HEAD5*5,HEAD4*4,CINP*3
	CHARACTER*60 FILEOUT,FILEIN
	character*1 hmol_scal
	REAL SECL(64),WDNSTY,WMXRAT,WDRAIR(MXLAY)
	REAL ZBNDOUT(MXFSC)
	INTEGER IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,        
     1       NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,       
     2       NLTEFL,LNFIL4,LNGTH4,IPTHRK,IPATHL,M
	character*8 XID,HMOLID,YID      
	COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,NOZERO,NOP,H1F,H2F, 
     1       ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,        
     2       XVBAR, HMINF,PHIF,IERRF,HSPACE                     
 	COMMON /BNDRY/ ZBND(MXFSC),PBND(MXFSC),TBND(MXFSC),
     1       ALORNZ(MXFSC),ADOPP(MXFSC),AVOIGT(MXFSC)       
	COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,
     1       IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,             
     2       IPHMID,IPDIM,KDIM,KMXNOM,KMAX                      
	COMMON /PATHD/ P,T,WKL,WBRODL,DVL,WTOTL,ALBL,
     1       ADBL,AVBL,H2OSL,IPTH,ITYL,SECNTA,HT1,HT2,ALTZ,         
     2       PZ,TZ                          
	common /profil_scal/ nmol_scal,hmol_scal(64),xmol_scal(64)

        COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(mx_xs),XAMNT(mx_xs,MXLAY)

	COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,  
     1       NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,       
     2       NLTEFL,LNFIL4,LNGTH4                                 
	COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,
     1       AVBAR,  
     1       AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      
     2       DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,      
     3       ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,     
     4       EXTID(10)    
	COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP            
	COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),  
     1       WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   
     2       EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    
        COMMON /CNTSCL/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL

	DATA CDOL / '$'/,CPRCNT / '%'/
	DATA CONE / '1'/,CTWO / '2'/,CTHREE / '3'/,CFOUR / '4'/,       
     1       CA / 'A'/,CB / 'B'/,CC / 'C'/                   
	character*6 idcntl(15)
	DATA IDCNTL / ' HIRAC',' LBLF4',' CNTNM',' AERSL',' EMISS', 
     *                ' SCNFN',' FILTR','  PLOT','  TEST','  IATM',    
     *                'CMRG_1','CMRG_2','  ILAS','  ISPD',' XSECT'/
	CHARACTER CEX*2,CEXST*2,CPRGID*60    
	DATA CEXST/'EX'/
	EQUIVALENCE (CXID,CXIDLINE)                    
	EQUIVALENCE (FSCDID(3),IXSCNT)
				
 20	READ (IRD,905,END=80,ERR=6000) CXIDLINE         	!---record 1.1
	IF (CXID.EQ.CPRCNT) THEN
	   WRITE(*,*) '-END OF FILE:',FILEIN
	   RETURN 
	ENDIF                       
	IF (CXID.NE.CDOL) GO TO 20                            
	READ (CXIDLINE,'(1x,10A8)') (XID(I),I=1,10)     

	READ(IRD,925,END=80,ERR=6000) IHIRAC,ILBLF4,    	!---record 1.2           
     1       ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,    
     2       ITEST,IATM, CMRG,ILAS, IOD,IXSECT,
     3       MPTS,NPTS,ISPD  

	IXSCNT = IXSECT*10 + ICNTNM

	WRITE (IPR,935) (IDCNTL(I),I=1,15)  
 935	FORMAT (15(A6,3X)) 

 	Write(ipr,940)                 IHIRAC,ILBLF4,          	!---record 1.2           
     1       ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,    
     2       ITEST,IATM,CMRG(1),CMRG(2),ILAS,
     3       ISPD,IXSECT
	
 940	FORMAT (1X,I4,9I9,2(8x,a1),3I9)

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
	   PRINT *, 'MONORTM DOES NOT HANDLE SCANNING/INTERPOL/FFT'
	   PRINT *, 'PLEASE CHECK YOUR FILE:',FILEIN
	   STOP
	ENDIF
	IF (IFILTR.NE.0) THEN
	   PRINT *, 'MONORTM DOES NOT HANDLE ANY FILTERING'
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


      IF (ICNTNM.EQ.0) THEN
         XSELF = 0.0
         XFRGN = 0.0
         XCO2C = 0.0
         XO3CN = 0.0
         XO2CN = 0.0
         XN2CN = 0.0
         XRAYL = 0.0
      ELSEIF (ICNTNM.EQ.1) THEN
         XSELF = 1.0
         XFRGN = 1.0
         XCO2C = 1.0
         XO3CN = 1.0
         XO2CN = 1.0
         XN2CN = 1.0
         XRAYL = 1.0
      ELSEIF (ICNTNM.EQ.2) THEN
         XSELF = 0.0
         XFRGN = 1.0
         XCO2C = 1.0
         XO3CN = 1.0
         XO2CN = 1.0
         XN2CN = 1.0
         XRAYL = 1.0
         ICNTNM = 1
      ELSEIF (ICNTNM.EQ.3) THEN
         XSELF = 1.0
         XFRGN = 0.0
         XCO2C = 1.0
         XO3CN = 1.0
         XO2CN = 1.0
         XN2CN = 1.0
         XRAYL = 1.0
         ICNTNM = 1
      ELSEIF (ICNTNM.EQ.4) THEN
         XSELF = 0.0
         XFRGN = 0.0
         XCO2C = 1.0
         XO3CN = 1.0
         XO2CN = 1.0
         XN2CN = 1.0
         XRAYL = 1.0
         ICNTNM = 1
      ELSEIF (ICNTNM.EQ.5) THEN
         XSELF = 1.0
         XFRGN = 1.0
         XCO2C = 1.0
         XO3CN = 1.0
         XO2CN = 1.0
         XN2CN = 1.0
         XRAYL = 0.0
         ICNTNM = 1
      ELSEIF (ICNTNM.EQ.6) THEN
         READ(IRD,*) XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, XRAYL   !---record 1.2a
         ICNTNM = 1
      ENDIF   

			
	IF (IEMIT.EQ.2) THEN
	   READ(IRD,1010,ERR=6000) INFLAG,IOTFLG,JULDAT         	!---record 1.2.1
	ENDIF
	IF (IEMIT.EQ.3) THEN
	   WRITE(*,*) 'CURRENTLY MONORTM DOES NOT HANDLE DERIVATIVES'
	   WRITE(*,*) 'CHECK YOUR INPUT FILE: ',FILEIN
	   STOP
	ENDIF	
			
	IF ((IHIRAC+IAERSL+IEMIT+IATM+ILAS).GT.0) THEN                	!---record 1.3
	   READ (IRD,970,END=80,ERR=6000) V1,V2,SAMPLE,DVSET, 
     1        ALFAL0,AVMASS,DPTMIN,DPTFAC,ILNFLG,DVOUT,nmol_scal

c       read in the profile scaling parameters
c
	   if (nmol_scal .gt. 0 ) then
	      if (nmol_scal .gt. 38) stop ' nmol_scal .gt. 38 '
	      read (ird,9701) (hmol_scal(m),m=1,nmol_scal)
	      read (ird,9702) (xmol_scal(m),m=1,nmol_scal)
 9701	      FORMAT (64a1)
 9702	      FORMAT (7e15.7,/,(8e15.7,/))
	   endif

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
     1          (V2.GT.0.)) THEN
	      WRITE(*,*) 'MONORTM REQUIRES POSITIVE DVSET,'
	      WRITE(*,*) 'OR (V1=V2 AND DVSET=0) IF YOU WANT ONLY '
	      WRITE(*,*) 'ONE FREQ PROCESSED'
	      WRITE(*,*) 'CHECK YOUR INPUT FILE: ',FILEIN,V1,V2,DVSET
	      STOP
	   ENDIF


	   IF ((V1.LT.0.).OR.(V2.LT.0.)) THEN                       !---record 1.3.1
	      READ(IRD,'(I8)',ERR=6000)NWN
	      IF (NWN.GT.NWNMX) THEN
		 WRITE(*,*) 'STOP: NUMBER OF WAVENUMBERS ',
     1		      'EXCEEDS LIMIT. ',NWN,NWNMX
		 WRITE(*,*) 'FIX: EXTEND NWNMX IN DECLAR.INCL'
		 STOP
	      ENDIF
	     
	      DO IWN=1,NWN
		 READ(IRD,'(E19.7)',ERR=6000) WN(IWN)               !---record 1.3.2
	      ENDDO
              dvset=0.
	   ELSE
	      IF (DVSET.NE.0.) THEN 
		 NWN=NINT(((V2-V1)/DVSET)+1.)
		 IF (NWN.GT.NWNMX) THEN
		    WRITE(*,*) 'STOP: NUMBER OF WAVENUMBERS ',
     1                   'EXCEEDS LIMIT. ',NWN,NWNMX
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
     1          (BNDEMI(IBND),IBND=1,3),            
     1          (BNDRFL(IBND),IBND=1,3)                    
	   WRITE (IPR,985)                TMPBND,
     *          (BNDEMI(IBND),IBND=1,3), 
     *          (BNDRFL(IBND),IBND=1,3)          ! surf_refl
 985  FORMAT (5(/),'0*********** BOUNDARY PROPERTIES ***********',/,      A07510
     *        '0 TBOUND   = ',F12.4,5X,'BOUNDARY EMISSIVITY   = ',        A07530
     *        3(1PE11.3),/,'0',29X,'BOUNDARY REFLECTIVITY = ',            A07540
     *        3(1PE11.3),/,'0',29X,' SURFACE REFLECTIVITY = ', A1)
C
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
	IF (IATM.EQ.0)  return
c
C  IPASSATM=0 means this is the first call to RDLBLINP
C  For the first call to this routine, we don't need to go any further
	IF (IPASSATM.EQ.0) THEN 
	    IPASSATM = 1
	    RETURN
	ENDIF

	IF (IATM.EQ.1) then

c          clw is not returned from lblatm	   
	   do l=1,nlayrs
	      clw(l) = 0.
	   enddo
c	   
	   CALL LBLATM
c
	ENDIF

	!---assignment of output variables
	IBMAXOUT=IBMAX
	H1Fout=H1F
	H2Fout=H2F
	do i=1,ibmax
	   ZBNDOUT(i)=ZBND(i)
	enddo
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
 901	FORMAT (1X,I1,I3,I5,F10.6,A20,F8.2,A4,F8.2,A5,F8.3,A7)
 905	FORMAT (A80)                                                   
 910	FORMAT (E15.7,F10.4,F10.6,A3,I2,1X,2(F7.2,F8.3,F7.2),E15.7)
 911	FORMAT (2F10.4,f10.6,A3,I2,1X,2(F7.2,F8.3,F7.2),E15.7)           
 915	FORMAT (E15.7,F10.4,F10.6,A3,I2,23X,(F7.2,F8.3,F7.2),E15.7)
 916	FORMAT (2F10.4,f10.6,A3,I2,23X,(F7.2,F8.3,F7.2),E15.7)         
 920	FORMAT (I3)                                                 
 925	FORMAT (10(4X,I1),3X,2A1,3(4X,I1),1X,I4,1X,I4,6X,I4)    
 9255	FORMAT (8E15.7)
 927	FORMAT (8E10.3)                                          
 930	FORMAT (I1)                                               
 945	FORMAT (A55,1X,I4)
 946	FORMAT (A55)
 970	FORMAT (8E10.3,4X,I1,5x,e10.3,i5)                           
 990	FORMAT (F20.8)                                   
 1000	FORMAT ('Layer',I2,': Changing molecule ',I2,' from ',E10.3,
     1       ' to 1.000E+20.')
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
     1             ' emissivity input.'
	      WRITE(*,*) ' VI = ',VI,' V1EMIS = ',V1EMIS,' V2EMIS = ',
     1             V2EMIS
	      STOP 'ERROR IN EMISFN'
	   ENDIF
	   V1A = V1EMIS+DVEMIS*NELMNT
	   V1B = V1EMIS+DVEMIS*(NELMNT+1)
	   CALL LINTCO(V1A,ZEMIS(NELMNT),V1B,ZEMIS(NELMNT+1),VI,ZINT,
     1	        ZDEL)
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
     1             ' reflectivity input.'
	      WRITE(*,*) ' VI = ',VI,' V1RFLT = ',V1RFLT,' V2RFLT = ',
     1             V2RFLT
	      STOP 'ERROR IN REFLFN'
	   ENDIF
	   V1A = V1RFLT+DVRFLT*NELMNT
	   V1B = V1RFLT+DVRFLT*(NELMNT+1)
	   CALL LINTCO(V1A,ZRFLT(NELMNT),V1B,ZRFLT(NELMNT+1),VI,ZINT,
     1          ZDEL)
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
	real*8            wn(*)
c
	INTEGER NWN,J
	REAL REFLFN,EMISFN,reflc(*),emiss(*)
	DO J=1,NWN
	   REFLC(J)=REFLFN(WN(J))
	   EMISS(J)=EMISFN(WN(J))
	ENDDO
	RETURN 
	END


	SUBROUTINE STOREOUT(NWN,WN,WKL,WBRODL,RAD,TB,TRTOT,
     1       NPR,WVCOLMN,CLWCOLMN,TMPSFC,REFLC,EMISS,
     4       NLAY,NMOL,ANGLE,IOT,FILEOUT)
	include "declar.incl"

	INTEGER I,J,NWN,NLAY,NMOL,NPR
	REAL FREQ,ANGLE
	REAL CLWCOLMN,TMPSFC,WVCOLMN
 	CHARACTER FILEOUT*60
        character*8 hmolc(mxmol),cmol(mxmol)
        logical giga

	REAL OTOT,otot_by_mol(mxmol),wk_tot(mxmol),odxtot
        integer id_mol(mxmol)

	COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     1       RADCN1,RADCN2 

        save cmol, id_mol, kount

      DATA HMOLC / '  H2O   ' , '  CO2   ' , '   O3   ' , '  N2O   ' ,   FA12260
     *             '   CO   ' , '  CH4   ' , '   O2   ' , '   NO   ' ,   FA12270
     *             '  SO2   ' , '  NO2   ' , '  NH3   ' , ' HNO3   ' ,   FA12280
     *             '   OH   ' , '   HF   ' , '  HCL   ' , '  HBR   ' ,   FA12290
     *             '   HI   ' , '  CLO   ' , '  OCS   ' , ' H2CO   ' ,   FA12300
     *             ' HOCL   ' , '   N2   ' , '  HCN   ' , ' CH3CL  ' ,   FA12310
     *             ' H2O2   ' , ' C2H2   ' , ' C2H6   ' , '  PH3   ' ,   FA12320
     *             ' COF2   ' , '  SF6   ' , '  H2S   ' , ' HCOOH  ' ,   FA12330
     *             '  HO2   ' , '   O+   ' , ' ClONO2 ' , '   NO+  ' ,
     *             '  HOBr  ' , ' C2H4   ' , ' CH3OH  '/

    ! set up headers: assumes same molecules used in all profiles!!!

        if (npr.eq.1) then
           if (nmol.lt.22) wkl(22,:) = wbrodl
           wk_tot = sum(wkl,2)
           kount = 0
           do im=1,mxmol
              if(wk_tot(im).gt.0) then
                kount=kount+1
                id_mol(kount) = im
                cmol(kount) = hmolc(im)
              end if
           end do
        end if

	!---Write header in output file
	WRITE(IOT,'(a)') 'MONORTM RESULTS:'
	WRITE(IOT,'(a)') '----------------' 
	WRITE(IOT,'(a5,I8,90x,a42)') 'NWN :',NWN ,
     &           'Molecular Optical Depths -->'

           
        write(iot,11)'PROF','FREQ','BT(K) ','  RAD(W/cm2_ster_cm-1)',
     &                    'TRANS','PWV',
     &                   'CLW','SFCT','EMIS','REFL','ANGLE',
     &                    'TOTAL_OD',cmol(1:kount), 'XSEC_OD' 
       ! Convert to GHz for small wavenumbers
        if (wn(1).lt.100) giga = .true.
	DO Iw=1,NWN
	   if (giga) then
              FREQ=WN(Iw)*CLIGHT/1.E9
           else
              FREQ=WN(Iw)
           end if
             
	   !----Computation of the integrated optical depths
	   OTOT=sum(o(iw,:))
	   odxtot = sum(odxsec(iw,:))
           do ik=1,kount
              otot_by_mol(ik) = sum(o_by_mol(iw,id_mol(ik),:))+
     &                          sum(oc(iw,id_mol(ik),:))
           end do

          WRITE(IOT,21) NPR,FREQ,TB(Iw),RAD(Iw),TRTOT(Iw),
     2          WVCOLMN,CLWCOLMN,TMPSFC,EMISS(Iw),REFLC(Iw),
     3          ANGLE,OTOT,otot_by_mol(1:kount), odxtot
	ENDDO
 11	format (a5,a9,a11,a22,a8,2a8,3a8,a9,36a12)
 21	format (i5,f9.3,f11.5,1p,E21.9,0p,f9.5,2f8.4,3f8.2,f9.3,
     1                                                  1p,36E12.4)
	RETURN
 1000	WRITE(*,*) 'ERROR OPENING FILE:',FILEOUT
	STOP
	END


	SUBROUTINE CORR_OPTDEPTH(INP,NLAY,SECNTA,NWN,ANGLE,  IRT)
	include "declar.incl"
	INTEGER INP,J,NLAY,NWN,I,IRT
	REAL PI,SECNT,ALPHA,ANGLE
	COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     1       RADCN1,RADCN2 
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
	REAL W(*),clw(*)
	REAL CLWCOLMN,WVCOLMN
	INTEGER I,NLAYRS
	WVCOLMN=0.
	CLWCOLMN=0.
	DO I=1,NLAYRS
	   WVCOLMN=WVCOLMN+W(I)
	   CLWCOLMN=CLWCOLMN+CLW(I)
	ENDDO 
c       value from vpayne 2006/07/23
	WVCOLMN=WVCOLMN*(2.99150e-23)
	RETURN
	END




	SUBROUTINE START(nprof,IATM)
	CHARACTER*15 HVRREL
	COMMON /CVRREL  / HVRREL

	
	WRITE(*,'(a40)') '**********************************'
	WRITE(*,'(a40)') '*        M O N O R T M           *'
	WRITE(*,'(a40)') '*        '//HVRREL//'         *'
	WRITE(*,'(a40)') '**********************************'
	WRITE(*,*) 'NUMBER OF PROFILES:',nprof
	IF (IATM.EQ.1) WRITE(*,*) 'INPUTS FROM MONORTM.IN'
	IF (IATM.EQ.0) WRITE(*,*) 'INPUTS FROM MONORTM_PROF.IN'
	WRITE(*,*)
	WRITE(*,*)
	WRITE(*,*)
	RETURN
	END


	SUBROUTINE GETPROFNUMBER(IATM,FILEIN,fileARMlist,fileprof,
     1       NPROF,filearmTAB)
	include "declar.incl"
	CHARACTER FILEIN*60,filearm*90,CXID*1
	CHARACTER CDOL*1,CPRCNT*1,fileprof*80
	CHARACTER fileARMlist*64
	CHARACTER*8      HMOD                           
	DATA CDOL / '$'/,CPRCNT / '%'/
	NPROF=0

	!---Get first the IATM info
	OPEN (90,FILE=FILEIN,STATUS='OLD',ERR=1000) 
 40	READ (90,'(a1)',END=80) CXID
	IF (CXID.NE.CDOL) THEN
	   GO TO 40     
	ENDIF
	READ (90,'(49X,I1,19x,I1)',END=80,ERR=6000) IATM,ixsect 
	close(90)

	IF (IATM.EQ.1) THEN
	   OPEN (90,FILE=FILEIN,STATUS='OLD',ERR=1000) 
 20	   READ (90,'(a1)',END=80) CXID  
	   IF (CXID.EQ.CDOL) NPROF=NPROF+1
	   GO TO 20
	ENDIF
	IF (IATM.EQ.0) THEN
	   open(90,file=fileprof,status='old',
     1          form='formatted',err=1000)
	   DO WHILE (.true.) 
	       READ (90,'(a1)') CXID 
	       IF (CXID.EQ.CPRCNT) THEN
		   GO TO 70 
	       ELSE
		   BACKSPACE(UNIT=90) 
		   READ (90,972,end=70,ERR=22) IFORM,
     1               LMAX,NMOL,SECNT0,HMOD,
     1	             HMOD,H1,H2,ANGLE,LEN  
		   NPROF=NPROF+1
  22		   CONTINUE
	       ENDIF
	   ENDDO
C  If there are cross-sections, there will be two records consistent with Format
C statement 972 for each profile.
  70	   continue
	   if (ixsect. eq. 1) nprof = nprof/2
	ENDIF
	
 80	IF (NPROF.EQ.0) THEN 
	   WRITE(*,*) 'NO PROFILE FOUND IN GETPROFNUMBER'
	   STOP
	ELSE
	   CLOSE(90)
	   RETURN
	ENDIF
 2000	WRITE(*,*) 'ERROR OPENING ARM FILE'
 1000	WRITE(*,*) 'ERROR OPENING OR READING FILE in GETPROFNUMBER'
 924	FORMAT (1X,I1,I3,I5,F10.6,A24) 
 972	FORMAT(1X,I1,I3,I5,F10.6,2A8,4X,F8.2,4X,F8.2,5X,F8.3,5X,I2) 
 6000	WRITE(*,*) ' EXIT; ERROR READING :',FILEIN          
	STOP
	END

	SUBROUTINE CHECKINPUTS(NWN,NPROF,NWNMX,NPROFMX)
	INTEGER NWN,NPROF,NWNMX,NPROFMX
	IF (NWN.GT.NWNMX) THEN
	   WRITE(*,*) 'Number of wavenumbers too big:',NWN
	   WRITE(*,*) 'Maximum Allowed:',NWNMX
	   WRITE(*,*) 'Must extend NWNMX in declar.incl AND RECOMPILE'
	   STOP
	ENDIF
	RETURN
	END

c******************************************************************************
c___________________________________________________________________
c___________________________________________________________________

	subroutine profil_scal_sub(nlayrs)

	include "declar.incl"

	REAL*8 V1,V2,SECANT,XALTZ 
	character*1 hmol_scal
	character*10 holn2
	character*8 XID,HMOLID,YID,HDATE,HTIME


	COMMON /PATHD/ P,T,WKL,WBRODL,DVL,WTOTL,ALBL,ADBL,AVBL,
     &     H2OSL,IPTH,ITYL,SECNTA,HT1,HT2,ALTZ,PZ,TZ
	common /profil_scal/ nmol_scal,hmol_scal(64),xmol_scal(64)

	COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4), 
     1    WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   
     2    EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    
	COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,  
     1    NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,       
     2    NLTEFL,LNFIL4,LNGTH4                                 

	DIMENSION WMT(64)

c  *** It should be noted that no attempt has been made to keep the 
c      mass in a given layer constant, i.e. level pressure nor retained ***

c      obtain accumulated amounts by molecule

	do m = 1, nmol
	   wmt(m) = 0.
	   do l = 1, nlayrs
	      wmt(m) = wmt(m) + wkl(m,l)
	   enddo
	enddo

	wsum_brod = 0.
	do l = 1, nlayrs
	   wsum_brod = wsum_brod + wbrodl(l)
	enddo

c        obtain dry air sum
c             check to see if nitrogen is included in the selected molecules

	if (nmol.ge.22) then
	   wsum_drair = 0.
	else
	   wsum_drair = wsum_brod
	endif

	do m = 2, nmol
	   wsum_drair = wsum_drair + wmt(m)
	enddo

	write (ipr,*)
	write (ipr,*) '   ',
     1          '*****************************************************'
	write (ipr,*)
	write (ipr,*) '               Profile Scaling          '  

	write (ipr,956) 
 956	format (/,4x,' molecule',2x, 
     1           'hmol_scale',3x, ' xmol_param ',3x, 'scale factor',/)

	do m = 1, nmol_scal

	   xmol_scal_m = xmol_scal(m)
	   if (hmol_scal(m).eq.' ') xmol_scal(m) = 1.
	   if (hmol_scal(m).eq.'0') xmol_scal(m) = 0.
	   if (hmol_scal(m).eq.'1') xmol_scal(m) = xmol_scal_m              ! Scale factor

	   if (hmol_scal(m).eq.'C' .or. hmol_scal(m).eq.'c')                ! Column Amount (molec/cm^2)
     1              xmol_scal(m) = xmol_scal_m/wmt(m)     

	   if (hmol_scal(m).eq.'M' .or. hmol_scal(m).eq.'m') then           ! Mixing ratio (molec/molec(dry_air))
	      if (wsum_drair.gt.0.) then
		 xmol_scal(m) = xmol_scal_m/(wmt(m)/wsum_drair)      
	      else
		 stop 'mixing ratio failure: wsum_drair = 0.'
	      endif
	   endif

	   if (hmol_scal(m).eq.'P' .or .hmol_scal(m).eq.'p') then           ! PWV for water vapor (cm)
	      if (m.eq.1) then 
c                value from vpayne 2006/07/23
		 xmol_scal(1) = (xmol_scal_m/2.99150e-23)/wmt(1)
	      else
		 write (ipr,*) 'm = ', m
		 stop ' (hmol_scal(m).eq."P" .and. m.ne.1) '
	      endif
	   endif

	   if (hmol_scal(m).eq.'D' .or. hmol_scal(m).eq.'d') ! Dobson Units (du)
     1              xmol_scal(m) =  (xmol_scal_m*2.68678e16)/wmt(m)

	   write (ipr,957) m, hmol_scal(m), xmol_scal_m, xmol_scal(m)
 957	   format (5x,i5,9x,a1,5x,1p, 4e15.7)

c                scale the amounts and recalculate the total

	   wmt(m) = 0.
	   do l = 1, nlayrs
	      wkl(m,l) = wkl(m,l) * xmol_scal(m)
	      wmt(m)   = wmt(m) +wkl(m,l)
	   enddo
	enddo

	write (ipr,*)
	write (ipr,*) '   ',
     1          '*****************************************************'
	write (ipr,*)

c       write  modified column amounts to ipr in lblatm format
C
	WRITE (IPR,970)    
C
C     --------------------------------------------------------------
C
C     Write out column densities for molecules to TAPE6
C
	iform = 1

	data holn2/'  OTHER  '/
C
	IF (IFORM.EQ.1) THEN                                             
	   WRITE (IPR,974) (HMOLID(I),I=1,7),HOLN2                       
	   DO L = 1, NLAYRS
	      WRITE (IPR,980) L,P(L),T(L),
     *                 (WKL(M,L),M=1,7),WBRODL(L)                            
	   enddo
	   IF (NLAYRS.GT.1) THEN                                         
	      WRITE (IPR,985)                                            
	      L = NLAYRS                                                 
	      WRITE (IPR,980) L,PWTD,TWTD,
     *                 (WMT(M),M=1,7),SUMN2    
	   ENDIF                                                       
	ELSE
	   WRITE (IPR,975) (HMOLID(I),I=1,7),HOLN2
	   DO  L = 1, NLAYRS
	      WRITE (IPR,982) L,P(L),T(L),
     *                 (WKL(M,L),M=1,7),WBRODL(L)
	   enddo
	   IF (NLAYRS.GT.1) THEN
	      WRITE (IPR,985)
	      L = NLAYRS
	      WRITE (IPR,991) L,PWTD,TWTD,
     *                  (WMT(M),M=1,7),SUMN2
	   ENDIF
	ENDIF
C
	IF (NMOL.GT.7) THEN                                            
	   DO MLO = 8, NMOL, 8                                         
	      MHI = MLO+7                                             
	      MHI = MIN(MHI,NMOL)                              
	      WRITE (IPR,970)                                  
	      IF (IFORM.EQ.1) THEN
		 WRITE (IPR,974) (HMOLID(I),I=MLO,MHI)            
		 DO L = 1, NLAYRS                                 
		    WRITE (IPR,980) L,P(L),T(L),
     *                       (WKL(M,L),M=MLO,MHI)                       
		 enddo
		 IF (NLAYRS.GT.1) THEN
		    WRITE (IPR,985)
		    L = NLAYRS
		    WRITE (IPR,990) L,PWTD,TWTD,
     *                       (WMT(M),M=MLO,MHI)
		 ENDIF
	      ELSE
		 WRITE (IPR,975) (HMOLID(I),I=MLO,MHI)              
		 DO L = 1, NLAYRS                                   
		    WRITE (IPR,982) L,P(L),T(L),
     *                       (WKL(M,L),M=MLO,MHI) 
		 enddo
		 IF (NLAYRS.GT.1) THEN                              
		    WRITE (IPR,985)
		    L = NLAYRS
		    WRITE (IPR,991) L,PWTD,TWTD,
     *                            (WMT(M),M=MLO,MHI)
		 ENDIF
	      ENDIF                                                 
C       
	   enddo                                                 
	ENDIF                                                       
C
C     --------------------------------------------------------------
C
C     Write out mixing ratios for molecules to TAPE6 in either
C     15.7 format (IFORM = 1) or 10.4 format (IFORM = 0).
C
C           Reset WDRAIR(L) for each layer
C           (WKL(M,L) now in column density)
C
C
	IF (IFORM.EQ.1) THEN
	   WRITE (IPR,976) (HMOLID(I),I=1,7),HOLN2
	   DO L = 1, NLAYRS
	      WDRAIR_l = WBRODL(L)
	      DO M = 2,NMOL
		 WDRAIR_l = WDRAIR_l + WKL(M,L)
	      enddo
	      IF (WDRAIR_l.EQ.0.0) THEN
		 WRITE(IPR,979)
	      ELSE
		 WRITE (IPR,980) L,P(L),T(L),
     *              (WKL(M,L)/WDRAIR_l,M=1,7),WBRODL(L)
	      ENDIF
	   enddo
	ELSE
	   WRITE (IPR,977) (HMOLID(I),I=1,7),HOLN2
	   DO L = 1, NLAYRS
	      WDRAIR_l = WBRODL(L)
	      DO M = 2,NMOL
		 WDRAIR_l = WDRAIR_l + WKL(M,L)
	      enddo
	      IF (WDRAIR_l.EQ.0.0) THEN
		 WRITE(IPR,979)
	      ELSE
		 WRITE (IPR,982) L,P(L),T(L),
     *                 (WKL(M,L)/WDRAIR_l,M=1,7),WBRODL(L)
	      ENDIF
	   enddo
	ENDIF
C
C
	IF (NMOL.GT.7) THEN
	   DO MLO = 8, NMOL, 8
	      MHI = MLO+7
	      MHI = MIN(MHI,NMOL)
	      IF (NLAYRS.LT.5) THEN
		 WRITE (IPR,970)
c       ELSE
c       WRITE (IPR,945) XID,(YID(M),M=1,2)
	      ENDIF
	      IF (IFORM.EQ.1) THEN
		 WRITE (IPR,976) (HMOLID(I),I=MLO,MHI)
		 DO L = 1, NLAYRS
		    IF (WDRAIR_l.EQ.0.0) THEN
		       WRITE(IPR,979)
		    ELSE
		       WRITE (IPR,980) L,P(L),T(L),
     *                       (WKL(M,L)/WDRAIR_l,M=MLO,MHI)
		    ENDIF
		 enddo
	      ELSE
		 WRITE (IPR,977) (HMOLID(I),I=MLO,MHI)
		 DO L = 1, NLAYRS
		    IF (WDRAIR_l.EQ.0.0) THEN
		       WRITE(IPR,979)
		    ELSE
		       WRITE (IPR,982) L,P(L),T(L),
     *                    (WKL(M,L)/WDRAIR_l,M=MLO,MHI)
		    ENDIF
		 enddo
	      ENDIF
	   enddo
	ENDIF
C
  970 FORMAT (////)                                                       A24360
  974 FORMAT ('0',53X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',/,13X,
     *        'P(MB)',6X,'T(K)',5X,8(A10,5X))
  975 FORMAT ('0',53X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',/,10X,    A24370
     *        'P(MB)',6X,'T(K)',3X,8(1X,A6,3X))
  976 FORMAT (/,'1',54X,'----------------------------------',
     *         /,'0',60X,'MIXING RATIOS BY LAYER ',/,10X,
     *        'P(MB)',6X,'T(K)',5X,8(A10,5X))                             A24380
  977 FORMAT (/,'1',54X,'----------------------------------',
     *         /,'0',60X,'MIXING RATIOS BY LAYER ',/,210X,
     *        'P(MB)',6X,'T(K)',3X,8(1X,A6,3X))
  979 FORMAT (/,'0','  MIXING RATIO IS UNDEFINED. DRYAIR DENSITY=0.0')
  980 FORMAT ('0',I3,F15.7,F9.2,2X,1P,8E15.7,0P)
  982 FORMAT ('0',I3,F12.5,F9.2,2X,1P,8E10.3,0P)                          A24390
  985 FORMAT ('0',54X,'ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH')     A24400
  990 FORMAT ('0',I3,F15.7,F9.2,2X,1P,8E15.7,0P,/,
     *         55X,1P,8E15.7,0P)
  991 FORMAT ('0',I3,F12.5,F9.2,2X,1P,8E10.3,0P)                          A24410
C     --------------------------------------------------------------

	return

	end
c___________________________________________________________________
c___________________________________________________________________

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
	character*8      hblnk
	data hblnk /'        '/
	hdate = hblnk
	RETURN
	END

        SUBROUTINE FTIME (HTIME)
	CHARACTER*8      HTIME   
	character*8      hblnk
	data hblnk /'        '/
	htime = hblnk
	RETURN
	END


	SUBROUTINE XSREAD (ipf,XV1,XV2)   	
C                                                                         E00020
      IMPLICIT REAL*8           (V)                                     ! E00030
C                                                                         E00040
C**********************************************************************   E00050
C     THIS SUBROUTINE READS IN THE DESIRED "CROSS-SECTION"                E00060
C     MOLECULES WHICH ARE THEN MATCHED TO THE DATA CONTAINED              E00070
C     ON INPUT FILE FSCDXS.                                               E00080
C**********************************************************************   E00090
      include "declar.incl"

C                                                                         E00100
C     IFIL CARRIES FILE INFORMATION                                       E00110
C                                                                         E00120
C
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         E00130
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        E00140
     *              NLTEFL,LNFIL4,LNGTH4                                  E00150
C                                                                         E00160
C     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE     E00170
C     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES       E00180
C     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR     E00190
C     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD         E00200
C     MOLECULES.                                                          E00210
C                                                                         E00220
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(mx_xs),XAMNT(mx_xs,MXLAY)              E00230
C                                                                         E00240
C     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES         E00250
C     FOR THE CROSS-SECTION MOLECULES.                                    E00260
C     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          E00270
C                                                                         E00280
c%%%%%LINUX_PGI90 (-i8)%%%%%      integer*4 iostat
      CHARACTER*10 XSFILE,XSNAME,ALIAS,XNAME,XFILS(6),BLANK               E00290
      COMMON /XSECTF/ XSFILE(6,5,mx_xs),XSNAME(mx_xs),ALIAS(4,mx_xs)
      COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS),     
     *                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),      
     *                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS), 
     *                NUMXS,IXSBIN                                   
C                                                                         E00330
      DIMENSION IXFLG(mx_xs)                                                 E00340
C                                                                         E00350
      CHARACTER*120 XSREC                                                 E00360
      CHARACTER*1 CFLG,CASTSK,CPRCNT,CFRM,CN,CF                           E00370
      EQUIVALENCE (CFLG,XSREC)                                            E00380
C                                                                         E00390
      DATA CASTSK / '*'/,CPRCNT / '%'/,CN / 'N'/,CF / 'F'/                E00400
      DATA BLANK / '          '/                                          E00410
C                                                                         E00411
C     T296 IS TEMPERATURE FOR INITAL CALCULATIN OF DOPPLER WIDTHS         E00412
C                                                                         E00413
      DATA T296 / 296.0 /                                                 E00414
C                                                                         E00420
      IXMAX = mx_xs                                                          E00430
      DO 10 I = 1, IXMAX                                                  E00440
         XSNAME(I) = BLANK                                                E00450
   10 CONTINUE                                                            E00460
C                                                                         E00470
C     READ IN THE NAMES OF THE MOLECULES                                  E00480
C                                                                         E00490
      IF (IXMOLS.GT.7) THEN                                               E00500
         READ (Ipf,'(7A10)') (XSNAME(I),I=1,7)                            E00510
         READ (Ipf,'(8A10)') (XSNAME(I),I=8,IXMOLS)                       E00520
      ELSE                                                                E00530
         READ (Ipf,'(7A10)') (XSNAME(I),I=1,IXMOLS)                       E00540
      ENDIF   

C                                                                         E00560
C     Left-justify all inputed names                                      E00570
C                                                                         E00580
      DO 15 I=1,IXMOLS                                                    E00582
         CALL CLJUST (XSNAME(I),10)                                       E00590
 15   CONTINUE
C                                                                         E00600
CPRT  WRITE(IPR,'(/,''  THE FOLLOWING MOLECULES ARE REQUESTED:'',//,      E00610
CPRT 1    (5X,I5,2X,A1))') (I,XSNAME(I),I=1,IXMOLS)                        E00620
C                                                                         E00630
C     MATCH THE NAMES READ IN AGAINST THE NAMES STORED IN ALIAS           E00640
C     AND DETERMINE THE INDEX VALUE.  STOP IF NO MATCH IS FOUND.          E00650
C     NAME MUST BE ALL IN CAPS.                                           E00660
C                                                                         E00670
      DO 40 I = 1, IXMOLS                                                 E00680
         DO 20 J = 1, IXMAX                                               E00690

            IF ((XSNAME(I).EQ.ALIAS(1,J)) .OR.                            E00700
     *          (XSNAME(I).EQ.ALIAS(2,J)) .OR.                            E00710
     *          (XSNAME(I).EQ.ALIAS(3,J)) .OR.                            E00720
     *          (XSNAME(I).EQ.ALIAS(4,J))) THEN                           E00730
               IXINDX(I) = J                                              E00740
               GO TO 30                                                   E00750
            ENDIF                                                         E00760
   20    CONTINUE                                                         E00770
C                                                                         E00780
C         NO MATCH FOUND                                                  E00790
C                                                                         E00800
         WRITE (IPR,900) XSNAME(I)                                        E00810
         STOP 'STOPPED IN XSREAD'                                         E00820
C                                                                         E00830
   30    CONTINUE                                                         E00840
         IXFLG(I) = 0                                                     E00850
   40 CONTINUE                                                            E00860
C                                                                         E00870
C     READ IN "CROSS SECTION" MASTER FILE FSCDXS                          E00880
C                                                                         E00890
      IXFIL = 8                                                           E00900
      OPEN (IXFIL,FILE='FSCDXS',STATUS='OLD',FORM='FORMATTED',
     *       IOSTAT=iostat)
        if (IOSTAT.gt.0) stop 'FSCDXS does not exist - XSREAD'
      REWIND IXFIL                                                        E00920
      READ (IXFIL,905)                                                    E00930
C                                                                         E00940
   50 READ (IXFIL,910,END=80) XSREC                                       E00950
C                                                                         E00960

      IF (CFLG.EQ.CASTSK) GO TO 50                                        E00970
      IF (CFLG.EQ.CPRCNT) GO TO 80                                        E00980
C                                                                         E00990
      READ (XSREC,915) XNAME,V1X,V2X,DVX,NTEMP,IFRM,CFRM,                 E01000
     *                 (XFILS(I),I=1,NTEMP)                               E01010
C                                                                         E01020
C     LEFT-JUSTIFY INPUTED NAME                                           E01030
C                                                                         E01040
      CALL CLJUST (XNAME,10)                                              E01050
C                                                                         E01060
C     CHECK MASTER FILE FOR CROSS SECTION AND STORE DATA                  E01070
C                                                                         E01080
      NUMXS = IXMOLS                                                      E01090
      DO 70 I = 1, IXMOLS                                                 E01100
         IF ((XNAME.EQ.ALIAS(1,IXINDX(I))) .OR.                           E01110
     *       (XNAME.EQ.ALIAS(2,IXINDX(I))) .OR.                           E01120
     *       (XNAME.EQ.ALIAS(3,IXINDX(I))) .OR.                           E01130
     *       (XNAME.EQ.ALIAS(4,IXINDX(I)))) THEN                          E01140
            IXFLG(I) = 1                                                  E01150
            IF (V2X.GT.XV1.AND.V1X.LT.XV2) THEN                           E01160
               NSPECR(I) = NSPECR(I)+1                                    E01170
               IF (NSPECR(I).GT.6) THEN                                   E01180
                  WRITE (IPR,920) I,XSNAME(I),NSPECR(I)                   E01190
                  STOP ' XSREAD - NSPECR .GT. 6'                          E01200
               ENDIF                                                      E01210
               IXFORM(NSPECR(I),I) = 91                                   E01220
               IF (IFRM.EQ.86) IXFORM(NSPECR(I),I) = IFRM                 E01230
               IF (CFRM.NE.CN)                                            E01240
     *             IXFORM(NSPECR(I),I) = IXFORM(NSPECR(I),I)+100          E01250
               IF (CFRM.EQ.CF)                                            E01260
     *             IXFORM(NSPECR(I),I) = -IXFORM(NSPECR(I),I)             E01270
               NTEMPF(NSPECR(I),I) = NTEMP                                E01280
               V1FX(NSPECR(I),I) = V1X                                    E01290
               V2FX(NSPECR(I),I) = V2X                                    E01300
C                                                                         E01301
C     3.58115E-07 = SQRT( 2.* LOG(2.)*AVOGAD*BOLTZ/(CLIGHT*CLIGHT) ) 
C                                                                         E01303
               XDOPLR(NSPECR(I),I)=3.58115E-07*(0.5*(V1X+V2X))*           E01304
     *                             SQRT(T296/XSMASS(IXINDX(I)))           E01305
C                                                                         E01306
               DO 60 J = 1, NTEMP                                         E01310
                  XSFILE(J,NSPECR(I),I) = XFILS(J)                        E01320
   60          CONTINUE                                                   E01330
            ENDIF                                                         E01340
         ENDIF                                                            E01350
   70 CONTINUE                                                            E01360
C                                                                         E01370
      GO TO 50                                                            E01380
C                                                                         E01390
   80 IXFLAG = 0                                                          E01400
      DO 90 I = 1, IXMOLS                                                 E01410
         IF (IXFLG(I).EQ.0) THEN                                          E01420
            WRITE (IPR,925) XSNAME(I)                                     E01430
            IXFLAG = 1                                                    E01440
         ENDIF                                                            E01450
   90 CONTINUE                                                            E01460
      IF (IXFLAG.EQ.1) STOP ' IXFLAG - XSREAD '                           E01470
C                                                                         E01480
      RETURN                                                              E01490
C                                                                         E01500
  900 FORMAT (/,'  THE NAME: ',A10, ' IS NOT ONE OF THE ',                E01510
     *        'CROSS SECTION MOLECULES. CHECK THE SPELLING.')             E01520
  905 FORMAT (/)                                                          E01530
  910 FORMAT (A120)                                                       E01540
  915 FORMAT (A10,2F10.4,F10.8,I5,5X,I5,A1,4X,6A10)                       E01550
  920 FORMAT (/,'******* ERROR IN XSREAD ** MOLECULE SECLECTED -',A10,    E01560
     *        '- HAS ',I2,' SPECTRAL REGIONS ON FILE FSCDXS, BUT THE',    E01570
     *        ' MAXIMUM ALLOWED IS 6 *******',/)                          E01580
  925 FORMAT (/,'******* MOLECULE SELECTED -',A10,'- IS NOT FOUND ON',    E01590
     *        ' FILE FSCDXS *******',/)                                   E01600
C                                                                         E01610

	RETURN
	END

      BLOCK DATA BXSECT                                                   E01630
C                                                                         E01640
      IMPLICIT REAL*8           (V)                                     ! E01650

	include "declar.incl"
c      parameter (mx_xs=38)
C                                                                         E01660
C**   XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          E01670
C**            (NOTE: ALL NAMES ARE LEFT-JUSTIFIED)                       E01680
C                                                                         E01690
      CHARACTER*10 XSFILE,XSNAME,ALIAS                                    E01700
      COMMON /XSECTI/ XSMAX(6,5,mx_xs),XSTEMP(6,5,mx_xs),   
     *                NPTSFX(5,mx_xs),NFILEX(5,mx_xs),NLIMX 
      COMMON /XSECTF/ XSFILE(6,5,mx_xs),XSNAME(mx_xs),ALIAS(4,mx_xs) 
      COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS),     
     *                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),      
     *                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS), 
     *                NUMXS,IXSBIN                                   
      COMMON /XSECTS/ JINPUT,NMODES,NPANEL,NDUM,V1XS,V2XS,DVXS,NPTSXS     E02870
C                                                                         E01750
      DATA NMODES / 1 /,NPANEL / 0 /,V1XS / 0.0 /,V2XS / 0.0 /,           E02990
     *     DVXS / 0.0 /,NPTSXS / 0 /                                      E03000
      DATA XSMAX / 1140*0.0 /                                             E03010
      DATA (ALIAS(1,I),I=1,mx_xs)/                                           E01760
     *    'CLONO2    ', 'HNO4      ', 'CHCL2F    ', 'CCL4      ',         E01770
     *    'CCL3F     ', 'CCL2F2    ', 'C2CL2F4   ', 'C2CL3F3   ',         E01780
     *    'N2O5      ', 'HNO3      ', 'CF4       ', 'CHCLF2    ',         E01790
     *    'CCLF3     ', 'C2CLF5    ', 'NO2       ',
     *	  23*' ZZZZZZZZ ' / 
      DATA (ALIAS(2,I),I=1,mx_xs)/                                           E01810
     *    'CLNO3     ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ',         E01820
     *    'CFCL3     ', 'CF2CL2    ', 'C2F4CL2   ', 'C2F3CL3   ',         E01830
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CHF2CL    ',         E01840
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 24*' ZZZZZZZZ ' /                   E01850
      DATA (ALIAS(3,I),I=1,mx_xs)/                                           E01860
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ',         E01870
     *    'CFC11     ', 'CFC12     ', 'CFC114    ', 'CFC113    ',         E01880
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC14     ', 'CFC22     ',         E01890
     *    'CFC13     ', 'CFC115    ', 24*' ZZZZZZZZ ' /                   E01900
      DATA (ALIAS(4,I),I=1,mx_xs)/                                           E01910
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F21       ', ' ZZZZZZZZ ',         E01920
     *    'F11       ', 'F12       ', 'F114      ', 'F113      ',         E01930
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F14       ', 'F22       ',         E01940
     *    'F13       ', 'F115      ', 24*' ZZZZZZZZ ' /                   E01950
C                                                                         E01960
C     XSMASS IS MASS OF EACH CROSS-SECTION                                E01961
C                                                                         E01962
c     note 23 =  mx_xs - 15 = 38 - 15

      DATA XSMASS/                                                        E01963
     1      97.46     ,   79.01     ,  102.92     ,  153.82     ,         E01964
     2     137.37     ,  120.91     ,  170.92     ,  187.38     ,         E01965
     3     108.01     ,   63.01     ,   88.00     ,   86.47     ,         E01966
     4     104.46     ,  154.47     ,   45.99     , 23*0.00 /             E01967
C                                                                         E01969
      DATA V1FX / 190*0.0 /,V2FX / 190*0.0 /,DVFX / 190*0.0 /,            E01970
     *     WXM / mx_xs*0.0 /                                                 E01980
      DATA NTEMPF / 190*0 /,NSPECR / mx_xs*0 /,IXFORM / 190*0 /,             E01990
     *     NUMXS / 0 /                                                    E02000
C                                                                         E02010
      END                                                                 E02020

      SUBROUTINE CLJUST (CNAME,NCHAR)                                     E02030
C                                                                         E02040
C     THIS SUBROUTINE LEFT-JUSTIFIES THE CHARACTER CNAME                  E02050
C                                                                         E02060
      CHARACTER*(*) CNAME                                                 E02070
      CHARACTER*25 CTEMP                                                  E02070
      CHARACTER*1  CTEMP1(25),BLANK                                       E02080
      EQUIVALENCE (CTEMP,CTEMP1(1))                                       E02090
C                                                                         E02100
      DATA BLANK / ' '/                                                   E02110
C                                                                         E02120
         CTEMP = CNAME                                                    E02140
         JJ=0                                                             E02145
         DO 10 J = 1, NCHAR                                               E02150
            IF (CTEMP1(J).NE.BLANK) THEN                                  E02160
               JJ = J                                                     E02170
               IF (JJ.EQ.1) GO TO 50                                      E02180
               GO TO 20                                                   E02190
            ENDIF                                                         E02200
   10    CONTINUE                                                         E02210
         IF (JJ .EQ. 0) GO TO 50                                          E02215
C                                                                         E02220
   20    KCNT = 0                                                         E02230
         DO 30 K = JJ, NCHAR                                              E02240
            KCNT = KCNT+1                                                 E02250
            CTEMP1(KCNT) = CTEMP1(K)                                      E02260
   30    CONTINUE                                                         E02270
C                                                                         E02280
         KK = NCHAR-JJ+2                                                  E02290
         DO 40 L = KK,NCHAR                                               E02300
            CTEMP1(L) = BLANK                                             E02310
   40    CONTINUE                                                         E02320
         CNAME = CTEMP                                                    E02330
   50 CONTINUE                                                            E02340
C                                                                         E02350
      RETURN                                                              E02360
C                                                                         E02370
      END                                                                 E02380

	FUNCTION MAXIMUM(ARRAY, NSIZE)
	DIMENSION ARRAY(NSIZE)
	MAXIMUM = ARRAY(1)
	DO 10 J=2,NSIZE
	    IF (MAXIMUM.LT.ARRAY(J)) MAXIMUM=ARRAY(J)
   10	CONTINUE
	RETURN
	END

	FUNCTION MINIMUM(ARRAY, NSIZE)
	DIMENSION ARRAY(NSIZE)
	MINIMUM = ARRAY(1)
	DO 10 J=2,NSIZE
	    IF (MINIMUM.GT.ARRAY(J)) MINIMUM=ARRAY(J)
   10	CONTINUE
	RETURN
	END

        SUBROUTINE MONORTM_XSEC_SUB(wn,nwn,p,t,nlay)

C----------------------------------------------------------------
C Authors: Eli Mlawer and Vivienne Payne, AER Inc.
C
C Created: July 2008
C
C Description: Calculates cross-section optical depths for MonoRTM.
C              Method follows that of LBLRTM. 
C              Results are consistent with LBLRTM. 
C              Note that only the total optical depth for all 
C              cross-section molecules is returned.  The "odxsec" array 
C              does not distinguish between different xsec molecules.     
C
C Input:
C      wn       array of wavenumbers
C      nwn      number of wavenumbers 
C      p        pressure of layer
C      t        temperature of layer
C      nlay     number of layers
C
C Output:  
C      odxsec   total optical depth for all cross-sections in array odxsec
C               (stored through declar.incl)
C    
C----------------------------------------------------------------

      IMPLICIT REAL*8           (V) ! for consistency with LBLRTM routines

      Include "declar.incl"


C                                                                         
C     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES         
C     FOR THE CROSS-SECTION MOLECULES.                                    
C     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          

C     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE    
C     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES       
C     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR     
C     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD         
C     MOLECULES.                                                          
C                                                                         
C     NUMXS IS THE NUMBER OF 'CROSS SECTION' MOLECULES TO BE USED         
C                                                                         
C     XSFILE(ITEMP,ISPEC,NXS) IS THE NAME OF THE FILE CONTAINING THE      
C                             'CROSS SECTION' DATA.  THE THREE INDICES    
C                             ARE DEFINED AS FOLLOWS:                     
C                                                                         
C                             ITEMP - DENOTES THE TEMPERATURE FOR WHICH   
C                                     THE 'CROSS SECTION' IS VALID        
C                                     (IMPLEMENTED FOR HITRAN 91 TAPE)    
C                             ISPEC - DENOTES THE SECTRAL REGION FOR      
C                                     WHICH THE FILE PERTAINS             
C                             NXS   - IS THE INCREMENT FOR THE 'CROSS     
C                                     SECTION' INDEX                      
C                                                                         
C     NTEMPF(ISPEC,NXS) IS THE NUMBER OF TEMPERATURE FILES TO BE USED     
C                       FOR EACH SPECTRAL REGION OF EACH MOLECULE         
C                                                                         
C     NSPECR(NXS) IS THE NUMBER OF SPECTRAL REGIONS FOR THE MOLECULE NX   
C**********************************************************************   
C                        

      CHARACTER*10 XSFILE,XSNAME,ALIAS,XNAME,XFILS(6),BLANK,xxfile,ctorr               
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     $     RADCN1,RADCN2 
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(mx_xs),XAMNT(mx_xs,MXLAY)        
      COMMON /XSECTF/ XSFILE(6,5,mx_xs),XSNAME(mx_xs),ALIAS(4,mx_xs)
      COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS),     
     *                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),      
     *                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS), 
     *                NUMXS,IXSBIN                                   
C                                                                         
      DIMENSION IXFLG(mx_xs),tx(6,5,mx_xs),pdx(6,5,mx_xs)
      dimension xsdat(150000,6),xspd(150000)
      dimension xspave(nwnmx),xsmoltot(nwnmx,mxlay),xstot(nwnmx,mxlay)
C                                                                         
      CHARACTER*120 XSREC
      character*10 source(3)

      data dvbuf /1.0/         ! used to check if a particular xs file needs to be processed
      DATA CTORR / '      TORR'/

C     DEFINE PRESSURE CONVERSIONS                                         E07160
C                                                                         E07170
C        PTORMB = 1013. MB / 760. TORR  (TORR TO MILLIBARS)               E07180
C        PATMMB = 1013. MB / 1.0  ATM   (ATMOPHERES TO MILLIBARS)         E07190
C                                                                         E07200
      PTORMB = 1013./760.                                                 E07210
      PATMMB = 1013.   

C Initialize variables.
      xstot(1:nwn,1:nlay) = 0.

c loop over cross-section molecules
      DO 6000 IXMOL = 1,IXMOLS
          xsmoltot(1:nwn,1:nlay) = 0.

C loop over the spectral regions for this molecule
          DO 5000 IXSR = 1,NSPECR(IXMOL)
C Check to see if this spectral region needs to be processed.
              ipt = 1
 2100         continue
              if (wn(ipt) .ge. v1fx(ixsr,ixmol)-dvbuf .and.
     &            wn(ipt) .le. v2fx(ixsr,ixmol)+dvbuf) goto 2200
              ipt = ipt + 1
              if (ipt .gt. nwn) goto 5000
              goto 2100

 2200         continue

C Getting here means that this spectral region for this molecule needs to be handled.
C Loop over the temperatures for which this molecule has stored values in this region.

C The filenames are arranged in FSCDXS in order of ascending temperature
              DO 3000 IXTEMP=1,NTEMPF(IXSR,IXMOL)
                  ifile = 100 + ixtemp
                  xxfile = xsfile(ixtemp,ixsr,ixmol)
                  open(ifile,file=xxfile,form='formatted')
                  READ (ifile,910) AMOL,V1x,V2x,NPTSx,
     &                TX(ixtemp,ixsr,ixmol),PRES,    
     &                              SMAX,SOURCE
                  if (source(3) .eq. ctorr) then
                      pdx(ixtemp,ixsr,ixmol) = pres * ptormb
                  else
                      pdx(ixtemp,ixsr,ixmol) = pres
                  endif
                  read(ifile,*)(xsdat(i,ixtemp),i=1,nptsx)
 3000         continue

C Interpolate absorption coefficients to given temperature.
              
              do 4000 il = 1, nlay
                  pave = p(il)
                  tave = t(il)
                  coef1 = 1.
                  ind2 = 1
                  coef2 = 0.
                  it = 1
C The temperatures at which the xsecs are stored are in ascending order (as in LBLRTM).
C it=1 is the lowest temperature
                  if (ntempf(ixsr,ixmol).eq.1 .or. 
     &                tave.le.tx(it,ixsr,ixmol)) then
                      ind1 = 1
                  else
 3100                 continue
                      it = it + 1
                      if (it. gt. ntempf(ixsr,ixmol)) then
                          ind1 = ntempf(ixsr,ixmol)
                          ind2 = ntempf(ixsr,ixmol)
                      elseif (tave. le. tx(it,ixsr,ixmol)) then                      
                          ind1 = it-1 ! ind1 is a lower temperature than ind2
                          ind2 = it
                          coef1 = (tave - tx(it,ixsr,ixmol))/
     &                        (tx(it-1,ixsr,ixmol) - tx(it,ixsr,ixmol))
                          coef2 = 1. - coef1
                      else
                          goto 3100
                      endif
                  endif

                  pd = coef1 * pdx(ind1,ixsr,ixmol)+ 
     &                coef2* pdx(ind2,ixsr,ixmol)
                  xkt1 = tx(ind1,ixsr,ixmol)/radcn2
                  xkt2 = tx(ind2,ixsr,ixmol)/radcn2
                  
                  delvx =  (v2x-v1x) / float(nptsx-1)
                  do 3300 i = 1, nptsx
                      vv = v1x + float(i-1) * delvx
                      xspd(i) = coef1 * xsdat(i,ind1) / radfn(vv,xkt1)                  
     &                    + coef2 * xsdat(i,ind2) / radfn(vv,xkt2)
 3300             continue

C Perform convolution from appropriate pressure for stored values to layer pressure.
C (296K is the temperature used in the calculation of the Doppler width)
                  hwdop = xdoplr(ixsr,ixmol) * sqrt(tave/296.)
                  call convolve(xspd,v1x,v2x,delvx,pd,hwdop,
     &                tave,pave,x,wn,nwn,xspave)
                  xsmoltot(1:nwn,il) = xsmoltot(1:nwn,il) + 
     &                xspave(1:nwn)
 4000         continue
C Finished processing this spectral range of this molecule.

 5000     continue

C Multiply by absorber amount.
          do 5500 il = 1, nlay
              xstot(1:nwn,il) = xstot(1:nwn,il) + 
     &            xamnt(ixmol,il) * xsmoltot(1:nwn,il)
 5500     continue
C Finished prcoessing this molecule.

 6000 continue
C Finished processing all molecules.
C Put the radiation field back in.
      do 6500 il = 1, nlay
          xkt = t(il)/radcn2
          do 6300 iwn = 1, nwn              
              odxsec(iwn,il) = xstot(iwn,il) * radfn(wn(iwn),xkt)
 6300     continue
 6500 continue

       RETURN  
  910 FORMAT (A10,2F10.4,I10,3G10.3,3A10)              

       END

      subroutine convolve(xspd,v1x,v2x,delvx,pd,hwdop,
     &    tave,pave,x,wn,nwn,xspave)

c performs pressure convolution for cross sections

      Include "declar.incl"
      dimension xspd(150000),xspave(nwnmx),xspd_int(0:10000000)
      data p0 /1013./

C     Set up the halfwidths needed.
      hwpave = 0.1 * (pave/p0) * (273.15/tave)
      hwd = 0.1 * (pd/p0) * (273.15/tave)
      hwd = max (hwd,hwdop)
C Don't want hwd to be bigger than hwpave
      if (hwd.gt.hwpave) then hwpave = 1.001*hwd
      hwb = hwpave - hwd

C Step size is at least ~4 pts per halfwidth.
      ratio = 0.25
      step = ratio * hwb
      if (step .gt. delvx) step = delvx
      npts = int((v2x-v1x)/step)
      step = (v2x - v1x) / float(npts)
      ratio = step / hwb
C Linearly interpolate incoming values to desired grid.
      do 500 i = 0, npts
          vv = v1x + float(i) * step
          delvv = vv - v1x
          ind = int(delvv / delvx)
          coef = (delvv - float(ind) * delvx) / delvx
          xspd_int(i) = (1.-coef)*xspd(ind+1) + coef*xspd(ind+2)
  500 continue

      hwb2 = hwb**2
      do 4000 iwn = 1, nwn
          if (wn(iwn) .lt. v1x .or. wn(iwn) .gt. v2x) then
              xspave(iwn) = 0.
              goto 4000
          endif
          if (hwb/hwd .gt. 0.1) then
C Perform convolution using Lorentzian with width hwb.
              wn_v1x = wn(iwn) - v1x
              ind = int(wn_v1x / step)
              dvlo = wn(iwn) - (v1x + float(ind) * step)
              dvhi = wn(iwn) - (v1x + float(ind+1) * step)
              answer = (hwb/(hwb2 + dvlo**2)) * xspd_int(ind) +
     &            (hwb/(hwb2 + dvhi**2)) * xspd_int(ind+1)
              j = 1
 1000         continue
              vlo = v1x + float(ind-j) * step
              if (vlo .gt. v1x) then
                  dvlo = wn(iwn) - vlo
                  contlo = (hwb/(hwb2 + dvlo**2)) * xspd_int(ind-j)
              else
                  contlo = 0.
              endif
              vhi = v1x + float(ind+j+1) * step
              if (vhi .lt. v2x) then
                  dvhi = wn(iwn) - vhi
                  conthi = (hwb/(hwb2 + dvhi**2)) * xspd_int(ind+j+1)
              else
                  conthi = 0.
              endif
              xincr = contlo + conthi
              if ((xincr/answer) .lt. ratio*1e-6) goto 1800
              answer = answer + xincr
              j = j + 1
              goto 1000
 1800         continue
C Done with convolution.  Multiply by step size and Lorentzian normalization.
              xspave(iwn) = answer * step / 3.14159
          else
C     Use linearly interpolated values.
              wn_v1x = wn(iwn) - v1x
              ind = int(wn_v1x / delvx)
              coef = (wn_v1x - float(ind) * delvx) / delvx
              xspave(iwn) = (1.-coef)*xspd(ind) + coef*xspd(ind+1)
          endif
 4000 continue

      return
      
      end

      subroutine calctmr(nlayrs, nwn, wn, T, tz, tmr)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CAuthor: Dave Turner, January 2008
C   
CThis routine computes the mean radiating temperature from the optical depth
C and temperature profiles.  The logic was provided by Vivienne Payne, AER,
C in an email on 10 Jan 2008.
C
C This routine is currently not connected to the body of the code, but can 
C easily be called from the subroutine rtm.
C
C Results have been checked against results from a modified version of rtm.f 
C that Tony Clough had supplied to Jim Liljegren and Maria Caddedu at ANL 
c pre-2006 for the purposes of calculations for Jim's statistical retrievals.
C
C Note that this routine is currently only applicable to downwelling 
C calculations
C
C Inputs:
C    nlayrs:	The number of layers in the model atmosphere
C    nwn: 	The number of spectral channels
C    wn:	The wavenumber array, in cm-1 (NWN)
C    T:		The layer-averaged temperature profile, in K (NLAYRS)
C    O:		The optical depth data, in nepers (NWN x NLAYRS) 
C Output:
C    Tmr:	The mean radiating temperature spectrum, in K (NWN)
C
C  Vivienne Payne, AER Inc, 2008
C------------------------------------------------------------------------------
	include "declar.incl"
	integer nlayrs, nwn
        integer ifr, ilay
	real    sumtau, sumexp
        real    bbvec(MXLAY),bbavec(0:MXLAY),odtot(NWNMX)
        real    odt, odvi, vv, beta, beta_a
        real    radtmr, x
	bb_fn(v,fbeta)  = radcn1*(v**3)/(exp(v*fbeta)-1.)
        real    tmr(*)
        COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     1       RADCN1,RADCN2 


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
C calculate the "effective emissivity" of the layer using "linear in tau"
C (see Clough et al 1992)
                pade=0.193*odvi+0.013*odvi**2
                beff = (bb + pade*bba)/(1.+pade)
                sumexp = sumexp + beff*trtot(ifr)*(1-tri)
            enddo

C this bit is based on Han & Westwater (2000) eq 14
            radtmr = sumexp / (1. - exp(-1*odtot(ifr)))
            x=radcn1*(wn(ifr)**3)/radtmr+1.
            tmr(ifr) = radcn2*wn(ifr) / log(x)

        enddo
        
        return
       end



