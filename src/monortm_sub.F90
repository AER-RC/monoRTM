  SUBROUTINE READEM(ICOEF)
!Reads in emission function values directly from file "EMISSION"
  IMPLICIT REAL*8           (V)
  PARAMETER (NMAXCO=4040)
  COMMON /EMSFIN/ V1EMIS,V2EMIS,DVEMIS,NLIMEM,ZEMIS(NMAXCO)
  READ (ICOEF,900,END=20,ERR=20) V1EMIS,V2EMIS,DVEMIS,NLIMEM
  DO 100 NGNU = 1,NLIMEM
     READ (ICOEF,910,END=20,ERR=20) ZEMIS(NGNU)
 100 CONTINUE
  RETURN
 900 FORMAT (3E10.3,5X,I5)
 910 FORMAT (E15.7)
 20  WRITE(*,*) 'INCONSISTENT DATA OR ERROR OPENING IN READEM'
  STOP
  END

  SUBROUTINE READRF(ICOEF)
  IMPLICIT REAL*8           (V)
  PARAMETER (NMAXCO=4040)
  COMMON /RFLTIN/ V1RFLT,V2RFLT,DVRFLT,NLIMRF,ZRFLT(NMAXCO)
  READ (ICOEF,900,END=20,ERR=20) V1RFLT,V2RFLT,DVRFLT,NLIMRF
  DO 100 NGNU = 1,NLIMRF
     READ (ICOEF,910,END=20,ERR=20) ZRFLT(NGNU)
 100 CONTINUE
  RETURN
 900 FORMAT (3E10.3,5X,I5)
 910 FORMAT (E15.7)
 20  WRITE(*,*) 'INCONSISTENT DATA OR ERROR OPENING IN READRF'
  END



  SUBROUTINE RDLBLINP(IATM,IPLOT,IOD,IRT,NWN,WN, &
       FILEIN,cntnmScaleFac,IXSECT,IBMAXOUT,TMPBND,ZBNDOUT, &
       H1fout,H2fout,ISPD,IPASSATM,IBRD)
!-------------------------------------------------------------------------------
!
!     SUBROUTINE:  RDLBLINP
!     -----------

!     AUTHOR: Sid-Ahmed Boukabara 
!     ------
!
!     AFFILIATION: ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!     -----------
!
!     DATE OF CREATION : April 2001 
!     Modified on Jan 7th 2002.
!
!
!     AIM: This subroutine reads the control parameters
!     ---- from the input file MONORTM.IN (former TAPE5) 
!          and makes available some of them via common 
!          blocks. This subroutine has been added to make 
!          MONORTM compatible with LBLRTM inputs.
!          S.A. Boukabara AER INC. 1999
!          IRT=1  
           !=1->Space
           !=2->limb
           !=3->ground
!     UPDATES:
!     --------
!     Extension of the MONORTM.IN option for MonoRTM
!     If V1 or V2 is negative, then we expect a 
!     finite number of wavenumbers to be included
!     in the input file. These wavenumbers do not have 
!     to be equally spaced. This is a special option
!     specific for MonoRTM. This is particularly
!     useful for those interested in simulating 
!     ARM monochromatic frequencies for instance (23.8 
!     and 31.4 GHz)
!
!     Sid Ahmed Boukabara, April 2001
!
!-------------------------------------------------------------------------------
  USE CntnmFactors, ONLY: CntnmFactors_t,applyCntnmCombo
  USE lblparams, ONLY: MXLAY,MXFSC,MX_XS
  USE RTMmono, ONLY: NWNMX
  !include "declar.incl"
  REAL*8           V1,V2,SECANT,XALTZ,WN(NWNMX) 
  character*4 ht1,ht2
  CHARACTER*1 CMRG(2),CONE,CTWO,CTHREE,CFOUR,CXIDLINE*80
  CHARACTER CDOL*1,CPRCNT*1,CXID*1,CA*1,CB*1,CC*1
  INTEGER IRD,IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,L, &
       ISCAN,IFILTR,IPLOT,ITEST,IATM,ILAS,ILNFLG, &
       IOD,IXSECT,MPTS,NPTS,INFLAG,IOTFLG,JULDAT
  INTEGER IPASSATM
  REAL SAMPLE,DVSET,ALFAL0,AVMASS,DPTMIN,DPTFAC,DVOUT 
  REAL TMPBND,XVMID,EMITST,REFTST
  INTEGER IBND,ICOEF,IMRG,LAYTOT,IFORM,NLAYRS,NMOL
  REAL PATH1,PTHODT,SECNT0,ZH1,ZH2,ZANGLE,PAVE,TAVE,SECNTK
  CHARACTER HEAD20*20,HEAD7*7,HEAD5*5,HEAD4*4,CINP*3
  CHARACTER*60 FILEOUT,FILEIN
  character*1 hmol_scal
  REAL SECL(64),WDNSTY,WMXRAT,WDRAIR(MXLAY)
  REAL ZBNDOUT(MXFSC)
  REAL CLW(MXLAY)
  INTEGER IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, &
       NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL, &
       NLTEFL,LNFIL4,LNGTH4,IPTHRK,IPATHL,M
  character*8 XID,HMOLID,YID
  TYPE(CntnmFactors_t) :: cntnmScaleFac
  COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,NOZERO,NOP,H1F,H2F, &
       ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH, &
       XVBAR, HMINF,PHIF,IERRF,HSPACE                     
  COMMON /BNDRY/ ZBND(MXFSC),PBND(MXFSC),TBND(MXFSC), &
       ALORNZ(MXFSC),ADOPP(MXFSC),AVOIGT(MXFSC)       
  COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX, &
       IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX, & 
       IPHMID,IPDIM,KDIM,KMXNOM,KMAX                      
  common /profil_scal/ nmol_scal,hmol_scal(64),xmol_scal(64)

  COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(mx_xs),XAMNT(mx_xs,MXLAY)

  COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, &
       NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL, &   
       NLTEFL,LNFIL4,LNGTH4                                 
  COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
       AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS, &     
       DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1, &    
       ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,EXTID(10)    
  COMMON /BNDPRP/ BNDEMI(3),BNDRFL(3)            
  COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4), &
       WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
       EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF 
  COMMON /CNTSCL/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL

  DATA CDOL / '$'/,CPRCNT / '%'/
  DATA CONE / '1'/,CTWO / '2'/,CTHREE / '3'/,CFOUR / '4'/, &
       CA / 'A'/,CB / 'B'/,CC / 'C'/                   
  character*6 idcntl(9)
  DATA IDCNTL / ' HIRAC',' CNTNM',' EMISS', &
                '  PLOT','  IATM','   IOD', &    
                ' XSECT','  ISPD','  IBRD'/
  CHARACTER CEXST*2    
  DATA CEXST/'EX'/
  EQUIVALENCE (CXID,CXIDLINE)                    
  EQUIVALENCE (FSCDID(3),IXSCNT)
 20   READ (IRD,905,END=80,ERR=6000) CXIDLINE      !---record 1.1
      IF (CXID.EQ.CPRCNT) THEN
         WRITE(*,*) '-END OF FILE:',FILEIN
         RETURN 
      ENDIF                       
      IF (CXID.NE.CDOL) GO TO 20                            
      READ (CXIDLINE,'(1x,10A8)') (XID(I),I=1,10)     

      print *,'before header read'
      READ(IRD,925,END=80,ERR=6000) IHIRAC, &!---record 1.2
         ICNTNM,IEMIT,IPLOT, &    
         IATM, IOD,IXSECT, &
         ISPD,IBRD
      print *,'after header read'

      IXSCNT = IXSECT*10 + ICNTNM

      WRITE (IPR,935) (IDCNTL(I),I=1,9)  
 935  FORMAT (15(A6,3X)) 

      Write(ipr,940)                 IHIRAC,&!---record 1.2           
         ICNTNM,IEMIT,IPLOT, &   
         IATM,IOD,IXSECT,ISPD,IBRD

 940  FORMAT (1X,I4,9I9,2(8x,a1),3I9)

!----CHECKING THE INPUTS FROM RECORD 1.2
      IF (IEMIT.NE.1) THEN
         WRITE(*,*) '----------------------------------------'
         WRITE(*,*) 'WARNING: IEMIT IS IGNORED IN MONORTM'
         WRITE(*,*) 'IT IS SET INTERNALLY TO ONE'
         WRITE(*,*) 'IN ',FILEIN,' IEMIT=',IEMIT
      ENDIF
      IF (IPLOT.NE.1) THEN 
         WRITE(*,*) 'WARNING: IPLOT MUST BE SET TO 1 TO OUTPUT TBs'
         WRITE(*,*)  'IN ', FILEIN, ' IPLOT=', IPLOT
      ENDIF
      IF (IOD.EQ.1) THEN
         WRITE(*,*) '----------------------------------------'
         WRITE(*,*) 'IOD FLAG SET TO OUTPUT LAYER OPTICAL DEPTHS'
         WRITE(*,*) 'IN ',FILEIN,' IOD=',IOD
      ENDIF

!------END CHECKING RECORD 1.2


      IF (ICNTNM.EQ.6) THEN
         READ(IRD,*) cntnmScaleFac%XSELF, & ! each is type REAL   |---record 1.2a
                 cntnmScaleFac%XFRGN, & 
                 cntnmScaleFac%XCO2C, & 
                 cntnmScaleFac%XO3CN, &
                 cntnmScaleFac%XO2CN, &
                 cntnmScaleFac%XN2CN, &
                 cntnmScaleFac%XRAYL
      ELSE
         CALL applyCntnmCombo(ICNTNM,cntnmScaleFac)
      ENDIF   


      IF (IEMIT.EQ.2) THEN
         READ(IRD,1010,ERR=6000) INFLAG,IOTFLG,JULDAT !---record 1.2.1
      ENDIF
      IF (IEMIT.EQ.3) THEN
         WRITE(*,*) 'CURRENTLY MONORTM DOES NOT HANDLE DERIVATIVES'
         WRITE(*,*) 'CHECK YOUR INPUT FILE: ',FILEIN
         STOP
      ENDIF

      IF ((IHIRAC+IAERSL+IEMIT+IATM+ILAS).GT.0) THEN !---record 1.3
         READ (IRD,970,END=80,ERR=6000) V1,V2,SAMPLE,DVSET, &
              ALFAL0,AVMASS,DPTMIN,DPTFAC,ILNFLG,DVOUT,nmol_scal

!        read in the profile scaling parameters
!
         if (nmol_scal .gt. 0 ) then
            if (nmol_scal .gt. 38) stop ' nmol_scal .gt. 38 '
            read (ird,9701) (hmol_scal(m),m=1,nmol_scal)
            read (ird,9702) (xmol_scal(m),m=1,nmol_scal)
 9701       FORMAT (64a1)
 9702       FORMAT (7e15.7,/,(8e15.7,/))
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
         IF ((DVSET.LE.0.).AND.(V1.NE.V2).AND.(V1.GT.0.).AND. &
               (V2.GT.0.)) THEN
            WRITE(*,*) 'MONORTM REQUIRES POSITIVE DVSET,'
            WRITE(*,*) 'OR (V1=V2 AND DVSET=0) IF YOU WANT ONLY '
            WRITE(*,*) 'ONE FREQ PROCESSED'
            WRITE(*,*) 'CHECK YOUR INPUT FILE: ',FILEIN,V1,V2,DVSET
            STOP
         ENDIF


         IF ((V1.LT.0.).OR.(V2.LT.0.)) THEN   !---record 1.3.1
            READ(IRD,'(I8)',ERR=6000)NWN
            IF (NWN.GT.NWNMX) THEN
               WRITE(*,*) 'STOP: NUMBER OF WAVENUMBERS ', &
                   'EXCEEDS LIMIT. ',NWN,NWNMX
               WRITE(*,*) 'FIX: EXTEND NWNMX IN DECLAR.INCL'
               STOP
            ENDIF
       
            DO IWN=1,NWN
               READ(IRD,'(E19.7)',ERR=6000) WN(IWN) !---record 1.3.2
            ENDDO
            dvset=0.
         ELSE
            IF (DVSET.NE.0.) THEN 
               NWN=NINT(((V2-V1)/DVSET)+1.)
               IF (NWN.GT.NWNMX) THEN
                  WRITE(*,*) 'STOP: NUMBER OF WAVENUMBERS ', &
                       'EXCEEDS LIMIT. ',NWN,NWNMX
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
         READ (IRD,970,END=80,ERR=6000) TMPBND, &
            (BNDEMI(IBND),IBND=1,3), &
            (BNDRFL(IBND),IBND=1,3)                    
         WRITE (IPR,985) TMPBND,(BNDEMI(IBND),IBND=1,3), &
            (BNDRFL(IBND),IBND=1,3)          ! surf_refl
 985     FORMAT (5(/),'0*********** BOUNDARY PROPERTIES ***********',/, &
               '0 TBOUND   = ',F12.4,5X,'BOUNDARY EMISSIVITY   = ',   &
               3(1PE11.3),/,'0',29X,'BOUNDARY REFLECTIVITY = ',       &
               3(1PE11.3),/,'0',29X,' SURFACE REFLECTIVITY = ', A1)
!
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
!---record 2.1
      IF (IATM.EQ.0)  return
!
!  IPASSATM=0 means this is the first call to RDLBLINP
!  For the first call to this routine, we don''t need to go any further
      IF (IPASSATM.EQ.0) THEN 
         IPASSATM = 1
         RETURN
      ENDIF

      IF (IATM.EQ.1) then

!          clw is not returned from lblatm         
         do l=1,nlayrs
            clw(l) = 0.
         enddo
!   
         CALL LBLATM
!
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

 10   FORMAT (i4,5f19.13)
 901  FORMAT (1X,I1,I3,I5,F10.6,A20,F8.2,A4,F8.2,A5,F8.3,A7)
 905  FORMAT (A80)                                                   
 910  FORMAT (E15.7,F10.4,F10.6,A3,I2,1X,2(F7.2,F8.3,F7.2),E15.7)
 911  FORMAT (2F10.4,f10.6,A3,I2,1X,2(F7.2,F8.3,F7.2),E15.7)           
 915  FORMAT (E15.7,F10.4,F10.6,A3,I2,23X,(F7.2,F8.3,F7.2),E15.7)
 916  FORMAT (2F10.4,f10.6,A3,I2,23X,(F7.2,F8.3,F7.2),E15.7)         
 920  FORMAT (I3)                                                 
 925  FORMAT (4X,I1,9X,I1,9X,I1,14X,I1,9X,I1,14X,I1,4X,I1,16X,I4,I4)
 9255 FORMAT (8E15.7)
 927  FORMAT (8E10.3)                                          
 930  FORMAT (I1)                                               
 945  FORMAT (A55,1X,I4)
 946  FORMAT (A55)
 970  FORMAT (8E10.3,4X,I1,5x,e10.3,i5)                           
 990  FORMAT (F20.8)                                   
 1000 FORMAT ('Layer',I2,': Changing molecule ',I2,' from ',E10.3, &
       ' to 1.000E+20.')
 1010 FORMAT (2I5,2X,I3)
 1015 FORMAT (I5)
      RETURN
 80   WRITE(*,*) ' EXIT; EOF ON :',FILEIN             
      STOP
 4000 WRITE(*,*) ' EXIT; ERROR OPENING EMISSION FILE'          
      STOP
 5000 WRITE(*,*) ' EXIT; ERROR OPENING REFLECTION FILE'          
      STOP
 6000 WRITE(*,*) ' EXIT; ERROR READING :',FILEIN          
      STOP
      END


  FUNCTION EMISFN (VI)                         
  IMPLICIT REAL*8           (V)                        
        !  FUNCTION EMISFN CALCULATES BOUNDARY EMISSIVITY FOR WAVE NUMBER      
        !  VALUE CORRESPONDING TO       
  PARAMETER (NMAXCO=4040)
  REAL TMPBND
  COMMON /EMSFIN/ V1EMIS,V2EMIS,DVEMIS,NLIMEM,ZEMIS(NMAXCO)
  COMMON /BNDPRP/ BNDEMI(3),BNDRFL(3)       
  EQUIVALENCE (BNDEMI(1),A) , (BNDEMI(2),B) , (BNDEMI(3),C)     
!---Check for A < 0->use inputs in file "EMISSION"
  IF (A.LT.0.) THEN
     NELMNT = INT((VI-V1EMIS)/DVEMIS)
     IF ((NELMNT.LE.0).OR.(NELMNT.GE.NLIMEM)) THEN
        WRITE(*,*) 'Frequency range of calculation exceeded', &
             ' emissivity input.'
        WRITE(*,*) ' VI = ',VI,' V1EMIS = ',V1EMIS,' V2EMIS = ', &
             V2EMIS
        STOP 'ERROR IN EMISFN'
     ENDIF
     V1A = V1EMIS+DVEMIS*NELMNT
     V1B = V1EMIS+DVEMIS*(NELMNT+1)
     CALL LINTCO(V1A,ZEMIS(NELMNT),V1B,ZEMIS(NELMNT+1),VI,ZINT,ZDEL)
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
  REAL TMPBND
  COMMON /RFLTIN/ V1RFLT,V2RFLT,DVRFLT,NLIMRF,ZRFLT(NMAXCO)
  COMMON /BNDPRP/ BNDEMI(3),BNDRFL(3)         
  EQUIVALENCE (BNDRFL(1),A) , (BNDRFL(2),B) , (BNDRFL(3),C)     
!---Check for A < 0->use values in from file"REFLECTION"
  IF (A.LT.0.) THEN
     NELMNT = INT((VI-V1RFLT)/DVRFLT)
     IF ((NELMNT.LE.0).OR.(NELMNT.GE.NLIMRF)) THEN
        WRITE(*,*) 'Frequency range of calculation exceeded', &
             ' reflectivity input.'
        WRITE(*,*) ' VI = ',VI,' V1RFLT = ',V1RFLT,' V2RFLT = ', &
             V2RFLT
        STOP 'ERROR IN REFLFN'
     ENDIF
     V1A = V1RFLT+DVRFLT*NELMNT
     V1B = V1RFLT+DVRFLT*(NELMNT+1)
     CALL LINTCO(V1A,ZRFLT(NELMNT),V1B,ZRFLT(NELMNT+1),VI,ZINT, &
          ZDEL)
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
  ZDEL = (Z2-Z1)/(V2-V1) !--ZDEL is the slope of the line
  ZCEPT = Z1 - ZDEL*V1   !--ZCEPT is the intercept for V = 0.0
  ZINT = ZDEL*VINT + ZCEPT !--Calculate ZINT value at VINT
  RETURN
  END


  SUBROUTINE EMISS_REFLEC(NWN,EMISS,REFLC,WN)
  real*8            wn(*)
!
  INTEGER NWN,J
  REAL REFLFN,EMISFN,reflc(*),emiss(*)
  DO J=1,NWN
     REFLC(J)=REFLFN(WN(J))
     EMISS(J)=EMISFN(WN(J))
  ENDDO
  RETURN 
  END


       SUBROUTINE STOREOUT(NWN,WN,WKL,WBRODL,RAD,TB,TRTOT, NPR, &
            O,O_BY_MOL, OC, O_CLW, ODXSEC, TMR, &
            WVCOLMN,CLWCOLMN,TMPSFC,REFLC,EMISS, &
            NLAY,NMOL,ANGLE,IOT,IOD) 
       USE PhysConstants, ONLY: getPhysConst
       USE RTMmono, ONLY: NWNMX
       USE netcdf_helper_mod
       !include "declar.incl"
       USE lblparams, ONLY: MXLAY,MXMOL
#ifdef USENETCDF
       USE netcdf
#endif

       INTEGER*4 NWN,NLAY
       INTEGER NMOL,NPR,IOT,IOD
       REAL ANGLE
       REAL CLWCOLMN,TMPSFC,WVCOLMN
       REAL*8 WN(:)
       REAL TMR(:)
       REAL WBRODL(:)
       REAL, DIMENSION(:) :: RAD,TRTOT,TB,EMISS,REFLC
       REAL O(:,:),OC(:,:,:), &
           O_BY_MOL(:,:,:),O_CLW(:,:), &
           odxsec(:,:),WKL(:,:)
  

       INTEGER I,J,IOL

       CHARACTER FILEOD*22
       character wnunits*12
       character*8 hmolc(mxmol),cmol(mxmol)
       logical giga

       REAL OTOT,otot_by_mol(mxmol,nwn),wk_tot(mxmol),odxtot
       REAL, ALLOCATABLE :: O_BY_MOL_LAYER(:,:,:)
       integer id_mol(mxmol)
       REAL CLIGHT
       REAL FREQ

       integer rv
       integer*4 ncid
       integer*4 wn_dimid
       integer*4 mol_dimid
       integer*4 dimids(2)
       integer*4 cdimids(2)
       integer*4 odimids(2)
       integer*4 fdimids(3)
       integer*4 s_dimid
       integer*4 lay_dimid

       integer*4 varid(20)

       real freqa(NWN)
       real wvcolmna(NWN)
       real clwcolmna(NWN)
       real tmpsfca(NWN)
       real anglea(NWN)
       real otota(NWN)
       real odxtota(NWN)

        character (len = *), parameter :: UNITS = "units"


       character (len = 40) :: NETCDF_FILE_NAME

       integer*4, save :: kount
       save cmol, id_mol

       DATA HMOLC / &
             '  H2O   ' , '  CO2   ' , '   O3   ' , '  N2O   ' ,&
             '   CO   ' , '  CH4   ' , '   O2   ' , '   NO   ' ,&
             '  SO2   ' , '  NO2   ' , '  NH3   ' , ' HNO3   ' ,&
             '   OH   ' , '   HF   ' , '  HCL   ' , '  HBR   ' ,&
             '   HI   ' , '  CLO   ' , '  OCS   ' , ' H2CO   ' ,&
             ' HOCL   ' , '   N2   ' , '  HCN   ' , ' CH3CL  ' ,&
             ' H2O2   ' , ' C2H2   ' , ' C2H6   ' , '  PH3   ' ,&
             ' COF2   ' , '  SF6   ' , '  H2S   ' , ' HCOOH  ' ,&
             '  HO2   ' , '   O+   ' , ' ClONO2 ' , '   NO+  ' ,&
             '  HOBr  ' , ' C2H4   ' , ' CH3OH  '/


       call getPhysConst(CLIGHT=CLIGHT)

    ! set up headers: assumes same molecules used in all profiles!!!

       if (npr.eq.1) then
          if (nmol.lt.22) wkl(22,:) = wbrodl
          wk_tot = sum(wkl,2)
          kount = 0
          do im=1,mxmol
             if (wk_tot(im).gt.0) then
                kount=kount+1
                id_mol(kount) = im
                cmol(kount) = hmolc(im)
             end if
          end do
       end if

!---Write header in output file
       WRITE(IOT,'(a)') 'MONORTM RESULTS:'
       WRITE(IOT,'(a)') '----------------' 
       WRITE(IOT,'(a5,I8,101x,a42)') 'NWN :',NWN , &
           'Molecular Optical Depths -->'

       ! Use GHz for small wavenumbers
       if (wn(1).lt.100) giga = .true.

       if (giga) then
          wnunits = 'FREQ(GHz)'
       else
          wnunits = 'FREQ(cm-1)'
       end if

       write(iot,11)'PROF ',wnunits,'BT(K) ','TMR(K)','  RAD(W/cm2_ster_cm-1)', &
                    'TRANS','PWV','CLW','TBOUND','EMIS','REFL','ANGLE',  &
                    'TOTAL_OD',cmol(1:kount), 'XSEC_OD'

       otot_by_mol(:,:) = 0.
       ! Convert to GHz for small wavenumbers
       DO Iw=1,NWN
          if (giga) then
             FREQ=WN(Iw)*CLIGHT/1.E9
          else
             FREQ=WN(Iw)
          end if
             
  !----Computation of the integrated optical depths
          otot  = 0.
          odxtot = 0.


          DO J = 1,NLAY
             OTOT = OTOT + O(IW,J)
             ODXTOT = ODXTOT + ODXSEC(IW,J)
             DO IK = 1,KOUNT
                OTOT_BY_MOL(IK,IW) =  OTOT_BY_MOL(IK,IW) + &
                   O_BY_MOL(IW,ID_MOL(IK),J) + OC(IW,ID_MOL(IK),J)
             ENDDO
          ENDDO

#ifdef USENETCDF

          ! Fill the arrays for netcdf output
          freqa(iw) = FREQ
          wvcolmna(iw) = wvcolmn
          clwcolmna(iw) = clwcolmn
          tmpsfca(iw) = tmpsfc
          anglea(iw) = angle
          otota(iw) = otot
          odxtota(iw) = odxtot


#endif          

          WRITE(IOT,21) NPR,FREQ,TB(Iw),TMR(iw),RAD(Iw),TRTOT(Iw), &
              WVCOLMN,CLWCOLMN,TMPSFC,EMISS(Iw),REFLC(Iw), &
              ANGLE,OTOT,otot_by_mol(1:kount,iw), odxtot
       ENDDO

       IF (IOD.EQ.1) THEN ! write out layer optical depths to ascii files
          IOL = 11
          DO J=1,NLAY
             WRITE(FILEOD,31) 'ODmono_prf', NPR, '_lay', J
             OPEN(IOL,FILE=FILEOD,STATUS='UNKNOWN')
             WRITE(IOL,'(a5,I8)') 'NWN :',NWN 
             WRITE(IOL,'(2a10)') wnunits, ' LAYER_OD'
             DO IW=1,NWN
                if (giga) then
                   FREQ=WN(Iw)*CLIGHT/1.E9
                else
                   FREQ=WN(Iw)
                end if
                WRITE(IOL,'(f10.3,e12.4)') FREQ, O(IW,J)
             ENDDO
             CLOSE(IOL)
          ENDDO
       ENDIF
    
! Conditional compilation of NETCDF code
! Requires preprocessing
#ifdef USENETCDF

       ! allocate the array to hold the sum of O_BY_MOL and OC
       allocate(O_BY_MOL_LAYER(NWNMX,MXMOL,MXLAY))
       O_BY_MOL_LAYER = O_BY_MOL + OC

       ! Build the NETCDF file name
       write(NETCDF_FILE_NAME,'(A,i5.5,A)') 'MONORTM.', NPR, '.nc'

! Write out netcdf file
       call checkm(nf90_create(NETCDF_FILE_NAME, NF90_CLOBBER, ncid))
       call checkm(nf90_def_dim(ncid, "FREQUENCY", NWN, wn_dimid))
       call checkm(nf90_def_dim(ncid, "MOLECULE", kount, mol_dimid))
       call checkm(nf90_def_dim(ncid, "LAYERS", nlay, lay_dimid))
       call checkm(nf90_def_dim(ncid, "STRING_LENGTH", int(8,4), s_dimid))

! Setup arrays for multi-dimensional output       
       dimids(1) = mol_dimid
       dimids(2) = wn_dimid

       cdimids(1) = s_dimid
       cdimids(2) = mol_dimid

       odimids(1) = wn_dimid
       odimids(2) = lay_dimid

       fdimids(1) = wn_dimid
       fdimids(2) = mol_dimid
       fdimids(3) = lay_dimid

! Define the variables
       call checkm(nf90_def_var(ncid, "FREQUENCY", nf90_xtype(freqa(1)), wn_dimid, varid(1)))
       call checkm(nf90_def_var(ncid, "BT", nf90_xtype(tb(1)), wn_dimid, varid(2)))
       call checkm(nf90_def_var(ncid, "RAD", nf90_xtype(rad(1)), wn_dimid, varid(3)))
       call checkm(nf90_def_var(ncid, "TRANS", nf90_xtype(trtot(1)), wn_dimid, varid(4)))
       call checkm(nf90_def_var(ncid, "PWV", nf90_xtype(wvcolmna(1)), wn_dimid, varid(5)))
       call checkm(nf90_def_var(ncid, "CLW", nf90_xtype(clwcolmna(1)), wn_dimid, varid(6)))
       call checkm(nf90_def_var(ncid, "SFCT", nf90_xtype(tmpsfca(1)), wn_dimid, varid(7)))
       call checkm(nf90_def_var(ncid, "EMIS", nf90_xtype(emiss(1)), wn_dimid, varid(8)))
       call checkm(nf90_def_var(ncid, "REFL", nf90_xtype(reflc(1)), wn_dimid, varid(9)))
       call checkm(nf90_def_var(ncid, "ANGLE", nf90_xtype(anglea(1)), wn_dimid, varid(10)))
       call checkm(nf90_def_var(ncid, "TMR", nf90_xtype(tmr(1)), wn_dimid, varid(11)))
       call checkm(nf90_def_var(ncid, "TOTAL_OD", nf90_xtype(otota(1)), wn_dimid, varid(12)))
       call checkm(nf90_def_var(ncid, "TOTAL_OD_BY_MOLECULE", nf90_xtype(otot_by_mol(1,1)), dimids, varid(13)))
       call checkm(nf90_def_var(ncid, "XSEC_OD", nf90_xtype(odxtota(1)), wn_dimid, varid(14)))
       call checkm(nf90_def_var(ncid, "MOLECULE", NF90_CHAR, cdimids, varid(15)))
       call checkm(nf90_def_var(ncid, "LAYER_OPTICAL_DEPTH", nf90_xtype(o(1,1)), odimids, varid(16)))
       call checkm(nf90_def_var(ncid, "LAYER_OPTICAL_DEPTH_BY_MOLECULE", NF90_FLOAT, fdimids, varid(17)))

! Define the dimensions 
       call NF_PUT_ATT_TEXT  (ncid, varid(1),"units", 11 ,wnunits)
       ! call NF_PUT_ATT_TEXT  (ncid, varid(3),"units", 18 ,"(W/cm2_ster_cm-1)")
       ! call NF_PUT_ATT_TEXT  (ncid, varid(10),"units", 7 ,"degrees")

! We need to end definition mode       
       call checkm(nf90_enddef(ncid))

! Write out the arrays to the variables in the netcdf file     
       call checkm(nf90_put_var(ncid, varid(1), freqa(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(2), tb(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(3), rad(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(4), trtot(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(5), wvcolmna(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(6), clwcolmna(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(7), tmpsfca(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(8), emiss(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(9), reflc(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(10), anglea(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(11), tmr(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(12), otota(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(13), otot_by_mol(1:kount, 1:NWN)))
       call checkm(nf90_put_var(ncid, varid(14), odxtota(1:NWN)))
       call checkm(nf90_put_var(ncid, varid(15), cmol(1:kount)))
       call checkm(nf90_put_var(ncid, varid(16), o(1:NWN, 1:nlay)))
       call checkm(nf90_put_var(ncid, varid(17), O_BY_MOL_LAYER(1:NWN, 1:kount, 1:nlay)))

! Close the netcdf file       
       call checkm(nf90_close(ncid))
! Clean up dynamic memory usage       
       deallocate(O_BY_MOL_LAYER)
#endif
 
 11    format (a5,a10,2a11,a22,a8,2a8,3a8,a9,36a12)
 21    format (i5,f10.3,2f11.5,1p,E21.9,0p,f9.5,2f8.4,3f8.2,f9.3,&
                                                  1p,36E12.4)
 31    format (a10,i4.4,a4,i4.4) 
       RETURN
 1000  WRITE(*,*) 'ERROR OPENING FILE:',FILEOUT
       STOP
       END

#ifdef USENETCDF
  subroutine checkm(status)
! This interface is used to check the status of netcdf calls.
! It display an error message and exit the program if an error occurs.
    USE netcdf
    integer*4, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, nf90_strerror(status)
      stop "Stopped"
    end if
    end subroutine checkm
#endif

  SUBROUTINE CORR_OPTDEPTH(INP,NLAY,SECNTA,NWN,ANGLE,O,IRT)
  USE PhysConstants, ONLY: getPhysConst
  USE RTMmono, ONLY: NWNMX
  !include "declar.incl"
  USE lblparams, ONLY: MXLAY
  INTEGER INP,J,NLAY,NWN,I,IRT
  REAL SECNT,ALPHA,ANGLE
  REAL PI
  REAL O(:,:),SECNTA(MXLAY)
  CALL getPhysConst(PI=PI)
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
!       value from vpayne 2006/07/23
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


      SUBROUTINE GETPROFNUMBER(IATM,FILEIN,fileprof, &
           NPROF)
      !include "declar.incl"
      CHARACTER FILEIN*60,CXID*1
      CHARACTER CDOL*1,CPRCNT*1,fileprof*80
      CHARACTER fileARMlist*64
      CHARACTER*8      HMOD                           
      DATA CDOL / '$'/,CPRCNT / '%'/
      NPROF=0

!---Get first the IATM info
      OPEN (90,FILE=FILEIN,STATUS='OLD',ERR=1000) 
 40   READ (90,'(a1)',END=80) CXID
      IF (CXID.NE.CDOL) THEN
         GO TO 40     
      ENDIF
      READ (90,'(49X,I1,19x,I1)',END=80,ERR=6000) IATM,ixsect 
      close(90)

      IF (IATM.EQ.1) THEN
         OPEN (90,FILE=FILEIN,STATUS='OLD',ERR=1000) 
 20      READ (90,'(a1)',END=80) CXID  
         IF (CXID.EQ.CDOL) NPROF=NPROF+1
         GO TO 20
      ENDIF
      IF (IATM.EQ.0) THEN
         open(90,file=fileprof,status='old',form='formatted',err=1000)
         DO WHILE (.true.) 
            READ (90,972,end=70,ERR=22) IFORM,LMAX,NMOL,SECNT0,HMOD, &
                  HMOD,H1,H2,ANGLE,LEN  
            NPROF=NPROF+1
 22         CONTINUE
         ENDDO
!  If there are cross-sections, there will be two records consistent with Format
!  statement 972 for each profile.
 70      continue
         if (ixsect .eq. 1) nprof = nprof/2
      ENDIF

 80   IF (NPROF.EQ.0) THEN 
         WRITE(*,*) 'NO PROFILE FOUND IN GETPROFNUMBER'
         STOP
      ELSE
         CLOSE(90)
         RETURN
      ENDIF
 2000 WRITE(*,*) 'ERROR OPENING ARM FILE'
 1000 WRITE(*,*) 'ERROR OPENING OR READING FILE in GETPROFNUMBER'
 924  FORMAT (1X,I1,I3,I5,F10.6,A24) 
 972  FORMAT(1X,I1,I3,I5,F10.6,2A8,4X,F8.2,4X,F8.2,5X,F8.3,5X,I2) 
 6000 WRITE(*,*) ' EXIT; ERROR READING :',FILEIN          
      STOP
      END

  SUBROUTINE CHECKINPUTS(NWN,NPROF,NWNMX)
  INTEGER NWN,NPROF,NWNMX
  IF (NWN.GT.NWNMX) THEN
     WRITE(*,*) 'Number of wavenumbers too big:',NWN
     WRITE(*,*) 'Maximum Allowed:',NWNMX
     WRITE(*,*) 'Must extend NWNMX in declar.incl AND RECOMPILE'
     STOP
  ENDIF
  RETURN
  END

!******************************************************************************
!___________________________________________________________________
!___________________________________________________________________

     subroutine profil_scal_sub(nlayrs)

     !include "declar.incl"
     USE lblparams, ONLY: MXMOL,MXLAY
     REAL WKL(MXMOL,MXLAY)
     REAL, DIMENSION(MXLAY) :: P,T,WBRODL
     REAL,    DIMENSION(MXLAY)   :: DVL,WTOTL,ALBL,ADBL,AVBL,H2OSL,SECNTA
     INTEGER, DIMENSION(MXLAY)   :: IPATH,ITYL
     REAL,    DIMENSION(0:MXLAY) :: ALTZ,PZ,TZ
     CHARACTER*4                 :: HT1,HT2
     REAL*8 V1,V2,SECANT,XALTZ 
     character*1 hmol_scal
     character*10 holn2
     character*8 XID,HMOLID,YID,HDATE,HTIME
     

     COMMON /PATHD/ P,T,WKL,WBRODL,DVL,WTOTL,ALBL,ADBL,AVBL, &
          H2OSL,IPTH,ITYL,SECNTA,HT1,HT2,ALTZ,PZ,TZ
     common /profil_scal/ nmol_scal,hmol_scal(64),xmol_scal(64)

     COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4), &
          WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
          EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF   
     COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, &
          NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL, &
          NLTEFL,LNFIL4,LNGTH4                                 

     DIMENSION WMT(64)

!  *** It should be noted that no attempt has been made to keep the 
!      mass in a given layer constant, i.e. level pressure nor retained ***

!      obtain accumulated amounts by molecule

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

!        obtain dry air sum
!             check to see if nitrogen is included in the selected molecules

     if (nmol.ge.22) then
        wsum_drair = 0.
     else
        wsum_drair = wsum_brod
     endif

     do m = 2, nmol
        wsum_drair = wsum_drair + wmt(m)
     enddo

     write (ipr,*)
     write (ipr,*) '   ', &
          '*****************************************************'
     write (ipr,*)
     write (ipr,*) '               Profile Scaling          '  

     write (ipr,956) 
 956 format (/,4x,' molecule',2x, &
           'hmol_scale',3x, ' xmol_param ',3x, 'scale factor',/)

     do m = 1, nmol_scal

        xmol_scal_m = xmol_scal(m)
        if (hmol_scal(m).eq.' ') xmol_scal(m) = 1.
        if (hmol_scal(m).eq.'0') xmol_scal(m) = 0.
        if (hmol_scal(m).eq.'1') xmol_scal(m) = xmol_scal_m     ! Scale factor
        if (hmol_scal(m).eq.'C' .or. hmol_scal(m).eq.'c') &     ! Column Amount (molec/cm^2)
              xmol_scal(m) = xmol_scal_m/wmt(m)     

        if (hmol_scal(m).eq.'M' .or. hmol_scal(m).eq.'m') then  ! Mixing ratio (molec/molec(dry_air))
           if (wsum_drair.gt.0.) then
              xmol_scal(m) = xmol_scal_m/(wmt(m)/wsum_drair)      
           else
              stop 'mixing ratio failure: wsum_drair = 0.'
           endif
        endif

        if (hmol_scal(m) .eq. 'P' .or. hmol_scal(m) .eq.'p') then  ! PWV for water vapor (cm)
           if (m.eq.1) then 
!                value from vpayne 2006/07/23
              xmol_scal(1) = (xmol_scal_m/2.99150e-23)/wmt(1)
           else
              write (ipr,*) 'm = ', m
              stop ' (hmol_scal(m).eq."P" .and. m.ne.1) '
           endif
        endif

        if (hmol_scal(m).eq.'D' .or. hmol_scal(m).eq.'d') & ! Dobson Units (du)
              xmol_scal(m) =  (xmol_scal_m*2.68678e16)/wmt(m)

        write (ipr,957) m, hmol_scal(m), xmol_scal_m, xmol_scal(m)
 957    format (5x,i5,9x,a1,5x,1p, 4e15.7)

!                scale the amounts and recalculate the total

        wmt(m) = 0.
        do l = 1, nlayrs
           wkl(m,l) = wkl(m,l) * xmol_scal(m)
           wmt(m)   = wmt(m) +wkl(m,l)
        enddo
     enddo

     write (ipr,*)
     write (ipr,*) '   ', &
          '*****************************************************'
     write (ipr,*)

!       write  modified column amounts to ipr in lblatm format
!
     WRITE (IPR,970)    
!
!     --------------------------------------------------------------
!
!     Write out column densities for molecules to TAPE6
!
     iform = 1

     data holn2/'  OTHER  '/
!
     IF (IFORM.EQ.1) THEN                                             
        WRITE (IPR,974) (HMOLID(I),I=1,7),HOLN2                       
        DO L = 1, NLAYRS
           WRITE (IPR,980) L,P(L),T(L),(WKL(M,L),M=1,7),WBRODL(L)
        enddo
        IF (NLAYRS.GT.1) THEN                                         
           WRITE (IPR,985)                                            
           L = NLAYRS                                                 
           WRITE (IPR,980) L,PWTD,TWTD,(WMT(M),M=1,7),wsum_brod
        ENDIF                                                       
     ELSE
        WRITE (IPR,975) (HMOLID(I),I=1,7),HOLN2
        DO  L = 1, NLAYRS
            WRITE (IPR,982) L,P(L),T(L),(WKL(M,L),M=1,7),WBRODL(L)
        enddo
        IF (NLAYRS.GT.1) THEN
           WRITE (IPR,985)
           L = NLAYRS
           WRITE (IPR,991) L,PWTD,TWTD,(WMT(M),M=1,7),wsum_brod
        ENDIF
     ENDIF

     IF (NMOL.GT.7) THEN                                            
        DO MLO = 8, NMOL, 8                                         
           MHI = MLO+7                                             
           MHI = MIN(MHI,NMOL)                              
           WRITE (IPR,970)                                  
           IF (IFORM.EQ.1) THEN
               WRITE (IPR,974) (HMOLID(I),I=MLO,MHI)            
               DO L = 1, NLAYRS                                 
                  WRITE (IPR,980) L,P(L),T(L),(WKL(M,L),M=MLO,MHI)
               enddo
               IF (NLAYRS.GT.1) THEN
                  WRITE (IPR,985)
                  L = NLAYRS
                  WRITE (IPR,990) L,PWTD,TWTD,(WMT(M),M=MLO,MHI)
              ENDIF
           ELSE
              WRITE (IPR,975) (HMOLID(I),I=MLO,MHI)              
              DO L = 1, NLAYRS                                   
                 WRITE (IPR,982) L,P(L),T(L),(WKL(M,L),M=MLO,MHI)
              enddo
              IF (NLAYRS.GT.1) THEN                              
                 WRITE (IPR,985)
                 L = NLAYRS
                 WRITE (IPR,991) L,PWTD,TWTD,(WMT(M),M=MLO,MHI)
              ENDIF
           ENDIF
        enddo                                                 
     ENDIF                                                       

!     --------------------------------------------------------------
!
!    Write out mixing ratios for molecules to TAPE6 in either
!    15.7 format (IFORM = 1) or 10.4 format (IFORM = 0).
!
!           Reset WDRAIR(L) for each layer
!           (WKL(M,L) now in column density)
!

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
              WRITE (IPR,980) L,P(L),T(L), &
                  (WKL(M,L)/WDRAIR_l,M=1,7),WBRODL(L)
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
               WRITE (IPR,982) L,P(L),T(L), &
                   (WKL(M,L)/WDRAIR_l,M=1,7),WBRODL(L)
           ENDIF
        enddo
     ENDIF

     IF (NMOL.GT.7) THEN
        DO MLO = 8, NMOL, 8
           MHI = MLO+7
           MHI = MIN(MHI,NMOL)
           IF (NLAYRS.LT.5) THEN
              WRITE (IPR,970)
!       ELSE
!       WRITE (IPR,945) XID,(YID(M),M=1,2)
          ENDIF
          IF (IFORM.EQ.1) THEN
             WRITE (IPR,976) (HMOLID(I),I=MLO,MHI)
             DO L = 1, NLAYRS
                IF (WDRAIR_l.EQ.0.0) THEN
                   WRITE(IPR,979)
                ELSE
                   WRITE (IPR,980) L,P(L),T(L), &
                       (WKL(M,L)/WDRAIR_l,M=MLO,MHI)
                ENDIF
             enddo
          ELSE
             WRITE (IPR,977) (HMOLID(I),I=MLO,MHI)
             DO L = 1, NLAYRS
                IF (WDRAIR_l.EQ.0.0) THEN
                   WRITE(IPR,979)
                ELSE
                   WRITE (IPR,982) L,P(L),T(L), &
                       (WKL(M,L)/WDRAIR_l,M=MLO,MHI)
                ENDIF
             enddo
          ENDIF
       enddo
    ENDIF

  970 FORMAT (////)
  974 FORMAT ('0',53X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',/,13X, &
        'P(MB)',6X,'T(K)',5X,8(A10,5X))
  975 FORMAT ('0',53X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',/,10X, &
        'P(MB)',6X,'T(K)',3X,8(1X,A6,3X))
  976 FORMAT (/,'1',54X,'----------------------------------', &
        /,'0',60X,'MIXING RATIOS BY LAYER ',/,10X, &
        'P(MB)',6X,'T(K)',5X,8(A10,5X))
  977 FORMAT (/,'1',54X,'----------------------------------', &
        /,'0',60X,'MIXING RATIOS BY LAYER ',/,210X, &
        'P(MB)',6X,'T(K)',3X,8(1X,A6,3X))
  979 FORMAT (/,'0','  MIXING RATIO IS UNDEFINED. DRYAIR DENSITY=0.0')
  980 FORMAT ('0',I3,F15.7,F9.2,2X,1P,8E15.7,0P)
  982 FORMAT ('0',I3,F12.5,F9.2,2X,1P,8E10.3,0P)
  985 FORMAT ('0',54X,'ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH')
  990 FORMAT ('0',I3,F15.7,F9.2,2X,1P,8E15.7,0P,/, &
        55X,1P,8E15.7,0P)
  991 FORMAT ('0',I3,F12.5,F9.2,2X,1P,8E10.3,0P)
!     --------------------------------------------------------------

return

end
!___________________________________________________________________
!___________________________________________________________________

  SUBROUTINE EXPINT (X,X1,X2,A)                                   
!********************************************************************
!       THIS SUBROUTINE EXPONENTIALLY INTERPOLATES X1 AND X2 TO X BY  
!       THE FACTOR A             . NEEDED BY LBLATM                   
!********************************************************************
  IF (X1.EQ.0.0.OR.X2.EQ.0.0) GO TO 10                          
  X = X1*(X2/X1)**A                                        
  RETURN                                                   
 10 X = X1+(X2-X1)*A                                         
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

      USE lblparams, ONLY: MX_XS,MXLAY
      IMPLICIT REAL*8           (V)

!********************************************************************
!     THIS SUBROUTINE READS IN THE DESIRED "CROSS-SECTION" 
!     MOLECULES WHICH ARE THEN MATCHED TO THE DATA CONTAINED
!     ON INPUT FILE FSCDXS.
!********************************************************************
      !include "declar.incl"
!
!     IFIL CARRIES FILE INFORMATION
!
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, &
              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL, &
              NLTEFL,LNFIL4,LNGTH4
!     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE 
!     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES 
!     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR 
!     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD 
!     MOLECULES. 
!
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(mx_xs),XAMNT(mx_xs,MXLAY)
!
!     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES
!     FOR THE CROSS-SECTION MOLECULES.
!     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES
!
!%%%%%LINUX_PGI90 (-i8)%%%%%      integer*4 iostat
      CHARACTER*10 XSFILE,XSNAME,ALIAS,XNAME,XFILS(6),BLANK
      COMMON /XSECTF/ XSFILE(6,5,mx_xs),XSNAME(mx_xs),ALIAS(4,mx_xs)
      COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS), &    
                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),  &   
                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS), & 
                NUMXS,IXSBIN                                   
      DIMENSION IXFLG(mx_xs)
      CHARACTER*120 XSREC
      CHARACTER*1 CFLG,CASTSK,CPRCNT,CFRM,CN,CF
      EQUIVALENCE (CFLG,XSREC)
      DATA CASTSK / '*'/,CPRCNT / '%'/,CN / 'N'/,CF / 'F'/
      DATA BLANK / '          '/
!
!     T296 IS TEMPERATURE FOR INITAL CALCULATIN OF DOPPLER WIDTHS
!
      DATA T296 / 296.0 /
!
      IXMAX = mx_xs
      DO 10 I = 1, IXMAX
         XSNAME(I) = BLANK
   10 CONTINUE
!
!     READ IN THE NAMES OF THE MOLECULES
!
      IF (IXMOLS.GT.7) THEN
         READ (Ipf,'(7A10)') (XSNAME(I),I=1,7)
         READ (Ipf,'(8A10)') (XSNAME(I),I=8,IXMOLS)
      ELSE
         READ (Ipf,'(7A10)') (XSNAME(I),I=1,IXMOLS)
      ENDIF

!
!     Left-justify all input names
!
      DO 15 I=1,IXMOLS
         CALL CLJUST (XSNAME(I),10)
 15   CONTINUE
!PRT  WRITE(IPR,'(/,''  THE FOLLOWING MOLECULES ARE REQUESTED:'',//,
!PRT 1    (5X,I5,2X,A1))') (I,XSNAME(I),I=1,IXMOLS)
!
!     MATCH THE NAMES READ IN AGAINST THE NAMES STORED IN ALIAS
!     AND DETERMINE THE INDEX VALUE.  STOP IF NO MATCH IS FOUND.
!     NAME MUST BE ALL IN CAPS.
      DO 40 I = 1, IXMOLS
         DO 20 J = 1, IXMAX

            IF ((XSNAME(I).EQ.ALIAS(1,J)) .OR. &
                (XSNAME(I).EQ.ALIAS(2,J)) .OR. &
                (XSNAME(I).EQ.ALIAS(3,J)) .OR. &
                (XSNAME(I).EQ.ALIAS(4,J))) THEN
               IXINDX(I) = J
               GO TO 30
            ENDIF
 20      CONTINUE
!
!         NO MATCH FOUND
!
         WRITE (IPR,900) XSNAME(I)
         STOP 'STOPPED IN XSREAD'
 30      CONTINUE
         IXFLG(I) = 0
 40   CONTINUE
!
!     READ IN "CROSS SECTION" MASTER FILE FSCDXS
!
      IXFIL = 8
      OPEN (IXFIL,FILE='FSCDXS',STATUS='OLD',FORM='FORMATTED', &
            IOSTAT=iostat)
      if (IOSTAT.gt.0) stop 'FSCDXS does not exist - XSREAD'
      REWIND IXFIL
      READ (IXFIL,905)
!
 50   READ (IXFIL,910,END=80) XSREC
!

      IF (CFLG.EQ.CASTSK) GO TO 50
      IF (CFLG.EQ.CPRCNT) GO TO 80
!
      READ (XSREC,915) XNAME,V1X,V2X,DVX,NTEMP,IFRM,CFRM, &
                 (XFILS(I),I=1,NTEMP)
!
!     LEFT-JUSTIFY INPUTED NAME
!
      CALL CLJUST (XNAME,10)
!     CHECK MASTER FILE FOR CROSS SECTION AND STORE DATA
      NUMXS = IXMOLS
      DO 70 I = 1, IXMOLS
         IF ((XNAME.EQ.ALIAS(1,IXINDX(I))) .OR. &
             (XNAME.EQ.ALIAS(2,IXINDX(I))) .OR. &
             (XNAME.EQ.ALIAS(3,IXINDX(I))) .OR. &
             (XNAME.EQ.ALIAS(4,IXINDX(I)))) THEN
            IXFLG(I) = 1
            IF (V2X.GT.XV1.AND.V1X.LT.XV2) THEN
               NSPECR(I) = NSPECR(I)+1
               IF (NSPECR(I).GT.6) THEN  
                  WRITE (IPR,920) I,XSNAME(I),NSPECR(I)
                  STOP ' XSREAD - NSPECR .GT. 6'
               ENDIF
               IXFORM(NSPECR(I),I) = 91
               IF (IFRM.EQ.86) IXFORM(NSPECR(I),I) = IFRM
               IF (CFRM.NE.CN) &
                   IXFORM(NSPECR(I),I) = IXFORM(NSPECR(I),I)+100
               IF (CFRM.EQ.CF) &
                   IXFORM(NSPECR(I),I) = -IXFORM(NSPECR(I),I)
               NTEMPF(NSPECR(I),I) = NTEMP
               V1FX(NSPECR(I),I) = V1X
               V2FX(NSPECR(I),I) = V2X
!
!     3.58115E-07 = SQRT( 2.* LOG(2.)*AVOGAD*BOLTZ/(CLIGHT*CLIGHT) ) 
!
               XDOPLR(NSPECR(I),I)=3.58115E-07*(0.5*(V1X+V2X))* &
                             SQRT(T296/XSMASS(IXINDX(I)))
!
               DO 60 J = 1, NTEMP
                  XSFILE(J,NSPECR(I),I) = XFILS(J)
 60            CONTINUE
            ENDIF
         ENDIF
 70   CONTINUE
!
      GO TO 50
!
 80   IXFLAG = 0
      DO 90 I = 1, IXMOLS
         IF (IXFLG(I).EQ.0) THEN
            WRITE (IPR,925) XSNAME(I)
            IXFLAG = 1
         ENDIF
 90   CONTINUE
      IF (IXFLAG.EQ.1) STOP ' IXFLAG - XSREAD '
!
      RETURN
!
  900 FORMAT (/,'  THE NAME: ',A10, ' IS NOT ONE OF THE ', &
        'CROSS SECTION MOLECULES. CHECK THE SPELLING.')
  905 FORMAT (/)
  910 FORMAT (A120)
  915 FORMAT (A10,2F10.4,F10.8,I5,5X,I5,A1,4X,6A10)
  920 FORMAT (/,'******* ERROR IN XSREAD ** MOLECULE SECLECTED -',A10, &
        '- HAS ',I2,' SPECTRAL REGIONS ON FILE FSCDXS, BUT THE', &
        ' MAXIMUM ALLOWED IS 6 *******',/)
  925 FORMAT (/,'******* MOLECULE SELECTED -',A10,'- IS NOT FOUND ON', &
        ' FILE FSCDXS *******',/)
!

      RETURN
      END

      BLOCK DATA BXSECT
      USE lblparams, ONLY: MX_XS
      IMPLICIT REAL*8           (V)

!include "declar.incl"
!      parameter (mx_xs=38)
!
!**   XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES
!**            (NOTE: ALL NAMES ARE LEFT-JUSTIFIED)
!
      CHARACTER*10 XSFILE,XSNAME,ALIAS
      COMMON /XSECTI/ XSMAX(6,5,mx_xs),XSTEMP(6,5,mx_xs), &
                NPTSFX(5,mx_xs),NFILEX(5,mx_xs),NLIMX 
      COMMON /XSECTF/ XSFILE(6,5,mx_xs),XSNAME(mx_xs),ALIAS(4,mx_xs) 
      COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS), &
                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),  &   
                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS), & 
                NUMXS,IXSBIN                                   
      COMMON /XSECTS/ JINPUT,NMODES,NPANEL,NDUM,V1XS,V2XS,DVXS,NPTSXS
!
      DATA NMODES / 1 /,NPANEL / 0 /,V1XS / 0.0 /,V2XS / 0.0 /, &
           DVXS / 0.0 /,NPTSXS / 0 /
      DATA XSMAX / 1140*0.0 /
      DATA (ALIAS(1,I),I=1,mx_xs)/ &
         'CLONO2    ', 'HNO4      ', 'CHCL2F    ', 'CCL4      ', &
         'CCL3F     ', 'CCL2F2    ', 'C2CL2F4   ', 'C2CL3F3   ', &
         'N2O5      ', 'HNO3      ', 'CF4       ', 'CHCLF2    ', &
         'CCLF3     ', 'C2CLF5    ', 'NO2       ', 23*' ZZZZZZZZ ' / 
      DATA (ALIAS(2,I),I=1,mx_xs)/ &
         'CLNO3     ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ', &
         'CFCL3     ', 'CF2CL2    ', 'C2F4CL2   ', 'C2F3CL3   ', &
         ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CHF2CL    ', &
         ' ZZZZZZZZ ', ' ZZZZZZZZ ', 24*' ZZZZZZZZ ' /
      DATA (ALIAS(3,I),I=1,mx_xs)/ &
         ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ', &
         'CFC11     ', 'CFC12     ', 'CFC114    ', 'CFC113    ', &
         ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC14     ', 'CFC22     ', &
         'CFC13     ', 'CFC115    ', 24*' ZZZZZZZZ ' /
      DATA (ALIAS(4,I),I=1,mx_xs)/ &
         ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F21       ', ' ZZZZZZZZ ', &
         'F11       ', 'F12       ', 'F114      ', 'F113      ', &
         ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F14       ', 'F22       ', &
         'F13       ', 'F115      ', 24*' ZZZZZZZZ ' /
!
!     XSMASS IS MASS OF EACH CROSS-SECTION
!
!     note 23 =  mx_xs - 15 = 38 - 15

      DATA XSMASS/ &
           97.46     ,   79.01     ,  102.92     ,  153.82     , &
          137.37     ,  120.91     ,  170.92     ,  187.38     , &
          108.01     ,   63.01     ,   88.00     ,   86.47     , &
          104.46     ,  154.47     ,   45.99     , 23*0.00 /
!
      DATA V1FX / 190*0.0 /,V2FX / 190*0.0 /,DVFX / 190*0.0 /, &
          WXM / mx_xs*0.0 /
      DATA NTEMPF / 190*0 /,NSPECR / mx_xs*0 /,IXFORM / 190*0 /, &
          NUMXS / 0 /
!
      END

      SUBROUTINE CLJUST (CNAME,NCHAR)
!
!     THIS SUBROUTINE LEFT-JUSTIFIES THE CHARACTER CNAME
!
      CHARACTER*(*) CNAME
      CHARACTER*25 CTEMP
      CHARACTER*1  CTEMP1(25),BLANK
      EQUIVALENCE (CTEMP,CTEMP1(1))
!
      DATA BLANK / ' '/
!
         CTEMP = CNAME
         JJ=0
         DO 10 J = 1, NCHAR
            IF (CTEMP1(J).NE.BLANK) THEN
               JJ = J
               IF (JJ.EQ.1) GO TO 50
               GO TO 20
            ENDIF
 10      CONTINUE
         IF (JJ .EQ. 0) GO TO 50
!
 20      KCNT = 0
         DO 30 K = JJ, NCHAR
            KCNT = KCNT+1
            CTEMP1(KCNT) = CTEMP1(K)
 30      CONTINUE
!
         KK = NCHAR-JJ+2
         DO 40 L = KK,NCHAR
            CTEMP1(L) = BLANK
 40      CONTINUE
         CNAME = CTEMP
 50   CONTINUE
!
      RETURN
!
      END

  FUNCTION MAXIMUM(ARRAY, NSIZE)
  DIMENSION ARRAY(NSIZE)
  MAXIMUM = ARRAY(1)
  DO 10 J=2,NSIZE
     IF (MAXIMUM.LT.ARRAY(J)) MAXIMUM=ARRAY(J)
 10 CONTINUE
  RETURN
  END

  FUNCTION MINIMUM(ARRAY, NSIZE)
  DIMENSION ARRAY(NSIZE)
  MINIMUM = ARRAY(1)
  DO 10 J=2,NSIZE
     IF (MINIMUM.GT.ARRAY(J)) MINIMUM=ARRAY(J)
 10 CONTINUE
  RETURN
  END

        SUBROUTINE MONORTM_XSEC_SUB(wn,nwn,p,t,nlay, odxsec)

!----------------------------------------------------------------
! Authors: Eli Mlawer and Vivienne Payne, AER Inc.
!
! Created: July 2008
!
! Description: Calculates cross-section optical depths for MonoRTM.
!              Method follows that of LBLRTM. 
!              Results are consistent with LBLRTM. 
!              Note that only the total optical depth for all 
!              cross-section molecules is returned.  The "odxsec" array 
!              does not distinguish between different xsec molecules.     
!
! Input:
!      wn       array of wavenumbers
!      nwn      number of wavenumbers 
!      p        pressure of layer
!      t        temperature of layer
!      nlay     number of layers
!
! Output:  
!      odxsec   total optical depth for all cross-sections in array odxsec
!               (stored through declar.incl)
!    
!----------------------------------------------------------------

      USE PhysConstants, ONLY: getPhysConst
      USE RTMmono, ONLY: NWNMX
      USE lblparams, ONLY: MX_XS,MXLAY
      IMPLICIT REAL*8           (V) ! for consistency with LBLRTM routines
      REAL*8 WN(NWNMX),P(MXLAY),T(MXLAY)
      !Include "declar.incl"


!                                                                         
!     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES         
!     FOR THE CROSS-SECTION MOLECULES.                                    
!     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          

!     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE    
!     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES       
!     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR     
!     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD         
!     MOLECULES.                                                          
!                                                                         
!     NUMXS IS THE NUMBER OF 'CROSS SECTION' MOLECULES TO BE USED         
!                                                                         
!     XSFILE(ITEMP,ISPEC,NXS) IS THE NAME OF THE FILE CONTAINING THE      
!                             'CROSS SECTION' DATA.  THE THREE INDICES    
!                             ARE DEFINED AS FOLLOWS:                     
!                                                                         
!                             ITEMP - DENOTES THE TEMPERATURE FOR WHICH   
!                                     THE 'CROSS SECTION' IS VALID        
!                                     (IMPLEMENTED FOR HITRAN 91 TAPE)    
!                             ISPEC - DENOTES THE SECTRAL REGION FOR      
!                                     WHICH THE FILE PERTAINS             
!                             NXS   - IS THE INCREMENT FOR THE 'CROSS     
!                                     SECTION' INDEX                      
!                                                                         
!     NTEMPF(ISPEC,NXS) IS THE NUMBER OF TEMPERATURE FILES TO BE USED     
!                       FOR EACH SPECTRAL REGION OF EACH MOLECULE         
!                                                                         
!     NSPECR(NXS) IS THE NUMBER OF SPECTRAL REGIONS FOR THE MOLECULE NX   
!**********************************************************************   
!                        
      CHARACTER*10 XSFILE,XSNAME,ALIAS,XNAME,XFILS(6),BLANK,xxfile,ctorr               
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(mx_xs),XAMNT(mx_xs,MXLAY)        
      COMMON /XSECTF/ XSFILE(6,5,mx_xs),XSNAME(mx_xs),ALIAS(4,mx_xs)
      COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS), &
                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),  &
                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS), &
                NUMXS,IXSBIN    
      REAL ODXSEC(NWNMX,MXLAY)
!
      REAL RADCN2
!                                                                         
      DIMENSION IXFLG(mx_xs),tx(6,5,mx_xs),pdx(6,5,mx_xs)
      dimension xsdat(150000,6),xspd(150000)
      dimension xspave(nwnmx),xsmoltot(nwnmx,mxlay),xstot(nwnmx,mxlay)
!                                                                         
      CHARACTER*120 XSREC
      character*10 source(3)

      data dvbuf /1.0/         ! used to check if a particular xs file needs to be processed
      DATA CTORR / '      TORR'/
!
      call getPhysConst(RADCN2=RADCN2)

!     DEFINE PRESSURE CONVERSIONS
!
!        PTORMB = 1013. MB / 760. TORR  (TORR TO MILLIBARS)
!        PATMMB = 1013. MB / 1.0  ATM   (ATMOPHERES TO MILLIBARS)
!
      PTORMB = 1013./760
      PATMMB = 1013.

! Initialize variables.
      xstot(1:nwn,1:nlay) = 0.

! loop over cross-section molecules
      DO 6000 IXMOL = 1,IXMOLS
          xsmoltot(1:nwn,1:nlay) = 0.

! loop over the spectral regions for this molecule
          DO 5000 IXSR = 1,NSPECR(IXMOL)
! Check to see if this spectral region needs to be processed.
              ipt = 1
 2100         continue
              if (wn(ipt) .ge. v1fx(ixsr,ixmol)-dvbuf .and. &
                   wn(ipt) .le. v2fx(ixsr,ixmol)+dvbuf) goto 2200
              ipt = ipt + 1
              if (ipt .gt. nwn) goto 5000
              goto 2100

 2200         continue

! Getting here means that this spectral region for this molecule needs to be handled.
! Loop over the temperatures for which this molecule has stored values in this region.

! The filenames are arranged in FSCDXS in order of ascending temperature
              DO 3000 IXTEMP=1,NTEMPF(IXSR,IXMOL)
                  ifile = 100 + ixtemp
                  xxfile = xsfile(ixtemp,ixsr,ixmol)
                  open(ifile,file=xxfile,form='formatted')
                  READ (ifile,910) AMOL,V1x,V2x,NPTSx, &
                        TX(ixtemp,ixsr,ixmol),PRES,    &
                        SMAX,SOURCE
                  if (source(3) .eq. ctorr) then
                      pdx(ixtemp,ixsr,ixmol) = pres * ptormb
                  else
                      pdx(ixtemp,ixsr,ixmol) = pres
                  endif
                  read(ifile,*)(xsdat(i,ixtemp),i=1,nptsx)
 3000         continue

! Interpolate absorption coefficients to given temperature.
              
              do 4000 il = 1, nlay
                  pave = p(il)
                  tave = t(il)
                  coef1 = 1.
                  ind2 = 1
                  coef2 = 0.
                  it = 1
! The temperatures at which the xsecs are stored are in ascending order (as in LBLRTM).
! it=1 is the lowest temperature
                  if (ntempf(ixsr,ixmol).eq.1 .or. &
                       tave.le.tx(it,ixsr,ixmol)) then
                     ind1 = 1
                  else
3100                 continue
                     it = it + 1
                     if (it .gt. ntempf(ixsr,ixmol)) then
                        ind1 = ntempf(ixsr,ixmol)
                        ind2 = ntempf(ixsr,ixmol)
                     elseif (tave .le. tx(it,ixsr,ixmol)) then
                        ind1 = it-1 ! ind1 is a lower temperature than ind2
                        ind2 = it
                        coef1 = (tave - tx(it,ixsr,ixmol))/ &
                             (tx(it-1,ixsr,ixmol) - tx(it,ixsr,ixmol))
                        coef2 = 1. - coef1
                     else
                        goto 3100
                     endif
                  endif

                  pd = coef1 * pdx(ind1,ixsr,ixmol)+ &
                      coef2* pdx(ind2,ixsr,ixmol)
                  xkt1 = tx(ind1,ixsr,ixmol)/radcn2
                  xkt2 = tx(ind2,ixsr,ixmol)/radcn2
                  
                  delvx =  (v2x-v1x) / float(nptsx-1)
                  do 3300 i = 1, nptsx
                      vv = v1x + float(i-1) * delvx
                      xspd(i) = coef1 * xsdat(i,ind1) / radfn(vv,xkt1) &
                         + coef2 * xsdat(i,ind2) / radfn(vv,xkt2)
 3300             continue

! Perform convolution from appropriate pressure for stored values to layer pressure.
! (296K is the temperature used in the calculation of the Doppler width)
                  hwdop = xdoplr(ixsr,ixmol) * sqrt(tave/296.)
                  call convolve(xspd,v1x,v2x,delvx,pd,hwdop, &
                     tave,pave,x,wn,nwn,xspave)
                  xsmoltot(1:nwn,il) = xsmoltot(1:nwn,il) + &
                     xspave(1:nwn)
 4000         continue
! Finished processing this spectral range of this molecule.

 5000     continue

! Multiply by absorber amount.
          do 5500 il = 1, nlay
              xstot(1:nwn,il) = xstot(1:nwn,il) + &
                 xamnt(ixmol,il) * xsmoltot(1:nwn,il)
 5500     continue
! Finished prcoessing this molecule.

 6000 continue
! Finished processing all molecules.
! Put the radiation field back in.
      do 6500 il = 1, nlay
          xkt = t(il)/radcn2
          do 6300 iwn = 1, nwn              
              odxsec(iwn,il) = xstot(iwn,il) * radfn(wn(iwn),xkt)
 6300     continue
 6500 continue

       RETURN  
  910 FORMAT (A10,2F10.4,I10,3G10.3,3A10)              

       END

      subroutine convolve(xspd,v1x,v2x,delvx,pd,hwdop, &
         tave,pave,x,wn,nwn,xspave)

! performs pressure convolution for cross sections

      USE RTMmono, ONLY: NWNMX
      !Include "declar.incl"
      dimension xspd(150000),xspave(nwnmx),xspd_int(0:10000000)
      real*8 wn(nwnmx)
      data p0 /1013./

!     Set up the halfwidths needed.
      hwpave = 0.1 * (pave/p0) * (273.15/tave)
      hwd = 0.1 * (pd/p0) * (273.15/tave)
      hwd = max (hwd,hwdop)
! Don''t want hwd to be bigger than hwpave
      if (hwd.gt.hwpave) hwpave = 1.001*hwd
      hwb = hwpave - hwd

! Step size is at least ~4 pts per halfwidth.
      ratio = 0.25
      step = ratio * hwb
      if (step .gt. delvx) step = delvx
      npts = int((v2x-v1x)/step)
      step = (v2x - v1x) / float(npts)
      ratio = step / hwb
! Linearly interpolate incoming values to desired grid.
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
! Perform convolution using Lorentzian with width hwb.
              wn_v1x = wn(iwn) - v1x
              ind = int(wn_v1x / step)
              dvlo = wn(iwn) - (v1x + float(ind) * step)
              dvhi = wn(iwn) - (v1x + float(ind+1) * step)
              answer = (hwb/(hwb2 + dvlo**2)) * xspd_int(ind) + &
                 (hwb/(hwb2 + dvhi**2)) * xspd_int(ind+1)
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
! Done with convolution.  Multiply by step size and Lorentzian normalization.
              xspave(iwn) = answer * step / 3.14159
          else
!     Use linearly interpolated values.
              wn_v1x = wn(iwn) - v1x
              ind = int(wn_v1x / delvx)
              coef = (wn_v1x - float(ind) * delvx) / delvx
              xspave(iwn) = (1.-coef)*xspd(ind) + coef*xspd(ind+1)
          endif
 4000 continue

      return
      
      end
