MODULE LNFL_MOD

   USE struct_types

   PARAMETER (MXMOL=39,IIM=  75000)

   TYPE(LINE_DATA)    :: BUFR

   INTEGER :: NBLM(mxmol),ISO(mxmol,IIM)
   REAL*8 :: XNU0(mxmol,IIM)
   REAL, dimension(mxmol,iim) ::  DELTNU,E,ALPS,ALPF,X,XG,S0,RMOL,SDEP
   integer*4, dimension(mxbrdmol,IIM) ::  brd_mol_flg(mxmol,mxbrdmol,IIM)
   REAL, dimension(mxmol,mxbrdmol,IIM) ::  brd_mol_tmp,brd_mol_hw,brd_mol_shft
   PRIVATE

   PUBLIC :: GET_LNFL
   PUBLIC :: NBLM,ISO,XNU0,DELTNU,E,ALPS,ALPF,X,XG,S0,Rmol, &
             brd_mol_flg,brd_mol_tmp,brd_mol_hw,brd_mol_shft,sdep

CONTAINS

      SUBROUTINE GET_LNFL(IPR,ICP,HFILE,ISPD,v1,v2)

      use struct_types, ONLY: line_data
      real*8  ::    v1,v2

      INTEGER, INTENT(IN) :: IPR
      REAL amol
      REAL*4 xmol
      INTEGER*4 linfil
      Character Q1*9,Q2*9,HFILE*80

      DATA I_1/1/, I_100/100/, I_1000/1000/

! Initialize
      nblm(1:mxmol) = 0
      ieof = 0

      linfil = 9
      OPEN(linfil,FILE=HFILE,form='unformatted',status='old',ERR=20)
      call prlnhd (linfil,ipr)

      do while (ieof.eq.0)
         call rdlnfl (v1,ipr,linfil,ieof,ilo,ihi)
         do ik=ilo,ihi
            select case (bufr%iflg(ik))
               case (0:100)
		  mo = mod(bufr%mol(ik),i_100)
               case (-3:-1)
		   mo = mod(bufr%mol(ik-1),i_100)
               case (-5)
        ! Special code for handling cases with both foreign and self LC: assumes two lines with LC info
                   if (bufr%iflg(ik-1).ge.0) then
        ! First line (foreign)
		      mo = mod(bufr%mol(ik-1),i_100)
                      mo_prev = mo
                   else
        ! Second line (foreign)
		      mo =  mo_prev
		  endif
               case default
                  print *,'LC flag not recongnized: ', bufr%iflg(ik),'. Must be 1, 3 or 5.'
                  stop
            end select
            nblm(mo) = nblm(mo)+1
            ii = nblm(mo)
            iso(mo,ii)    =  mod(bufr%mol(ik),i_1000)/100
            xnu0(mo,ii)   =  bufr%vnu(ik)
            s0(mo,ii)     =  bufr%sp(ik)
            alpf(mo,ii)   =  bufr%alfa(ik)
            alps(mo,ii)   =  bufr%hwhm(ik)
            e(mo,ii)      =  bufr%epp(ik)
            X(mo,ii)      =  bufr%tmpalf(ik)
            deltnu(mo,ii) =  bufr%pshift(ik)
            if (bufr%iflg(ik).ge.0) then
               XG(mo,ii) =  -1*bufr%iflg(ik)
            else
               XG(mo,ii) =  bufr%iflg(ik)
            endif
            xmol         =  transfer(bufr%mol(ik),xmol)     !int*4 to real*4
            amol         = xmol                             ! real*4 to real*8 (if compiled as r8)
            rmol(mo,ii) = transfer(amol,rmol(mo,ii))        ! real*8 

            brd_mol_flg(mo,:,ii)  = bufr%brd_mol_flg(:,ik)
            brd_mol_hw(mo,:,ii)   = bufr%brd_mol_hw(:,ik)
            brd_mol_tmp(mo,:,ii) =  bufr%brd_mol_tmp(:,ik)
            brd_mol_shft(mo,:,ii) = bufr%brd_mol_shft(:,ik)

            sdep(mo,ii) = bufr%speed_dep(ik)

            print *,mo,xnu0(mo,ii),sdep(mo,ii)
         end do
         if (bufr%vnu(ihi).GT.(v2+25.)) ieof=1
      end do


! keep this as a reminder that we need to deal with ISPD

     ! IF (ISPD .eq. 1) then
     !    IF (iflg.ne.mo) then 
     !       IF ((iref3.EQ.-1).OR.(iref3.EQ.-3)) THEN 
     !          READ(9,2,END=3,ERR=30)JF,XF,DF,SF,EF,AS,AF,XS,XF
     !       ENDIF
     !       goto 22
     !    ENDIF
     ! ENDIF

! Keep this as a reminder that LC can be ignored in MONORTM
!      IF (((iref3.EQ.-1).OR.(iref3.EQ.-3)).AND.(ICP.EQ.0)) THEN !Lines Cplng ignored
!         READ(9,2,END=3,ERR=30)JF,XF,DF,SF,EF,AS,AF,XS,XF
!         XG(JJ,II)=0.
!      ENDIF
      !IF (ISPD.EQ.1) THEN
      !   WRITE(*,*) '****************************************'
      !   WRITE(*,*) '*            W A R N I N G             *'
      !   WRITE(*,*) '****************************************'
      !   WRITE(*,*) 'FAST VERSION IS RUNNING.'
      !   WRITE(*,*) 'CURRENTLY VALID IN THE MICROWAVE ONLY.'
      !ENDIF
      WRITE(*,*)
      WRITE(*,*) '****************************************'
      CLOSE(LINFIL)
      RETURN
 20   PRINT *, 'ERROR OPENING HITRAN FILE:',HFILE
      STOP
      END SUBROUTINE GET_LNFL


       SUBROUTINE RDLNFL (VLO,IPR,LINFIL,IEOF,ILO,IHI) 
!                                                                       
      REAL*8  ::          VLO
      REAl*4  ::          XMOL, AMOL
!                                                                       
      integer*4 leof,npnlhd,linfil 
!                                                                       
      TYPE(INPUT_HEADER) :: RDLNPNL
      TYPE(INPUT_BLOCK)  :: RDLNBUF, DUMBUF

      real*4 dum(2)
      integer*4 i_1 
!
      DATA I_1/1/  
!                                                                       
      IPASS = 1 
      IF (ILO.GT.0) IPASS = 2 
!                                                                       
      ILO = 1 
      IHI = 0 
!                                                                       
      npnlhd = 6 

   10 CALL BUFIN_sgl(linfil,LEOF,RDLNPNL,npnlhd)

!                                                                       
      IF (LEOF.EQ.0) GO TO 30
      IF (RDLNPNL%VMAX.LT.VLO) THEN
         CALL BUFIN_sgl(linfil,LEOF,DUMBUF,i_1)
         GO TO 10
      ELSE
         CALL BUFIN_sgl(linfil,LEOF,rdlnbuf,RDLNPNL%NWDS)
      ENDIF
!                                                                       
      IF ((IPASS.EQ.1).AND.(RDLNBUF%VNU(1).GT.VLO)) WRITE (ipr,900)
!                                                                       
      IJ = 0
!                                                                       
!     precision conversion occurs here:                                 
!     incoming on right: vlin is real*8;  others are real*4 and integer*4
!                                                                       
      do 15 i=1,RDLNPNL%NREC
         BUFR%IFLG(i)  = RDLNBUF%IFLG(i)
         BUFR%VNU(i)   = RDLNBUF%VNU(i)
         BUFR%SP(i)    = RDLNBUF%SP(i)
         BUFR%ALFA(i)   = RDLNBUF%ALFA(i)
         BUFR%EPP(i)   = RDLNBUF%EPP(i)
	 BUFR%MOL(i)  = RDLNBUF%MOL(i)   ! int*4 to int*8
         BUFR%HWHM(i) = RDLNBUF%HWHM(i)
         BUFR%TMPALF(i)= RDLNBUF%TMPALF(i)
         BUFR%PSHIFT(i)= RDLNBUF%PSHIFT(i)
         do j=1,mxbrdmol
	    BUFR%BRD_MOL_FLG(j,i) = RDLNBUF%BRD_MOL_FLG_IN(j,i)
         end do
         j=1
         do j1=1,mxbrdmol
            BUFR%BRD_MOL_HW(j1,i)   = RDLNBUF%BRD_MOL_DAT(j,i)
            BUFR%BRD_MOL_TMP(j1,i)  = RDLNBUF%BRD_MOL_DAT(j+1,i)
            BUFR%BRD_MOL_SHFT(j1,i) = RDLNBUF%BRD_MOL_DAT(j+2,i)
            j=j+3
        end do
        BUFR%SPEED_DEP(i) = RDLNBUF%SPEED_DEP(i)
   15 continue
!                                                                       
      IHI = RDLNPNL%NREC
      RETURN
   30 WRITE (ipr,905)
      IEOF = 1
      RETURN
!                                                                       
  900 FORMAT ('0 FIRST LINE USED IN RDLNFL--- CHECK THE LINEFIL  ')
  905 FORMAT ('0 EOF ON LINFIL IN RDLNFL -- CHECK THE LINFIL ')
!                                                                       
      END SUBROUTINE RDLNFL

     SUBROUTINE PRLNHD (LINFIL,IPR) 
!                                                                       
      USE lblparams, ONLY: MXMOL
      IMPLICIT REAL*8           (V)
!                                                                       
!     PRLNHD PRINTS OUT LINE FILE HEADER                                
!                                                                       
!     -------------------------                                         
!                                                                       
!                                                                       
      integer *4 linmol,                                         &
                 lincnt,ilinlc,ilinnl,irec,irectl,               &
                 linfil
!                                                                       
      character*8 HLINID(10),BMOLID(64),HID1(2),HMOLID(60)
      real*4 SUMSTR(64)
      integer*4 MOLCNT(64),MCNTLC(64),MCNTNL(64)

      integer*4 n_negepp(64),n_resetepp(64)
      real*4 xspace(4096)
         
      real *4 flinlo,flinhi
      integer *4 lnfil
      integer *4 negepp_flag
!                                                                       
!     LSTWD (LAST WORD) IS DUMMY, DOES NOT NEED TO BE COUNTED           
!                                                                       
!                                                                       
      CHARACTER CHID10*8,CHARID*5,CHARDT*2,CHARI*1,CHTST*1
      CHARACTER*1 CNEGEPP(8)
      CHARACTER*6 CDUM,SPCRT
!                                                                       
!                                                                       
      DATA CHARI / 'I'/
!                                                                       
      REWIND LINFIL

      lnfil = linfil
      negepp_flag = 0

      read (lnfil,end=777)    HLINID,BMOLID,MOLCNT,MCNTLC,              &
     &                MCNTNL,SUMSTR,LINMOL,FLINLO,FLINHI,               &
     &                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1


!     Test for negative values of ENERGY identified in lnfl             
!     and read in second header for line information, if needed         

      READ (HLINID(7),950) CNEGEPP
      IF (CNEGEPP(8).eq.'^') THEN
         negepp_flag = 1
         read (lnfil) n_negepp,n_resetepp,xspace
         endif
!                                                                       
         go to 5
!                                                                       
  777    STOP 'LAYER; TAPE3 DOES NOT EXIST'
!                                                                       
    5    continue
!                                                                       
         DO 10 M = 1, LINMOL
            HMOLID(M) = BMOLID(M)
   10    END DO
         WRITE (IPR,900)
         WRITE (IPR,905) HLINID,HID1

!     Output header information regarding lines; if negative values of  
!     ENERGY were identified in lnfl, output extra header information   

         if (CNEGEPP(8).eq.'^') THEN
         WRITE (IPR,960)
         WRITE (IPR,965) (BMOLID(I),MOLCNT(I),MCNTLC(I),MCNTNL(I),      &
         N_NEGEPP(I),N_RESETEPP(I), SUMSTR(I),I=1,LINMOL)
         else
         WRITE (IPR,910)
         WRITE (IPR,915) (BMOLID(I),MOLCNT(I),MCNTLC(I),MCNTNL(I),      &
         SUMSTR(I),I=1,LINMOL)
         endif
!                                                                       
         WRITE (IPR,920) FLINLO,FLINHI,LINCNT
!                                                                       
!     When calculating derivative, check make sure the                  
!     appropriate molecule is included in the linefile.                 
!     If not, then stop and issue message.                              
!                                                                       
!     CHECK HEADER FOR FLAG INDICATING COMPATIBILITY WITH ISOTOPES      
!                                                                       
   30    WRITE (CHID10,925) HLINID(10)
         READ (CHID10,930) CHARID,CHARDT,CHTST
         IF (CHTST.NE.CHARI) THEN
            WRITE (IPR,935) CHARID,CHARDT,CHTST
            STOP ' PRLNHD - NO ISOTOPE INFO ON LINFIL '
         ENDIF
!                                                                       
         RETURN

!                                                                       
  900 FORMAT ('0'/'0',20X,'   LINE FILE INFORMATION ') 
  905 FORMAT ('0',10A8,2X,2(1X,A8,1X)) 
  910 FORMAT ('0',/,23X,'COUPLED',4X,'NLTE',3X,'SUM LBLRTM ',/,7X,      &
     &        'MOL',5X,'LINES',4X,'LINES',4X,'LINES',4X,'STRENGTHS',/)  
  915 FORMAT (' ',4X,A6,' = ',I6,3X,I6,3X,I6,2X,1PE12.4,0P) 
  920 FORMAT (/,'0 LOWEST LINE = ',F10.3,5X,'  HIGHEST LINE = ',F10.3,  &
     &        5X,' TOTAL NUMBER OF LINES =',I8)        
  925 FORMAT (A8) 
  930 FORMAT (A5,A2,A1) 
  935 FORMAT (3(/),10X,'LINEFILE PROGRAM: ',A5,3X,'VERSION: ',A2,A1,    &
     &        3(/),3X,52('*'),/,3X,'* THE LINEFILE (TAPE3) IS NOT ',    &
     &        'COMPATIBLE WITH THIS *',/,3X,'* VERSION OF LBLRTM .',    &
     &        '  ISOTOPIC INFORMATION (FROM  *',/,3X,'* HITRAN) ',      &
     &        'MUST BE PRESERVED ON TAPE3.  USE A TAPE3 *',/,3X,        &
     &        '* CREATED WITH THE 91I OR LATER VERSION OF LNFL.   *',   &
     &        /,3X,52('*'))        
  940 FORMAT (' Molecule to be retrieved: ',A6,' not in linefile.',/,   &
     &        ' Molecules in linefile: ')        
  945 FORMAT (24X,A6) 
  950 FORMAT (8a1) 
  960 FORMAT ('0',/,23X,'COUPLED',4X,'NLTE',3X,'NEGATIVE',3X,           &
     &        'RESET',4X,'SUM LBLRTM',/,7X,'MOL',5X,'LINES',4X,         &
     &        'LINES',4X,'LINES',6X,'EPP',6X,'EPP',                     &
     &        6X,'STRENGTHS',/)        
  965 FORMAT (' ',4X,A6,' = ',I6,                                       &
     &        3X,I6,3X,I6,3X,I6,3X,i6,3X,1PE12.4)  
!                             
      END SUBROUTINE PRLNHD

      SUBROUTINE BUFIN_sgl (IFILE,IEOF,IARRAY,IWORDS) 
!                                                                       
!     THIS SUBROUTINE BUFFERS IN (READS) IWORDS INTO  IARRAY STARTING   
!     AT LOCATION IARRAY                                                
!                                                                       
!     IFILE IS THE FILE DESIGNATION                                     
!                                                                       
      implicit integer*4 (i-n) 
      implicit real*4    (a-h,o-z) 
                                                                        
      DATA i_4 / 4 / 
                                                                        
      DIMENSION IARRAY(IWORDS) 

      IEOF = 1 
!                                                                       
!#    BUFFER IN (IFILE,1) (IARRAY(ILO),IARRAY(IHI))                     
!#    IF (UNIT(IFILE).EQ.0.) GO TO 10                                   
!                                                                       
      READ (IFILE,END=10) IARRAY 
      ITEST = MIN(IWORDS,i_4) 
      IF (IARRAY(ITEST).EQ.-99) IEOF = -99 
!                                                                       
      RETURN 
!                                                                       
   10 IEOF = 0 
!                                                                       
      RETURN 
!                                                                       
      END  SUBROUTINE BUFIN_sgl                                         
!_______________________________________________________________________

END MODULE LNFL_MOD
   
