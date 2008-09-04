      SUBROUTINE XINT (V1A,V2A,DVA,A,AFACT,VFT,DVR3,R3,N1R3,N2R3)         B17520
C                                                                         B17530
      IMPLICIT REAL*8           (V)                                     ! B17540
C                                                                         B17550
C     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED                     B17560
C     FROM V1A TO V2A IN INCREMENTS OF DVA USING A MULTIPLICATIVE         B17570
C     FACTOR AFACT, INTO THE R3 ARRAY FROM LOCATION N1R3 TO N2R3 IN       B17580
C     INCREMENTS OF DVR3                                                  B17590
C                                                                         B17600
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           B17610
      DIMENSION A(*),R3(*)                                                B17620
C                                                                         B17630
      RECDVA = 1./DVA                                                     B17640
      ILO = (V1A+DVA-VFT)/DVR3+1.+ONEMI                                   B17650
      ILO = MAX(ILO,N1R3)                                                 B17660
      IHI = (V2A-DVA-VFT)/DVR3+ONEMI                                      B17670
      IHI = MIN(IHI,N2R3)                                                 B17680
C                                                                         B17690
      DO 10 I = ILO, IHI                                                  B17700
         VI = VFT+DVR3* REAL(I-1)                                         B17710
         J = (VI-V1A)*RECDVA+ONEPL                                        B17720
         VJ = V1A+DVA* REAL(J-1)                                          B17730
         P = RECDVA*(VI-VJ)                                               B17740
         C = (3.-2.*P)*P*P                                                B17750
         B = 0.5*P*(1.-P)                                                 B17760
         B1 = B*(1.-P)                                                    B17770
         B2 = B*P                                                         B17780
         CONTI = -A(J-1)*B1+A(J)*(1.-C+B2)+A(J+1)*(C+B1)-A(J+2)*B2        B17790
         R3(I) = R3(I)+CONTI*AFACT                                        B17800
   10 CONTINUE                                                            B17810
C                                                                         B17820
      RETURN                                                              B17830
C                                                                         B17840
      END                                             

      FUNCTION RADFN (VI,XKT)                                             B17860
C                                                                         B17870
      IMPLICIT REAL*8           (V)                                     ! B17880
C                                                                         B17890
C     FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE     B17900
C                                                                         B17910
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B17920
C                                                                         B17930
C               LAST MODIFICATION:    12 AUGUST 1991                      B17940
C                                                                         B17950
C                  IMPLEMENTATION:    R.D. WORSHAM                        B17960
C                                                                         B17970
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         B17980
C                                     R.D. WORSHAM                        B17990
C                                     J.L. MONCET                         B18000
C                                                                         B18010
C                                                                         B18020
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B18030
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B18040
C                                                                         B18050
C----------------------------------------------------------------------   B18060
C                                                                         B18070
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     B18080
C                                     OFFICE OF ENERGY RESEARCH           B18090
C                                     DEPARTMENT OF ENERGY                B18100
C                                                                         B18110
C                                                                         B18120
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B18130
C                                                                         B18140
C                                             FASCOD3                     B18150
C                                                                         B18160
C                                                                         B18170
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B18180
C                                                                         B18190
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           B18200
C                                                                         B18210
C      IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED                         B18220
C                                                                         B18230
      XVI = VI                                                            B18240
C                                                                         B18250
      IF (XKT.GT.0.0) THEN                                                B18260
C                                                                         B18270
         XVIOKT = XVI/XKT                                                 B18280
C                                                                         B18290
         IF (XVIOKT.LE.0.01) THEN                                         B18300
            RADFN = 0.5*XVIOKT*XVI                                        B18310
C                                                                         B18320
         ELSEIF (XVIOKT.LE.10.0) THEN                                     B18330
            EXPVKT = EXP(-XVIOKT)                                         B18340
            RADFN = XVI*(1.-EXPVKT)/(1.+EXPVKT)                           B18350
C                                                                         B18360
         ELSE                                                             B18370
            RADFN = XVI                                                   B18380
         ENDIF                                                            B18390
C                                                                         B18400
      ELSE                                                                B18410
         RADFN = XVI                                                      B18420
      ENDIF
C                                                                         B18440
      RETURN                                                              B18450
C                                                                         B18460
      END         
