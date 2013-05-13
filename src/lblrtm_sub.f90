      SUBROUTINE XINT (V1A,V2A,DVA,A,AFACT,VFT,DVR3,R3,N1R3,N2R3)
!
      IMPLICIT REAL*8           (V)
!
!     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED
!     FROM V1A TO V2A IN INCREMENTS OF DVA USING A MULTIPLICATIVE
!     FACTOR AFACT, INTO THE R3 ARRAY FROM LOCATION N1R3 TO N2R3 IN
!     INCREMENTS OF DVR3
!
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
      DIMENSION A(*),R3(*)
!
      RECDVA = 1./DVA
      ILO = (V1A+DVA-VFT)/DVR3+1.+ONEMI
      ILO = MAX(ILO,N1R3)
      IHI = (V2A-DVA-VFT)/DVR3+ONEMI
      IHI = MIN(IHI,N2R3)
!
      DO 10 I = ILO, IHI
         VI = VFT+DVR3* REAL(I-1)
         J = (VI-V1A)*RECDVA+ONEPL
         VJ = V1A+DVA* REAL(J-1)
         P = RECDVA*(VI-VJ)
         C = (3.-2.*P)*P*P
         B = 0.5*P*(1.-P)
         B1 = B*(1.-P)
         B2 = B*P
         CONTI = -A(J-1)*B1+A(J)*(1.-C+B2)+A(J+1)*(C+B1)-A(J+2)*B2
         R3(I) = R3(I)+CONTI*AFACT
   10 CONTINUE
!
      RETURN
!
      END                                             

      FUNCTION RADFN (VI,XKT)
!
      IMPLICIT REAL*8           (V)
!
!     FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    12 AUGUST 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
!
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
!
!      IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED
!
      XVI = VI
!
      IF (XKT.GT.0.0) THEN
!
         XVIOKT = XVI/XKT
!
         IF (XVIOKT.LE.0.01) THEN
            RADFN = 0.5*XVIOKT*XVI
!
         ELSEIF (XVIOKT.LE.10.0) THEN
            EXPVKT = EXP(-XVIOKT)
            RADFN = XVI*(1.-EXPVKT)/(1.+EXPVKT)
!
         ELSE
            RADFN = XVI
         ENDIF
!
      ELSE
         RADFN = XVI
      ENDIF
!
      RETURN
!
      END

