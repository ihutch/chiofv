c Plasma susceptibility for complex frequency.
c********************************************************************
c
c The plasma susceptibility for a maxwellian distribution of electrons
c can be considered to be 1/(k\lambda_{De})^2 times a function 
c that can be written 
c               \hat{\chi}(z) = (-i\sqrt{\pi}/2) w'(z), 
c with z = (\omega/k)/sqrt{2T_e/m_e}
c and w(z)=\exp(-z^2)\mbox{erfc}(-iz) is the Faddeeva function.
c To calculate w' one can use the recurrence relation
c  w'(z) = -2zw(z) + 2i/\sqrt{\pi} 
c so that
c  \hat{\chi}(z) = 1 +  i z \sqrt{\pi} w(z)
c
c The value of w for complex z is obtained from the netlib-distributed
c Algorithm 680 below.

      complex function chihat(z)
      complex z
      double precision xi,yi,u,v
      logical flag
      xi=real(z)
      yi=imag(z)
      call wofz(xi,yi,u,v,flag)
      if(flag)then
         write(*,*)'chihat overflow'
         chihat=0.
      else
         chihat=1. + complex(0.,sqrt(3.1415926))*z*complex(u,v)
      endif
c      write(*,*)xi,yi,chihat
      end
c**************************************************************************
C      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 16, NO. 1, PP. 47.
      SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
C
C  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
C  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
C  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
C  MEANS SQRT(-1).
C  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
C  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
C  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
C  OF THE FUNCTION.
C  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
C
C
C  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
C     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
C                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
C                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
C                FLOATING-POINT ARITHMETIC
C     RMAXEXP  = LN(RMAX) - LN(2)
C     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
C                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
C  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
C  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
C
C
C  PARAMETER LIST
C     XI     = REAL      PART OF Z
C     YI     = IMAGINARY PART OF Z
C     U      = REAL      PART OF W(Z)
C     V      = IMAGINARY PART OF W(Z)
C     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C              OCCUR OR NOT; TYPE LOGICAL;
C              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C              MEANING :
C              FLAG=.FALSE. : NO ERROR CONDITION
C              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C                             BECOMES INACTIVE
C  XI, YI      ARE THE INPUT-PARAMETERS
C  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
C
C  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
C
C  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
C  PUT TO 0 UPON UNDERFLOW;
C
C  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
C  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE (1990)
C  VOL. 16, NO. 1, PP. 47.
C
*
*
*
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      LOGICAL A, B, FLAG
      PARAMETER (FACTOR   = 1.12837916709551257388D0,
     *           RMAXREAL = 0.5D+154,
     *           RMAXEXP  = 708.503061461606D0,
     *           RMAXGONI = 3.53711887601422D+15)
c IHH silence uninitialization warnings:
      data u2/0./v2/0./h2/0./qlambda/0./
*
      FLAG = .FALSE.
*
      XABS = DABS(XI)
      YABS = DABS(YI)
      X    = XABS/6.3
      Y    = YABS/4.4
*
C
C     THE FOLLOWING IF-STATEMENT PROTECTS
C     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
C
      IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100
*
      QRHO = X**2 + Y**2
*
      XABSQ = XABS**2
      XQUAD = XABSQ - YABS**2
      YQUAD = 2*XABS*YABS
*
      A     = QRHO.LT.0.085264D0
*
      IF (A) THEN
C
C  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
C  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
C  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
        QRHO  = (1-0.85*Y)*DSQRT(QRHO)
        N     = IDNINT(6 + 72*QRHO)
        J     = 2*N+1
        XSUM  = 1.0/J
        YSUM  = 0.0D0
        DO 10 I=N, 1, -1
          J    = J - 2
          XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
          YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
          XSUM = XAUX + 1.0/J
 10     CONTINUE
        U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
        V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
        DAUX =  DEXP(-XQUAD)
        U2   =  DAUX*DCOS(YQUAD)
        V2   = -DAUX*DSIN(YQUAD)
*
        U    = U1*U2 - V1*V2
        V    = U1*V2 + V1*U2
*
      ELSE
C
C  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
C  CONTINUED FRACTION
C  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
C  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
C  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
C  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
C  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
C  TO OBTAIN THE REQUIRED ACCURACY
C  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
C  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
C
*
        IF (QRHO.GT.1.0) THEN
          H    = 0.0D0
          KAPN = 0
          QRHO = DSQRT(QRHO)
          NU   = IDINT(3 + (1442/(26*QRHO+77)))
        ELSE
          QRHO = (1-Y)*DSQRT(1-QRHO)
          H    = 1.88*QRHO
          H2   = 2*H
          KAPN = IDNINT(7  + 34*QRHO)
          NU   = IDNINT(16 + 26*QRHO)
        ENDIF
*
        B = (H.GT.0.0)
*
        IF (B) QLAMBDA = H2**KAPN
*
        RX = 0.0
        RY = 0.0
        SX = 0.0
        SY = 0.0
*
        DO 11 N=NU, 0, -1
          NP1 = N + 1
          TX  = YABS + H + NP1*RX
          TY  = XABS - NP1*RY
          C   = 0.5/(TX**2 + TY**2)
          RX  = C*TX
          RY  = C*TY
          IF ((B).AND.(N.LE.KAPN)) THEN
            TX = QLAMBDA + SX
            SX = RX*TX - RY*SY
            SY = RY*TX + RX*SY
            QLAMBDA = QLAMBDA/H2
          ENDIF
 11     CONTINUE
*
        IF (H.EQ.0.0) THEN
          U = FACTOR*RX
          V = FACTOR*RY
        ELSE
          U = FACTOR*SX
          V = FACTOR*SY
        END IF
*
        IF (YABS.EQ.0.0) U = DEXP(-XABS**2)
*
      END IF
*
*
C
C  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
C
*
      IF (YI.LT.0.0) THEN
*
        IF (A) THEN
          U2    = 2*U2
          V2    = 2*V2
        ELSE
          XQUAD =  -XQUAD
*
C
C         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
C         AGAINST OVERFLOW
C
          IF ((YQUAD.GT.RMAXGONI).OR.
     *        (XQUAD.GT.RMAXEXP)) GOTO 100
*
          W1 =  2*DEXP(XQUAD)
          U2  =  W1*DCOS(YQUAD)
          V2  = -W1*DSIN(YQUAD)
        END IF
*
        U = U2 - U
        V = V2 - V
        IF (XI.GT.0.0) V = -V
      ELSE
        IF (XI.LT.0.0) V = -V
      END IF
*
      RETURN
*
  100 FLAG = .TRUE.
      RETURN
*
      END
c************************************************************************

