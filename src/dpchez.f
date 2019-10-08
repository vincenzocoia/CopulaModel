      SUBROUTINE PCHEZ(N,X,F,D,SPLINE,WK,LWK,IERR)
C***BEGIN PROLOGUE  PCHEZ
C***DATE WRITTEN   870821   (YYMMDD)
C***REVISION DATE  870908   (YYMMDD)
C***CATEGORY NO.  E1B
C***KEYWORDS  CUBIC HERMITE MONOTONE INTERPOLATION, SPLINE
C             INTERPOLATION, EASY TO USE PIECEWISE CUBIC INTERPOLATION
C***AUTHOR  KAHANER, D.K., (NBS)
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             GAITHERSBURG, MARYLAND 20899
C             (301) 975-3808
C***PURPOSE  Easy to use spline or cubic Hermite interpolation.
C***DESCRIPTION
C
C          PCHEZ:  Piecewise Cubic Interpolation, Easy to Use.
C
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C
C     Sets derivatives for spline (two continuous derivatives) or
C     Hermite cubic (one continuous derivative) interpolation.
C     Spline interpolation is smoother, but may not "look" right if the
C     data contains both "steep" and "flat" sections.  Hermite cubics
C     can produce a "visually pleasing" and monotone interpolant to
C     monotone data. This is an easy to use driver for the routines 
C     by F. N. Fritsch in reference (4) below. Various boundary 
C     conditions are set to default values by PCHEZ. Many other choices 
C     are available in the subroutines PCHIC, PCHIM and PCHSP.
C
C     Use PCHEV to evaluate the resulting function and its derivative.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:   CALL  PCHEZ (N, X, F, D, SPLINE, WK, LWK, IERR) 
C
C     INTEGER  N, IERR,  LWK
C     REAL  X(N), F(N), D(N), WK(*)
C     LOGICAL SPLINE
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C           If N=2, simply does linear interpolation.
C
C     X -- (input) real array of independent variable values.  The 
C           elements of X must be strictly increasing: 
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.) 
C
C     F -- (input) real array of dependent variable values to be inter- 
C           polated.  F(I) is value corresponding to X(I). 
C
C     D -- (output) real array of derivative values at the data points. 
C
C     SPLINE -- (input) logical variable to specify if the interpolant 
C           is to be a spline with two continuous derivaties 
C           (set SPLINE=.TRUE.) or a Hermite cubic interpolant with one
C           continuous derivative (set SPLINE=.FALSE.).
C        Note: If SPLINE=.TRUE. the interpolating spline satisfies the
C           default "not-a-knot" boundary condition, with a continuous
C           third derivative at X(2) and X(N-1). See reference (3).
C              If SPLINE=.FALSE. the interpolating Hermite cubic will be
C           monotone if the input data is monotone. Boundary conditions are
C           computed from the derivative of a local quadratic unless this
C           alters monotonicity.
C
C     WK -- (scratch) real work array, which must be declared by the calling
C           program to be at least 2*N if SPLINE is .TRUE. and not used
C           otherwise.
C
C     LWK -- (input) length of work array WK. (Error return if 
C           LWK.LT.2*N and SPLINE is .TRUE., not checked otherwise.) 
C
C     IERR -- (output) error flag. 
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  (can only occur when SPLINE=.FALSE.) means that 
C                 IERR switches in the direction of monotonicity were detected.
C                 When SPLINE=.FALSE.,  PCHEZ guarantees that if the input
C                 data is monotone, the interpolant will be too. This warning
C                 is to alert you to the fact that the input data was not 
C                 monotone.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 . 
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -7  if LWK is less than 2*N and SPLINE is .TRUE.
C             (The D-array has not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C ----------------------------------------------------------------------
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING
C                 LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,' LLNL
C                 PREPRINT UCRL-87559 (APRIL 1982).
C               3. CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
C                 VERLAG (NEW YORK, 1978).  (ESP. CHAPTER IV, PP.49-62.)
C               4. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION 
C                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
C                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194,
C                 AUGUST 1982.
C***ROUTINES CALLED  PCHIM,PCHSP
C***END PROLOGUE  PCHEZ
      implicit none
      INTEGER  N, LWK, IERR 
      double precision  X(N), F(N), D(N), WK(LWK)
      LOGICAL SPLINE
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER IC(2), INCFD
      double precision  VC(2)
      DATA IC(1) /0/
      DATA IC(2) /0/
      DATA INCFD /1/
C
C
C***FIRST EXECUTABLE STATEMENT  PCHEZ
C
      IF ( SPLINE ) THEN
        CALL  PCHSP (IC, VC, N, X, F, D, INCFD, WK, LWK, IERR)
      ELSE
        CALL  PCHIM (N, X, F, D, INCFD, IERR)
      ENDIF
C
C  ERROR CONDITIONS ALREADY CHECKED IN PCHSP OR PCHIM

      RETURN
C------------- LAST LINE OF PCHEZ FOLLOWS ------------------------------
      END 
      SUBROUTINE PCHIM(N,X,F,D,INCFD,IERR)
C***BEGIN PROLOGUE  PCHIM
C     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE
C     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C***END PROLOGUE  PCHIM
C
C  DECLARE ARGUMENTS.
C
      implicit none
      INTEGER  N, INCFD, IERR
      double precision  X(N), F(INCFD,N), D(INCFD,N)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, NLESS1
      double precision  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE,
     *      H1, H2, HSUM, HSUMT3, THREE, W1, W2, ZERO
      double precision  PCHST
      DATA  ZERO /0.d0/,  THREE /3.d0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHIM
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
      IERR = 0
      NLESS1 = N - 1
      H1 = X(2) - X(1)
      DEL1 = (F(1,2) - F(1,1))/H1
      DSAVE = DEL1
C
C  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
C
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
C
C  NORMAL CASE  (N .GE. 3).
C
   10 CONTINUE
      H2 = X(3) - X(2)
      DEL2 = (F(1,3) - F(1,2))/H2
C
C  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
C     SHAPE-PRESERVING.
C
      HSUM = H1 + H2
      W1 = (H1 + HSUM)/HSUM
      W2 = -H1/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( PCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN
C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL1
         IF (dABS(D(1,1)) .GT. dABS(DMAX))  D(1,1) = DMAX
      ENDIF
C
C  LOOP THROUGH INTERIOR POINTS.
C
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
C
         H1 = H2
         H2 = X(I+1) - X(I)
         HSUM = H1 + H2
         DEL1 = DEL2
         DEL2 = (F(1,I+1) - F(1,I))/H2
   40    CONTINUE
C
C        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
C
         D(1,I) = ZERO
         IF ( PCHST(DEL1,DEL2) )  42, 41, 45
C
C        COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
C
   41    CONTINUE
         IF (DEL2 .EQ. ZERO)  GO TO 50
         IF ( PCHST(DSAVE,DEL2) .LT. ZERO)  IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
C
   42    CONTINUE
         IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
C
C        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
C
   45    CONTINUE
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H1)/HSUMT3
         W2 = (HSUM + H2)/HSUMT3
         DMAX = dMAX1( dABS(DEL1), dABS(DEL2) )
         DMIN = dMIN1( dABS(DEL1), dABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
C
   50 CONTINUE
C
C  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
C     SHAPE-PRESERVING.
C
      W1 = -H2/HSUM
      W2 = (H2 + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( PCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN
C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL2
         IF (dABS(D(1,N)) .GT. dABS(DMAX))  D(1,N) = DMAX
      ENDIF
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
c     CALL XERROR ('PCHIM -- NUMBER OF DATA POINTS LESS THAN TWO'
c    *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
c     CALL XERROR ('PCHIM -- INCREMENT LESS THAN ONE'
c    *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
c     CALL XERROR ('PCHIM -- X-ARRAY NOT STRICTLY INCREASING'
c    *           , 40, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHIM FOLLOWS ------------------------------
      END
      double precision FUNCTION PCHST(ARG1,ARG2)
C***BEGIN PROLOGUE  PCHST
C***REFER TO  PCHCE,PCHCI,PCHCS,PCHIM
C***END PROLOGUE  PCHST
      implicit none
      double precision  ARG1, ARG2
C
C  DECLARE LOCAL VARIABLES.
C
      double precision  ONE, ZERO
      DATA  ZERO /0.d0/,  ONE /1.d0/
C
C  PERFORM THE TEST.
C
C***FIRST EXECUTABLE STATEMENT  PCHST
      PCHST = dSIGN(ONE,ARG1) * dSIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  PCHST = ZERO
C
      RETURN
C------------- LAST LINE OF PCHST FOLLOWS ------------------------------
      END
      SUBROUTINE PCHSP(IC,VC,N,X,F,D,INCFD,WK,NWK,IERR)
C***BEGIN PROLOGUE  PCHSP
C     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE
C     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C***END PROLOGUE  PCHSP
C
C  DECLARE ARGUMENTS.
C
      implicit none
      INTEGER  IC(2), N, INCFD, NWK, IERR
      double precision  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(2,N)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  IBEG, IEND, INDEX, J, NM1
      double precision G,HALF,ONE,STEMP(3),THREE,TWO,XTEMP(4),ZERO
      double precision PCHDF
C
      DATA  ZERO/0.d0/, HALF/0.5d0/, ONE/1.d0/, TWO/2.d0/, THREE/3.d0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHSP
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  J = 2, N
         IF ( X(J).LE.X(J-1) )  GO TO 5003
    1 CONTINUE
C
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
      IF ( (IBEG.LT.0).OR.(IBEG.GT.4) )  IERR = IERR - 1
      IF ( (IEND.LT.0).OR.(IEND.GT.4) )  IERR = IERR - 2
      IF ( IERR.LT.0 )  GO TO 5004
C
C  FUNCTION DEFINITION IS OK -- GO ON.
C
      IF ( NWK .LT. 2*N )  GO TO 5007
C
C  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
C  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
      DO 5  J=2,N
         WK(1,J) = X(J) - X(J-1)
         WK(2,J) = (F(1,J) - F(1,J-1))/WK(1,J)
    5 CONTINUE
C
C  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
C
      IF ( IBEG.GT.N )  IBEG = 0
      IF ( IEND.GT.N )  IEND = 0
C
C  SET UP FOR BOUNDARY CONDITIONS.
C
      IF ( (IBEG.EQ.1).OR.(IBEG.EQ.2) )  THEN
         D(1,1) = VC(1)
      ELSE IF (IBEG .GT. 2)  THEN
C        PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
         DO 10  J = 1, IBEG
            INDEX = IBEG-J+1
C           INDEX RUNS FROM IBEG DOWN TO 1.
            XTEMP(J) = X(INDEX)
            IF (J .LT. IBEG)  STEMP(J) = WK(2,INDEX)
   10    CONTINUE
C                 --------------------------------
         D(1,1) = PCHDF (IBEG, XTEMP, STEMP, IERR)
C                 --------------------------------
         IF (IERR .NE. 0)  GO TO 5009
         IBEG = 1
      ENDIF
C
      IF ( (IEND.EQ.1).OR.(IEND.EQ.2) )  THEN
         D(1,N) = VC(2)
      ELSE IF (IEND .GT. 2)  THEN
C        PICK UP LAST IEND POINTS.
         DO 15  J = 1, IEND
            INDEX = N-IEND+J
C           INDEX RUNS FROM N+1-IEND UP TO N.
            XTEMP(J) = X(INDEX)
            IF (J .LT. IEND)  STEMP(J) = WK(2,INDEX+1)
   15    CONTINUE
C                 --------------------------------
         D(1,N) = PCHDF (IEND, XTEMP, STEMP, IERR)
C                 --------------------------------
         IF (IERR .NE. 0)  GO TO 5009
         IEND = 1
      ENDIF
C
C --------------------( BEGIN CODING FROM CUBSPL )--------------------
C
C  **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
C  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
C  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
C     WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
C
C  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
C             WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
C
      IF (IBEG .EQ. 0)  THEN
         IF (N .EQ. 2)  THEN
C           NO CONDITION AT LEFT END AND N = 2.
            WK(2,1) = ONE
            WK(1,1) = ONE
            D(1,1) = TWO*WK(2,2)
         ELSE
C           NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
            WK(2,1) = WK(1,3)
            WK(1,1) = WK(1,2) + WK(1,3)
            D(1,1) =((WK(1,2) + TWO*WK(1,1))*WK(2,2)*WK(1,3)
     *                        + WK(1,2)**2*WK(2,3)) / WK(1,1)
         ENDIF
      ELSE IF (IBEG .EQ. 1)  THEN
C        SLOPE PRESCRIBED AT LEFT END.
         WK(2,1) = ONE
         WK(1,1) = ZERO
      ELSE
C        SECOND DERIVATIVE PRESCRIBED AT LEFT END.
         WK(2,1) = TWO
         WK(1,1) = ONE
         D(1,1) = THREE*WK(2,2) - HALF*WK(1,2)*D(1,1)
      ENDIF
C
C  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
C  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
C  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
C
      NM1 = N-1
      IF (NM1 .GT. 1)  THEN
         DO 20 J=2,NM1
            IF (WK(2,J-1) .EQ. ZERO)  GO TO 5008
            G = -WK(1,J+1)/WK(2,J-1)
            D(1,J) = G*D(1,J-1)
     *                  + THREE*(WK(1,J)*WK(2,J+1) + WK(1,J+1)*WK(2,J))
            WK(2,J) = G*WK(1,J-1) + TWO*(WK(1,J) + WK(1,J+1))
   20    CONTINUE
      ENDIF
C
C  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
C           (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
C
C     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
C     SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
C     AT THIS POINT.
      IF (IEND .EQ. 1)  GO TO 30
C
      IF (IEND .EQ. 0)  THEN
         IF (N.EQ.2 .AND. IBEG.EQ.0)  THEN
C           NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
            D(1,2) = WK(2,2)
            GO TO 30
         ELSE IF ((N.EQ.2) .OR. (N.EQ.3 .AND. IBEG.EQ.0))  THEN
C           EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
C           NOT-A-KNOT AT LEFT END POINT).
            D(1,N) = TWO*WK(2,N)
            WK(2,N) = ONE
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -ONE/WK(2,N-1)
         ELSE
C           NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
C           KNOT AT LEFT END POINT.
            G = WK(1,N-1) + WK(1,N)
C           DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
            D(1,N) = ((WK(1,N)+TWO*G)*WK(2,N)*WK(1,N-1)
     *                  + WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -G/WK(2,N-1)
            WK(2,N) = WK(1,N-1)
         ENDIF
      ELSE
C        SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
         D(1,N) = THREE*WK(2,N) + HALF*WK(1,N)*D(1,N)
         WK(2,N) = TWO
         IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
         G = -ONE/WK(2,N-1)
      ENDIF
C
C  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
C
      WK(2,N) = G*WK(1,N-1) + WK(2,N)
      IF (WK(2,N) .EQ. ZERO)   GO TO 5008
      D(1,N) = (G*D(1,N-1) + D(1,N))/WK(2,N)
C
C  CARRY OUT BACK SUBSTITUTION
C
   30 CONTINUE
      DO 40 J=NM1,1,-1
         IF (WK(2,J) .EQ. ZERO)  GO TO 5008
         D(1,J) = (D(1,J) - WK(1,J)*D(1,J+1))/WK(2,J)
   40 CONTINUE
C --------------------(  END  CODING FROM CUBSPL )--------------------
C
C  NORMAL RETURN.
C
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
c     CALL XERROR ('PCHSP -- NUMBER OF DATA POINTS LESS THAN TWO'
c    *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
c     CALL XERROR ('PCHSP -- INCREMENT LESS THAN ONE'
c    *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
c     CALL XERROR ('PCHSP -- X-ARRAY NOT STRICTLY INCREASING'
c    *           , 40, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     IC OUT OF RANGE RETURN.
      IERR = IERR - 3
c     CALL XERROR ('PCHSP -- IC OUT OF RANGE'
c    *           , 24, IERR, 1)
      RETURN
C
 5007 CONTINUE
C     NWK TOO SMALL RETURN.
      IERR = -7
c     CALL XERROR ('PCHSP -- WORK ARRAY TOO SMALL'
c    *           , 29, IERR, 1)
      RETURN
C
 5008 CONTINUE
C     SINGULAR SYSTEM.
C   *** THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES   ***
C   *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
      IERR = -8
c     CALL XERROR ('PCHSP -- SINGULAR LINEAR SYSTEM'
c    *           , 31, IERR, 1)
      RETURN
C
 5009 CONTINUE
C     ERROR RETURN FROM PCHDF.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
c     CALL XERROR ('PCHSP -- ERROR RETURN FROM PCHDF'
c    *           , 32, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHSP FOLLOWS ------------------------------
      END
      double precision FUNCTION PCHDF(K,X,S,IERR)
C***BEGIN PROLOGUE  PCHDF
C***REFER TO  PCHCE,PCHSP
C***END PROLOGUE  PCHDF
      implicit none
      INTEGER  K, IERR
      double precision  X(K), S(K)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, J
      double precision  VALUE, ZERO
      DATA  ZERO /0.d0/
C
C  CHECK FOR LEGAL VALUE OF K.
C
C***FIRST EXECUTABLE STATEMENT  PCHDF
      IF (K .LT. 3)  GO TO 5001
C
C  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
C
      DO 10  J = 2, K-1
         DO 9  I = 1, K-J
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
C
C  EVALUATE DERIVATIVE AT X(K).
C
      VALUE = S(1)
      DO 20  I = 2, K-1
         VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
C
C  NORMAL RETURN.
C
      IERR = 0
      PCHDF = VALUE
      RETURN
C
C  ERROR RETURN.
C
 5001 CONTINUE
C     K.LT.3 RETURN.
      IERR = -1
c     CALL XERROR ('PCHDF -- K LESS THAN THREE'
c    *           , 26, IERR, 1)
      PCHDF = ZERO
      RETURN
C------------- LAST LINE OF PCHDF FOLLOWS ------------------------------
      END
