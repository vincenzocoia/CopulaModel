      SUBROUTINE PCHEV(N,X,F,D,NVAL,XVAL,FVAL,DVAL,IERR) 
C***BEGIN PROLOGUE  PCHEV
C***DATE WRITTEN   870828   (YYMMDD)
C***REVISION DATE  870828   (YYMMDD)
C***CATEGORY NO.  E3,H1
C***KEYWORDS  CUBIC HERMITE OR SPLINE DIFFERENTIATION,CUBIC HERMITE 
C             EVALUATION,EASY TO USE SPLINE OR CUBIC HERMITE EVALUATOR
C***AUTHOR  KAHANER, D.K., (NBS)
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             ROOM A161, TECHNOLOGY BUILDING
C             GAITHERSBURG, MARYLAND 20899
C             (301) 975-3808
C***PURPOSE  Evaluates the function and first derivative of a piecewise
C            cubic Hermite or spline function at an array of points XVAL,
C            easy to use.
C***DESCRIPTION
C
C          PCHEV:  Piecewise Cubic Hermite or Spline Derivative Evaluator,
C                  Easy to Use.
C
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C                 Prentice Hall 1988
C
C     Evaluates the function and first derivative of the cubic Hermite
C     or spline function defined by  N, X, F, D, at the array of points XVAL.
C
C     This is an easy to use driver for the routines by F.N. Fritsch
C     described in reference (2) below. Those also have other capabilities.
C
C ----------------------------------------------------------------------
C
C  Calling sequence: CALL  PCHEV (N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR)
C
C     INTEGER  N, NVAL, IERR
C     REAL  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The 
C           elements of X must be strictly increasing: 
C             X(I-1) .LT. X(I),  I = 2(1)N. (Error return if not.) 
C
C     F -- (input) real array of function values.  F(I) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(I) is
C           the value corresponding to X(I).
C
C  NVAL -- (input) number of points at which the functions are to be
C           evaluated. ( Error return if NVAL.LT.1 )
C
C  XVAL -- (input) real array of points at which the functions are to
C           be evaluated.
C
C          NOTES:
C           1. The evaluation will be most efficient if the elements
C              of XVAL are increasing relative to X;
C              that is,   XVAL(J) .GE. X(I)
C              implies    XVAL(K) .GE. X(I),  all K.GE.J .
C           2. If any of the XVAL are outside the interval [X(1),X(N)], 
C              values are extrapolated from the nearest extreme cubic,
C              and a warning error is returned.
C
C  FVAL -- (output) real array of values of the cubic Hermite function
C           defined by  N, X, F, D  at the points  XVAL.
C
C  DVAL -- (output) real array of values of the first derivative of
C           the same function at the points  XVAL.
C
C  IERR -- (output) error flag. 
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  means that extrapolation was performed at
C                 IERR points.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 . 
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if NVAL.LT.1 .
C           (Output arrays have not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C              IERR = -5  if an error has occurred in the lower-level
C                         routine CHFDV.  NB: this should never happen.
C                         Notify the author **IMMEDIATELY** if it does.
C
C ----------------------------------------------------------------------
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION 
C                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
C                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194,
C                 AUGUST 1982.
C***ROUTINES CALLED  PCHFD
C***END PROLOGUE  PCHEV
      implicit none
      INTEGER  N, NVAL, IERR
      double precision  X(N),F(N),D(N),XVAL(NVAL),FVAL(NVAL),DVAL(NVAL)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER INCFD
      LOGICAL SKIP
      DATA SKIP /.TRUE./
      DATA INCFD /1/

C
C
C***FIRST EXECUTABLE STATEMENT  PCHEV
C
      CALL PCHFD(N,X,F,D,INCFD,SKIP,NVAL,XVAL,FVAL,DVAL,IERR)   
C
C
C5000 CONTINUE
      RETURN
C
C------------- LAST LINE OF PCHEV FOLLOWS ------------------------------
      END 
      SUBROUTINE PCHFD(N,X,F,D,INCFD,SKIP,NE,XE,FE,DE,IERR)
C***BEGIN PROLOGUE  PCHFD
C     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE
C     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C***END PROLOGUE  PCHFD
C
C  DECLARE ARGUMENTS.
C
      implicit none
      INTEGER  N, INCFD, NE, IERR
      double precision  X(N),F(INCFD,N),D(INCFD,N),XE(NE),FE(NE),DE(NE)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IERC, IR, J, JFIRST, NEXT(2), NJ
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHFD
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      IF ( NE.LT.1 )  GO TO 5004
      IERR = 0
      SKIP = .TRUE.
C
C  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
C                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
      JFIRST = 1
      IR = 2
   10 CONTINUE
C
C     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
C
         IF (JFIRST .GT. NE)  GO TO 5000
C
C     LOCATE ALL POINTS IN INTERVAL.
C
         DO 20  J = JFIRST, NE
            IF (XE(J) .GE. X(IR))  GO TO 30
   20    CONTINUE
         J = NE + 1
         GO TO 40
C
C     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
C
   30    CONTINUE
         IF (IR .EQ. N)  J = NE + 1
C
   40    CONTINUE
         NJ = J - JFIRST
C
C     SKIP EVALUATION IF NO POINTS IN INTERVAL.
C
         IF (NJ .EQ. 0)  GO TO 50
C
C     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
C
C       ----------------------------------------------------------------
        CALL CHFDV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR),
     *              NJ, XE(JFIRST), FE(JFIRST), DE(JFIRST), NEXT, IERC)
C       ----------------------------------------------------------------
         IF (IERC .LT. 0)  GO TO 5005
C
         IF (NEXT(2) .EQ. 0)  GO TO 42
C        IF (NEXT(2) .GT. 0)  THEN
C           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
C           RIGHT OF X(IR).
C
            IF (IR .LT. N)  GO TO 41
C           IF (IR .EQ. N)  THEN
C              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(2)
               GO TO 42
   41       CONTINUE
C           ELSE
C              WE SHOULD NEVER HAVE GOTTEN HERE.
               GO TO 5005
C           ENDIF
C        ENDIF
   42    CONTINUE
C
         IF (NEXT(1) .EQ. 0)  GO TO 49
C        IF (NEXT(1) .GT. 0)  THEN
C           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
C           LEFT OF X(IR-1).
C
            IF (IR .GT. 2)  GO TO 43
C           IF (IR .EQ. 2)  THEN
C              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(1)
               GO TO 49
   43       CONTINUE
C           ELSE
C              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
C              EVALUATION INTERVAL.
C
C              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
               DO 44  I = JFIRST, J-1
                  IF (XE(I) .LT. X(IR-1))  GO TO 45
   44          CONTINUE
C              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
C                     IN CHFDV.
               GO TO 5005
C
   45          CONTINUE
C              RESET J.  (THIS WILL BE THE NEW JFIRST.)
               J = I
C
C              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
               DO 46  I = 1, IR-1
                  IF (XE(J) .LT. X(I)) GO TO 47
   46          CONTINUE
C              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
C
   47          CONTINUE
C              AT THIS POINT, EITHER  XE(J) .LT. X(1)
C                 OR      X(I-1) .LE. XE(J) .LT. X(I) .
C              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
C              CYCLING.
               IR = MAX0(1, I-1)
C           ENDIF
C        ENDIF
   49    CONTINUE
C
         JFIRST = J
C
C     END OF IR-LOOP.
C
   50 CONTINUE
      IR = IR + 1
      IF (IR .LE. N)  GO TO 10
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
c     CALL XERROR ('PCHFD -- NUMBER OF DATA POINTS LESS THAN TWO'
c    *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
c     CALL XERROR ('PCHFD -- INCREMENT LESS THAN ONE'
c    *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
c     CALL XERROR ('PCHFD -- X-ARRAY NOT STRICTLY INCREASING'
c    *           , 40, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -4
c     CALL XERROR ('PCHFD -- NUMBER OF EVALUATION POINTS LESS THAN ONE'
c    *           , 50, IERR, 1)
      RETURN
C
 5005 CONTINUE
C     ERROR RETURN FROM CHFDV.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -5
c     CALL XERROR ('PCHFD -- ERROR RETURN FROM CHFDV -- FATAL'
c    *           , 41, IERR, 2)
      RETURN
C------------- LAST LINE OF PCHFD FOLLOWS ------------------------------
      END
      SUBROUTINE CHFDV(X1,X2,F1,F2,D1,D2,NE,XE,FE,DE,NEXT,IERR)
C***BEGIN PROLOGUE  CHFDV
C     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE
C     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C***END PROLOGUE  CHFDV
C
C  DECLARE ARGUMENTS.
C
      implicit none
      INTEGER  NE, NEXT(2), IERR
      double precision  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I
      double precision C2,C2T2,C3,C3T3,DEL1,DEL2,DELTA,H,X,XMI,XMA,ZERO
      DATA  ZERO /0.d0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  CHFDV
      IF (NE .LT. 1)  GO TO 5001
      H = X2 - X1
      IF (H .EQ. ZERO)  GO TO 5002
C
C  INITIALIZE.
C
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = dMIN1(ZERO, H)
      XMA = dMAX1(ZERO, H)
C
C  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
C
      DELTA = (F2 - F1)/H
      DEL1 = (D1 - DELTA)/H
      DEL2 = (D2 - DELTA)/H
C                                           (DELTA IS NO LONGER NEEDED.)
      C2 = -(DEL1+DEL1 + DEL2)
      C2T2 = C2 + C2
      C3 = (DEL1 + DEL2)/H
C                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
      C3T3 = C3+C3+C3
C
C  EVALUATION LOOP.
C
      DO 500  I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
         DE(I) = D1 + X*(C2T2 + X*C3T3)
C          COUNT EXTRAPOLATION POINTS.
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
C        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 CONTINUE
C
C  NORMAL RETURN.
C
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -1
c     CALL XERROR ('CHFDV -- NUMBER OF EVALUATION POINTS LESS THAN ONE'
c    *           , 50, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     X1.EQ.X2 RETURN.
      IERR = -2
c     CALL XERROR ('CHFDV -- INTERVAL ENDPOINTS EQUAL'
c    *           , 33, IERR, 1)
      RETURN
C------------- LAST LINE OF CHFDV FOLLOWS ------------------------------
      END
