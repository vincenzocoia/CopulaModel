c  this fortran77 code is extracted from the mvtnorm R package source

      double precision function mvbvt( nu, lower, upper, infin, correl )
*
*     A function for computing bivariate normal and t probabilities.
*
*  Parameters
*
*     NU     INTEGER degrees of freedom parameter; NU < 1 gives normal case.
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, correlation coefficient.
*
      DOUBLE PRECISION LOWER(*), UPPER(*), CORREL, MVBVN, MVBVTL
      INTEGER NU, INFIN(*)
      IF ( NU .LT. 1 ) THEN
            MVBVT =  MVBVN ( LOWER, UPPER, INFIN, CORREL )
      ELSE
         IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 2 ) THEN
            MVBVT =  MVBVTL ( NU, UPPER(1), UPPER(2), CORREL )
     +           - MVBVTL ( NU, UPPER(1), LOWER(2), CORREL )
     +           - MVBVTL ( NU, LOWER(1), UPPER(2), CORREL )
     +           + MVBVTL ( NU, LOWER(1), LOWER(2), CORREL )
         ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 1 ) THEN
            MVBVT =  MVBVTL ( NU, -LOWER(1), -LOWER(2), CORREL )
     +           - MVBVTL ( NU, -UPPER(1), -LOWER(2), CORREL )
         ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 2 ) THEN
            MVBVT =  MVBVTL ( NU, -LOWER(1), -LOWER(2), CORREL )
     +           - MVBVTL ( NU, -LOWER(1), -UPPER(2), CORREL )
         ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 0 ) THEN
            MVBVT =  MVBVTL ( NU, UPPER(1), UPPER(2), CORREL )
     +           - MVBVTL ( NU, LOWER(1), UPPER(2), CORREL )
         ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 2 ) THEN
            MVBVT =  MVBVTL ( NU, UPPER(1), UPPER(2), CORREL )
     +           - MVBVTL ( NU, UPPER(1), LOWER(2), CORREL )
         ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 0 ) THEN
            MVBVT =  MVBVTL ( NU, -LOWER(1), UPPER(2), -CORREL )
         ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 1 ) THEN
            MVBVT =  MVBVTL ( NU, UPPER(1), -LOWER(2), -CORREL )
         ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 1 ) THEN
            MVBVT =  MVBVTL ( NU, -LOWER(1), -LOWER(2), CORREL )
         ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 0 ) THEN
            MVBVT =  MVBVTL ( NU, UPPER(1), UPPER(2), CORREL )
         ELSE
            MVBVT = 1
         END IF
      END IF
      END
*
      DOUBLE PRECISION FUNCTION MVBVTC( NU, L, U, INFIN, RHO )      
*
*     A function for computing complementary bivariate normal and t 
*       probabilities.
*
*  Parameters
*
*     NU     INTEGER degrees of freedom parameter.
*     L      REAL, array of lower integration limits.
*     U      REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(1) INFIN(2),        then MVBVTC computes
*                 0         0              P( X>U(1), Y>U(2) )
*                 1         0              P( X<L(1), Y>U(2) )
*                 0         1              P( X>U(1), Y<L(2) )
*                 1         1              P( X<L(1), Y<L(2) )
*                 2         0      P( X>U(1), Y>U(2) ) + P( X<L(1), Y>U(2) )
*                 2         1      P( X>U(1), Y<L(2) ) + P( X<L(1), Y<L(2) )
*                 0         2      P( X>U(1), Y>U(2) ) + P( X>U(1), Y<L(2) )
*                 1         2      P( X<L(1), Y>U(2) ) + P( X<L(1), Y<L(2) )
*                 2         2      P( X>U(1), Y<L(2) ) + P( X<L(1), Y<L(2) )
*                               +  P( X>U(1), Y>U(2) ) + P( X<L(1), Y>U(2) )
*
*     RHO    REAL, correlation coefficient.
*
      DOUBLE PRECISION L(*), U(*), LW(2), UP(2), B, RHO, MVBVT
      INTEGER I, NU, INFIN(*), INF(2)
*
      DO I = 1, 2
         IF ( MOD( INFIN(I), 2 ) .EQ. 0 ) THEN
            INF(I) = 1
            LW(I) = U(I) 
         ELSE
            INF(I) = 0
            UP(I) = L(I) 
         END IF
      END DO
      B = MVBVT( NU, LW, UP, INF, RHO )
      DO I = 1, 2
         IF ( INFIN(I) .EQ. 2 ) THEN
            INF(I) = 0
            UP(I) = L(I) 
            B = B + MVBVT( NU, LW, UP, INF, RHO )
         END IF
      END DO
      IF ( INFIN(1) .EQ. 2 .AND. INFIN(2) .EQ. 2 ) THEN
         INF(1) = 1
         LW(1) = U(1) 
         B = B + MVBVT( NU, LW, UP, INF, RHO )
      END IF
      MVBVTC = B
      END
*
      double precision function mvbvtl( nu, dh, dk, r )
*
*     a function for computing bivariate t probabilities.
*
*       Alan Genz
*       Department of Mathematics
*       Washington State University
*       Pullman, Wa 99164-3113
*       Email : alangenz@wsu.edu
*
*    this function is based on the method described by 
*        Dunnett, C.W. and M. Sobel, (1954),
*        A bivariate generalization of Student's t-distribution
*        with tables for certain special cases,
*        Biometrika 41, pp. 153-169.
*
* mvbvtl - calculate the probability that x < dh and y < dk. 
*
* parameters
*
*   nu number of degrees of freedom
*   dh 1st lower integration limit
*   dk 2nd lower integration limit
*   r   correlation coefficient
*
      integer nu, j, hs, ks
      double precision dh, dk, r
      double precision tpi, pi, ors, hrk, krh, bvt, snu 
      double precision gmph, gmpk, xnkh, xnhk, qhrk, hkn, hpk, hkrn
      double precision btnckh, btnchk, btpdkh, btpdhk, one
      parameter ( pi = 3.14159265358979323844d0, tpi = 2*pi, one = 1 )
      snu = sqrt( dble(nu) )
      ors = 1 - r*r  
      hrk = dh - r*dk  
      krh = dk - r*dh  
      if ( abs(hrk) + ors .gt. 0 ) then
         xnhk = hrk**2/( hrk**2 + ors*( nu + dk**2 ) ) 
         xnkh = krh**2/( krh**2 + ors*( nu + dh**2 ) ) 
      else
         xnhk = 0
         xnkh = 0  
      end if
      hs = sign( one, dh - r*dk )  
      ks = sign( one, dk - r*dh ) 
      if ( mod( nu, 2 ) .eq. 0 ) then
         bvt = atan2( sqrt(ors), -r )/tpi 
         gmph = dh/sqrt( 16*( nu + dh**2 ) )  
         gmpk = dk/sqrt( 16*( nu + dk**2 ) )  
         btnckh = 2*atan2( sqrt( xnkh ), sqrt( 1 - xnkh ) )/pi  
         btpdkh = 2*sqrt( xnkh*( 1 - xnkh ) )/pi 
         btnchk = 2*atan2( sqrt( xnhk ), sqrt( 1 - xnhk ) )/pi  
         btpdhk = 2*sqrt( xnhk*( 1 - xnhk ) )/pi 
         do j = 1, nu/2
            bvt = bvt + gmph*( 1 + ks*btnckh ) 
            bvt = bvt + gmpk*( 1 + hs*btnchk ) 
            btnckh = btnckh + btpdkh  
            btpdkh = 2*j*btpdkh*( 1 - xnkh )/( 2*j + 1 )  
            btnchk = btnchk + btpdhk  
            btpdhk = 2*j*btpdhk*( 1 - xnhk )/( 2*j + 1 )  
            gmph = gmph*( 2*j - 1 )/( 2*j*( 1 + dh**2/nu ) ) 
            gmpk = gmpk*( 2*j - 1 )/( 2*j*( 1 + dk**2/nu ) ) 
         end do
      else
         qhrk = sqrt( dh**2 + dk**2 - 2*r*dh*dk + nu*ors )  
         hkrn = dh*dk + r*nu  
         hkn = dh*dk - nu  
         hpk = dh + dk 
         bvt = atan2(-snu*(hkn*qhrk+hpk*hkrn),hkn*hkrn-nu*hpk*qhrk)/tpi
         if ( bvt .lt. -1d-15 ) bvt = bvt + 1
         gmph = dh/( tpi*snu*( 1 + dh**2/nu ) )  
         gmpk = dk/( tpi*snu*( 1 + dk**2/nu ) )  
         btnckh = sqrt( xnkh )  
         btpdkh = btnckh 
         btnchk = sqrt( xnhk )  
         btpdhk = btnchk  
         do j = 1, ( nu - 1 )/2
            bvt = bvt + gmph*( 1 + ks*btnckh ) 
            bvt = bvt + gmpk*( 1 + hs*btnchk ) 
            btpdkh = ( 2*j - 1 )*btpdkh*( 1 - xnkh )/( 2*j )  
            btnckh = btnckh + btpdkh  
            btpdhk = ( 2*j - 1 )*btpdhk*( 1 - xnhk )/( 2*j )  
            btnchk = btnchk + btpdhk  
            gmph = 2*j*gmph/( ( 2*j + 1 )*( 1 + dh**2/nu ) ) 
            gmpk = 2*j*gmpk/( ( 2*j + 1 )*( 1 + dk**2/nu ) ) 
         end do
      end if
      mvbvtl = bvt 
*
*     end mvbvtl
*
      end
      DOUBLE PRECISION FUNCTION MVBVN( LOWER, UPPER, INFIN, CORREL )
*
*     A function for computing bivariate normal probabilities.
*
*  Parameters
*
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, correlation coefficient.
*
      DOUBLE PRECISION LOWER(*), UPPER(*), CORREL, MVBVU
      INTEGER INFIN(*)
      IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 2 ) THEN
         MVBVN =  MVBVU ( LOWER(1), LOWER(2), CORREL )
     +           - MVBVU ( UPPER(1), LOWER(2), CORREL )
     +           - MVBVU ( LOWER(1), UPPER(2), CORREL )
     +           + MVBVU ( UPPER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 1 ) THEN
         MVBVN =  MVBVU ( LOWER(1), LOWER(2), CORREL )
     +           - MVBVU ( UPPER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 2 ) THEN
         MVBVN =  MVBVU ( LOWER(1), LOWER(2), CORREL )
     +           - MVBVU ( LOWER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 0 ) THEN
         MVBVN =  MVBVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - MVBVU ( -LOWER(1), -UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 2 ) THEN
         MVBVN =  MVBVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - MVBVU ( -UPPER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 0 ) THEN
         MVBVN =  MVBVU ( LOWER(1), -UPPER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 1 ) THEN
         MVBVN =  MVBVU ( -UPPER(1), LOWER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 1 ) THEN
         MVBVN =  MVBVU ( LOWER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 0 ) THEN
         MVBVN =  MVBVU ( -UPPER(1), -UPPER(2), CORREL )
      ELSE
         MVBVN = 1
      END IF
      END 
      DOUBLE PRECISION FUNCTION MVBVU( SH, SK, R )
*
*     A function for computing bivariate normal probabilities;
*       developed using 
*         Drezner, Z. and Wesolowsky, G. O. (1989),
*         On the Computation of the Bivariate Normal Integral,
*         J. Stat. Comput. Simul.. 35 pp. 101-107.
*       with extensive modications for double precisions by    
*         Alan Genz and Yihong Ge
*         Department of Mathematics
*         Washington State University
*         Pullman, WA 99164-3113
*         Email : alangenz@wsu.edu
*
* BVN - calculate the probability that X is larger than SH and Y is
*       larger than SK.
*
* Parameters
*
*   SH  REAL, integration limit
*   SK  REAL, integration limit
*   R   REAL, correlation coefficient
*   LG  INTEGER, number of Gauss Rule Points and Weights
*
      DOUBLE PRECISION BVN, SH, SK, R, ZERO, TWOPI 
      INTEGER I, LG, NG
      PARAMETER ( ZERO = 0, TWOPI = 6.283185307179586D0 ) 
      DOUBLE PRECISION X(10,3), W(10,3), AS, A, B, C, D, RS, XS
      DOUBLE PRECISION MVPHI, SN, ASR, H, K, BS, HS, HK
      SAVE X, W
*     Gauss Legendre Points and Weights, N =  6
      DATA ( W(I,1), X(I,1), I = 1,3) /
     *  0.1713244923791705D+00,-0.9324695142031522D+00,
     *  0.3607615730481384D+00,-0.6612093864662647D+00,
     *  0.4679139345726904D+00,-0.2386191860831970D+00/
*     Gauss Legendre Points and Weights, N = 12
      DATA ( W(I,2), X(I,2), I = 1,6) /
     *  0.4717533638651177D-01,-0.9815606342467191D+00,
     *  0.1069393259953183D+00,-0.9041172563704750D+00,
     *  0.1600783285433464D+00,-0.7699026741943050D+00,
     *  0.2031674267230659D+00,-0.5873179542866171D+00,
     *  0.2334925365383547D+00,-0.3678314989981802D+00,
     *  0.2491470458134029D+00,-0.1252334085114692D+00/
*     Gauss Legendre Points and Weights, N = 20
      DATA ( W(I,3), X(I,3), I = 1,10) /
     *  0.1761400713915212D-01,-0.9931285991850949D+00,
     *  0.4060142980038694D-01,-0.9639719272779138D+00,
     *  0.6267204833410906D-01,-0.9122344282513259D+00,
     *  0.8327674157670475D-01,-0.8391169718222188D+00,
     *  0.1019301198172404D+00,-0.7463319064601508D+00,
     *  0.1181945319615184D+00,-0.6360536807265150D+00,
     *  0.1316886384491766D+00,-0.5108670019508271D+00,
     *  0.1420961093183821D+00,-0.3737060887154196D+00,
     *  0.1491729864726037D+00,-0.2277858511416451D+00,
     *  0.1527533871307259D+00,-0.7652652113349733D-01/
      IF ( ABS(R) .LT. 0.3 ) THEN
         NG = 1
         LG = 3
      ELSE IF ( ABS(R) .LT. 0.75 ) THEN
         NG = 2
         LG = 6
      ELSE 
         NG = 3
         LG = 10
      ENDIF
      H = SH
      K = SK 
      HK = H*K
      BVN = 0
      IF ( ABS(R) .LT. 0.925 ) THEN
         HS = ( H*H + K*K )/2
         ASR = ASIN(R)
         DO I = 1, LG
            SN = SIN(ASR*( X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
            SN = SIN(ASR*(-X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
         END DO
         BVN = BVN*ASR/(2*TWOPI) + MVPHI(-H)*MVPHI(-K) 
      ELSE
         IF ( R .LT. 0 ) THEN
            K = -K
            HK = -HK
         ENDIF
         IF ( ABS(R) .LT. 1 ) THEN
            AS = ( 1 - R )*( 1 + R )
            A = SQRT(AS)
            BS = ( H - K )**2
            C = ( 4 - HK )/8 
            D = ( 12 - HK )/16
            BVN = A*EXP( -(BS/AS + HK)/2 )
     +             *( 1 - C*(BS - AS)*(1 - D*BS/5)/3 + C*D*AS*AS/5 )
            IF ( HK .GT. -160 ) THEN
               B = SQRT(BS)
               BVN = BVN - EXP(-HK/2)*SQRT(TWOPI)*MVPHI(-B/A)*B
     +                    *( 1 - C*BS*( 1 - D*BS/5 )/3 ) 
            ENDIF
            A = A/2
            DO I = 1, LG
               XS = ( A*(X(I,NG)+1) )**2
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I,NG)*
     +              ( EXP( -BS/(2*XS) - HK/(1+RS) )/RS 
     +              - EXP( -(BS/XS+HK)/2 )*( 1 + C*XS*( 1 + D*XS ) ) )
               XS = AS*(-X(I,NG)+1)**2/4
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I,NG)*EXP( -(BS/XS + HK)/2 )
     +                    *( EXP( -HK*(1-RS)/(2*(1+RS)) )/RS 
     +                       - ( 1 + C*XS*( 1 + D*XS ) ) )
            END DO
            BVN = -BVN/TWOPI
         ENDIF
         IF ( R .GT. 0 ) BVN =  BVN + MVPHI( -MAX( H, K ) )
         IF ( R .LT. 0 ) BVN = -BVN + MAX( ZERO, MVPHI(-H) - MVPHI(-K) )
      ENDIF
      MVBVU = BVN
      END
*
      DOUBLE PRECISION FUNCTION MVPHI( Z )
*     
*     Normal distribution probabilities accurate to 1.e-15.
*     Z = no. of standard deviations from the mean.
*     
*     Based upon algorithm 5666 for the error function, from:
*     Hart, J.F. et al, 'Computer Approximations', Wiley 1968
*     
*     Programmer: Alan Miller
*     
*     Latest revision - 30 March 1986
*     
      DOUBLE PRECISION P0, P1, P2, P3, P4, P5, P6, 
     *     Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     *     Z, P, EXPNTL, CUTOFF, ROOTPI, ZABS
      PARAMETER(
     *     P0 = 220.20 68679 12376 1D0,
     *     P1 = 221.21 35961 69931 1D0, 
     *     P2 = 112.07 92914 97870 9D0,
     *     P3 = 33.912 86607 83830 0D0,
     *     P4 = 6.3739 62203 53165 0D0,
     *     P5 = .70038 30644 43688 1D0, 
     *     P6 = .035262 49659 98910 9D0 )
      PARAMETER(
     *     Q0 = 440.41 37358 24752 2D0,
     *     Q1 = 793.82 65125 19948 4D0, 
     *     Q2 = 637.33 36333 78831 1D0,
     *     Q3 = 296.56 42487 79673 7D0, 
     *     Q4 = 86.780 73220 29460 8D0,
     *     Q5 = 16.064 17757 92069 5D0, 
     *     Q6 = 1.7556 67163 18264 2D0,
     *     Q7 = .088388 34764 83184 4D0 )
      PARAMETER( ROOTPI = 2.5066 28274 63100 1D0 )
      PARAMETER( CUTOFF = 7.0710 67811 86547 5D0 )
*     
      ZABS = ABS(Z)
*     
*     |Z| > 37
*     
      IF ( ZABS .GT. 37 ) THEN
         P = 0
      ELSE
*     
*     |Z| <= 37
*     
         EXPNTL = EXP( -ZABS**2/2 )
*     
*     |Z| < CUTOFF = 10/SQRT(2)
*     
         IF ( ZABS .LT. CUTOFF ) THEN
            P = EXPNTL*( (((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS
     *           + P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS
     *           + Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS
     *           + Q0 )
*     
*     |Z| >= CUTOFF.
*     
         ELSE
            P = EXPNTL/( ZABS + 1/( ZABS + 2/( ZABS + 3/( ZABS 
     *                        + 4/( ZABS + 0.65D0 ) ) ) ) )/ROOTPI
         END IF
      END IF
      IF ( Z .GT. 0 ) P = 1 - P
      MVPHI = P
      END
*
      DOUBLE PRECISION FUNCTION MVPHNV(P)
*
*	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
*
*	Produces the normal deviate Z corresponding to a given lower
*	tail area of P.
*
*	The hash sums below are the sums of the mantissas of the
*	coefficients.   They are included for use in checking
*	transcription.
*
      DOUBLE PRECISION SPLIT1, SPLIT2, CONST1, CONST2, 
     *     A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7, 
     *     C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7, 
     *     E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, 
     *     P, Q, R
      PARAMETER ( SPLIT1 = 0.425, SPLIT2 = 5,
     *            CONST1 = 0.180625D0, CONST2 = 1.6D0 )
*     
*     Coefficients for P close to 0.5
*     
      PARAMETER (
     *     A0 = 3.38713 28727 96366 6080D0,
     *     A1 = 1.33141 66789 17843 7745D+2,
     *     A2 = 1.97159 09503 06551 4427D+3,
     *     A3 = 1.37316 93765 50946 1125D+4,
     *     A4 = 4.59219 53931 54987 1457D+4,
     *     A5 = 6.72657 70927 00870 0853D+4,
     *     A6 = 3.34305 75583 58812 8105D+4,
     *     A7 = 2.50908 09287 30122 6727D+3,
     *     B1 = 4.23133 30701 60091 1252D+1,
     *     B2 = 6.87187 00749 20579 0830D+2,
     *     B3 = 5.39419 60214 24751 1077D+3,
     *     B4 = 2.12137 94301 58659 5867D+4,
     *     B5 = 3.93078 95800 09271 0610D+4,
     *     B6 = 2.87290 85735 72194 2674D+4,
     *     B7 = 5.22649 52788 52854 5610D+3 )
*     HASH SUM AB    55.88319 28806 14901 4439
*     
*     Coefficients for P not close to 0, 0.5 or 1.
*     
      PARAMETER (
     *     C0 = 1.42343 71107 49683 57734D0,
     *     C1 = 4.63033 78461 56545 29590D0,
     *     C2 = 5.76949 72214 60691 40550D0,
     *     C3 = 3.64784 83247 63204 60504D0,
     *     C4 = 1.27045 82524 52368 38258D0,
     *     C5 = 2.41780 72517 74506 11770D-1,
     *     C6 = 2.27238 44989 26918 45833D-2,
     *     C7 = 7.74545 01427 83414 07640D-4,
     *     D1 = 2.05319 16266 37758 82187D0,
     *     D2 = 1.67638 48301 83803 84940D0,
     *     D3 = 6.89767 33498 51000 04550D-1,
     *     D4 = 1.48103 97642 74800 74590D-1,
     *     D5 = 1.51986 66563 61645 71966D-2,
     *     D6 = 5.47593 80849 95344 94600D-4,
     *     D7 = 1.05075 00716 44416 84324D-9 )
*     HASH SUM CD    49.33206 50330 16102 89036
*
*	Coefficients for P near 0 or 1.
*
      PARAMETER (
     *     E0 = 6.65790 46435 01103 77720D0,
     *     E1 = 5.46378 49111 64114 36990D0,
     *     E2 = 1.78482 65399 17291 33580D0,
     *     E3 = 2.96560 57182 85048 91230D-1,
     *     E4 = 2.65321 89526 57612 30930D-2,
     *     E5 = 1.24266 09473 88078 43860D-3,
     *     E6 = 2.71155 55687 43487 57815D-5,
     *     E7 = 2.01033 43992 92288 13265D-7,
     *     F1 = 5.99832 20655 58879 37690D-1,
     *     F2 = 1.36929 88092 27358 05310D-1,
     *     F3 = 1.48753 61290 85061 48525D-2,
     *     F4 = 7.86869 13114 56132 59100D-4,
     *     F5 = 1.84631 83175 10054 68180D-5,
     *     F6 = 1.42151 17583 16445 88870D-7,
     *     F7 = 2.04426 31033 89939 78564D-15 )
*     HASH SUM EF    47.52583 31754 92896 71629
*     
      Q = ( 2*P - 1 )/2
      IF ( ABS(Q) .LE. SPLIT1 ) THEN
         R = CONST1 - Q*Q
         MVPHNV = Q*( ( ( ((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     *                  *R + A2 )*R + A1 )*R + A0 )
     *            /( ( ( ((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     *                  *R + B2 )*R + B1 )*R + 1 )
      ELSE
         R = MIN( P, 1 - P )
         IF ( R .GT. 0 ) THEN
            R = SQRT( -LOG(R) )
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               MVPHNV = ( ( ( ((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     *                      *R + C2 )*R + C1 )*R + C0 ) 
     *                /( ( ( ((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     *                      *R + D2 )*R + D1 )*R + 1 )
            ELSE
               R = R - SPLIT2
               MVPHNV = ( ( ( ((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     *                      *R + E2 )*R + E1 )*R + E0 )
     *                /( ( ( ((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     *                      *R + F2 )*R + F1 )*R + 1 )
            END IF
         ELSE
            MVPHNV = 9
         END IF
         IF ( Q .LT. 0 ) MVPHNV = - MVPHNV
      END IF
      END
