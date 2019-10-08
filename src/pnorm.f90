! extracted from file mvt.f in the src of R package mvtnorm
! input: z is real-valued
! output: standard normal cdf at z
!DOUBLE PRECISION FUNCTION MVPHI(Z)
DOUBLE PRECISION FUNCTION pnorms(Z)
!     Normal distribution probabilities accurate to 1d-15.
!     Reference: J.L. Schonfelder, Math Comp 32(1978), pp 1232-1240. 
  implicit none
  INTEGER I, IM
  DOUBLE PRECISION A(0:43), BM, B, BP, P, RTWO, T, XA, Z
  PARAMETER( RTWO = 1.414213562373095048801688724209D0, IM = 24 )
  SAVE A
  DATA ( A(I), I = 0, 43 )/                                         &
       6.10143081923200417926465815756D-1,                          & 
      -4.34841272712577471828182820888D-1,                          &
       1.76351193643605501125840298123D-1,                          &
      -6.0710795609249414860051215825D-2,                           &
       1.7712068995694114486147141191D-2,                           &
      -4.321119385567293818599864968D-3,                            &
       8.54216676887098678819832055D-4,                             &
      -1.27155090609162742628893940D-4,                             &
       1.1248167243671189468847072D-5, 3.13063885421820972630152D-7,&      
      -2.70988068537762022009086D-7, 3.0737622701407688440959D-8,   &
       2.515620384817622937314D-9, -1.028929921320319127590D-9,     &
       2.9944052119949939363D-11, 2.6051789687266936290D-11,        &
      -2.634839924171969386D-12, -6.43404509890636443D-13,          &
       1.12457401801663447D-13, 1.7281533389986098D-14,             &
      -4.264101694942375D-15, -5.45371977880191D-16,                &
       1.58697607761671D-16, 2.0899837844334D-17,                   &
      -5.900526869409D-18, -9.41893387554D-19, 2.14977356470D-19,   &
       4.6660985008D-20, -7.243011862D-21, -2.387966824D-21,        &
       1.91177535D-22, 1.20482568D-22, -6.72377D-25, -5.747997D-24, &
      -4.28493D-25, 2.44856D-25, 4.3793D-26, -8.151D-27, -3.089D-27,& 
       9.3D-29, 1.74D-28, 1.6D-29, -8.0D-30, -2.0D-30 /             
     
  XA = ABS(Z)/RTWO
  IF ( XA .GT. 100 ) THEN
     P = 0
  ELSE
     T = ( 8*XA - 30 ) / ( 4*XA + 15 )
     BM = 0
     B  = 0
     DO I = IM, 0, -1 
        BP = B
        B  = BM
        BM = T*B - BP  + A(I)
     END DO
     P = EXP( -XA*XA )*( BM - BP )/4
  END IF
  IF ( Z .GT. 0 ) P = 1 - P
  !MVPHI = P
  pnorms = P
END


! input: z is real-valued
! output: standard normal pdf at z
double precision function dnorms(z)
  implicit none
  double precision z,sqtwpi
  parameter ( sqtwpi = 2.506628274631001D0 )
  dnorms=exp(-z*z/2.d0)/sqtwpi
  return
  end
