#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random t(df)  
   gcc -DMAIN -DNOTR -o rt rt.c -lm
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double df,x;
  long seed;
  double s1,s2;
  int isim,nsim;
  double rt(double);
  scanf("%lf %d %ld", &df,&nsim,&seed);
  while(df>0.)
  { srand(seed);
    printf("\ndf=%f nsim=%d seed=%d\n", df,nsim,seed);
    for(isim=1,s1=0.,s2=0.;isim<=nsim;isim++)
    { x=rt(df);
      s1+=x; s2+=x*x; 
    }
    s1/=nsim; s2/=nsim;  s2-=s1*s1;
    printf("mean=%f, var=%f\n", s1,s2);
    printf("theorvar=%f\n", df/(df-2.));
    scanf("%lf %d %ld", &df,&nsim,&seed);
  }
}
double urand()
{ return(rand()/2147483648.); }

#endif

/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Calls rchisq and rnorm to do the real work */

/* input
     df = degree of freedom >0
   output
     random variate from t(df)
*/
double rt(double df)
/*{ return snorm() / sqrt(rchisq(df) / df); }*/
{ double snorm();
  double rgamma(double,double);
  return snorm() / sqrt(rgamma(df/2.,2.) / df); 
}

/*double rchisq(double df)
{ return rgamma(df / 2.0, 2.0); } */

/*
 *	References:
 *
 *	[1] Shape parameter a >= 1.  Algorithm GD in:
 *
 *		Ahrens, J.H. and Dieter, U. (1982).
 *		Generating gamma variates by a modified
 *		rejection technique.
 *		Comm. ACM, 25, 47-54.
 *
 *
 *	[2] Shape parameter 0 < a < 1. Algorithm GS in:
 *
 *              Ahrens, J.H. and Dieter, U. (1974).
 *		Computer methods for sampling from gamma, beta,
 *		poisson and binomial distributions.
 *		Computing, 12, 223-246.                       C
 *
 *	Input: a = parameter (mean) of the standard gamma distribution.
 *	Output: a variate from the gamma(a)-distribution
 *
 *	Coefficients q(k) - for q0 = sum(q(k)*a**(-k))
 *	Coefficients a(k) - for q = q0+(t*t/2)*sum(a(k)*v**k)
 *	Coefficients e(k) - for exp(q)-1 = sum(e(k)*q**k)
 */

static double a1 = 0.3333333;
static double a2 = -0.250003;
static double a3 = 0.2000062;
static double a4 = -0.1662921;
static double a5 = 0.1423657;
static double a6 = -0.1367177;
static double a7 = 0.1233795;
static double e1 = 1.0;
static double e2 = 0.4999897;
static double e3 = 0.166829;
static double e4 = 0.0407753;
static double e5 = 0.010293;
static double q1 = 0.04166669;
static double q2 = 0.02083148;
static double q3 = 0.00801191;
static double q4 = 0.00144121;
static double q5 = -7.388e-5;
static double q6 = 2.4511e-4;
static double q7 = 2.424e-4;
static double sqrt32 = 5.656854;
static double aa = 0.;
static double aaa = 0.;

#define repeat for(;;)

double rgamma(double a, double scale)
{ static double b, c, d, e, p, q, r, s, t, u, v, w, x;
  static double q0, s2, si;
  double retval;
  double snorm(),urand();
  
  if (a < 1.0)
  { /* alternate method for parameters a below 1 */
    /* 0.36787944117144232159 = exp(-1) */
    aa = 0.0;
    b = 1.0 + 0.36787944117144232159 * a;
    repeat
    { //p = b * sunif();
      p = b *urand();
      if (p >= 1.0)
      { retval = -log((b - p) / a);
        //if (sexp() >= (1.0 - a) * log(retval)) break;
        if (-log(urand()) >= (1.0 - a) * log(retval)) break;
      }
      else
      { retval = exp(log(p) / a);
        //if (sexp() >= retval) break;
        if (-log(urand()) >= retval) break;
      }
    }
    return scale * retval;
  }
  /* Step 1: Recalculations of s2, s, d if a has changed */
  if (a != aa)
  { aa = a;
    s2 = a - 0.5;
    s = sqrt(s2);
    d = sqrt32 - s * 12.0;
  }
  /* Step 2: t = standard normal deviate, */
  /* x = (s,1/2)-normal deviate. */
  /* immediate acceptance (i) */
  
  t = snorm();
  x = s + 0.5 * t;
  retval = x * x;
  if (t >= 0.0) return scale * retval;
  
  /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
  //u = sunif();
  u=urand();
  if (d * u <= t * t * t)
  { return scale * retval; }
  /* Step 4: recalculations of q0, b, si, c if necessary */
  
  if (a != aaa)
  { aaa = a;
    r = 1.0 / a;
    q0 = ((((((q7*r + q6) *r + q5) *r + q4) *r + q3) *r + q2) *r + q1) * r;
    
    /* Approximation depending on size of parameter a */
    /* The constants in the expressions for b, si and */
    /* c were established by numerical experiments */
    
    if (a <= 3.686)
    { b = 0.463 + s + 0.178 * s2;
      si = 1.235;
      c = 0.195 / s - 0.079 + 0.16 * s;
    }
    else if (a <= 13.022)
    { b = 1.654 + 0.0076 * s2;
      si = 1.68 / s + 0.275;
      c = 0.062 / s + 0.024;
    }
    else
    { b = 1.77;
      si = 0.75;
      c = 0.1515 / s;
    }
  }
  /* Step 5: no quotient test if x not positive */
  
  if (x > 0.0)
  { /* Step 6: calculation of v and quotient q */
    v = t / (s + s);
    if (fabs(v) <= 0.25)
        q = q0 + 0.5 * t * t * ((((((a7 * v + a6)
        * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
    else
        q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
    
    /* Step 7: quotient acceptance (q) */
    
    if (log(1.0 - u) <= q) return scale * retval;
  }
  /* Step 8: e = standard exponential deviate */
  /* u= 0,1 -uniform deviate */
  /* t=(b,si)-double exponential (laplace) sample */
  
  repeat
  { // e = sexp();
    e=-log(urand());
    //u = sunif();
    u=urand();
    u = u + u - 1.0;
    if (u < 0.0) t = b - si * e;
    else t = b + si * e;
    /* Step  9:  rejection if t < tau(1) = -0.71874483771719 */
    if (t >= -0.71874483771719)
    { /* Step 10:  calculation of v and quotient q */
      v = t / (s + s);
      if (fabs(v) <= 0.25)
      {   q = q0 + 0.5 * t * t * ((((((a7 * v + a6)
            * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
      }
      else
      {   q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v); }
      /* Step 11:  hat acceptance (h) */
      /* (if q not positive go to step 8) */
      if (q > 0.0)
      { if (q <= 0.5) w = ((((e5 * q + e4) * q + e3) * q + e2) * q + e1) * q;
        else w = exp(q) - 1.0;
        /* if t is rejected */
        /* sample again at step 8 */
        if (c * fabs(u) <= w * exp(e - 0.5 * t * t)) break;
      }
    }
  }
  x = s + 0.5 * t;
  return scale * x * x;
}


//#define repeat for(;;)

/*  Kinderman A. J. and Ramage J. G. (1976).
   *  Computer generation of normal random variables.
   *  JASA 71, 893-896.
   */

#define C1  0.398942280401433
#define C2  0.180025191068563
#define g(x)  (C1*exp(-x*x/2.0)-C2*(a-fabs(x)))

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

static double a =  2.216035867166471;

double snorm()
{ double t, u1, u2, u3;
  //double sunif();
  double urand();
  
  u1 = urand(); // sunif();
  if( u1<0.884070402298758 )
  { u2 = urand(); // sunif();
    return a*(1.13113163544180*u1+u2-1);
  }
  
  if( u1>=0.973310954173898 )
  { tail:
    u2 = urand(); // sunif();
    u3 = urand(); // sunif();
    t = (a*a-2*log(u3));
    if( u2*u2<(a*a)/t ) return (u1<0.986655477086949) ? sqrt(t) : -sqrt(t) ;
    goto tail;
  }
  
  if( u1>=0.958720824790463 )
  { region3:
    u2 = urand(); // sunif();
    u3 = urand(); // sunif();
    t = a-0.630834801921960*min(u2,u3);
    if( max(u2,u3)<=0.755591531667601 ) return (u2<u3) ? t : -t ;
    if( 0.034240503750111*fabs(u2-u3)<=g(t) ) return (u2<u3) ? t : -t ;
    goto region3;
  }
  
  if( u1>=0.911312780288703 )
  { region2:
    u2 = urand(); // sunif();
    u3 = urand(); // sunif();
    t = 0.479727404222441+1.105473661022070*min(u2,u3);
    if( max(u2,u3)<=0.872834976671790 ) return (u2<u3) ? t : -t ;
    if( 0.049264496373128*fabs(u2-u3)<=g(t) ) return (u2<u3) ? t : -t ;
    goto region2;
  }
  
  region1:
  u2 = urand(); // sunif();
  u3 = urand(); // sunif();
  t = 0.479727404222441-0.595507138015940*min(u2,u3);
  if( max(u2,u3)<=0.805577924423817 ) return (u2<u3) ? t : -t ;
  goto region1;
}

