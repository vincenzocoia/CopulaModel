//#include "Mathlib.h"
#include <math.h>
#include <stdio.h>
#include <float.h>
/* random t(df)  
   gcc -DMAIN -o qt qt.c -lm
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double df,p,x;
  double qt(double,double);
  int i,n;
  scanf("%lf %d", &df,&n);
  while(df>0.)
  { printf("\ndf=%f\n", df);
    for(i=1;i<=n;i++)
    { p=i/(n+1.);
      x=qt(p,df);
      printf("%f %f\n", p,x,"\n");
    }
    scanf("%lf %d", &df,&n);
  }
}
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
 *
 *  Reference:
 *  Algorithm 396: Student's t-quantiles by G.W. Hill
 *  Comm. A.C.M., vol.13(10), 619-620, October 1970
 */

#define M_PI_half 1.57079632679489661923  /* pi/2 */
static double eps = 1.e-12;
double qt(double p, double ndf)
{ double a, b, c, d, prob, P, q, x, y;
  double qnorms(double);
  int neg;
  
  //if (ndf < 1 || p >= 1 || p <= 0) DOMAIN_ERROR;
  if (ndf > 1e20) return qnorms(p);
  if(p > 0.5) { neg = 0; P = 2 * (1 - p); }
  else { neg = 1; P = 2 * p; }
  
  if (fabs(ndf - 2) < eps)
  { /* df ~= 2 */
    q = sqrt(2 / (P * (2 - P)) - 2);
  }
  else if (ndf < 1 + eps)
  { /* df ~= 1 */
    prob = P * M_PI_half;
    q = cos(prob) / sin(prob);
  }
  else
  { /*-- normal case;  including, e.g.,  df = 1.1 */
    a = 1 / (ndf - 0.5);
    b = 48 / (a * a);
    c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
    d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * M_PI_half) * ndf;
    y = pow(d * P, 2 / ndf);
    
    if (y > 0.05 + a)
    { /* Asymptotic inverse expansion about normal */
      x = qnorms(0.5 * P);
      y = x * x;
      if (ndf < 5) c = c + 0.3 * (ndf - 4.5) * (x + 0.6);
      c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
      y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x;
      y = a * y * y;
      if (y > 0.002) y = exp(y) - 1;
      else  y = 0.5 * y * y + y; /* Taylor of  e^y -1 : */
    }
    else
    { y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822)
          * (ndf + 2) * 3) + 0.5 / (ndf + 4))
          * y - 1) * (ndf + 1) / (ndf + 2) + 1 / y;
    }
    q = sqrt(ndf * y);
  }
  if(neg) q = -q;
  return q;
}



/*#include "Mathlib.h"*/
#define M_1_SQRT_2PI  0.398942280401432677939946059934

static double a0 = 2.50662823884;
static double a1 = -18.61500062529;
static double a2 = 41.39119773534;
static double a3 = -25.44106049637;
static double b1 = -8.47351093090;
static double b2 = 23.08336743743;
static double b3 = -21.06224101826;
static double b4 = 3.13082909833;
static double c0 = -2.78718931138;
static double c1 = -2.29796479134;
static double c2 = 4.85014127135;
static double c3 = 2.32121276858;
static double d1 = 3.54388924762;
static double d2 = 1.63706781897;
static double zero = 0.0;
static double half = 0.5;
static double one = 1.0;
static double split = 0.42;

/*double qnorm(double p, double mean, double sd)*/
double qnorms(double p)
{ double q, r, val;
  double pnorms(double),dnorms(double);
  
  /* if (p <= 0.0 || p >= 1.0) DOMAIN_ERROR;*/
  q = p - half;
  if (fabs(q) <= split)
  { /* 0.08 < p < 0.92 */
    r = q * q;
    val = q * (((a3 * r + a2) * r + a1) * r + a0)
        / ((((b4 * r + b3) * r + b2) * r + b1) * r + one);
  }
  else
  { /* p < 0.08 or p > 0.92, set r = min(p,1-p) */
    r = p;
    if (q > zero) r = one - p;
    r = sqrt(-log(r));
    val = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + one);
    if (q < zero) val = -val;
  }
  /* val = val - (pnorm(val, 0.0, 1.0) - p) / dnorm(val, 0.0, 1.0);*/
  val = val - (pnorms(val) - p) / dnorms(val);
  /* return mean + sd * val;*/
  return val;
}

/*double dnorm(double x, double mean, double sd)*/
double dnorms(double x)
{ /* if (sd <= 0.0) DOMAIN_ERROR;
     x = (x - mean) / sd;
     return M_1_SQRT_2PI * exp(-0.5 * x * x) / sd;*/
  return M_1_SQRT_2PI * exp(-0.5 * x * x);
  
}

/* Reference:
 * Cody, W.D. (1993). ALGORITHM 715: SPECFUN - A Portable FORTRAN
 * Package of Special Function Routines and Test Drivers"
 * ACM Transactions on Mathematical Software. 19, 22-32.
 *
 * This function evaluates the normal distribution function:
 * The main computation evaluates near-minimax approximations
 * derived from those in "Rational Chebyshev approximations for
 * the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
 * This transportable program uses rational functions that
 * theoretically approximate the normal distribution function to
 * at least 18 significant decimal digits.  The accuracy achieved
 * depends on the arithmetic system, the compiler, the intrinsic
 * functions, and proper selection of the machine-dependent
 * constants.
 *
 * Mathematical Constants:
 * sqrpi = 1 / sqrt(2*pi),
 * root32 = sqrt(32),
 * thrsh = the argument for which pnorm(thrsh,0,1) = 0.75.
 */


/*float.h:#define DBL_EPSILON 2.2204460492503131e-16*/
/*float.h:#define DBL_MIN 2.2250738585072014e-308*/

/*#include "Mathlib.h"*/

static double c[9] = {
	0.39894151208813466764, 8.8831497943883759412,
	93.506656132177855979, 597.27027639480026226,
	2494.5375852903726711, 6848.1904505362823326,
	11602.651437647350124, 9842.7148383839780218,
	1.0765576773720192317e-8 };

static double d[8] = {
	22.266688044328115691, 235.38790178262499861,
	1519.377599407554805, 6485.558298266760755,
	18615.571640885098091, 34900.952721145977266,
	38912.003286093271411, 19685.429676859990727 };

static double p[6] = {
	0.21589853405795699, 0.1274011611602473639,
	0.022235277870649807, 0.001421619193227893466,
	2.9112874951168792e-5, 0.02307344176494017303 };

static double q[5] = {
	1.28426009614491121, 0.468238212480865118,
	0.0659881378689285515, 0.00378239633202758244,
	7.29751555083966205e-5 };

static double a[5] = {
	2.2352520354606839287, 161.02823106855587881,
	1067.6894854603709582, 18154.981253343561249,
	0.065682337918207449113 };

static double b[4] = {
	47.20258190468824187, 976.09855173777669322,
	10260.932208618978205, 45507.789335026729956} ;

//static double one = 1.0;
//static double half = 0.5;
//static double zero = 0.0;
static double sixten = 1.6;
static double sqrpi = 0.39894228040143267794;
static double thrsh = 0.66291;
static double root32 = 5.656854248;

/*double pnorm(double x, double mean, double sd)*/
double pnorms(double x)
{ static double xden, temp, xnum, result, ccum;
  static double del, min, eps, xsq;
  static double y;
  static int i;
  double fint(double);
  
  eps = DBL_EPSILON * .5;
  min = DBL_MIN;
  y = fabs(x);
  if (y <= thrsh)
  { /* Evaluate pnorm for |z| <= 0.66291 */
    xsq = zero;
    if (y > eps) { xsq = x * x; }
    xnum = a[4] * xsq; xden = xsq;
    for (i = 1; i <= 3; ++i)
    { xnum = (xnum + a[i - 1]) * xsq;
      xden = (xden + b[i - 1]) * xsq;
    }
    result = x * (xnum + a[3]) / (xden + b[3]);
    temp = result; result = half + temp; ccum = half - temp;
  }
  else if (y <= root32)
  { /* Evaluate pnorm for 0.66291 <= |z| <= sqrt(32) */
    xnum = c[8] * y; xden = y;
    for (i = 1; i <= 7; ++i)
    { xnum = (xnum + c[i - 1]) * y;
      xden = (xden + d[i - 1]) * y;
    }
    result = (xnum + c[7]) / (xden + d[7]);
    xsq = fint(y * sixten) / sixten;
    del = (y - xsq) * (y + xsq);
    result = exp(-xsq * xsq * half) * exp(-del * half) * result;
    ccum = one - result;
    if (x > zero) { temp = result; result = ccum; ccum = temp; }
  }
  else
  { /* Evaluate pnorm for |z| > sqrt(32) */
    result = zero;
    xsq = one / (x * x); xnum = p[5] * xsq; xden = xsq;
    for (i = 1; i <= 4; ++i)
    { xnum = (xnum + p[i - 1]) * xsq;
      xden = (xden + q[i - 1]) * xsq;
    }
    result = xsq * (xnum + p[4]) / (xden + q[4]);
    result = (sqrpi - result) / y;
    xsq = fint(x * sixten) / sixten;
    del = (x - xsq) * (x + xsq);
    result = exp(-xsq * xsq * half) * exp(-del * half) * result;
    ccum = one - result;
    if (x > zero)
    { temp = result; result = ccum; ccum = temp; }
  }
  if (result < min) { result = 0.0; }
  if (ccum < min) { ccum = 0.0; }
  return result;
}

double fint(double x)
{ return (x >= 0.0) ? floor(x) : -floor(-x); }

