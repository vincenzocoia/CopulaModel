#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define DBL_EPSILON 2.225074e-308  // don't know if this is correct
/* gcc -DMAIN -o ptbeta ptbeta.c -lm
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double pt_(double*,double*);
  double x,nu,pr;
  
  scanf("%lf", &x); 
  for(nu=1;nu<=10;nu+=1)
  { pr=pt_(&x,&nu);
    printf("%f %f %f\n", x,nu,pr);
  }
}
#endif

// from R source in 1997
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

/*
 *    Incomplete Beta Integral
 *
 *    Taken from Cephes Math Library, Release 2.3:  March, 1995
 *    Copyright 1984, 1995 by Stephen L. Moshier
 *    Changes for R : Copyright 1997 by Ross Ihaka
 *		      More comments + minor cosmetic: Martin Maechler, May 1997
 */

//#include "Mathlib.h"

static double incbcf(double, double, double), incbd(double, double, double), 
  pseries(double, double, double, double);

static double big = 4.503599627370496e15;
static double biginv = 2.22044604925031308085e-16;


double pbetab(double xx, double aa, double bb, double logbeta)
{
  /* logbeta == log(beta(aa,bb)) = log(beta(bb,aa)) */
  double a, b, t, x, xc, w, y;
  int swap_tail;
  
  //if (aa <= 0 || bb <= 0) goto domerr;
  if (xx <= 0 || xx >= 1)
  { if (xx == 0) return 0;
    if (xx == 1) return 1;
    //domerr: errno = EDOM;
    //return 0;
  }
  swap_tail = 0;
  if (bb * xx <= 1 && xx <= 0.95)
  { t = pseries(aa, bb, xx, logbeta);
    goto done;
  }
  w = 1 - xx;
  
  /* Reverse a and b if x is greater than the mean. */
  if (xx > (aa / (aa + bb)))
  { swap_tail = 1;
    a = bb; b = aa;
    xc = xx; x = w;
  }
  else
  { a = aa; b = bb;
    xc = w; x = xx;
  }
  if (swap_tail && (b * x) <= 1 && x <= 0.95)
  { t = pseries(a, b, x, logbeta);
    goto done;
  }
  
  /* Take __Continued Fraction__ expansion.
   * Choose the one with better convergence. 
   * Both are from  Abramowitz & Stegun, sec. 26.5
   */
  
  if (x * (a + b - 2) < a - 1) w = incbcf(a, b, x);
  else w = incbd(a, b, x) / xc; /* A & S, 26.5.9 */
  
  /* Multiply w by the factor
      a     b    _             _     _
     x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */
  y = a * log(x);
  t = b * log(xc);
  
  /* Resort to logarithms.  */
  y += t - logbeta; /* = lgamma(a + b) - lgamma(a) - lgamma(b); */
  y += log(w / a);
  t = exp(y);
  done: 
  if (swap_tail)  // what is DBL_EPSILON ?
  { if (t <= DBL_EPSILON) t = 1 - DBL_EPSILON;
    else t = 1 - t;
  }
  return t;
}


/* Continued fraction expansion #1.  Abramowitz & Stegun, 26.5.8  */

static double incbcf(double a, double b, double x)
{
  double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
  double k1, k2, k3, k4, k5, k6, k7, k8;
  double r, rel_err, ans, rel_tol = 3 * DBL_EPSILON;
  int n;
  
  k1 = a; k2 = a + b;
  k3 = a; k4 = a + 1;
  k5 = 1;  k6 = b - 1;
  k7 = a + 1; k8 = a + 2;
  pkm2 = 0; qkm2 = 1;
  pkm1 = 1; qkm1 = 1;
  ans = 1;
  r = 1;
  for(n=0; n < 300; n++)
  { xk = -(x * k1 * k2) / (k3 * k4);
    pk = pkm1 + pkm2 * xk;
    qk = qkm1 + qkm2 * xk;
    pkm2 = pkm1; qkm2 = qkm1;
    pkm1 = pk;   qkm1 = qk;
    
    xk = (x * k5 * k6) / (k7 * k8);
    pk = pkm1 + pkm2 * xk;
    qk = qkm1 + qkm2 * xk;
    pkm2 = pkm1; qkm2 = qkm1;
    pkm1 = pk;   qkm1 = qk;
    
    if (qk != 0) r = pk / qk;
    if (r != 0)
    { rel_err = fabs((ans - r) / r);
      ans = r;
    }
    else rel_err = 1;
    
    if (rel_err < rel_tol) break;
    k1 += 1; k2 += 1;
    k3 += 2; k4 += 2;
    k5 += 1; k6 -= 1;
    k7 += 2; k8 += 2;
    
    /* re-normalize  numerators and denominators if necessary */
    if ((fabs(qk) + fabs(pk)) > big)
    { pkm2 *= biginv; pkm1 *= biginv;
      qkm2 *= biginv; qkm1 *= biginv;
    }
    if ((fabs(qk) < biginv) || (fabs(pk) < biginv))
    { pkm2 *= big; pkm1 *= big;
      qkm2 *= big; qkm1 *= big;
    }
  }
  return ans;
}

/* Continued fraction expansion #2.  Abramowitz & Stegun, 26.5.9  */
static double incbd(double a, double b, double x)
{
  double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
  double k1, k2, k3, k4, k5, k6, k7, k8;
  double r, rel_err, ans, z, rel_tol = 3 * DBL_EPSILON;
  int n;
  
  k1 = a;  k2 = b - 1;
  k3 = a;  k4 = a + 1;
  k5 = 1;  k6 = a + b;
  k7 = a + 1; k8 = a + 2;
  pkm2 = 0; qkm2 = 1;
  pkm1 = 1; qkm1 = 1;
  z = x / (1 - x);
  ans = 1;
  r = 1;
  for(n=0; n < 300; n++)
  { xk = -(z * k1 * k2) / (k3 * k4);
    pk = pkm1 + pkm2 * xk;
    qk = qkm1 + qkm2 * xk;
    pkm2 = pkm1; qkm2 = qkm1;
    pkm1 = pk;   qkm1 = qk;
    
    xk = (z * k5 * k6) / (k7 * k8);
    pk = pkm1 + pkm2 * xk;
    qk = qkm1 + qkm2 * xk;
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    
    if (qk != 0) r = pk / qk;
    if (r != 0)
    { rel_err = fabs((ans - r) / r);
      ans = r;
      if (rel_err < rel_tol) break;
    }
    k1 += 1; k2 -= 1;
    k3 += 2; k4 += 2;
    k5 += 1; k6 += 1;
    k7 += 2; k8 += 2;
    
    /* re-normalize  numerators and denominators if necessary */
    if ((fabs(qk) + fabs(pk)) > big)
    { pkm2 *= biginv; pkm1 *= biginv;
      qkm2 *= biginv; qkm1 *= biginv;
    }
    if ((fabs(qk) < biginv) || (fabs(pk) < biginv))
    { pkm2 *= big; pkm1 *= big;
      qkm2 *= big; qkm1 *= big;
    }
  }
  return ans;
}


/* Power series for incomplete beta integral.
 *  Use when b*x is small and x not too close to 1, here  (b*x) <= 1 & x <= 0.95
   */
static double pseries(double a, double b, double x, double logbeta)
{
  /* logbeta == log(beta(a,b)) = log(beta(b,a)) */
  double s, t, u, v, n, t1, z, ai;
  
  ai = 1 / a;
  u = (1 - b) * x;
  v = u / (a + 1);
  t1 = v;
  t = u;
  n = 2;
  s = 0;
  z = DBL_EPSILON * ai;
  while (fabs(v) > z)
  { u = (n - b) * x / n;
    t *= u;
    v = t / (a + n);
    s += v;
    n += 1;
  }
  s += t1;
  s += ai;
  u = a * log(x);
  t = u + log(s) - logbeta;
  s = exp(t);
  return s;
}
//#include "Mathlib.h"

double lbeta(double a, double b)
{ return lgamma(a) + lgamma(b) - lgamma(a + b); }

double pbeta(double x, double a, double b)
{ return pbetab(x, a, b, lbeta(a,b)); }

double pt_(double *x0, double *nu0)
{ double val,x,nu;
  x=*x0; nu=*nu0;
  val = 0.5 * pbeta(nu / (nu + x * x), nu / 2.0, 0.5);
  //printf ("%f %f %f\n", x,nu,val);
  return (x > 0.0) ? 1 - val : val;
}
