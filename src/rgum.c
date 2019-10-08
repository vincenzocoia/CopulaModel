#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random pairs from bivariate Gumbel copula 
   gcc -DMAIN -DNOTR -o rgum rgum.c -lm
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double de,u1,u2,tau;
  long seed;
  double s1,s2,s12;
  int isim,nsim;
  void rgum(double, double *, double *);
  scanf("%lf %d %ld", &de,&nsim,&seed);
  while(de>1.)
  { srand(seed);
    tau=(de-1.)/de;
    printf("\ndelta=%f tau=%f nsim=%d seed=%d\n", de,tau,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rgum(de,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    scanf("%lf %d %ld", &de,&nsim,&seed);
  }
}

double urand()
{ return(rand()/2147483648.); }
#endif

#ifdef NOTR
/* input
     de = cpar = copula parameter >1 
   outputs
     u1, u2 = random pair in (0,1) with Gumbel copula 
*/
void rgum(double de, double *u1, double *u2)
{ double p,pp;
  void qcondgum(double *p, double *u, double *cpar, double *);
  void pcondgum(double *v, double *u, double *cpar, double *);
  double urand();
  p=urand();
  *u1=urand();
  qcondgum(&p,u1,&de,u2);
  pcondgum(u2,u1,&de,&pp); // backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else 
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     cpar = copula parameter >1
   outputs
     uvec, vvec = random sample of size nn from Gumbel copula 
*/
void rgum(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondgum(double *,double *,double *, double*);
  double u,v,qq,de;
  int i,n;

  n=*nn; de=*cpar;
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    qq=unif_rand();
    qcondgum(&qq,&u,&de,&v);
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* inputs
     v0, u0 = values in (0,1)
     cpar = copula parameter >1
   output
     p = conditional cdf
*/
void pcondgum(double *v0,double *u0,double *cpar, double *p)
{ double u,v,x,y,tem1,tem2,sm,tem,ccdf,de;
  de=*cpar; u=*u0; v=*v0;
  x= -log(u); y= -log(v);
  tem1=pow(x,de); tem2=pow(y,de); sm=tem1+tem2; tem=pow(sm,1./de);
  ccdf=exp(-tem);
  ccdf*=pow(1.+tem2/tem1,-1.+1./de);
  ccdf/=u;
  *p=ccdf;
}

/* G(x,y)=exp(-(x^de+y^de)^(1/de)), de>=1
   conditional P(Y>=y|x)=exp(x)*G(x,y)* (x^de+y^de)^(1/de-1) * x^(de-1)
                           = p
   take logs, and let z=(x^de+y^de)^(1/de) >=x
   g(z) = log(p) - x + z + (de-1)*log(z) -(de-1)*log(x) =0
   g'(z) = 1 - (de-1)/z
   solve for z with Newton-Raphson, then solve for y=y(z,x,p,de)
   y = (z^de - x^de)^(1/de)
*/
/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter >1
   output
     vv = inverse of conditional cdf
*/
void qcondgum(double *p0, double *u0, double *cpar, double *vv)
{ double p,u,de,z,g,gp,x,y,con,de1,dif,mxdif,eps;
  int iter,mxiter;

  de=*cpar; p=*p0; u=*u0; mxiter=20; eps=1.e-6;
  x=-log(u); de1=de-1.;
  con=log(p)-x-de1*log(x); 
  z=x*pow(2.,1./de); // this is for z 
  mxdif=1; iter=0; 
  dif=.1;  // needed in case first step leads to NaN
  while(mxdif>eps && iter<mxiter)
  { g=z+de1*log(z)+con;
    gp=1.+de1/z;
    if(isnan(g) || isnan(gp) || isnan(g/gp) ) { dif/=-2.; }  // added for de>50
    else dif=g/gp;
    z-=dif; iter++;
    while(z<=x) { dif/=2.; z+=dif; }
#ifdef DEBUG
    printf("%d %f %f\n", iter, dif, z);
#endif
    mxdif=fabs(dif);
  }
#ifdef MAIN
  if(iter>=mxiter) 
  { printf("***did not converge\n");
    printf("p=%f, x=%f, delta=%f, lastz=%f\n", p,x,de,z);
  }
#endif
/*
#ifndef NOTR
  if(iter>=mxiter) 
  { Rprintf("** did not converge **\n");
    Rprintf("p=%f, x=%f, delta=%f, lastz=%f\n", p,x,de,z);
  }
#endif
*/
  y=pow(pow(z,de)-pow(x,de),1./de);
  *vv=exp(-y);
}

/* reflected Gumbel */
/* input
     p0, u0 = values in (0,1)
     cpar = copula parameter >1
   output
     vv = inverse of conditional cdf
*/
void qcondgumr(double *p0, double *u0, double *cpar, double *vv)
{ void qcondgum(double *p, double *u, double *cpar, double *);
  double p1,u1;
  p1 = 1.-*p0; u1=1.-*u0;
  qcondgum(&p1,&u1,cpar,vv);
  *vv=1.- *vv;
}

