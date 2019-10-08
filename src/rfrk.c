#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random pairs from bivariate Frank copula 
   gcc -DMAIN -DNOTR -o rfrk rfrk.c -lm
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double cpar,u1,u2;
  long seed;
  double s1,s2,s12;
  int isim,nsim;
  void rfrk(double, double *, double *);
  scanf("%lf %d %ld", &cpar,&nsim,&seed);
  while(cpar>1.)
  { srand(seed);
    printf("\ncpar=%f nsim=%d seed=%d\n", cpar,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rfrk(cpar,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    scanf("%lf %d %ld", &cpar,&nsim,&seed);
  }
}

double urand()
{ return(rand()/2147483648.); }
#endif

#ifdef NOTR
/* input
     cpar = copula parameter >0 or <0
   output
     u1, u2 = random pair in (0,1) with Frank copula 
*/
void rfrk(double cpar, double *u1, double *u2)
{ double p,pp;
  void qcondfrk(double *p, double *u, double *cpar, double *);
  void pcondfrk(double *v, double *u, double *cpar, double *);
  double urand();
  p=urand();
  *u1=urand();
  qcondfrk(&p,u1,&cpar,u2);
  pcondfrk(u2,u1,&cpar,&pp); // backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* input
     nn = simulation sample size
     param = cpar = copula parameter >0 or <0
   output
     uvec, vvec = random sample of size nn from Frank copula 
*/
void rfrk(int *nn, double *param, double *uvec, double *vvec)
{ void qcondfrk(double *,double *,double *, double*);
  double u,v,qq,cpar;
  int i,n;
  n=*nn; cpar=*param;
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    qq=unif_rand();
    qcondfrk(&qq,&u,&cpar,&v);
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* input
     v0, u0 = values in (0,1)
     cpar = copula parameter >0 or <0
   output
     p = conditional cdf
*/
void pcondfrk(double *v0, double *u0, double *cpar, double *p)
{ double u,v,de,de1,tem,ccdf;
  de=*cpar; v=*v0; u=*u0; 
  de1=1.-exp(-de);
  tem=1.-exp(-de*u);
  ccdf=(1.-tem)/(de1/(1.-exp(-de*v))-tem);
  *p=ccdf;
}

/* input
     p0, u0 = values in (0,1)
     cpar = copula parameter >0 or <0
   output
     vv = inverse of conditional cdf
*/
void qcondfrk(double *p0, double *u0, double *cpar, double *vv)
{ double p,u,de,de1,tem;
  de=*cpar; p=*p0; u=*u0; 
  if(de==0.) { *vv=p; return; }
  de1=1.-exp(-de);
  tem=1.-de1/((1./p-1.)*exp(-de*u)+1.);
  *vv=(-log(tem))/de;
}
