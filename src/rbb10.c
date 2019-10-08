#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random bivariate U(0,1) from bb10 copula 
   gcc -DMAIN -o rbb10 rbb10.c -lm
   gcc -DDEBUG -DMAIN -o rbb10d rbb10.c -lm
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double th,pi,u1,u2;
  long seed;
  double s1,s2,s12;
  int isim,nsim;
  void rbb10(double, double, double *, double *);
  setbuf(stdout,NULL);
  scanf("%lf %lf %d %ld", &th,&pi,&nsim,&seed);
  // th>0, 0<pi<=1
  while(th>0.)
  { srand(seed);
    printf("\ntheta=%f pi=%f nsim=%d seed=%ld\n", th,pi,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rbb10(th,pi,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    printf("\n============================================================\n");
    scanf("%lf %lf %d %ld", &th,&pi,&nsim,&seed);
  }
}

/* inputs
     th = theta parameter >0 
     pi = pi parameter in (0,1]
   outputs
     u1, u2 = random pair in (0,1) with BB10 copula 
*/
void rbb10(double th, double pi, double *u1, double *u2)
{ double p,pp;
  void qcondbb10(double *p, double *u, double *, double *, double *);
  void pcondbb10(double *v, double *u, double *, double *, double *);
  *u1=rand()/2147483648.;
  p=rand()/2147483648.;
  qcondbb10(&p,u1,&th,&pi,u2);  // return u2
  pcondbb10(u2,u1,&th,&pi,&pp); // return pp, backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     cpar = copula parameter (th.pi). th>0, 0<pi<=1
   outputs
     uvec, vvec = random sample of size nn from BB10 copula 
*/
void rbb10(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondbb10(double *,double *,double *, double *, double*);
  double u,v,p,th,pi;
  int i,n;
  n=*nn; th=cpar[0]; pi=cpar[1];
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    p=unif_rand();
    qcondbb10(&p,&u,&th,&pi,&v); // return v
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* inputs
     v0, u0 = values in (0,1)
     cpar = copula parameter (th.pi). th>0, 0<pi<=1
   output
     p = conditional cdf
*/
void pcondbb10(double *v0,double *u0,double *th0, double *pi0, double *p)
{ double th,pi,u,v,ut,vt,tem,ttem,ccdf;
  th=*th0; pi=*pi0; u=*u0; v=*v0; 
  ut=pow(u,th); vt=pow(v,th);
  tem=1.-pi*(1.-ut)*(1.-vt);
  ttem=pow(tem,-1./th);
  ccdf=ttem/tem*v*(1.-pi+pi*vt);
  *p=ccdf;
}

/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter (th.pi). th>0, 0<pi<=1
   output
     vv = inverse of conditional cdf
*/
void qcondbb10(double *p0, double *u0, double *th0, double *pi0, double *vv)
{ double p,u,th,pi,th1,ut,vt,tem,ttem,ccdf,pdf,h,hp,v;
  double diff,mxdif,eps;
  int iter,mxiter;

  th=*th0; pi=*pi0; p=*p0; u=*u0; 
  //mxiter=30;
  mxiter=20; 
  eps=1.e-6;
  th1=1./th;
  ut=pow(u,th);
  mxdif=1.; iter=0; 
  diff=1.; 
  v=.7*u; // what is good starting point?
  while(mxdif>eps && iter<mxiter)
  { vt=pow(v,th);
    tem=1.-pi*(1.-ut)*(1.-vt);
    ttem=pow(tem,-th1);
    ccdf=ttem/tem*v*(1.-pi+pi*vt);
    pdf=ttem/tem/tem;
    pdf=pdf*(1.-pi+pi*(1.+th)*ut*vt-pi*(1.-pi)*(1.-ut)*(1.-vt));
    h=ccdf-p; hp=pdf;
    diff=h/hp;
    v=v-diff;
    iter++;
    while(v<=0. || v>=1.) { diff/=2.; v+=diff; }
#ifdef DEBUG
    printf("%d %f %f\n", iter, diff, v);
#endif
    mxdif=fabs(diff);
  }
#ifdef MAIN
  if(iter>=mxiter) 
  { printf("***did not converge\n");
    printf("p=%f, u=%f, theta=%f, pi=%f, lastv=%f\n", p,u,th,pi,v);
  }
#else
  if(iter>=mxiter) 
  { Rprintf("** did not converge **\n");
    Rprintf("p=%f, u=%f, theta=%f, pi=%f, lastv=%f\n", p,u,th,pi,v);
  }
#endif
  *vv=v;
#ifdef DEBUG
  printf("p=%f u=%f v=%f\n", p,u,*vv);
#endif
}

