#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random pairs from bivariate Huesler-Reiss copula 
   gcc -DMAINHR -o rhr rhr.c qt.c -lm
   gcc -DDEBUG -DMAINHR -o rhrd rhr.c qt.c -lm
*/
#ifdef MAINHR
int kount;
main(int argc, char *argv[])
{ double cpar,u1,u2,lm;
  long seed;
  double s1,s2,s12,s11,s22;
  int isim,nsim;
  void rhr(double, double *, double *);
  setbuf(stdout,NULL);
  scanf("%lf %d %ld", &cpar,&nsim,&seed);
  while(cpar>0.)
  { srand(seed);
    kount=0;
    lm=pow(2.,-1./cpar);
    printf("\ncpar=%f lm=%f nsim=%d seed=%ld\n", cpar,lm,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rhr(cpar,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
      s11+=u1*u1; s22+=u2*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    s11/=nsim; s22/=nsim; s11=sqrt(s11-s1*s1); s22=sqrt(s22-s2*s2);
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,s12/s11/s22);
    printf("3rd dec place problems: %d\n", kount);
    scanf("%lf %d %ld", &cpar,&nsim,&seed);
  }
}

/* input
     cpar = copula parameter >0 
   outputs
     u1, u2 = random pair in (0,1) with Huesler-Reiss copula 
*/
void rhr(double cpar, double *u1, double *u2)
{ double p,pp;
  void qcondhr(double *p, double *u, double *cpar, double *);
  void pcondhr(double *v, double *u, double *cpar, double *);
  p=rand()/2147483648.;
  *u1=rand()/2147483648.;
  //printf("p=%f, u=%f\n", p,*u1);
  qcondhr(&p,u1,&cpar,u2);
  pcondhr(u2,u1,&cpar,&pp); // backcheck
  if(fabs(p-pp)>1.e-3) 
  { kount++; 
    if(kount<10) printf("*** 3rddecpl %.8f %f %.8f %f\n", *u1,p,*u2,pp);
  }
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     param = cpar = copula parameter >0
   outputs
     uvec, vvec = random sample of size nn from Huesler-Reiss copula 
*/
void rhr(int *nn, double *param, double *uvec, double *vvec)
{ void qcondhr(double *,double *,double *, double*);
  double u,v,qq,cpar;
  int i,n;

  n=*nn; cpar=*param;
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    qq=unif_rand();
    qcondhr(&qq,&u,&cpar,&v);
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* inputs
     v0, u0 = values in (0,1)
     param = cpar = copula parameter >0
   output
     p = conditional cdf
*/
void pcondhr(double *v0,double *u0,double *param, double *p)
{ double u,v,x,y,z,lz,tem1,tem2,p1,p2,lcdf,cdf,ccdf,cpar,cpar1;
  double pnorms(double);
  cpar=*param; u=*u0; v=*v0;
  x= -log(u); y= -log(v);
  z=x/y; cpar1=1./cpar; lz=log(z);
  tem1=cpar1+.5*cpar*lz; tem2=cpar1-.5*cpar*lz;  
  p1=pnorms(tem1); p2=pnorms(tem2);
  lcdf=-x*p1-y*p2;  cdf=exp(lcdf);
  ccdf=cdf*p1/u;
  *p=ccdf;
}

/* inputs
     p0, u0 = values in (0,1)
     param = cpar = copula parameter >0
   output
     vv = inverse of conditional cdf
*/
void qcondhr(double *p0, double *u0, double *param, double *vv)
{ double p,u,cpar,x,g,gp,y,z,lz,tem1,tem2,con,cpar1,p1,p2,diff,mxdif,eps;
  double v1,v2,di,ccdf,vold,v;
  double pnorms(double),dnorms(double);
  int iter,mxiter;
  void pcondhr(double *v0,double *u0,double *param, double *p);
  cpar=*param; p=*p0; u=*u0; mxiter=30; 
  eps=1.e-5;
  x=-log(u);
  con=-log(p); cpar1=1./cpar;
  y=.5*x; 
  if(cpar>1.8)
  { // bisection
    v=u;
    v1=0.; v2=1.;
    diff=1.; vold=v;
    while(diff>eps)
    { pcondhr(&v,u0,param,&ccdf); di=ccdf-p;
      if(di<0) { v1=v; } else { v2=v; }
      diff=v2-v1; vold=v; v=(v1+v2)/2.;
#ifdef DEBUG
      printf("%f %f %f\n", vold,ccdf,di);
#endif
    }
    *vv=v;
    return;
  } 
  /* Newton-Raphson for cpar<=1.8 */
  mxdif=1; iter=0;
  while(mxdif>eps && iter<mxiter)
  { z=x/y; lz=log(z);
    tem1=cpar1+.5*cpar*lz; tem2=cpar1-.5*cpar*lz;
    p1=pnorms(tem1); p2=pnorms(tem2);
    g=x*(1.-p1)-y*p2+log(p1)+con;
    gp= -p2-0.5*dnorms(tem1)/p1/y;
    diff=g/gp;
    y=y-diff;
#ifdef DEBUG
    printf("%d %f %f\n", iter, diff, y);
#endif
    while(y<=0. || fabs(diff)>5.) { diff=diff/2.; y=y+diff; } 
    mxdif=fabs(diff);
    iter=iter+1;
  }
#ifdef DEBUG
  if(iter>=mxiter) printf("***did not converge u=%f p=%f\n",u,p);
#endif
  *vv=exp(-y);
  return;
}

