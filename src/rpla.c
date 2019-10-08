#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random pairs from bivariate Plackett copula 
   gcc -DMAINP -DNOTR -o rpla rpla.c -lm
*/
#ifdef MAINP
main(int argc, char *argv[])
{ double cpar,u1,u2,rho;
  long seed;
  double s1,s2,s12;
  int isim,nsim;
  void rpla(double, double *, double *);
  scanf("%lf %d %ld", &cpar,&nsim,&seed);
  while(cpar>1.)
  { srand(seed);
    rho=(cpar+1.)/(cpar-1.) - (2.*cpar*log(cpar))/(cpar-1.)/(cpar-1.);
    printf("\ncpar=%f rho=%f nsim=%d seed=%d\n", cpar,rho,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rpla(cpar,&u1,&u2);
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
     cpar = copula parameter >0 
   output
     u1, u2 = random pair in (0,1) with Plackett copula 
*/
void rpla(double cpar, double *u1, double *u2)
{ double p,pp;
  void qcondpla(double *p, double *u, double *cpar, double *);
  void pcondpla(double *v, double *u, double *cpar, double *);
  double urand();
  p=urand();
  *u1=urand();
  qcondpla(&p,u1,&cpar,u2);
  pcondpla(u2,u1,&cpar,&pp); // backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else 
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* input
     nn = simulation sample size
     cpar = copula parameter >0
   output
     uvec, vvec = random sample of size nn from Plackett copula 
*/
void rpla(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondpla(double *,double *,double *, double*);
  double u,v,qq,de;
  int i,n;

  n=*nn; de=*cpar;
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    qq=unif_rand();
    qcondpla(&qq,&u,&de,&v);
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif


/* Plackett: C_{2|1}^{-1}(p|u;cpar) */
/* input
     p0, u0 = values in (0,1)
     cpar0 = cpar = copula parameter >0
   output
     v0 = inverse of conditional cdf
*/
void qcondpla(double *p0, double *u0, double *cpar0, double *v0)
{ double p,u,cpar,v;
  double dif,eps,cpar1,tem,tem0,tem1,tem2,tem3,pdf,ccdf;
  int iter,mxiter;
  p=*p0; u=*u0; cpar=*cpar0; mxiter=30; eps=1.e-8;
  if(cpar==1.) { *v0=p; return; }
  cpar1=cpar-1.;
  iter=0; dif=1.;
  v=u;
  while(iter<mxiter && fabs(dif)>eps)
  { tem=1.+cpar1*(u+v);
    tem0=2.*cpar1*u*v;
    tem1=tem*tem-2.*cpar*tem0;  tem2=sqrt(tem1);
    tem3=tem-tem0;
    pdf=cpar*tem3/tem1/tem2;
    ccdf=.5-.5*(cpar1*u+1.-(cpar1+2.)*v)/tem2;
    ccdf-=p;
    dif=ccdf/pdf;
    v-=dif;
    while(v<0. || v>1.) { dif/=2.; v+=dif;}
    iter++;
#ifdef QU
    printf("%3d %8.4f %10.3e\n", iter, v,dif);
#endif
  }
#ifdef QU
  if(iter>=mxiter) printf("*** did not converge \n");
#endif
  *v0=v;
}


/* C_{2|1}(v|u;cpar) */
/* input
     v0, u0 = values in (0,1)
     cpar0 = cpar = copula parameter >0
   output
     p = conditional cdf
*/
void pcondpla(double *v0,double *u0,double *cpar0, double *p)
{ double v,u,cpar,cpar1,tem1,tem2,ccdf,tem;
  v=*v0; u=*u0; cpar=*cpar0;
  if(cpar==1.) { *p=v; return; }
  cpar1=cpar-1.;
  tem=1.+cpar1*(u+v); tem1=tem*tem-4.*cpar*cpar1*u*v;
  tem2=sqrt(tem1);
  ccdf=(cpar1*u+1.-(cpar1+2.)*v)/tem2;
  ccdf=.5*(1.-ccdf);
  *p=ccdf;
}


