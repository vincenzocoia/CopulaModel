#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random pairs from bivariate Joe/B5 copula */
// gcc -DMAIN -o rb5 rb5.c -lm
// gcc -DDEBUG -DMAIN -o rb5d rb5.c -lm
#ifdef MAIN
main(int argc, char *argv[])
{ double de,u1,u2,be,lm;
  long seed;
  double s1,s2,s12;
  int isim,nsim;
  void rjoeb5(double, double *, double *);
  setbuf(stdout,NULL);
  scanf("%lf %d %ld", &de,&nsim,&seed);
  // de>1 
  while(de>1.)  
  { srand(seed);
    lm=2.-pow(2.,1./de);
    be=2.*pow(.5,de)-pow(.25,de);
    be=1.-pow(be,1./de); be=4.*be-1.;
    printf("\ndelta=%f beta=%f lm=%f, nsim=%d seed=%d\n", de,be,lm,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rjoeb5(de,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    printf("============================================================\n");
    scanf("%lf %d %ld", &de,&nsim,&seed);
  }
}

/* input
     de = cpar = copula parameter >1 
   output
     u1, u2 = random pair in (0,1) with Joe/B5 copula 
*/
void rjoeb5(double de, double *u1, double *u2)
{ double p,pp;
  void qcondjoe(double *p, double *u, double *de, double *);
  void pcondjoe(double *v, double *u, double *de, double *);
  *u1=rand()/2147483648.;
  p=rand()/2147483648.;
  qcondjoe(&p,u1,&de,u2);  // return u2
  pcondjoe(u2,u1,&de,&pp); // return pp, backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* input
     nn = simulation sample size
     cpar = copula parameter >1
   output
     uvec, vvec = random sample of size nn from Joe/B5 copula 
*/
void rjoeb5(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondjoe(double *,double *,double *, double*);
  double u,v,p,de;
  int i,n;
  n=*nn; de=*cpar;
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    p=unif_rand();
    qcondjoe(&p,&u,&de,&v); // return v
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* input
     v0, u0 = values in (0,1)
     cpar = copula parameter >1
   output
     p = conditional cdf
*/
void pcondjoe(double *v0,double *u0,double *cpar, double *p)
{ double u,v,temv,temu,ccdf,de;
  de=*cpar; u=*u0; v=*v0;
  temv=pow(1.-v,de); temu=pow(1.-u,de);
  ccdf=1.+temv/temu-temv;
  ccdf=pow(ccdf,-1.+1./de);
  ccdf*=(1.-temv);
  *p=ccdf;
}

// can be problems for cpar>=30 
/* solve C_{2|1}(v|u)=p for v given u,p */
/* input
     p0, u0 = values in (0,1)
     cpar = copula parameter >1
   output
     vv = inverse of conditional cdf
*/
void qcondjoe(double *p0, double *u0, double *cpar, double *vv)
{ double ubar,ud,vbar,vd,sm,di,smd,eps;
  double p,u,de,c21,pdf;
  int iter,mxiter;
  double diff,v,de1,dtem,de1inv,tem;

  de=*cpar; p=*p0; u=*u0; mxiter=20; eps=1.e-6;
  ubar = 1.0-u; ud = pow(ubar,de);
  di = 1./de;
  de1=de-1;  // may need better modification for large delta
  dtem=-de1/(1.+de1); de1inv=-1./de1;

#ifdef DEBUG
  printf("p=%f u=%f\n", p,u);
#endif
  // v = 0.5 * (p+u); // old starting guess
  // Use a better starting point based on reflected MTCJ copula
  // A good starting point is crucial when delta is large because
  //    C_{2|1} will be steep
  // C_{R,2|1}(v|u)=1-C_{2|1}(1-v|1-u), 
  // C_{R,2|1}^{-1}(p|u)=1-C_{2|1}^{-1}(1-p|1-u) 
  tem=pow(1.-p,dtem)-1.;
  tem=tem*pow(1.-u,-de1)+1.;
  v=pow(tem,de1inv); v=1.-v;
  diff=1; iter=0;
  while(fabs(diff)>eps && iter<mxiter)
  { vbar = 1.-v;
    vd = pow(vbar,de);
    sm=ud+vd-ud*vd;
    smd = pow(sm,di);
    c21=1.+vd/ud-vd;
    c21=pow(c21,-1.+di);
    c21=c21*(1.-vd);
    pdf=smd*ud*vd*(de1+sm)/sm/sm/ubar/vbar;
    iter++;
    if(isnan(pdf) || isnan(c21) ) { diff/=-2.; }  // added for de>=30
    else diff=(c21-p)/pdf;
    v-=diff;
    while(v<=0 || v>=1 || fabs(diff)>0.25 ) { diff/=2.; v+=diff; }
#ifdef DEBUG
    printf("%d %f\n", iter,v);
#endif
  }
#ifdef MAIN
  if(iter>=mxiter) 
  { printf("** did not converge **\n");
    printf("p=%f, u=%f, delta=%f, lastv=%f, c21=%f\n", p,u,de,v,c21);
  }
#endif
/*
  if(iter>=mxiter) 
  { Rprintf("** did not converge **\n");
    Rprintf("p=%f, u=%f, delta=%f, lastv=%f, c21=%f\n", p,u,de,v,c21);
  }
*/
  *vv=v;
}

