#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random bivariate U(0,1) from BB8 copula 
   gcc -DMAIN -o rbb8 rbb8.c -lm
   gcc -DDEBUG -DMAIN -o rbb8d rbb8.c -lm
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double vth,de,u1,u2;
  long seed;
  double s1,s2,s12;
  int isim,nsim;
  void rbb8(double, double, double *, double *);
  setbuf(stdout,NULL);
  scanf("%lf %lf %d %ld", &vth,&de,&nsim,&seed);
  // vth>1, 0<de<=1
  while(vth>1.)
  { srand(seed);
    printf("\nvtheta=%f delta=%f nsim=%d seed=%d\n", vth,de,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rbb8(vth,de,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    printf("\n============================================================\n");
    scanf("%lf %lf %d %ld", &vth,&de,&nsim,&seed);
  }
}

/* inputs
     vth = theta parameter >1 
     de = delta parameter in (0,1) 
   outputs
     u1, u2 = random pair in (0,1) with BB8 copula 
*/
void rbb8(double vth, double de, double *u1, double *u2)
{ double p,pp;
  void qcondbb8(double *p, double *u, double *, double *, double *);
  void pcondbb8(double *v, double *u, double *, double *, double *);
  *u1=rand()/2147483648.;
  p=rand()/2147483648.;
  qcondbb8(&p,u1,&vth,&de,u2);  // return u2
  pcondbb8(u2,u1,&vth,&de,&pp); // return pp, backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     cpar = copula parameter (vth.de). vth>1, 0<de<1
   outputs
     uvec, vvec = random sample of size nn from BB8 copula 
*/
void rbb8(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondbb8(double *,double *,double *, double *, double*);
  double u,v,p,vth,de;
  int i,n;
  n=*nn; vth=cpar[0]; de=cpar[1];
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    p=unif_rand();
    qcondbb8(&p,&u,&vth,&de,&v); // return v
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* inputs
     v0, u0 = values in (0,1)
     cpar = copula parameter (vth.de). vth>1, 0<de<1
   output
     p = conditional cdf
*/
void pcondbb8(double *v0,double *u0,double *vth0, double *de0, double *p)
{ double u,v,vth,de,ga1,ut,x,y,tem,den,ccdf;
  vth=*vth0; de=*de0; u=*u0; v=*v0; 
  ga1=1./(1.-pow(1.-de,vth)); // reciprocal of ga
  ut=pow(1.-de*u,vth);
  x=1.-ut;
  y=1.-pow(1.-de*v,vth);
  tem=pow(1.-ga1*x*y,1./vth-1.);
  den=(1.-de*u)/ut;
  ccdf=ga1*y*tem/den;
  *p=ccdf;
}

/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter (vth.de). vth>1, 0<de<1
   output
     vv = inverse of conditional cdf
*/
void qcondbb8(double *p0, double *u0, double *vth0, double *de0, double *vv)
{ double p,u,vth,vth1,de,ga,ga1,ut,tem,h,hp,x,y,v;
  double con,diff,mxdif,eps;
  int iter,mxiter;

  vth=*vth0; de=*de0; p=*p0; u=*u0; 
  vth1=1./vth;
  ga=1.-pow(1.-de,vth); ga1=1/ga;
  mxiter=30; eps=1.e-6;
  ut=pow(1.-de*u,vth); x=1.-ut;
  con=log(ga1)+(1.-vth1)*log(1.-x)-log(p);
  // starting guess, need y=.9*x for strong dependence
  y=.5*x; // what is good starting point?
  mxdif=1; iter=0; 
  diff=.1;  
  while(mxdif>eps && iter<mxiter)
  { tem=1.-ga1*x*y;
    h=log(y)+(vth1-1.)*log(tem)+con;
    hp=1./y+(1.-vth1)*ga1*x/tem;
    //if(isnan(h) || isnan(hp) || isnan(h/hp) ) { diff/=-2.; }  // added for th>50
    diff=h/hp;
    y-=diff; iter++;
    while(y<=0. || y>=ga) { diff/=2.; y+=diff; }
#ifdef DEBUG
    printf("%d %f %f\n", iter, diff, y);
#endif
    mxdif=fabs(diff);
  }
#ifdef MAIN
  if(iter>=mxiter) 
  { printf("***did not converge\n");
    printf("p=%f, x=%f, vtheta=%f, de=%f, lasty=%f\n", p,x,vth,de,y);
  }
#else
  if(iter>=mxiter) 
  { Rprintf("** did not converge **\n");
    Rprintf("p=%f, x=%f, vtheta=%f, de=%f, lasty=%f\n", p,x,vth,de,y);
  }
#endif
  v=pow(1.-y,vth1); v=(1.-v)/de;
  *vv=v;
#ifdef DEBUG
  printf("p=%f u=%f v=%f\n", p,u,*vv);
#endif
}

