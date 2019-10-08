#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random bivariate U(0,1) from bb4 copula 
   gcc -DMAIN -o rbb4 rbb4.c -lm
   gcc -DDEBUG -DMAIN -o rbb4d rbb4.c -lm
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double th,de,u1,u2,tau,be;
  long seed;
  double s1,s2,s12;
  int isim,nsim;
  void rbb4(double, double, double *, double *);
  setbuf(stdout,NULL);
  scanf("%lf %lf %d %ld", &th,&de,&nsim,&seed);
  // th>0, de>=0
  while(th>0.)
  { srand(seed);
    printf("\ntheta=%f delta=%f nsim=%d seed=%ld\n", th,de,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rbb4(th,de,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    printf("\n============================================================\n");
    scanf("%lf %lf %d %ld", &th,&de,&nsim,&seed);
  }
}

/* inputs
     th = theta parameter >0 
     de = delta parameter >0 
   outputs
     u1, u2 = random pair in (0,1) with BB4 copula 
*/
void rbb4(double th, double de, double *u1, double *u2)
{ double p,pp;
  void qcondbb4(double *p, double *u, double *, double *, double *);
  void pcondbb4(double *v, double *u, double *, double *, double *);
  *u1=rand()/2147483648.;
  p=rand()/2147483648.;
  qcondbb4(&p,u1,&th,&de,u2);  // return u2
  pcondbb4(u2,u1,&th,&de,&pp); // return pp, backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     cpar = copula parameter (th.de). th>0, de>0
   outputs
     uvec, vvec = random sample of size nn from BB4 copula 
*/
void rbb4(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondbb4(double *,double *,double *, double *, double*);
  double u,v,p,th,de;
  int i,n;
  n=*nn; th=cpar[0]; de=cpar[1];
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    p=unif_rand();
    qcondbb4(&p,&u,&th,&de,&v); // return v
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* inputs
     v0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>0, de>0
   output
     p = conditional cdf
*/
void pcondbb4(double *v0,double *u0,double *th0, double *de0, double *p)
{ double th,de,u,v,ut,vt,x,y,xy,tem,xtem,ccdf;
  th=*th0; de=*de0; u=*u0; v=*v0; 
  ut=pow(u,-th)-1.; vt=pow(v,-th)-1.;
  x=pow(ut,-de); y=pow(vt,-de); xy=x+y;
  tem=pow(xy,-1./de);
  ccdf=ut+vt+1.-tem;
  xtem=ut/x-tem/xy;
  ccdf=pow(ccdf,-1./th-1.)*xtem*x/ut*(1.+ut)/u;
  *p=ccdf;
}

/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>0, de>0
   output
     vv = inverse of conditional cdf
*/
void qcondbb4(double *p0, double *u0, double *th0, double *de0, double *vv)
{ double p,u,th,de,th1,de1,ut,vt,h,hp,x,y,xtem,ytem,xy,tem,sm;
  double con,diff,mxdif,eps;
  int iter,mxiter;

  th=*th0; de=*de0; p=*p0; u=*u0; 
  mxiter=30;
  //mxiter=20; 
  eps=1.e-6;
  de1=1./de; th1=1./th;
  ut=pow(u,-th)-1.; 
  x=pow(ut,-de); 
  con=(de1+1.)*log(x)-(th+1.)*log(u)-log(p);
  mxdif=1.; iter=0; 
  diff=1.; 
  y=.5*x; // what is good starting point?
  // use y=x when dependence is larger?
  if(de>1.4 && th>1.4) y=.9*x;
  // also boundary cases is de<0.1, th<0.1
  while(mxdif>eps && iter<mxiter)
  { xy=x+y; vt=pow(y,-de1);
    tem=pow(xy,-de1);
    sm=ut+vt+1.-tem;
    xtem=ut/x-tem/xy;
    ytem=vt/y-tem/xy;
    h=-(th1+1.)*log(sm)+log(xtem)+con;
    hp=(th1+1.)*ytem/sm/de + (1.+de1)*tem/xtem/xy/xy;
    diff=h/hp;
    y=y-diff;
    iter++;
    while(y<=0.) { diff/=2.; y+=diff; }
#ifdef DEBUG
    printf("%d %f %f\n", iter, diff, y);
#endif
    mxdif=fabs(diff);
  }
#ifdef MAIN
  if(iter>=mxiter) 
  { printf("***did not converge\n");
    printf("p=%f, x=%f, theta=%f, delta=%f, lasty=%f\n", p,x,th,de,y);
  }
#else
  if(iter>=mxiter) 
  { Rprintf("** did not converge **\n");
    Rprintf("p=%f, x=%f, theta=%f, delta=%f, lasty=%f\n", p,x,th,de,y);
  }
#endif
  *vv=pow(1.+pow(y,-de1),-th1);
#ifdef DEBUG
  printf("p=%f u=%f v=%f\n", p,u,*vv);
#endif
}

