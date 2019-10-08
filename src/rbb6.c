#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random bivariate U(0,1) from BB6 copula 
   gcc -DMAIN -o rbb6 rbb6.c -lm
   gcc -DDEBUG -DMAIN -o rbb6d rbb6.c -lm
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double th,de,u1,u2,tau,be;
  long seed;
  double s1,s2,s12;
  int isim,nsim;
  void rbb6(double, double, double *, double *);
  setbuf(stdout,NULL);
  scanf("%lf %lf %d %ld", &th,&de,&nsim,&seed);
  // th>1, de>=1
  while(th>1.)
  { srand(seed);
    printf("\ntheta=%f delta=%f nsim=%d seed=%ld\n", th,de,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rbb6(th,de,&u1,&u2);
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
     u1, u2 = random pair in (0,1) with BB6 copula 
*/
void rbb6(double th, double de, double *u1, double *u2)
{ double p,pp;
  void qcondbb6(double *p, double *u, double *, double *, double *);
  void pcondbb6(double *v, double *u, double *, double *, double *);
  *u1=rand()/2147483648.;
  p=rand()/2147483648.;
  qcondbb6(&p,u1,&th,&de,u2);  // return u2
  pcondbb6(u2,u1,&th,&de,&pp); // return pp, backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     cpar = copula parameter (th.de). th>1, de>1
   outputs
     uvec, vvec = random sample of size nn from BB6 copula 
*/
void rbb6(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondbb6(double *,double *,double *, double *, double*);
  double u,v,p,th,de;
  int i,n;
  n=*nn; th=cpar[0]; de=cpar[1];
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    p=unif_rand();
    qcondbb6(&p,&u,&th,&de,&v); // return v
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* inputs
     v0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>1, de>1
   output
     p = conditional cdf
*/
void pcondbb6(double *v0,double *u0,double *th0, double *de0, double *p)
{ double th,de,u,v,ubar,vbar,x,y,zu,xd,yd,sm,tem,w,ccdf;
  th=*th0; de=*de0; u=*u0; v=*v0; 
  ubar=1.-u; vbar=1.-v;
  zu=1.-pow(ubar,th);
  x=-log(zu); y=-log(1.-pow(vbar,th));
  xd=pow(x,de); yd=pow(y,de); sm=xd+yd; tem=pow(sm,1./de);
  w=exp(-tem); 
  ccdf=pow((1.-w)/(1.-zu),1./th-1.) *(w/zu) *(tem/sm) * (xd/x);
  *p=ccdf;
}

/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>1, de>1
   output
     vv = inverse of conditional cdf
*/
void qcondbb6(double *p0, double *u0, double *th0, double *de0, double *vv)
{ double p,u,th,de,th1,de1,ubar,zu,h,hp,x,y,xd,yd,con,diff,mxdif,eps;
  double sm,tem,w;
  int iter,mxiter;

  th=*th0; de=*de0; p=*p0; u=*u0; 
  mxiter=30; eps=1.e-6;
  de1=1./de; th1=1./th;
  ubar=1.-u; 
  zu=1.-pow(ubar,th);
  x=-log(zu); 
  xd=pow(x,de); 
  con=(de-1)*log(x)-(th1-1)*log(1-zu)+x-log(p);
  mxdif=1.; iter=0; 
  diff=1.; 
  y=.5*x; // what is good starting point?
  while(mxdif>eps && iter<mxiter)
  { yd=pow(y,de); sm=xd+yd; tem=pow(sm,1./de);
    w=exp(-tem); 
    h=(th1-1.)*log(1.-w)-tem+(de1-1.)*log(sm)+con;
    hp=(th1-1.)*w*tem*yd/(1.-w)/sm/y - tem*yd/sm/y +(1.-de)*yd/y/sm;
    diff=h/hp;
    y=y-diff;
    //if(isnan(h) || isnan(hp) || isnan(h/hp) ) { diff/=-2.; }  
    //else diff=h/hp;
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
  *vv=1.-pow(1.-exp(-y),th1);
#ifdef DEBUG
  printf("p=%f u=%f v=%f\n", p,u,*vv);
#endif
}

