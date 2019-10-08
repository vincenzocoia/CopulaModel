#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random bivariate U(0,1) from BB9 copula 
   gcc -DMAIN -o rbb9 rbb9.c -lm
   gcc -DDEBUG -DMAIN -o rbb9d rbb9.c -lm
   change from alp -> 1/ga to get increasing in concordance
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double th,ga,u1,u2,tau,be,ga1;
  long seed;
  double s1,s2,s12;
  int isim,nsim;
  void rbb9(double, double, double *, double *);
  setbuf(stdout,NULL);
  scanf("%lf %lf %d %ld", &th,&ga,&nsim,&seed);
  // th>1, ga>=0
  while(th>1.)
  { srand(seed);
    ga1=1./ga;
    be=2.*pow(ga1+log(2),th)-pow(ga1,th);
    be=pow(be,1./th)-ga1; be=exp(-be);
    be=4.*be-1.;
    printf("\ntheta=%f gamma=%f be=%f nsim=%d seed=%d\n", th,ga,be,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rbb9(th,ga,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    printf("\n============================================================\n");
    scanf("%lf %lf %d %ld", &th,&ga,&nsim,&seed);
  }
}

/* inputs
     th = theta parameter >1 
     ga = gamma parameter >0 
   outputs
     u1, u2 = random pair in (0,1) with BB9 copula 
*/
void rbb9(double th, double ga, double *u1, double *u2)
{ double p,pp;
  void qcondbb9(double *p, double *u, double *, double *, double *);
  void pcondbb9(double *v, double *u, double *, double *, double *);
  *u1=rand()/2147483648.;
  p=rand()/2147483648.;
  qcondbb9(&p,u1,&th,&ga,u2);  // return u2
  pcondbb9(u2,u1,&th,&ga,&pp); // return pp, backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     cpar = copula parameter (th.ga). th>1, ga>0
   outputs
     uvec, vvec = random sample of size nn from BB9 copula 
*/
void rbb9(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondbb9(double *,double *,double *, double *, double*);
  double u,v,p,th,ga;
  int i,n;
  n=*nn; th=cpar[0]; ga=cpar[1];
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    p=unif_rand();
    qcondbb9(&p,&u,&th,&ga,&v); // return v
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* inputs
     v0, u0 = values in (0,1)
     cpar = copula parameter (th.ga). th>1, ga>0
   output
     p = conditional cdf
*/
void pcondbb9(double *v0,double *u0,double *th0, double *ga0, double *p)
{ double u,v,x,y,temx,temy,sm,smt,ccdf,th,ga,ga1;
  th=*th0; ga=*ga0; u=*u0; v=*v0; ga1=1./ga;
  x= ga1-log(u); y= ga1-log(v);
  temx=pow(x,th); temy=pow(y,th); sm=temx+temy-pow(ga1,th); smt=pow(sm,1./th);
  ccdf=exp(-smt+ga1);
  ccdf=ccdf*smt*temx/sm/x/u;
  *p=ccdf;
}

// x=ga1-log(u), y=ga1-log(v)
/* G(x,y)=exp(-(x^th+y^th-ga1^th)^(1/th)), th>=1
   conditional P(Y>=y|x)=
     exp(x-ga1)*G(x,y)* (x^th+y^th-ga1^th)^(1/th-1) * x^(th-1) = p
   take logs, and let z=(x^th+y^th-ga1^th)^(1/th) >= x
   g(z) = z + (th-1)*log(z) + log(p) - x -(th-1)*log(x) =0
   g'(z) = 1 + (th-1)/z
   solve for z with NR, then solve for y=y(z,x,p,th,ga)
   y = (z^th - x^th + ga1^th)^(1/th)
*/
/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter (th.ga). th>1, ga>0
   output
     vv = inverse of conditional cdf
*/
void qcondbb9(double *p0, double *u0, double *th0, double *ga0, double *vv)
{ double p,u,th,ga,ga1,z,h,hp,x,y,con,th1,diff,mxdif,eps;
  int iter,mxiter;

  th=*th0; ga=*ga0; p=*p0; u=*u0; ga1=1./ga;
  mxiter=20; eps=1.e-6;
  x=ga1-log(u); th1=th-1.;
  con=log(p)-x-th1*log(x); 
  // starting guess
  z=pow(2.*pow(x,th)-pow(ga1,th), 1./th); 
  mxdif=1; iter=0; 
  diff=.1;  // needed in case first step leads to NaN??
  while(mxdif>eps && iter<mxiter)
  { h=z+th1*log(z)+con;
    hp=1.+th1/z;
    if(isnan(h) || isnan(hp) || isnan(h/hp) ) { diff/=-2.; }  // added for th>50
    else diff=h/hp;
    z-=diff; iter++;
    while(z<=x) { diff/=2.; z+=diff; }
#ifdef DEBUG
    printf("%d %f %f\n", iter, diff, z);
#endif
    mxdif=fabs(diff);
  }
#ifdef MAIN
  if(iter>=mxiter) 
  { printf("***did not converge\n");
    printf("p=%f, x=%f, theta=%f, ga=%f, lastz=%f\n", p,x,th,ga,z);
  }
#else
  if(iter>=mxiter) 
  { Rprintf("** did not converge **\n");
    Rprintf("p=%f, x=%f, theta=%f, ga=%f, lastz=%f\n", p,x,th,ga,z);
  }
#endif
  y=pow(pow(z,th)-pow(x,th)+pow(ga1,th),1./th);
  *vv=exp(-y+ga1);
#ifdef DEBUG
  printf("p=%f u=%f v=%f\n", p,u,*vv);
#endif
}

