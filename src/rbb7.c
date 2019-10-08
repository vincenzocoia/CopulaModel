#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random bivariate uniforms from BB7 copula */
// gcc -DMAIN -o rbb7 rbb7.c -lm
// gcc -DMAIN -DDEBUG -o rbb7d rbb7.c -lm
// also pcondbb7, qcondbb7 with pointers, could be called from R
#ifdef MAIN
main(int argc, char *argv[])
{ double th,de,u1,u2,be,lml,lmu;
  void rbb7(double, double, double *, double *);
  int nsim,isim;
  long seed;
  double s1,s2,s12;
  
  // th>1, de>0 
  setbuf(stdout,NULL);
  scanf("%lf %lf %d %ld", &th,&de,&nsim,&seed);
  while(th>1.)
  { srand(seed);
    lmu=2-pow(2.,1./th); lml=pow(2,-1./de);
    be=2*pow(1-pow(.5,th),-de)-1;
    be=1.-pow(be,-1./de); be=1-pow(be,1./th);
    be=4*be-1.;
    printf("\nth=%f de=%f, beta=%f, lmL=%f lmU=%f\n", th,de,be,lml,lmu);
    printf("nsim=%d seed=%d\n", nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rbb7(th,de,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    printf("============================================================\n");
    scanf("%lf %lf %d %ld", &th,&de,&nsim,&seed);
  }
}

/* inputs
     th = theta parameter >1 
     de = delta parameter >0 
   outputs
     u1, u2 = random pair in (0,1) with BB7 copula 
*/
void rbb7(double th, double de, double *u1, double *u2)
{ double p,pp;
  void qcondbb7(double *p, double *u1, double *th, double *de, double *);
  void pcondbb7(double *p, double *u1, double *th, double *de, double *);
  *u1=rand()/2147483648.;
  p=rand()/2147483648.;
  qcondbb7(&p,u1,&th,&de,u2);  // return u2
  pcondbb7(u2,u1,&th,&de,&pp); // return pp, backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     cpar = copula parameter (th.de). th>1, de>0
   outputs
     uvec, vvec = random sample of size nn from BB7 copula 
*/
void rbb7(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondbb7(double*,double*,double*,double*, double*);
  double u,v,p,th,de;
  int i,n;
  n=*nn; th=cpar[0]; de=cpar[1];
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    p=unif_rand();
    qcondbb7(&p,&u,&th,&de,&v); // return v
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* inputs
     v0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>1, de>0
   output
     p = conditional cdf
*/
void pcondbb7(double *v0, double *u0, double *th0, double *de0, double *p)
{ double de1,th1,ut,vt,x,y,sm,smd,tem,ccdf;
  double v,u,th,de;
  th=*th0; de=*de0; u=*u0; v=*v0;
  de1=1./de; th1=1./th;
  ut=1.-pow(1.-u,th); vt=1.-pow(1.-v,th); 
  x=pow(ut,-de)-1.; y=pow(vt,-de)-1.;
  sm=x+y+1; smd=pow(sm,-de1);
  tem=pow(1.-smd,th1-1.);
  ccdf=tem*smd*(x+1.)*(1.-ut)/sm/ut/(1.-u);
  *p=ccdf;
}

// code using transformed variables, without bisection method
/* $x=(1-[1-u]^{\th})^{-\de}-1$ and $y=(1-[1-v]^{\th})^{-\de}-1$
   $u=1-[1-(x+1)^{-1/\de}]^{1/\th}$ and $v=1-[1-(y+1)^{-1/\de}]^{1/\th}$
   $C(u,v)=G(x,y)=1-[1-(x+y+1)^{-1/\de}]^{1/\th}$
   C_{2|1}(v|u)= \{ 1- (x+y+1)^{-1/\de} \}^{1/\th-1} 
     (x+y+1)^{-1/\de-1}  (x+1)^{1+1/\de} [1-u]^{\th-1} $
*/
/* G_{2|1}(y|x)=p for y given x,p */
/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>1, de>0
   output
     vv = inverse of conditional cdf
*/
void qcondbb7(double *p0, double *u0, double *th0, double *de0, double *vv)
{ double de1,th1,ut,x,y,den,pden,rhs,diff,sm,smd,eps,v,epsx,epsy,tem;
  double G21,gpdf;
  int iter,mxiter;
  double p,u,th,de;

  th=*th0; de=*de0; u=*u0; p=*p0;
  mxiter=30; eps=1.e-6;
  de1=1./de; th1=1./th;
  ut=1.-pow(1.-u,th); x=pow(ut,-de)-1.; 
  // if ut rounds to 1, x rounds to 0 and den becomes oo
  // for random inputs, this might occur once in 1.e6
  if(x<=0.) 
  { v=.999999; if(u>v) v=u;
    *vv=v; return;
  }
  den=pow(1.-ut,th1-1)*ut/(x+1.); // density of univariate margin
  pden=p*den; // term 1/(th*de) cancels
  // starting guess for y 
  rhs=pden*pow(1.-pow(2*x+1.,-de1),1.-th1);
  y=pow(rhs,-de/(de+1))-1-x;
  if(y<=0. || isnan(y)) y=0.1; // need better starting point if x<eps
  if(x<1.e-5) // empirically based
  { epsx=de*(1.-ut);   //x
    tem=p*(1-(1+de1)*epsx);
    tem=pow(tem,-th/(th-1.))-1.;
    epsy=tem*epsx;
    if(epsy>1.e-5) epsy=1.e-5;
    y=epsy;
  }
  
  diff=1; iter=0;
  while(fabs(diff/y)> eps && iter<mxiter)
  { sm=x+y+1; smd=pow(sm,-de1);
    G21=pow(1.-smd,th1-1.) * smd/sm; 
    gpdf=-G21; 
    gpdf=gpdf/(1.-smd)/sm/de/th;
    gpdf=gpdf*(th*(de+1.)-(th*de+1.)*smd);
    iter=iter+1;
    diff=(G21-pden)/gpdf;
    y=y-diff;
    while(y<=0.) { diff=diff/2.; y=y+diff; }
#ifdef DEBUG
    printf("%d %f %f\n",  iter, y, diff);
#endif
  }
  v=1.-pow(1.-pow(y+1.,-de1),th1);
#ifdef MAIN
  if(iter>=mxiter && y>0.001) 
  { printf("** did not converge **\n");
    printf("(p,u,th,de)= %f %f %.2f %.2f lasty=%f v=%f\n", p,u,th,de,y,v);
  }
#endif
#ifdef DEBUG
  printf("v= %f  #iter=%d\n",v,iter);
#endif
  *vv=v;
}

