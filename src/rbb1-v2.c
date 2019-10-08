#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random bivariate pairs from BB1 copula */
// also pcondbb1, qcondbb1 with pointers, could be called from R
#ifdef MAINB
main(int argc, char *argv[])
{ double th,de,u1,u2,be,lml,lmu,tau;
  void rbb1(double, double, double *, double *);
  int nsim,isim;
  long seed;
  double s1,s2,s12;
  
  // th>0, de>1
  setbuf(stdout,NULL);
  scanf("%lf %lf %d %ld", &th,&de,&nsim,&seed);
  while(th>0)
  { srand(seed);
    lml=pow(2.,-1./(de*th));
    lmu=2.-pow(2,1./de);
    be=1+(2.-lmu)*(pow(2.,th)-1);
    be=pow(be,-1./th);
    be=4*be-1.;
    tau=1-2/(de*(th+2.)); 
    printf("\nth=%f de=%f, tau=%f beta=%f, lmL=%f lmU=%f\n", th,de,tau,be,lml,lmu);
    printf("nsim=%d seed=%d\n", nsim,seed); 
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rbb1(th,de,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    printf("============================================================\n");
    scanf("%lf %lf %d %ld", &th,&de,&nsim,&seed);
  }
}

double urand()
{ return(rand()/2147483648.); }
#endif

#ifdef NOTR
/* inputs
     th = theta parameter >0 
     de = delta parameter >1 
   outputs
     u1, u2 = random pair in (0,1) with BB1 copula 
*/
void rbb1(double th, double de, double *u1, double *u2)
{ double p,pp;
  double cpar[2];
  void qcondbb1(double *p, double *u1, double *cpar, double *);
  void pcondbb1(double *p, double *u1, double *cpar, double *);
  double urand();
  *u1=urand();
  p=urand();
  cpar[0]=th; cpar[1]=de;
  qcondbb1(&p,u1,cpar,u2);  // return u2
  pcondbb1(u2,u1,cpar,&pp); // return pp, backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     cpar = copula parameter (th.de). th>0, de>1
   outputs
     uvec, vvec = random sample of size nn from BB1 copula 
*/
void rbb1(int *nn, double *cpar, double *uvec, double *vvec)
{ 
  void qcondbb1(double *p, double *u1, double *cpar, double *);
  double u,v,p,th,de;
  int i,n;
  n=*nn; th=cpar[0]; de=cpar[1];
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    p=unif_rand();
    qcondbb1(&p,&u,cpar,&v); // return v
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* inputs
     v0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>0, de>1
   output
     p = conditional cdf
*/
void pcondbb1(double *v0, double *u0, double *cpar, double *p)
{ double de1,th1,ut,vt,x,y,sm,smd,tem,ccdf;
  double v,u,th,de;
  //th=*th0; de=*de0; u=*u0; v=*v0;
  th=cpar[0]; de=cpar[1]; u=*u0; v=*v0;
  de1=1./de; th1=1./th;
  ut=pow(u,-th)-1.; vt=pow(v,-th)-1.;
  x=pow(ut,de); y=pow(vt,de);
  sm=x+y; smd=pow(sm,de1);
  tem=pow(1.+smd,-th1-1.);
  ccdf=tem*smd*x*(ut+1.)/sm/ut/u;
  *p=ccdf;
}


// code using transformed variables, without bisection method
/* $x=(u^{-\th}-1)^{\de}$ and $y=(v^{-\th}-1)^{\de}$
   $u=(1+x^{1/\de})^{-1/\th}$ and $v=(1+y^{1/\de})^{-1/\th}$.
   $C(u,v)=G(x,y)=[1+(x+y)^{1/\de}]^{-1/\th}$
   C_{2|1}(v|u)= \{ 1+ (x+y)^{1/\de} \}^{-1/\th-1} 
      (x+y)^{1/\de-1} x^{1-1/\de} u^{-\th-1} $
*/
/* G_{2|1}(y|x)=p for y given x,p */
/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>0, de>1
   output
     vv = inverse of conditional cdf
*/
void qcondbb1(double *p0, double *u0, double *cpar, double *vv)
{ double de1,th1,ut,x,y,den,pden,diff,sm,smd,eps,v,r;
  double G21,gpdf,tem,thr;
  int iter,mxiter;
  double p,u,th,de;
  void qcondgum(double *p0, double *u0, double *cpar, double *vv);
  //void qcondmtcj(double *p0, double *u0, double *cpar, double *vv);

  //th=*th0; de=*de0; u=*u0; p=*p0;
  th=cpar[0]; de=cpar[1]; u=*u0; p=*p0;
  mxiter=30; eps=1.e-6;
  de1=1./de; th1=1./th;
  ut=pow(u,-th)-1.; x=pow(ut,de); 
  den=pow(1.+ut,-th1-1)*ut/x; // density of univariate margin
  pden=p*den; // term 1/(th*de) cancels
  y=pow(pden,-1./(th1*de1+1))-x; // starting guess for y
  // in cases where x and ystart large pcond(qcond(q))-q might be .001
  if(x<1.e-5) // empirically based
  { y=x*(pow(p,-de/(de-1))-1.);  // take x=~ pow(th*(1-u),de) with u near 1?
    // need bound to prevent y large 
    if(y>1.e-5) y=1.e-5;
#ifdef DEBUG
    printf("\nsmall x case, x=%f, y=%f\n",x,y);
#endif
  }
  if(x>1.e5) // empirically based
  { r=pow(p,-de*th/(1.+de*th))-1.;
    y=r*x;
#ifdef DEBUG
    printf("\nlarge x case, x=%f y=%f\n",x,y);
#endif
    eps*=.0001; // good enough
    // some cases left where pcond(qcond(q))-q might be .001
  }
#ifdef DEBUG
  printf("\np=%f, u=%f\n", p,u);
  printf(" x=%f and starting y=%f\n", x, y);
#endif
  // previous tolerance of 0.01
  // if th<0.1 or de<1.1 use boundary cases of Gumbel and MTCJ as starting points
  // *** biggest problem for de>5 and th small
  // v=(y^(1/de)+1)^(-1/th) so y=(v^(-th)-1)^de
  if(de<=1.1)   // MTCJ boundary
  { //qcondmtcj(&p,&u,&th,&v); y=pow(pow(v,-th)-1.,de); 
    thr=-th/(1.+th); tem=pow(p,thr)-1.;
    tem=tem*pow(u,-th)+1.; y=pow(tem-1.,de);
  }
  // Gumbel boundary if th<.2;
  else if(th<0.2) { qcondgum(&p,&u,&de,&v); y=pow(pow(v,-th)-1.,de); }
  
  diff=1; iter=0;
  while(fabs(diff/y)> eps && iter<mxiter)
  { sm=x+y; smd=pow(sm,de1);
    G21=pow(1.+smd,-th1-1.) * smd/sm; 
    gpdf=-G21; 
    gpdf=gpdf/(1.+smd)/sm/de/th;
    gpdf=gpdf*(th*(de-1.)+(th*de+1.)*smd);
    iter=iter+1;
    diff=(G21-pden)/gpdf;
    y=y-diff;
#ifdef DEBUG
    printf("%d %f %f\n",  iter, y, diff);
#endif
    while(y<=0.) { diff=diff/2.; y=y+diff; }
  }
  v=pow(pow(y,de1)+1.,-th1);
#ifdef MAIN
  if(iter>=mxiter) 
  { printf("** did not converge **\n");
    printf("(p,u,th,de)= %f %f %.2f %.2f lasty=%f v=%f\n", p,u,th,de,y,v);
  }
#endif
#ifdef DEBUG
  printf("v= %f  #iter=%d\n",v,iter);
#endif
  *vv=v;
}

