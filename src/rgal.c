#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random pairs from bivariate Galambos copula 
   gcc -DMAIN -o rgal rgal.c -lm
*/
#ifdef MAIN
main(int argc, char *argv[])
{ double de,u1,u2,lm;
  long seed;
  double s1,s2,s12;
  int isim,nsim;
  void rgal(double, double *, double *);
  setbuf(stdout,NULL);
  scanf("%lf %d %ld", &de,&nsim,&seed);
  while(de>0.)
  { srand(seed);
    lm=pow(2.,-1./de);
    printf("\ndelta=%f lm=%f nsim=%d seed=%d\n", de,lm,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rgal(de,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    scanf("%lf %d %ld", &de,&nsim,&seed);
  }
}

/* input
     de = cpar = copula parameter >0 
   outputs
     u1, u2 = random pair in (0,1) with Galambos copula 
*/
void rgal(double de, double *u1, double *u2)
{ double p,pp;
  void qcondgal(double *p, double *u, double *de, double *);
  void pcondgal(double *v, double *u, double *de, double *);
  p=rand()/2147483648.;
  *u1=rand()/2147483648.;
  qcondgal(&p,u1,&de,u2);
  pcondgal(u2,u1,&de,&pp); // backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %.8f %f %.8f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     cpar = copula parameter >0
   outputs
     uvec, vvec = random sample of size nn from Galambos copula 
*/
void rgal(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondgal(double *,double *,double *, double*);
  double u,v,qq,de;
  int i,n;

  n=*nn; de=*cpar;
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    qq=unif_rand();
    qcondgal(&qq,&u,&de,&v);
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif
 
/* inputs
     v0, u0 = values in (0,1)
     cpar = copula parameter >0
   output
     p = conditional cdf
*/
void pcondgal(double *v0,double *u0,double *cpar, double *p)
{ double u,v,x,y,tem1,tem2,sm,tem,ccdf,de;
  de=*cpar; u=*u0; v=*v0;
  x= -log(u); y= -log(v);
  tem1=pow(x,-de); tem2=pow(y,-de); sm=tem1+tem2; tem=pow(sm,-1./de);
  ccdf=exp(-(x+y-tem));
  ccdf*=1.-pow(1+tem2/tem1,-1.-1./de);
  ccdf/=u;
  *p=ccdf;
}


/* qcondgal Ok for de<.0001, minimal changes from de=.04 and smaller */
/* G(x,y)=exp(-x-y+(x^dn+y^dn)^(1/dn)), de>=0, dn=-de
   conditional P(Y>=y|x)= p
     = exp(x)*G(x,y)*[1-x^(dn-1)*(x^dn+y^dn)^(1/dn-1)]
     = exp(-y+y)*[1-x^(dn-1)*(x^dn+y^dn)^(1/dn-1)]
     = exp(-y+y)* x^(dn-1)* [x^(d+1)-(x^dn+y^dn)^(1/dn-1)]
   take logs, and let y=(x^dn+y^dn)^(1/dn)
   g(y) = log(p) +y -(x^dn+y^dn)^(1/dn)  + 
           - log[1-x^(dn-1)*(x^dn+y^dn)^(1/dn-1)]
   g'(y)=  1 -  y^{dn-1}(x^{dn}+y^{dn})^{1/dn-1}
  + {(1+de)x^{dn-1}y^{dn-1} (x^{dn}+y^{dn})^{1/dn-2}  \over 
    [1 -x^{dn-1}(x^{dn}+y^{dn})^{1/dn-1} ]}.

*/
/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter >0
   output
     vv = inverse of conditional cdf
*/
void qcondgal(double *p0, double *u0, double *cpar, double *vv)
{ double p,u,de,x,tem,g,gp,y,z,con,dn,dif,dn1,xd,yd,rat,den,mxdif,eps;
  double r,rd,r1,r1d,h,hp;
  int iter,mxiter;

  de=*cpar; p=*p0; u=*u0; mxiter=30; eps=1.e-6;
  x=-log(u);
  con=log(p); dn=-de; dn1=1./dn;
  xd=pow(x,dn);
  y=(x*de-con)/(1.+de);
  if(de>2.7) y=x;
  // different equation r=[(-log(u))^de]/[(-log(v))^de]; 
  // or y=x * r^(1/de) when de is large and u is close to 1
  if(xd>1.e6) // empirical
  { con=log(p);
    r=1.;
#ifdef DEBUG
    cat("entering loop solving for r=(y/x)^de\n")
#endif
    mxdif=1; iter=0;
    while(mxdif>eps && iter<mxiter)
    { rd=pow(r,1./de); r1=1.+1./r; r1d=pow(r1,dn1);
      h=con+rd*x - x*r1d- log(1.-r1d/r1);
      hp= rd*x/r/de - x*r1d/r1/de/r/r + (1+1/de)*r1d/r1/r1/r/r/(1-r1d/r1);
      dif=h/hp;
#ifdef DEBUG
      printf("%d %f %f\n", iter, dif, r);
#endif
      r=r-dif; iter=iter+1;
      while(r<=0.) { dif=dif/2.; r=r+dif; }
      mxdif=fabs(dif);
    }
    y=x* pow(r,1./de);
    *vv=exp(-y);
    return;
  }

  //if(isnan(y)) { printf("*** p=%f, u=%f, y=%f\n", p,u,y); }
  mxdif=1; iter=0;
  while(mxdif>eps && iter<mxiter)
  { yd=pow(y,dn); z=pow(xd+yd,dn1);
    tem=pow(z,de+1.);
    //den=1.-tem*xd/x; // this can evaluate as <0 when u=~, x large
    rat=yd/xd; den=1.-pow(1.+rat,-1.+dn1); 
    g=con-z+y-log(den);
    gp=1.-tem*yd/y+(1.+de)*tem*xd*yd/(den*x*y*(xd+yd));
    dif=g/gp;
#ifdef DEBUG
    printf("%d %f %f\n", iter, dif, y);
#endif
    y-=dif; iter++;
#ifdef DEBUG
    if(isnan(y)) 
    { printf("*** p=%f, u=%f, g=%f, gp=%f z=%f den=%f\n\n", p,u,g,gp,z,den); }
#endif
    // above loop with h, hp should be replace this
    if(isnan(dif) || fabs(dif)>1.e10)  // can happen with large de and u near 1
    { y=x; mxdif=1.; }
    else 
    { while(y<=0.) { dif/=2.; y+=dif; }
      mxdif=fabs(dif);
    }
  }
#ifdef MAIN
  if(iter>=mxiter)
  { printf("***did not converge\n");
    printf("p=%f, x=%f, delta=%f, lasty=%f\n", p,x,de,y);
  }
#endif
/*
  if(iter>=mxiter) 
  { Rprintf("** did not converge **\n");
    Rprintf("p=%f, x=%f, delta=%f, lasty=%f\n", p,x,de,y);
  }
*/
  *vv=exp(-y);
}

