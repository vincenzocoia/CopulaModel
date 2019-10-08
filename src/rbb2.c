#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random bivariate uniforms from BB2 copula */
// gcc -DMAIN -o rbb2 rbb2pkg.c -lm
// gcc -DDEBUG -DMAIN -o rbb2d rbb2.c -lm
// also pcondbb2, qcondbb2 with pointers, could be called from R
#ifdef MAIN
main(int argc, char *argv[])
{ double th,de,u1,u2,be;
  void rbb2(double, double, double *, double *);
  int isim,nsim;
  long seed;
  double s1,s2,s12;
  
  /* th>0, de>0 */
  setbuf(stdout,NULL);
  scanf("%lf %lf %d %ld", &th,&de,&nsim,&seed);
  while(th>0)
  { srand(seed);
    be=2.*exp(de*(pow(2.,th)-1))-1;
    be=1+log(be)/de;
    be=pow(be,-1./th); be=4*be-1.;
    printf("\nth=%f de=%f, beta=%f nsim=%d seed=%ld\n", th,de,be,nsim,seed);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rbb2(th,de,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;;
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
   output
     u1, u2 = random pair in (0,1) with BB2 copula 
*/
void rbb2(double th, double de, double *u1, double *u2)
{ double p,pp;
  void qcondbb2(double *p, double *u, double *, double *, double *);
  void pcondbb2(double *v, double *u, double *, double *, double *);
  *u1=rand()/2147483648.;
  p=rand()/2147483648.;
  qcondbb2(&p,u1,&th,&de,u2);  // return u2
  pcondbb2(u2,u1,&th,&de,&pp); // return pp, backcheck
  if(fabs(p-pp)>1.e-3) printf("*** 3rddecpl %f %f %f %f\n", *u1,p,*u2,pp);
}
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
/* inputs
     nn = simulation sample size
     param = copula parameter (th.de). th>0, de>0
   outputs
     uvec, vvec = random sample of size nn from BB2 copula 
*/
void rbb2(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondbb2(double *,double *,double *, double *, double*);
  double u,v,p,th,de;
  int i,n;
  n=*nn; th=cpar[0]; de=cpar[1];
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    p=unif_rand();
    qcondbb2(&p,&u,&th,&de,&v); // return v
    uvec[i]=u; vvec[i]=v;
  }
  PutRNGstate();
}
#endif

/* inputs (later change to BB1 format)
     v0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>0, de>0
   output
     p = conditional cdf
*/
void pcondbb2(double *v0,double *u0,double *th0, double *de0, double *p)
{ double u,v,th,de,th1,de1,ut,vt,x,y,sm,smd,tem,ccdf,lr,r;
  th=*th0; de=*de0; u=*u0; v=*v0;
  de1=1./de; th1=1./th;
  ut=pow(u,-th)-1.; vt=pow(v,-th)-1.;
  x=exp(ut*de)-1.; y=exp(vt*de)-1.;
  if(isinf(x) || isinf(y))
  { lr=de*(vt-ut); r=exp(lr);
    tem=pow(1.+ut+de1*log(1.+r),-th1-1);
    ccdf=tem*(ut+1.)/u/(1.+r);
  }
  else
  { sm=x+y+1.; smd=de1*log(sm);
    tem=pow(1.+smd,-th1-1.);
    ccdf=tem*(x+1.)*(ut+1.)/sm/u;
  }
  *p=ccdf;
}


// code using transformed variables, without bisection method
/* $x=\exp\{\de(u^{-\th}-1)\}-1$ $y=\exp\{\de(v^{-\th}-1)\}-1$
   $u=[1+\de^{-1}\log (x+1)]^{-1/\th}$ and $v=[1+\de^{-1}\log (y+1)]^{-1/\th}$
   $C(u,v)=G(x,y)=[1+\de^{-1}\log(x+y+1)]^{-1/\th}$
   $C_{2|1}(v|u)= \{ 1+ \de^{-1}\log(x+y+1) \}^{-1/\th-1} 
   (x+y+1)^{-1} (x+1) u^{-\th-1} $
*/
/* G_{2|1}(y|x)=p for y given x,p */
/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>0, de>0
   output
     vv = inverse of conditional cdf
*/
void qcondbb2(double *p0, double *u0, double *th0, double *de0, double *vv)
{ double p,u,th,de,de1,th1,ut,x,y,den,pden,rhs,diff,sm,smd,eps,v;
  double G21,gpdf,r,g,gp,lx,tha;
  int iter,mxiter;

  th=*th0; de=*de0; p=*p0; u=*u0; 
  mxiter=30; eps=1.e-6;
  de1=1./de; th1=1./th;
  ut=pow(u,-th)-1.; 
  x=exp(ut*de)-1; // log(x+1)=ut*de 
  den=pow(1.+ut,-th1-1.)/(x+1.); // density of univariate margin
  // empirical for tolerance or den here; seems to prevent getting y=inf
  if(isinf(x) || den<=1.e-100)
  { // solve a different equation based on y=rx as x->oo
    // correction tha=-th/(1+th)
    //  g(r)= 1+log(1+r) /lx -p^tha (1+r)^tha; 
    //    g'(r)=1/(1+r)/lx -p^tha * tha *(1+r)^(tha-1);
    // v=(1+de1*(log(r)+lx))^{-th1}
#ifdef DEBUG
    printf("\n** infinite x, p=%f u=%f th=%f de=%f ut=%f\n", p,u,th,de,ut);
#endif
    lx=ut*de; r=1.; tha=-th/(1.+th);
    diff=1; iter=0;
    while(fabs(diff)>eps  && iter<mxiter)
    { sm=1.+r; smd=de1*log(sm); rhs=pow(p*sm,tha);
      g=1.+log(sm)/lx-rhs;
      gp=1./sm/lx-tha*rhs/sm;;
      iter=iter+1;
      diff=g/gp;
      r=r-diff;
      while(r<=0.) { diff=diff/2.; r=r+diff; }
#ifdef DEBUG
      printf("%d r=%f %f\n",  iter, r, diff);
#endif
    }
#ifdef DEBUG
    if(iter>=mxiter) 
    { printf("** did not converge **\n");
      printf("(p,u,th,de)= %f %f %.2f %.2f lastr=%f\n", p,u,th,de,r);
    }
#endif
    v=pow(1.+de1*(log(r)+lx),-th1);
#ifdef DEBUG
    printf("** infinite x, p=%f u=%f th=%f de=%f ut=%f, r=%f, v=%f\n",
       p,u,th,de,ut,r,v);
#endif
    *vv=v;
    return;
  }
  pden=p*den; // term 1/(th*de) cancels
  // starting guess for y, 
  rhs=pow(1.+log(2*x+1),-1.-th1)/pden;
  y=rhs-1-x;  
  if(y<=0.) y=0.1;
  if(y>1000.) y=1000.;
#ifdef DEBUG
  printf("\np=%f, u=%f\n", p,u);
  printf(" x=%f, den=%f, rhs=%f and starting y=%f\n", x,den,rhs,y);
#endif
  
  diff=1; iter=0;
  while((fabs(diff/y)> eps) && iter<mxiter)
  { sm=x+y+1; smd=de1*log(sm);
    G21=pow(1.+smd,-th1-1.)/sm; 
    gpdf=-G21; 
    gpdf=gpdf/(1.+smd)/sm/de/th; 
    gpdf=gpdf*(1.+th+th*de*(1.+smd));
    iter=iter+1;
    diff=(G21-pden)/gpdf;
    y=y-diff;
    while(y<=0.) { diff=diff/2.; y=y+diff; }
#ifdef DEBUG
    printf("%d %f %f\n",  iter, y, diff);
#endif
  }
  v=pow(1.+de1*log(y+1.),-th1);
#ifdef MAIN
  if(iter>=mxiter || v<=1.e-10) 
  { printf("** did not converge **\n");
    printf("(p,u,th,de)= %f %f %.2f %.2f lasty=%f v=%f\n", p,u,th,de,y,v);
  }
#endif
#ifdef DEBUG
  printf("v= %f  #iter=%d\n",v,iter);
#endif
  *vv=v;
}

