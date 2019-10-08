#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* random bivariate uniforms from BB3 copula */
// gcc -DMAIN -o rbb3 rbb3.c -lm
// gcc -DDEBUG -DMAIN -o rbb3d rbb3.c -lm
#ifdef MAIN
main(int argc, char *argv[])
{ double th,de,u1,u2,be;
  void rbb3(double, double, double *, double *);
  int isim,nsim;
  long seed;
  double s1,s2,s12;
  
  setbuf(stdout,NULL);
  scanf("%lf %lf %d %ld", &th,&de,&nsim,&seed);
  while(th>1) // th>1, de>0
  { srand(seed);
    be=2.*exp(de*pow(log(2.),th))-1.;
    be=log(be)/de; be=exp(-pow(be,1./th));
    be=4.*be-1.;
    printf("\nth=%f, de=%f, beta=%f nsim=%d\n", th,de,be,nsim);
    for(isim=1,s1=0.,s2=0.,s12=0.;isim<=nsim;isim++)
    { rbb3(th,de,&u1,&u2);
      s1+=u1; s2+=u2; s12+=u1*u2;;
    }
    s1/=nsim; s2/=nsim; s12/=nsim; s12-=s1*s2;
    printf("average u = %f, average v = %f\n", s1,s2);
    printf("correlation = %f\n", 12*s12);
    printf("\n============================================================\n");
    scanf("%lf %lf %d %ld", &th,&de,&nsim,&seed);
  }
}

/* input
     th = theta parameter >1 
     de = delta parameter >0 
   output
     u1, u2 = random pair in (0,1) with BB3 copula 
*/
void rbb3(double th, double de, double *u1, double *u2)
{ double p,pp;
  void qcondbb3(double *p, double *u, double *, double *, double *);
  void pcondbb3(double *v, double *u, double *, double *, double *);
  *u1=rand()/2147483648.;
  p=rand()/2147483648.;
  qcondbb3(&p,u1,&th,&de,u2);  // return u2
  pcondbb3(u2,u1,&th,&de,&pp); // return pp, backcheck
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
     uvec, vvec = random sample of size nn from BB3 copula 
*/
void rbb3(int *nn, double *cpar, double *uvec, double *vvec)
{ void qcondbb3(double *,double *,double *, double *, double*);
  double u,v,p,th,de;
  int i,n;
  n=*nn; th=cpar[0]; de=cpar[1];
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { u=unif_rand();
    p=unif_rand();
    qcondbb3(&p,&u,&th,&de,&v); // return v
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
void pcondbb3(double *v0,double *u0,double *th0, double *de0, double *p)
{ double u,v,th,de;
  double de1,th1,dt,ul,vl,ut,vt,x,y,sm,sml,tem,cdf,ccdf,lr,r,lx,xx,yy;
  th=*th0; de=*de0; u=*u0; v=*v0;
  de1=1./de; th1=1./th; dt=pow(de,th1);
  ul=-log(u); vl=-log(v);
  ut=pow(ul,th); vt=pow(vl,th);
  x=exp(ut*de)-1.; y=exp(vt*de)-1.;
  if(isinf(x) || isinf(y) || x>1.e200 || y>1.e200) 
  { lr=de*(vt-ut); r=exp(lr); lx=de*ut;
    sm=r+1.; sml=log(sm)+lx; 
    tem=pow(sml,th1);
    cdf=exp(-tem/dt);
    ccdf=cdf*tem*lx/sml/sm/ul/u/dt;
  }
  else if(x<=1.e-10 || y<=1.e-10)   // empirical
  { xx=de*ut; yy=vt*de; r=vt/ut;
    sm=xx+yy; 
    tem=pow(1.+r,th1-1);
    cdf=exp(-pow(sm,th1)/dt);
    ccdf=cdf*tem*(1+xx)/(1+xx+yy)/u;
  }
  else
  { sm=x+y+1.; sml=log(sm);
    tem=pow(sml,th1);
    cdf=exp(-tem/dt);
    ccdf=cdf*tem*de*ut*(x+1)/sml/sm/ul/u/dt;
  }
  *p=ccdf;
}

// code using transformed variables, without bisection method
/* $x=\exp\{\de(-\log u)^{\th}\}-1$ and $y=\exp\{\de(-\log v)^{\th}\}-1$
   $u=[\de^{-1}\log (x+1)]^{1/\th}$ and $v=[\de^{-1}\log (y+1)]^{1/\th}$
   $C(u,v)=G(x,y)=\exp\{-[\de^{-1}\log(x+y+1)]^{1/\th}\}$
   $C_{2|1}(v|u)= G(x,y) [\log(x+y+1) ]^{1/\th-1} 
     (x+y+1)^{-1} \de^{1-1/\th} (x+1) (-\log u)^{\th-1}/u$
*/
/* log G_{2|1}(y|x)= log p for y given x,p */
/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter (th.de). th>1, de>0
   output
     vv = inverse of conditional cdf
*/
void qcondbb3(double *p0, double *u0, double *th0, double *de0, double *vv)
{ double p,u,th,de;
  double de1,th1,dt,ul,ut,x,y,diff,sm,sml,smlt,eps,v;
  double r,con,h,hp,lx,llx,lxt,xx,xxt,smt;
  int iter,mxiter;

  th=*th0; de=*de0; p=*p0; u=*u0;
  mxiter=30; eps=1.e-6;
  de1=1./de; th1=1./th; dt=pow(de,th1);
  ul=-log(u); //vl=-log(v);
  ut=pow(ul,th); //vt=pow(vl,th);
  x=exp(ut*de)-1.; //y=exp(vt*de)-1.;
  if(isinf(x) || x>1.e200)
  { 
#ifdef DEBUG
    printf("\n** infinite x, p=%f u=%f th=%f de=%f ut=%f\n", p,u,th,de,ut);
#endif
    lx=ut*de; // log(x+1)
    llx=log(lx); lxt=pow(lx,th1);
    con=log(p) -lx +(th1-1.)*llx -lxt/dt;
    r=1.; diff=1; iter=0;
    while(fabs(diff)>eps  && iter<mxiter)
    { sm=r+1.; sml=log(sm)+lx; smlt=pow(sml,th1); 
      h=sml +(1.-th1)*log(sml) +smlt/dt +con;
      hp= 1./sm +(1.-th1)/sml/sm + th1*smlt/sml/sm/dt;
      iter=iter+1;
      diff=h/hp;
      r=r-diff;
      while(r<=0.) { diff=diff/2.; r=r+diff; }
#ifdef DEBUG
      printf("%d r=%f %f\n",  iter, r, diff);
#endif
    }
    v=pow(log(r)+lx,th1)/dt; v=exp(-v);
#ifdef DEBUG
    printf("** infinite x, p=%f u=%f th=%f de=%f ut=%f, v=%f\n",
       p,u,th,de,ut,v);
#endif
    *vv=v; return;
  }

  if(x<=1.e-10)
  { 
#ifdef DEBUG
    printf("\n** x near 0, p=%f u=%f th=%f de=%f ut=%f\n", p,u,th,de,ut);
#endif
    xx=ut*de;  xxt=pow(xx,th1);
    con=log(p) -xxt/dt;
    r=1.; diff=1; iter=0;
    while(fabs(diff)>eps  && iter<mxiter)
    { sm=r+1.; sml=log(sm); smt=pow(sm,th1); 
      h=xx*r +(1.-th1)*sml +smt*xxt/dt +con;
      hp= xx +(1.-th1)/sm + th1*smt*xxt/sm/dt;
      iter=iter+1;
      diff=h/hp;
      r=r-diff;
      while(r<=0.) { diff=diff/2.; r=r+diff; }
#ifdef DEBUG
      printf("%d r=%f %f\n",  iter, r, diff);
#endif
    }
    v=pow(r*xx/de,th1); v=exp(-v);
#ifdef DEBUG
    printf("** x near 0, p=%f u=%f th=%f de=%f ut=%f, v=%f\n",
       p,u,th,de,ut,v);
#endif
    *vv=v; return;
  }

  lx=ut*de; // log(x+1)
  llx=log(lx);  lxt=pow(lx,th1);
  con=log(p) -lx +(th1-1.)*llx -lxt/dt;
  // starting guess for y
  y=x;
  if(y<=0.) y=0.1;
  if(y>1000.) y=1000.;
#ifdef DEBUG
  printf("\np=%f, u=%f\n", p,u);
  printf(" x=%f, con=%f, and starting y=%f\n", x,con,y);
#endif
  
  diff=1; iter=0;
  while((fabs(diff/y)> eps) && iter<mxiter)
  { sm=x+y+1; sml=log(sm); smlt=pow(sml,th1);
    h=sml +(1.-th1)*log(sml) +smlt/dt + con;
    hp= 1./sm +(1.-th1)/sml/sm + th1*smlt/sml/sm/dt;
    iter=iter+1;
    diff=h/hp;
    y=y-diff;
    while(y<=0.) { diff=diff/2.; y=y+diff; }
#ifdef DEBUG
    printf("%d %f %f\n",  iter, y, diff);
#endif
  }
  v=pow(de1*log(y+1.),th1); v=exp(-v);
  if(iter>=mxiter && y>1.e-2) 
  { 
#ifdef DEBUG
    printf("** did not converge **\n");
    printf("(p,u,th,de)= %f %f %.2f %.2f lasty=%f v=%f\n", p,u,th,de,y,v);
#endif
  }
#ifdef DEBUG
  printf("v= %f  #iter=%d\n",v,iter);
#endif
  *vv=v;
}

