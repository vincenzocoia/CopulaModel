#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// there is isinf function in math.h
// gcc -DMAINC -o rmfrk rmfrk.c -lm
#ifdef MAINC
main(int argc, char *argv[])
{ double cpar,*uu;
  long seed;
  double s1,s2,s12; 
  int i,j,n,d;
  void rmfrk(int *, int *, double*, double *);
  scanf("%lf %d %d %ld", &cpar,&n,&d,&seed);
  while(cpar>1.)
  { srand(seed);
    printf("\ncpar=%f nsim=%d d=%d, seed=%ld\n", cpar,n,d,seed);
    uu=(double *) malloc(n*d * sizeof(double));
    rmfrk(&n,&d,&cpar,uu);
    for(i=0;i<5;i++)
    { for(j=0;j<d;j++) printf("%f ", uu[j*n+i]); printf("\n"); }
    for(i=0,s1=0.,s2=0.,s12=0.;i<n;i++)
    { s1+=uu[i]; s2+=uu[n+i]; s12+=uu[i]*uu[n+i]; }
    s1/=n; s2/=n; s12/=n; s12-=s1*s2;
    printf("mean1=%f, mean2=%f, cov=%f, corr=%f\n", s1,s2,s12,12*s12);
    printf("\n");
    free(uu);
    scanf("%lf %d %d %ld", &cpar,&n,&d,&seed);
  }
}


double urand()
{ return(rand()/2147483648.); }
#else
// use rng in R for linking
#include <R.h>
#include <Rmath.h>
#endif

// double urand() { return(unif_rand()); }  // in rfact.c

// logseries random variable, Kemp's algorithm
// This works better when parameter is cpar
/* input
     cpar = copula parameter for Frank copula >0
   output
     random variate from log series distribution, based on Kemp's method
*/
double rkemp1(double cpar)
{ double alp,v,u,tem,xx,x;
  double urand();
  alp=1.-exp(-cpar);
  x=1;
  v=urand();
  if(v>=alp) return(x);
  u=urand();
  tem=exp(-cpar*u);
  xx= (1+log(v)/log(1.-tem));
  //printf("%f %f %e %f\n", cpar,u,tem,xx);
  if(abs(isinf(xx))==1) 
  { x=floor(1.-log(v)/tem); 
    //printf("v=%f, u=%f, tem=%e, xx=%f, x=%f\n", v,u,tem,xx,x);
  }
  else { x=floor(xx); }
  // occurs where u large enough, e.g., >0.8 for cpar=45
  return(x);
}

// logseries stochastic representation of Archimedean copula
/* input
     nn = simulation sample size
     dd = dimension
     param = cpar = copula parameter >0
   output
     uu = nnxdd matrix with random d-vectors from multivariate Frank copula 
*/
void rmfrk(int *nn, int *dd, double *param, double *uu)
{ double v,cpar,cpar0,alp,tem,tem1,r;
  int i,j,n,d;
  double rkemp1(double cpar);
  double urand();
#ifndef MAINC
  GetRNGstate();  // in R/R_ext/Random.h
#endif
  n=*nn; d=*dd; cpar=*param;
  cpar0=exp(-cpar); alp=1.-cpar0;
  for(i=0;i<n;i++)
  { r=rkemp1(cpar);
    for(j=0;j<d;j++)
    { v=urand(); tem=-log(v)/r;
      tem1= -log(1.-alp*exp(-tem))/cpar;
      if(abs(isinf(tem1))==1) tem1=-log(cpar0+alp*tem)/cpar;
      uu[j*n+i]=tem1; 
    }
  }
#ifndef MAINC
  PutRNGstate();  
#endif
}

