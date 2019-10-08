#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#define POSINF 1.e310 
//#define NEGINF -1.e310 
/* sample main program */
// gcc -DMAINB -DQ -o qcgum qcgum.c -lm
// gcc -DQ -fpic -c qcgum.c
// gcc -shared -o qcgum.so qcgum.o
#ifdef MAINB
main(int argc, char *argv[])
{ int n,i;
  double *pvec,*uvec,*parvec,*vvec;
  void qcgum(int *n, double *pvec, double *uvec, double *parvec, double *vvec);
  scanf("%d", &n);
  pvec=(double *) malloc(n * sizeof(double));
  uvec=(double *) malloc(n * sizeof(double));
  parvec=(double *) malloc(n * sizeof(double));
  vvec=(double *) malloc(n * sizeof(double));
  for(i=0;i<n;i++) scanf("%lf %lf %lf", &pvec[i],&uvec[i],&parvec[i]);
  qcgum(&n,pvec,uvec,parvec,vvec);
  for(i=0;i<n;i++) printf("%f %f %f %f\n", pvec[i],uvec[i],parvec[i],vvec[i]);
  free(pvec); free(uvec); free(parvec); free(vvec); 
}
#endif

/* vectorized version of qcondgum in rgum.c */
/* inputs
     n0 = n = number of qcondgum to compute
     pvec, uvec = n-vectors of values in (0,1)
     parvec = n-vector of copula parameters, each >1
   output
     vvec = vector of inverse of conditional cdf
*/
void qcgum(int *n0, double *pvec, double *uvec, double *parvec, double *vvec)
{ int n,i;
  double p,u,param,v;
  void qcondgum(double *p0, double *u0, double *param, double *vv);
  n=*n0;
  for(i=0;i<n;i++)
  { p=pvec[i]; u=uvec[i]; param=parvec[i];
    qcondgum(&p, &u, &param, &v);
    vvec[i]=v;
  }
  return;
}

#ifdef P
/* vectorized version of pcondgum in rgum.c */
/* inputs
     n0 = n = number of pcondgum to compute
     vvec, uvec = n-vectors of values in (0,1)
     parvec = n-vector of copula parameters, each >1
   output
     pvec = vector of conditional cdf
*/
void pcgum(int *n0, double *vvec, double *uvec, double *parvec, double *pvec)
{ int n,i;
  double v,u,param,p;
  void pcondgum(double *v0, double *u0, double *param, double *pp);
  n=*n0;
  for(i=0;i<n;i++)
  { v=vvec[i]; u=uvec[i]; param=parvec[i];
    pcondgum(&v, &u, &param, &p);
    pvec[i]=p;
  }
  return;
}
#endif


#ifdef Q
/* inputs
     p0, u0 = values in (0,1)
     cpar = copula parameter >1
   output
     vv = inverse of conditional cdf
*/
void qcondgum(double *p0, double *u0, double *cpar, double *vv)
{ double p,u,de,z,g,gp,x,y,con,de1,dif,mxdif,eps;
  int iter,mxiter;
  de=*cpar; p=*p0; u=*u0; mxiter=20; eps=1.e-6;
  x=-log(u); de1=de-1.;
  con=log(p)-x-de1*log(x); 
  z=x*pow(2.,1./de); // this is for z 
  mxdif=1; iter=0; 
  dif=.1;  // needed in case first step leads to NaN
  while(mxdif>eps && iter<mxiter)
  { g=z+de1*log(z)+con;
    gp=1.+de1/z;
    if(isnan(g) || isnan(gp) || isnan(g/gp) ) { dif/=-2.; }  // added for de>50
    else dif=g/gp;
    z-=dif; iter++;
    while(z<=x) { dif/=2.; z+=dif; }
#ifdef DEBUG
    printf("%d %f %f\n", iter, dif, z);
#endif
    mxdif=fabs(dif);
  }
#ifdef MAINB
  if(iter>=mxiter) 
  { printf("***did not converge\n");
    printf("p=%f, x=%f, delta=%f, lastz=%f\n", p,x,de,z);
  }
#endif
  y=pow(pow(z,de)-pow(x,de),1./de);
  *vv=exp(-y);
}
#endif
