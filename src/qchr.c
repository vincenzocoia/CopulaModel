#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* sample main program */
// gcc -DMAINB -o qchr qchr.c subrhr.c ../src/qt.c -lm
#ifdef MAINB
main(int argc, char *argv[])
{ int n,i;
  double *pvec,*uvec,*parvec,*vvec;
  void qchr(int *n, double *pvec, double *uvec, double *parvec, double *vvec);
  scanf("%d", &n);
  pvec=(double *) malloc(n * sizeof(double));
  uvec=(double *) malloc(n * sizeof(double));
  parvec=(double *) malloc(n * sizeof(double));
  vvec=(double *) malloc(n * sizeof(double));
  for(i=0;i<n;i++) scanf("%lf %lf %lf", &pvec[i],&uvec[i],&parvec[i]);
  qchr(&n,pvec,uvec,parvec,vvec);
  for(i=0;i<n;i++) printf("%f %f %f %f\n", pvec[i],uvec[i],parvec[i],vvec[i]);
  free(pvec); free(uvec); free(parvec); free(vvec); 
}
#endif

/* vectorized version of qcondhr in rhr.c */
/* inputs
     n0 = n = number of qcondhr to compute
     pvec, uvec = n-vectors of values in (0,1)
     parvec = n-vector of copula parameters, each >0
   output
     vvec = vector of inverse of conditional cdf
*/
void qchr(int *n0, double *pvec, double *uvec, double *parvec, double *vvec)
{ int n,i;
  double p,u,param,v;
  void qcondhr(double *p0, double *u0, double *param, double *vv);
  n=*n0;
  for(i=0;i<n;i++)
  { p=pvec[i]; u=uvec[i]; param=parvec[i];
    qcondhr(&p, &u, &param, &v);
    vvec[i]=v;
  }
  return;
}

#ifdef P
/* vectorized version of pcondhr in rhr.c */
/* inputs
     n0 = n = number of pcondhr to compute
     vvec, uvec = n-vectors of values in (0,1)
     parvec = n-vector of copula parameters, each >1
   output
     pvec = vector of conditional cdf
*/
void pchr(int *n0, double *vvec, double *uvec, double *parvec, double *pvec)
{ int n,i;
  double v,u,param,p;
  void pcondhr(double *v0, double *u0, double *param, double *pp);
  n=*n0;
  for(i=0;i<n;i++)
  { v=vvec[i]; u=uvec[i]; param=parvec[i];
    pcondhr(&v, &u, &param, &p);
    pvec[i]=p;
  }
  return;
}
#endif


