#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* vectorization of qcond for other copulas, as in qcgum.c */
/* vectorized version of qcondbb1 in rbb1_v2.c */
/* inputs
     n0 = n = number of qcondbb1 to compute
     pvec, uvec = n-vectors of values in (0,1)
     thvec = n-vector of theta parameters, each >0
     devec = n-vector of delta parameters, each >1
   output
     vvec = vector of inverse of conditional cdf
*/
void qcbb1(int *n0, double *pvec, double *uvec, double *thvec, double *devec, 
  double *vvec)
{ int n,i;
  //double p,u,th,de,v;
  double p,u,v,cpar[2];
  void qcondbb1(double *p0, double *u0, double *cpar, double *vv);
  n=*n0;
  for(i=0;i<n;i++)
  { p=pvec[i]; u=uvec[i]; cpar[0]=thvec[i]; cpar[1]=devec[i];
    qcondbb1(&p, &u, cpar, &v);
    //qcondbb1(&p, &u, &th, &de, &v);
    vvec[i]=v;
  }
  return;
}

