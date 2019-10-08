#include <stdio.h>
#include <math.h>
/*  routine to link from R to mulnor (in pmnorm.c) */
/* inputs
     ub = vector of upper limit
     lb = vector of lower limit
     sig0 = vector from lower triangle of correlation matrix, by rows
     eps0 = tolerance for accuracy
     n0 = dimension of ub and lb
     inf0(i) = 0 if ith range is (lb(i),oo)
     inf0(i) = 1 if ith range is (-oo,ub(i))
     inf0(i) = 2 if ith range is (lb(i),ub(i))
   outputs
     prob = rectangle probability
     bound = error bound
     ifault = error code
     ifault = error code
              4 is sig0 is not positive definite
              1,2,3,5 : see Schervish's paper
              0 OK
*/

void mvnscher(double ub[], double lb[], double sig0[], double *eps0, int *n0, 
   int inf0[], double *prob, double *bound, int *ifault)
{  double eps;
   int n;
 /*int i,j; */
   void mulnor(double [], double [], double [], double, int,
               int [], double *, double *, int *);
   eps= *eps0; n= *n0;

 /* for(i=0;i<n;i++) printf("%8.4f", lb[i]);  printf("\n");
    for(i=0;i<n;i++) printf("%8.4f", ub[i]);  printf("\n");
    for(j=0;j<(n*(n-1))/2;j++) printf("%8.4f", sig0[j]); printf("\n"); 
 */
   mulnor(ub,lb,sig0,eps,n,inf0,prob,bound,ifault);
 /* printf("%f %f %d\n", *prob,*bound,*ifault); */
}
