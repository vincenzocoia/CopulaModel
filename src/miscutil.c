#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* utility functions used by a number of other functions */

#ifdef DMAT
/* allocate a matrix contiguous in memory
   inputs
     nrow = number of rows
     ncol = number of columns
   output
     amat = pointer to beginning of array
     amat[ ] = vector of pointers to beginning of rows of amat
*/
double **dmatrix(int nrow, int ncol)
{ int i; double **amat,*avec;
  avec=(double *)malloc((unsigned) (nrow*ncol)*sizeof(double));
  amat=(double **)malloc((unsigned) nrow*sizeof(double*));
  for(i=0;i<nrow;i++) amat[i]=avec+i*ncol;
  return amat;
}
#endif

/* square of a real number */
double sqr(double x) { return(x*x); }

