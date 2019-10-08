#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* Kendall's tau, computed with O(n log(n)) algorithm
   Knight's method, JASA 1966, v 61, pp 436-439.
*/
// gcc -fpic -c fastktau.c
// gcc -shared -o  fastktau.so fastktau.o
/* define boolean type for C */
typedef unsigned int boolean;
#define FALSE 0
#define TRUE (!FALSE)

#ifdef MAIN
int main()
{ void ktau(double *x0, double *y0, int *n0, double *tau, 
    double *numer, double *den, int *tx, int *ty, int *txy, int *pflag);
  int n0, tx, ty, txy, i;
  //double x[] = {1, 2.5, 2.5, 4.5, 4.5, 6.5, 6.5, 8, 9.5, 9.5};
  //double y[] = {1, 2, 4.5, 4.5, 4.5, 4.5, 8, 8, 8, 10};
  // Adjusted for ties, tau should be 0.859.
  double tau;
  boolean pflag;
  double x[8], y[8], numer, den;
  n0 = 8; tau = 0.; numer = 0.; den = 0.;
  tx = 0; ty = 0; txy = 0; pflag = 2;
  for(i=0;i<n0;i++) x[i]=20-i;
  x[2]=14;
  x[7]=14;
  for(i=0;i<n0;i++) y[i]=10-i;
  y[5]=10;
  y[2]=3;
  //n0 = 10; tau = 0.; pflag = 2;
  ktau(x, y, &n0, &tau, &numer, &den, &tx, &ty, &txy, &pflag);
}
#endif

void ktau(double *x, double *y, int *n0, double *tau,
  double *numer, double *den, int *tx, int *ty, int *txy, int *pflag)
{ // Defining variables
  int k, ell, i, j, iend, jend;
  int ii, ix, iy, ixy, n;
  double n2, nexch;

  n= *n0;
  double *y2=(double *)malloc(n*sizeof(double));
  double *x2=(double *)malloc(n*sizeof(double));
  double *xptr,*yptr; // HJ addition for swapping
  boolean iflag, jflag, xflag;
  nexch = 0; *den = 0.; *tx = 0; *ty = 0; *txy = 0;
  
  // Print for checking
  if(*pflag > 1)
  { printf("Unsorted\n------------\n");
    for(ii=0; ii<n; ii++) printf("x = %6.3g, y = %6.3g\n",x[ii], y[ii]);
    printf("------------\n\n");
  }
  if(*pflag > 0) printf("sample size:  n = %d\n", n);
  
  /* 1.1 Sort x and y in x order */
  /* Break ties in x according to y */
  if(*pflag > 2) printf("Sorting by x\n");
  k=1;
  do
  { ell=0;
    do
    { i = ell;
      j = (i+k)<(n)?(i+k):(n);
      iend = j;
      jend = (j+k)<(n)?(j+k):(n);
      do
      { iflag = (i < iend);
        jflag = (j < jend);
        //xflag = ((x[i] > x[j]) | ((x[i] == x[j]) & (y[i] > y[j])));
        // replacement in VineCopula, via Ulf Schepsmeier
        // xflag should use x[i] and x[j] only if i<iend and j<jend
        // otherwise x[i] or x[j] has index out of bounds
        if (iflag & jflag) 
        { xflag = ((x[i] > x[j]) | ((x[i] == x[j]) & (y[i] > y[j]))); }
        else { xflag = FALSE; }
        if((iflag & !jflag) | (iflag & jflag & !xflag))
        { x2[ell] = x[i]; y2[ell] = y[i];
          i++; ell++;
        }
        //if((!iflag & jflag) | (iflag & jflag & xflag))
        if( ( (!iflag) & jflag) | (iflag & jflag & xflag))
        { x2[ell] = x[j]; y2[ell] = y[j];
          j++; ell++;
        }
      }
      while(iflag | jflag);
    }
    while(ell < n);
    // Print iterations
    if(*pflag > 2)
    { printf("\nk = %d\n---\n", k);
      for(ii=0; ii<n; ii++) printf("x = %6.3g, y = %6.3g\n",x[ii], y[ii]);
      printf("---\n\n");
    }
    // Swap lists
    xptr=x; x=x2; x2=xptr;
    yptr=y; y=y2; y2=yptr;
    k *= 2;
  }
  while (k < n);
  
  /* 1.2 Count pairs of tied x, tx */
  ix = 1;
  ixy = 1;
  for(ii= 1; ii< n; ii++)
  { if(x[ii] == x[ii-1])
    { ix++;
      //if(y[ii] == y[ii-1]) ixy++;
      if(y[ii] == y[ii-1]) ixy++; 
      else { *txy += ixy*(ixy-1)/2; ixy=1; }
    }
    else if(ix > 1)
    { *tx += ix * (ix - 1) / 2;
      if(ixy > 1) *txy += ixy * (ixy - 1) / 2;
      ixy = 1;
      ix = 1;
    }
  }
  *tx += ix*(ix-1)/2;
  //printf("tx=%d, ix=%d, ixy=%d, txy=%d\n", *tx,ix,ixy,*txy);
  *txy += ixy*(ixy-1)/2;
  // Print for checking
  if(*pflag > 1)
  { printf("Sorted by x\n");
    printf("------------\n");
    for(ii= 0; ii< n; ii++) printf("x = %6.3g, y = %6.3g\n", x[ii], y[ii]);
    printf("------------\n\n");
  }
  if(*pflag > 0)
  { printf("Pairs of tied x:  tx = %d\n", *tx);
    printf("Pairs of both x and y tied:  txy = %d\n", *txy);
  }
  
  /* 2.1 Sort y again and count exchanges, nexch */
  /* keep original relative order if tied */
  if(*pflag > 2) printf("Sorting by y\n");
  k=1;
  do
  { ell=0;
    do
    { i = ell;
      j = (i+k)<(n)?(i+k):(n);
      iend = j;
      jend = (j+k)<(n)?(j+k):(n);
      do
      { iflag = (i < iend);
        jflag = (j < jend);
        //xflag = (y[i] > y[j]);
        // replacement in VineCopula, via Ulf Schepsmeier
        if (iflag & jflag) { xflag = (y[i] > y[j]); }
        else { xflag = FALSE; }
        if((iflag & !jflag) | (iflag & jflag & !xflag))
        { x2[ell] = x[i]; y2[ell] = y[i];
          i++; ell++;
        }
        //if((!iflag & jflag) | (iflag & jflag & xflag))
        if( ( (!iflag) & jflag) | (iflag & jflag & xflag))
        { x2[ell] = x[j]; y2[ell] = y[j];
          nexch += iend - i;
          j++; ell++;
        }
      }
      while(iflag | jflag);
    }
    while(ell < n);
    // Print iterations
    if(*pflag > 2)
    { printf("\nk = %d\n---\n", k);
      for(ii= 0; ii< n; ii++) printf("x = %6.3g, y = %6.3g\n",x[ii], y[ii]);
      printf("---\n\n");
    }
    // Swap lists
    xptr=x; x=x2; x2=xptr;
    yptr=y; y=y2; y2=yptr;
    k *= 2;
  }
  while (k < n);
  
  /* 2.2 Count pairs of tied y, ty */
  iy=1;
  for(ii= 1; ii< n; ii++)
  { if(y[ii] == y[ii-1]) iy++;
    else if(iy> 1)
    { *ty += iy * (iy - 1) / 2;
      iy = 1;
    }
  }
  *ty += iy * (iy - 1) / 2;
  // Print for checking
  if(*pflag > 1)
  { printf("Sorted in y\n");
    printf("------------\n");
    for(ii=0; ii<n; i++) printf("x = %6.3g, y = %6.3g\n", x[ii], y[ii]);
    printf("------------\n\n");
  }
  if(*pflag > 0)
  { printf("Pairs of tied y:  ty = %d\n", *ty);
    printf("Exchange count:  nexch = %.0f\n", nexch);
  }
  
  /* 3. Calc. Kendall's Score and Denominator */
  n2 = 0.5 * n * (n - 1);
  *numer = n2 - (2.*nexch + *tx + *ty - *txy);
  //if(*tx > 0 | *ty > 0) // adjust for ties
  *den = sqrt((n2 - *tx) * (n2 - *ty));
  *tau = *numer / *den;
  // Print for checking // change to dblepr and intpr for R
  if(*pflag > 0)
  { printf("Kendall's numerator:  numer = %f\n", *numer);
    printf("Denominator:  denom = %f\n", *den);
    printf("Kendall's tau:  numer/denom = %6.3g\n\n", *tau);
  }
  
  free(y2); free(x2);
}
