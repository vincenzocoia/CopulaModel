#include <R.h>
#include <Rmath.h>

/* unif_rand() is the RNG in R: this function is in file urand.c */
//double urand() { return(unif_rand()); }


/* interface functions linking C functions with RNG and seed in R */
/* inputs
     n0 = simulation sample size
     d0 = dimension
     cpar = parameter vector for 1-factor copula
     copcode = copula code (see rcopgarch1fact.c for details)
   output
     udat = n0 x d0 matrix, random d-vectors from the 1-factor copula
*/
void r1fact(int *n0, int *d0, double *cpar, int *copcode, double *udat) 
{ void sim1fact(int *d, double *cpar, int *copcode, double *uvec);
  double *uvec;
  int i,j,n,d;
  n=*n0; d=*d0;
  //for(j=0;j<d;j++) Rprintf("%f ", cpar[j]);  Rprintf("\n");
  uvec=(double *) malloc(d * sizeof(double));
  //Rprintf("%d %d %d\n", n,d, *copcode);
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { //Rprintf("i=%d\n", i);
    sim1fact(&d,cpar,copcode,uvec);
    //for(j=0;j<d;j++) Rprintf("%f ", uvec[j]);  Rprintf("\n");
    for(j=0;j<d;j++) udat[i+j*n]=uvec[j];
  }
  PutRNGstate();
  free(uvec);
}

// random Gaussian/normal or t vector with 1-factor structure
//       df = cpar[d] for latter
/* inputs
     n0 = simulation sample size
     d0 = dimension
     cpar = parameter vector for mvt copula with 1-factor structure
     copcode = copula code (see rcopgarch1fact.c for details)
   output
     udat = n0 x d0 matrix, random d-vectors from mvt 1-factor copula
*/
void r1factmvt(int *n0, int *d0, double *cpar, int *copcode, double *zdat) 
{ void sim1factmvt(int *d, double *cpar, int *copcode, double *zvec);
  double *zvec;
  int i,j,n,d;
  n=*n0; d=*d0;
  zvec=(double *) malloc(d * sizeof(double));
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { sim1factmvt(&d,cpar,copcode,zvec);
    for(j=0;j<d;j++) zdat[i+j*n]=zvec[j];
  }
  PutRNGstate();
  free(zvec);
}

/* inputs
     n0 = simulation sample size
     d0 = dimension
     cpar = parameter vector for 2-factor copula
     copcode = copula code (see rcopgarch2fact.c for details)
   output
     udat = n0 x d0 matrix, random d-vectors from the 2-factor copula
*/
void r2fact(int *n0, int *d0, double *cpar, int *copcode, double *udat) 
{ void sim2fact(int *d, double *cpar, int *copcode, double *uvec);
  double *uvec;
  int i,j,n,d;
  n=*n0; d=*d0;
  //for(j=0;j<d;j++) Rprintf("%f ", cpar[j]);  Rprintf("\n");
  uvec=(double *) malloc(d * sizeof(double));
  //Rprintf("%d %d %d\n", n,d, *copcode);
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { sim2fact(&d,cpar,copcode,uvec);
    //for(j=0;j<d;j++) Rprintf("%f ", uvec[j]);  Rprintf("\n");
    for(j=0;j<d;j++) udat[i+j*n]=uvec[j];
  }
  PutRNGstate();
  free(uvec);
}

// random Gaussian/normal or t vector with 2-factor structure
//       df = cpar[d] for latter
/* inputs
     n0 = simulation sample size
     d0 = dimension
     cpar = parameter vector for 2-factor copula
     copcode = copula code (see rcopgarch2fact.c for details)
   output
     udat = n0 x d0 matrix, random d-vectors from the mvt 2-factor copula
*/
void r2factmvt(int *n0, int *d0, double *cpar, int *copcode, double *zdat) 
{ void sim2factmvt(int *d, double *cpar, int *copcode, double *zvec);
  double *zvec;
  int i,j,n,d;
  n=*n0; d=*d0;
  zvec=(double *) malloc(d * sizeof(double));
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { sim2factmvt(&d,cpar,copcode,zvec);
    for(j=0;j<d;j++) zdat[i+j*n]=zvec[j];
  }
  PutRNGstate();
  free(zvec);
}

// random Gaussian/normal or t vector with bi-factor structure
//       df = cpar[d] for latter
/* inputs
     n0 = simulation sample size
     d0 = dimension
     cpar = parameter vector for mvt copula with bi-factor structure 
     copcode = copula code (see rcopgarchbifact.c for details)
   output
     udat = n0 x d0 matrix, random d-vectors from the mvt bi-factor copula
*/
void rbifactmvt(int *n0, int *m0, int *grsize, double *cpar, int *copcode, double *udat) 
{ void simbifactmvt(int *m0, int *grsize, double *cpar, int *copcode, double *uvec);
  double *uvec;
  int i,j,n,d,m,jg;
  n=*n0; m=*m0; 
  for(jg=0,d=0;jg<m;jg++) d+=grsize[jg];
  uvec=(double *) malloc(d * sizeof(double));
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { simbifactmvt(&m,grsize,cpar,copcode,uvec);
    for(j=0;j<d;j++) udat[i+j*n]=uvec[j];
  }
  PutRNGstate();
  free(uvec);
}

/* inputs
     n0 = simulation sample size
     d0 = dimension
     cpar = parameter vector for bi-factor copula
     copcode = copula code (see rcopgarchbifact.c for details)
   output
     udat = n0 x d0 matrix, random d-vectors from the bi-factor copula
*/
void rbifact(int *n0, int *m0, int *grsize, double *cpar, int *copcode, double *udat) 
{ void simbifact(int *m0, int *grsize, double *cpar, int *copcode, double *uvec);
  double *uvec;
  int i,j,n,d,m,jg;
  n=*n0; m=*m0; 
  for(jg=0,d=0;jg<m;jg++) d+=grsize[jg];
  //for(j=0;j<d;j++) Rprintf("%f ", cpar[j]);  Rprintf("\n");
  uvec=(double *) malloc(d * sizeof(double));
  //Rprintf("%d %d %d\n", n,d, *copcode);
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { simbifact(&m,grsize,cpar,copcode,uvec);
    //for(j=0;j<d;j++) Rprintf("%f ", uvec[j]);  Rprintf("\n");
    for(j=0;j<d;j++) udat[i+j*n]=uvec[j];
  }
  PutRNGstate();
  free(uvec);
}

// random Gaussian/normal or t vector with nested-factor structure
//       df = cpar[d] for latter
/* inputs
     n0 = simulation sample size
     d0 = dimension
     cpar = parameter vector for mvt copula with nested-factor structure
     copcode = copula code (see rcopgarchnestfact.c for details)
   output
     udat = n0 x d0 matrix, random d-vectors from the nested-factor copula
*/
void rnestfactmvt(int *n0, int *m0, int *grsize, double *cpar, int *copcode, double *udat) 
{ void simnestfactmvt(int *m0, int *grsize, double *cpar, int *copcode, double *uvec);
  double *uvec;
  int i,j,n,d,m,jg;
  n=*n0; m=*m0; 
  for(jg=0,d=0;jg<m;jg++) d+=grsize[jg];
  uvec=(double *) malloc(d * sizeof(double));
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { simnestfactmvt(&m,grsize,cpar,copcode,uvec);
    for(j=0;j<d;j++) udat[i+j*n]=uvec[j];
  }
  PutRNGstate();
  free(uvec);
}

/* inputs
     n0 = simulation sample size
     d0 = dimension
     cpar = parameter vector for nested-factor copula
     copcode = copula code (see rcopgarchnestfact.c for details)
   output
     udat = n0 x d0 matrix, random d-vectors from the nested-factor copula
*/
void rnestfact(int *n0, int *m0, int *grsize, double *cpar, int *copcode, double *udat) 
{ void simnestfact(int *m0, int *grsize, double *cpar, int *copcode, double *uvec);
  double *uvec;
  int i,j,n,d,m,jg;
  n=*n0; m=*m0; 
  for(jg=0,d=0;jg<m;jg++) d+=grsize[jg];
  //for(j=0;j<d;j++) Rprintf("%f ", cpar[j]);  Rprintf("\n");
  uvec=(double *) malloc(d * sizeof(double));
  //Rprintf("%d %d %d\n", n,d, *copcode);
  GetRNGstate();  // in R/R_ext/Random.h
  for(i=0;i<n;i++)
  { simnestfact(&m,grsize,cpar,copcode,uvec);
    //for(j=0;j<d;j++) Rprintf("%f ", uvec[j]);  Rprintf("\n");
    for(j=0;j<d;j++) udat[i+j*n]=uvec[j];
  }
  PutRNGstate();
  free(uvec);
}

