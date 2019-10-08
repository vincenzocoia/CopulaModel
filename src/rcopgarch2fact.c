#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// gcc -DNOTR -o rcopgarch2fact rcopgarch2fact.c qt.c rgum.c rfrk.c rbb1-v2.c miscutil.c -lm
//  rcopgarch2fact 3 < garchpar.in
//  rcopgarch2fact 5 < garchpar.in
//  rcopgarch2fact 9 < garchpar.in
#ifdef MAIN2F
#define NN 501
main(int argc, char *argv[])
{ int d,n,i,j,copcode,np1;
  double *gmu,*gar1,*gom,*galp1,*gbe1,*gnu,*gsig0;
  double *lgret0, *portfret,*cpar,*uvec,df;
  //double **dmatrix(int,int);
  void sim2fact(int *d, double *cpar, int *copcode, double *uvec);
  void sim2factmvt(int *d, double *cpar, int *copcode, double *uvec);
  void rgarch2fact(int *d, int *n, double *gmu, double *gar1, double *gom,
     double *galp1, double *gbe1, double *gnu, double *gsig0,
     double *cpar, int *copcode, double  *lgret, double *portfret);

  d=5; n=3;
  /* check on sim2fact() */
  cpar=(double *) malloc((3*d) * sizeof(double));
  uvec=(double *) malloc(d * sizeof(double));
  srand(123);
  if(argc>=2) copcode=atoi(argv[1]); else copcode=3;
  printf("\ncopcode=%d\n", copcode); 
  // cpar has length d except it is 2*d for BB1
  for(j=0;j<d;j++) cpar[j]=1.1+j*0.1;
  if(copcode==9) 
  { for(j=1;j<2*d;j+=2) cpar[j]=1.1+j*0.1;
    for(j=0;j<2*d-1;j+=2) cpar[j]=0.2;
    cpar[0]=.003;  // Gumbel for variable 1
    cpar[3]=1.003;  // MTCJ for variable 2
    np1=2;
  }
  else np1=1; 
  for(j=np1*d;j<np1*d+d;j++) cpar[j]=1.1;
  df=4.; if(copcode==1) df=1000.;
  if(copcode==1 || copcode==2)
  { cpar[2*d]=df;
    for(j=0;j<d;j++) cpar[j]=0.9*(j+1.)/d;
    for(j=d;j<2*d;j++) cpar[j]=0.4;
    for(i=0;i<n;i++)
    { sim2factmvt(&d,cpar,&copcode,uvec);
      // random normal or t
      for(j=0;j<d;j++) printf("%f ", uvec[j]); printf("\n");
    }
  }
  else
  { for(i=0;i<n;i++)
    { sim2fact(&d,cpar,&copcode,uvec);
      for(j=0;j<d;j++) printf("%f ", uvec[j]); printf("\n");
    }
  }

  gmu=(double *) malloc(d * sizeof(double));
  gar1=(double *) malloc(d * sizeof(double));
  gom=(double *) malloc(d * sizeof(double));
  galp1=(double *) malloc(d * sizeof(double));
  gbe1=(double *) malloc(d * sizeof(double));
  gnu=(double *) malloc(d * sizeof(double));
  gsig0=(double *) malloc(d * sizeof(double));
  srand(123); n=4;
  //lgret=dmatrix(n,d);
  lgret0=(double *) malloc((d*n) * sizeof(double));
  portfret=(double *) malloc(n * sizeof(double));
  for(j=0;j<d;j++)
  { scanf("%lf %lf %lf %lf %lf %lf", 
      &gmu[j],&gar1[j],&gom[j],&galp1[j],&gbe1[j],&gnu[j]);  
  }
  for(j=0;j<d;j++) scanf("%lf", &gsig0[j]);
  //rgarch2fact(&d,&n,gmu,gar1,gom,galp1,gbe1,gnu,cpar,&copcode,lgret,portfret);
  rgarch2fact(&d,&n,gmu,gar1,gom,galp1,gbe1,gnu,gsig0,cpar,&copcode,lgret0,portfret);
  for(i=0;i<n;i++)
  { //for(j=0;j<d;j++) printf("%f ", lgret[i][j]);
    for(j=0;j<d;j++) printf("%f ", lgret0[i*d+j]);
    printf("   %f\n", portfret[i]);
  }
  // free memory
  free(uvec); free(cpar);
  free(gmu); free(gar1); free(gom); free(galp1); free(gbe1); free(gnu);
  free(gsig0);
  free(portfret); free(lgret0);
}
  
double urand()
{ return(rand()/2147483648.); }
#endif

// cop=1: normal, 3: gumbel, 5: frank
#define BVN 1
#define BVT 2
#define GUM 3
#define GUMR -3
#define FRK 5
#define BB1 9
// BVN is multivariate Gaussian/normal with 2-factor structure
// BVT is multivariate t with 2-factor structure

/* inputs
     d0 = d = dimension
     n0 = n = simulation sample size
       AR(1)-GARCH(1,1) parameters (mu,ar1,om,alp1,be1,nu) for each asset
     gmu = d-vector of mu or location parameters 
     gar1 = d-vector of AR1 parameters 
     gom = d-vector of omega parameters 
     galp1 = d-vector of alpha parameters 
     gbe1 = d-vector of beta parameters 
     gnu = d-vector of nu parameters for Student t innovations
     gsig0 = d-vector of starting conditional SD values (one for each asset)
     cpar = parameter vector for 2-factor copula (dimension at least 2*d)
     copcode = copula code (see the above #define)
   outputs
     lgret0 = vector of length n*d of logreturns based on copula GARCH model 
              with 2-factor copula,
            lgret is  a vector for linking to R
     portfret = n-vector of portfolio returns (average of d assets)
*/
void rgarch2fact(int *d0, int *n0, double *gmu, double *gar1, double *gom, 
   double *galp1, double *gbe1, double *gnu, double *gsig0,
   double *cpar, int *copcode, double *lgret0, double *portfret)
{ int d,n,i,j,j0;
  double *tscale,*sigma2,*mu1,*z,*eps,*uvec,u,sportf,tem,prev,df;
  double qt(double,double),sqr(double),urand();
  double pt_(double*,double*),pnorms(double);
  void sim2fact(int *d, double *cpar, int * copcode, double *uvec);
  void sim2factmvt(int *d, double *cpar, int * copcode, double *uvec);

  n=*n0; d=*d0;
  df=1000.;
  if(*copcode==BVT) df=cpar[2*d];
#ifdef MAIN2F
  for(j=0;j<d;j++)
  { printf("%d : %f %f %f %f %f %f %f\n", j,
      gmu[j],gar1[j],gom[j],galp1[j],gbe1[j],gnu[j],gsig0[j]);  
  }
  printf("copcode=%d\n", *copcode);
  if(*copcode==BB1)
  { for(j=0;j<2*d;j++) printf("%f ", cpar[j]); printf("\n"); }
  else
  { for(j=0;j<d;j++) printf("%f ", cpar[j]); printf("\n"); }
  if(*copcode==BVT) { printf("df=%f ", df); printf("\n"); }
  // df=cpar[d];
#endif

  tscale=(double *) malloc(d * sizeof(double));
  mu1=(double *) malloc(d * sizeof(double));
  sigma2=(double *) malloc(d * sizeof(double));
  eps=(double *) malloc(d * sizeof(double));
  uvec=(double *) malloc(d * sizeof(double));
  z=(double *) malloc(d * sizeof(double));
  // set up temporary variables and initialize first observation 
  for(j0=0,sportf=0.;j0<d;j0++)
  { tscale[j0]=sqrt(1.-2./gnu[j0]); // multiplier so that variance of z's is 1
    mu1[j0]=gmu[j0]/(1-gar1[j0]);
    sigma2[j0]=sqr(gsig0[j0]);
    if(*copcode==BVN || *copcode==BVT) 
    { sim2factmvt(&d,cpar,copcode,uvec);
      if(*copcode==BVN) { for(j=0;j<d;j++) uvec[j]=pnorms(uvec[j]); }
      else { for(j=0;j<d;j++) uvec[j]=pt_(&uvec[j],&df); }
    }
    else { sim2fact(&d,cpar,copcode,uvec); }
    z[j0]=qt(uvec[j0],gnu[j0])*tscale[j0];
    eps[j0]=z[j0]*gsig0[j0];
    u=urand();
#ifdef MAIN2F
    printf("%d : %f %f %f\n", j0,tscale[j0],mu1[j0],sigma2[j0]);
#endif
    //lgret[0][j0]=mu1[j0]+0.05*(u-0.5);
    lgret0[j0]=mu1[j0]+0.05*(u-0.5);
    sportf = sportf + lgret0[j0];
  }
  portfret[0]=sportf/d;

  // generate remaining observations
  for(i=1;i<n;i++)
  { //printf("i=%d\n", i);
    if(*copcode==BVN || *copcode==BVT) 
    { sim2factmvt(&d,cpar,copcode,uvec);
      if(*copcode==BVN) { for(j=0;j<d;j++) uvec[j]=pnorms(uvec[j]); }
      else { for(j=0;j<d;j++) uvec[j]=pt_(&uvec[j],&df); }
    }
    else { sim2fact(&d,cpar,copcode,uvec); }
    //for(j=0;j<d;j++) printf("%f ", uvec[j]); printf("\n");
    for(j=0,sportf=0.;j<d;j++)
    { prev=lgret0[(i-1)*d+j];
      sigma2[j] = gom[j] + galp1[j] * sqr(eps[j]) + gbe1[j] * sigma2[j];
      z[j]=qt(uvec[j],gnu[j])*tscale[j];
      eps[j]=z[j]*sqrt(sigma2[j]);
      tem = gmu[j] + gar1[j]*prev + eps[j];
      lgret0[i*d+j] = tem;
      sportf = sportf + tem;
    }
    portfret[i]=sportf/d;
  }
#ifdef MAIN2F
  for(i=0;i<n;i++) printf("%f ", portfret[i]); printf("\n");
#endif
  free(mu1); free(tscale); free(sigma2); free(eps); free(uvec); free(z);
}

/* inputs
     d0 = d = dimension
     cpar = parameter vector with partial correlations with latent variables 
            and df for mvt
     copcode = copula code (see the above #define)
   output
     zvec = random normal or t d-vector with 2-factor structure
*/
void sim2factmvt(int *d0, double *cpar, int *copcode, double *zvec)
{ int cop,j,d;
  double df,w1,w2,z,rho,rhe,alp2,denom;
  double snorm(),rgamma(double,double); // in file rt.c
  cop=*copcode; d=*d0; df=cpar[2*d];
  //printf("copcode=%d\n", cop);
  w1=snorm();  // latent variable
  w2=snorm();  // second latent variable
  //printf("%f %f ", w1,w2);
  for(j=0;j<d;j++)
  { z=snorm();
    //printf("%f ", z);
    rho=cpar[j]; alp2=cpar[d+j]*sqrt(1.-rho*rho); 
    rhe=sqrt(1.-rho*rho-alp2*alp2);
    zvec[j]=rho*w1+alp2*w2+rhe*z;
  }
  //printf("\n");
  if(cop==BVT && df<300. && df>0.)  // bivariate t
  { denom=rgamma(df/2.,2.) / df; denom=sqrt(denom);
    for(j=0;j<d;j++) zvec[j]/=denom;
  }
}

/* inputs
     d0 = d = dimension
     cpar = parameter vector for 2-factor copula
     copcode = copula code (see the above #define)
   output
     uvec = random d-vector with U(0,1) margin and 2-factor structure
*/
void sim2fact(int *d0, double *cpar, int *copcode, double *uvec)
{ int cop,j,np1,jp1,jp2,d;
  void (*qcond1)(double*, double*, double *, double *);
  void (*qcond2)(double*, double*, double *, double *);
  double urand();
  double v1,v2,p,tem;
  void qcondgum(double *p0, double *u0, double *cpar, double *vv);
  void qcondgumr(double *p0, double *u0, double *cpar, double *vv);
  void qcondfrk(double *p0, double *u0, double *cpar, double *vv);
  void qcondbb1(double *p0, double *u0, double *cpar, double *vv);
  // qcond has arguments p in (0,1), condval in (0,1), cpar
  cop=*copcode;
  d=*d0;
  qcond1=qcondgum; qcond2=qcondgum; np1=1;
  switch(cop)
  { case GUM: qcond1=qcondgum; qcond2=qcondgum; np1=1; break;
    case GUMR: qcond1=qcondgumr; qcond2=qcondgumr; np1=1; break; 
    case FRK: qcond1=qcondfrk; qcond2=qcondfrk; np1=1; break;
    case BB1: qcond1=qcondbb1; qcond2=qcondfrk; np1=2; break;
    // np1 #parameter for factor1 linking copula
    // np2 =1 for all models considered here
    // for t(nu) [to handle later]
    // for BB1 cpar=(theta,delta)
  }
  //Rprintf("%d %d\n", d, *copcode);
  v1=urand();  // latent variable
  v2=urand();  // second latent variable
  //printf("%f ", v1);
  for(j=0,jp1=0,jp2=d*np1;j<d;j++)
  { p=urand();
    //printf("%f ", p);
    qcond2(&p,&v2,&cpar[jp2],&tem);
    qcond1(&tem,&v1,&cpar[jp1],&uvec[j]);
    jp1=jp1+np1; jp2=jp2+1;
  }
  /*
  printf("\n");
  if(cop==BB1)
  { for(j=0;j<2*d;j++) printf("%f ", cpar[j]); printf("\n"); }
  else
  { for(j=0;j<d;j++) printf("%f ", cpar[j]); printf("\n"); }
  for(j=0;j<d;j++) printf("%f ", uvec[j]); printf("\n");
  */
}
