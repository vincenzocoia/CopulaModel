#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// gcc -DNOTR -o rcopgarch1fact rcopgarch1fact.c qt.c rgum.c rfrk.c rbb1-v2.c miscutil.c -lm
//  rcopgarch1fact 3 < garchpar.in
//  rcopgarch1fact -3 < garchpar.in
//  rcopgarch1fact 5 < garchpar.in
//  rcopgarch1fact 9 < garchpar.in
//  rcopgarch1fact 1 < garchpar.in
//  rcopgarch1fact 2 < garchpar.in
#ifdef MAIN1F
#define NN 501
main(int argc, char *argv[])
{ int d,n,i,j,copcode;
  double *gmu,*gar1,*gom,*galp1,*gbe1,*gnu,*gsig0;
  double *lgret0,**lgret, *portfret,*cpar,*uvec,df;
  //double **dmatrix(int,int);
  void sim1fact(int *d, double *cpar, int *copcode, double *uvec);
  void sim1factmvt(int *d, double *cpar, int *copcode, double *uvec);
  void rgarch1fact(int *d, int *n, double *gmu, double *gar1, double *gom,
     double *galp1, double *gbe1, double *gnu, double *gsig0,
     double *cpar, int *copcode, double  *lgret, double *portfret);

  d=5; n=3;
  /* check on sim1fact() */
  cpar=(double *) malloc((2*d) * sizeof(double));
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
  }
  df=4.; if(copcode==1) df=1000.;
  if(copcode==1 || copcode==2)
  { cpar[d]=df;
    for(j=0;j<d;j++) cpar[j]=0.9*(j+1.)/d;
    for(i=0;i<n;i++)
    { sim1factmvt(&d,cpar,&copcode,uvec);
      // random normal or t
      for(j=0;j<d;j++) printf("%f ", uvec[j]); printf("\n");
    }
  }
  else
  { for(i=0;i<n;i++)
    { sim1fact(&d,cpar,&copcode,uvec);
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
  //rgarch1fact(&d,&n,gmu,gar1,gom,galp1,gbe1,gnu,cpar,&copcode,lgret,portfret);
  rgarch1fact(&d,&n,gmu,gar1,gom,galp1,gbe1,gnu,gsig0,cpar,&copcode,lgret0,portfret);
  for(i=0;i<n;i++)
  { //for(j=0;j<d;j++) printf("%f ", lgret[i][j]);
    for(j=0;j<d;j++) printf("%f ", lgret0[i*d+j]);
    printf("   %f\n", portfret[i]);
  }
  // free memory
  free(uvec); free(cpar);
  free(gmu); free(gar1); free(gom); free(galp1); free(gbe1); free(gnu);
  free(gsig0);
  free(portfret); free(lgret0); //free(lgret[0]); free(lgret);
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
// BVN is multivariate Gaussian/normal with 1-factor structure
// BVT is multivariate t with 1-factor structure

/* check on generation from garch model 
     eps[t] = z[t]*sigma[t]
     lgret[t]= mu + ar1*lgret[t-1] + eps[t]
     sigma[t]^2 = om + alp1 * eps[t-1]^2 + be1 * sigma[t-1]^2
   alp1>0, be1>0, alp1+be1<1
*/

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
     cpar = parameter vector for 1-factor copula
      (dimension d for copcode=1,3,-3,5; 2*d for copcode=9, d+1 fo copcode=2)
     copcode = copula code (see the above #define)
   outputs
     lgret0 = vector of length n*d of logreturns based on copula GARCH model 
              with 1-factor copula,
            lgret is  a vector for linking to R
     portfret = n-vector of portfolio returns (average of d assets)
*/
void rgarch1fact(int *d0, int *n0, double *gmu, double *gar1, double *gom, 
   double *galp1, double *gbe1, double *gnu, double *gsig0,
   double *cpar, int *copcode, double *lgret0, double *portfret)
   //double **lgret, double *portfret)
{ int d,n,i,j,j0;
  double *tscale,*sigma2,*mu1,*z,*eps,*uvec,u,sportf,tem,prev,df;
  double qt(double,double),sqr(double),urand();
  double pt_(double*,double*),pnorms(double);
  void sim1fact(int *d, double *cpar, int * copcode, double *uvec);
  void sim1factmvt(int *d, double *cpar, int *copcode, double *uvec);

  n=*n0; d=*d0;
  //for(i=0;i<n;i++) lgret[i]=lgret0+i*d;  // set memory addresses
  df=1000.;
  if(*copcode==BVT) df=cpar[d];
#ifdef MAIN1F
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
    //sigma2[j0]=(gom[j0]+galp1[j0]*sqr(mu1[j0]))/(1.-gbe1[j0]);
    sigma2[j0]=sqr(gsig0[j0]);
    if(*copcode==BVN || *copcode==BVT) 
    { sim1factmvt(&d,cpar,copcode,uvec);
      if(*copcode==BVN) { for(j=0;j<d;j++) uvec[j]=pnorms(uvec[j]); }
      else { for(j=0;j<d;j++) uvec[j]=pt_(&uvec[j],&df); }
    }
    else { sim1fact(&d,cpar,copcode,uvec); }
    z[j0]=qt(uvec[j0],gnu[j0])*tscale[j0];
    eps[j0]=z[j0]*gsig0[j0];
    u=urand();
#ifdef MAIN1F
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
    { sim1factmvt(&d,cpar,copcode,uvec);
      if(*copcode==BVN) { for(j=0;j<d;j++) uvec[j]=pnorms(uvec[j]); }
      else { for(j=0;j<d;j++) uvec[j]=pt_(&uvec[j],&df); }
    }
    else { sim1fact(&d,cpar,copcode,uvec); }
    //for(j=0;j<d;j++) printf("%f ", uvec[j]); printf("\n");
    for(j=0,sportf=0.;j<d;j++)
    { //sigma2[j] = gom[j] + galp1[j] * sqr(lgret[i-1][j]) + gbe1[j] * sigma2[j];
      prev=lgret0[(i-1)*d+j];
      sigma2[j] = gom[j] + galp1[j] * sqr(eps[j]) + gbe1[j] * sigma2[j];
      z[j]=qt(uvec[j],gnu[j])*tscale[j];
      eps[j]=z[j]*sqrt(sigma2[j]);
      //tem = gmu[j] + gar1[j]*lgret[i-1][j] + z[j]*sqrt(sigma2[j]);
      tem = gmu[j] + gar1[j]*prev + eps[j];
      //lgret[i][j] = tem;
      lgret0[i*d+j] = tem;
      sportf = sportf + tem;
    }
    portfret[i]=sportf/d;
  }
  //for(j=0;j<d;j++) printf("%f ", lgret0[0+j]); printf("\n");
  //for(j=0;j<d;j++) printf("%f ", lgret0[d+j]); printf("\n");
  //for(j=0;j<d;j++) printf("%f ", lgret0[2*d+j]); printf("\n");
#ifdef MAIN1F
  for(i=0;i<n;i++) printf("%f ", portfret[i]); printf("\n");
#endif
  free(mu1); free(tscale); free(sigma2); free(eps); free(uvec); free(z);
}

/* inputs
     d0 = d = dimension
     cpar = parameter vector with correlations with latent and df for mvt
     copcode = copula code (see the above #define)
   output
     zvec = random normal or t d-vector with 1-factor structure
*/
void sim1factmvt(int *d0, double *cpar, int *copcode, double *zvec)
{ int cop,j,d;
  double df,w1,z,rho,rhe,denom;
  double snorm(),rgamma(double,double); // in file rt.c
  cop=*copcode; d=*d0; df=cpar[d];
  //printf("copcode=%d\n", cop);
  w1=snorm();  // latent variable
  //printf("%f ", w1);
  for(j=0;j<d;j++)
  { z=snorm();
    //printf("%f ", z);
    rho=cpar[j]; rhe=sqrt(1.-rho*rho);
    zvec[j]=rho*w1+rhe*z;
  }
  //printf("\n");
  if(cop==BVT && df<300. && df>0.)  // bivariate t
  { denom=rgamma(df/2.,2.) / df; denom=sqrt(denom);
    for(j=0;j<d;j++) zvec[j]/=denom;
  }
}

/* inputs
     d0 = d = dimension
     cpar = parameter vector for 1-factor copula
     copcode = copula code (see the above #define)
   output
     uvec = random d-vector with U(0,1) margin and 1-factor structure
*/
void sim1fact(int *d0, double *cpar, int *copcode, double *uvec)
{ int cop,j,np,jp,d;
  void (*qcond)(double*, double*, double *, double *);
  double urand();
  double v1,p;
  void qcondgum(double *p0, double *u0, double *cpar, double *vv);
  void qcondgumr(double *p0, double *u0, double *cpar, double *vv);
  void qcondfrk(double *p0, double *u0, double *cpar, double *vv);
  void qcondbb1(double *p0, double *u0, double *cpar, double *vv);
  // qcond has arguments p in (0,1), condval in (0,1), cpar
  cop=*copcode; d=*d0;
  qcond=qcondgum; np=1;
  switch(cop)
  { case GUM: qcond=qcondgum; np=1; break;
    case GUMR: qcond=qcondgumr; np=1; break; 
    case FRK: qcond=qcondfrk; np=1; break;
    case BB1: qcond=qcondbb1; np=2; break;
    // for BB1 cpar=(theta,delta)
  }
  //Rprintf("%d %d\n", d, *copcode);
  v1=urand();  // latent variable
  //printf("%f ", v1);
  for(j=0,jp=0;j<d;j++)
  { p=urand();
    //printf("%f ", p);
    qcond(&p,&v1,&cpar[jp],&uvec[j]);
    jp=jp+np;
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
