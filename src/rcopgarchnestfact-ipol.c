#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// gcc -DNOTR -DAGAIN -o rcopgarchnestfact-ipol rcopgarchnestfact-ipol.c qt.c rgum.c rfrk.c rbb1-v2.c miscutil.c -lm
//  rcopgarchnestfact-ipol 3 < garchpar6.in
//  rcopgarchnestfact-ipol 5 < garchpar6.in
//  rcopgarchnestfact-ipol 9 < garchpar6.in
#ifdef MAIN2F
#define NN 501
main(int argc, char *argv[])
{ int d,n,i,j,copcode,np1,grsize[10],m,jg;
  double *gmu,*gar1,*gom,*galp1,*gbe1,*gnu,*gsig0;
  double *lgret0, *portfret,*cpar,*uvec,df;
  //double **dmatrix(int,int);
  void simnestfact(int *m, int *grsize, double *cpar, int *copcode, double *uvec);
  void simnestfactmvt(int *m, int *grsize, double *cpar, int *copcode, double *uvec);
  void rgarchnestfact(int *m, int *grsize, int *n, 
     double *gmu, double *gar1, double *gom,
     double *galp1, double *gbe1, double *gnu, double *gsig0,
     double *cpar, int *copcode, double  *lgret, double *portfret);

  d=6; n=3; m=3; 
  for(jg=0;jg<m;jg++) grsize[jg]=2;
  /* check on simnestfact() */
  cpar=(double *) malloc((3*d) * sizeof(double));
  uvec=(double *) malloc(d * sizeof(double));
  srand(123);
  if(argc>=2) copcode=atoi(argv[1]); else copcode=3;
  printf("\ncopcode=%d\n", copcode); 
  // cpar has length d except it is 2*d for BB1
  for(j=m;j<m+d;j++) cpar[j]=1.1+j*0.1;
  if(copcode==9) 
  { for(j=m+1;j<m+2*d;j+=2) cpar[j]=1.1+j*0.1;
    for(j=m;j<m+2*d-1;j+=2) cpar[j]=0.2;
    cpar[0]=.003;  // Gumbel for variable 1
    cpar[3]=1.003;  // MTCJ for variable 2
    np1=2;
  }
  else np1=1; 
  for(j=0;j<m;j++) cpar[j]=1.1;
  df=4.; if(copcode==1) df=1000.;
  if(copcode==1 || copcode==2)
  { cpar[d+m]=df;
    for(j=0;j<m;j++) cpar[j]=0.4+j/10.;
    for(j=m;j<d+m;j++) cpar[j]=0.9*(j-m+1.)/d;
    for(i=0;i<n;i++)
    { simnestfactmvt(&m,grsize,cpar,&copcode,uvec);
      // random normal or t
      for(j=0;j<d;j++) printf("%f ", uvec[j]); printf("\n");
    }
  }
  else
  { for(i=0;i<n;i++)
    { simnestfact(&m,grsize,cpar,&copcode,uvec);
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
  rgarchnestfact(&m,grsize,&n,gmu,gar1,gom,galp1,gbe1,gnu,gsig0,cpar,&copcode,lgret0,portfret);
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

#define BVT 2

/* inputs
     m0 = mgrp = #groups
     grsize = vector of group sizes, length mgrp; sum(mgrp)=d
       assume grsize[0] in group 1, next grsize[1] in group 2 etc.
     n0 = n = simulation sample size
       AR(1)-GARCH(1,1) parameters (mu,ar1,om,alp1,be1,nu) for each asset
     gmu = d-vector of mu or location parameters 
     gar1 = d-vector of AR1 parameters 
     gom = d-vector of omega parameters 
     galp1 = d-vector of alpha parameters 
     gbe1 = d-vector of beta parameters 
     gnu = d-vector of nu parameters for Student t innovations
     gsig0 = d-vector of starting conditional SD values (one for each asset)
     cpar = parameter vector for nested-factor copula
         cpar order grp[1],...,grp[m], then var[1],...,var[d]
     nipol = # points used for the monotone interpolation
     qq = quantiles of t(df+1), df = mvt parameter
     pp = probabilities of t(df+1)
     pder = derivatives for the interpolation
   outputs
     lgret0 = vector of length n*d of logreturns based on copula GARCH model 
              with nested-factor copula,
            lgret is  a vector for linking to R
     portfret = n-vector of portfolio returns (average of d assets)
*/
void rgarchnestfactmvt(int *m0, int *grsize, int *n0, 
   double *gmu, double *gar1, double *gom, 
   double *galp1, double *gbe1, double *gnu, double *gsig0,
   double *cpar, 
   int *nipol, double *qq, double *pp, double *pder, 
   double *lgret0, double *portfret)
{ int d,n,i,j,m,jg,j0;
  int copcode,nnew,ierr;
  double *tscale,*sigma2,*mu1,*z,*eps,*uvec,u,sportf,tem,prev,df;
  double ptem,dertem;
  double qt(double,double),sqr(double),urand();
  double pt_(double*,double*),pnorms(double);
  void simnestfactmvt(int *m, int *grsize, double *cpar, int * copcode, double *uvec);
  void pchev_(int *nipol,double *qq, double *pp, double *pder, int *nnew,
     double *xnew, double *ynew, double *dernew, int *ierr);

  n=*n0; m=*m0;
  for(jg=0,d=0;jg<m;jg++) d+=grsize[jg];
  df=cpar[m+d];
  copcode=BVT; 
#ifdef MAIN2F
  for(j=0;j<d;j++)
  { printf("%d : %f %f %f %f %f %f %f\n", j,
      gmu[j],gar1[j],gom[j],galp1[j],gbe1[j],gnu[j],gsig0[j]);  
  }
  for(j=0;j<d;j++) printf("%f ", cpar[j]); printf("\n"); 
  printf("copcode=%d\n", *copcode);
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
    simnestfactmvt(&m,grsize,cpar,&copcode,uvec);
    { for(j=0;j<d;j++) uvec[j]=pt_(&uvec[j],&df); }
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

  nnew=1;
  // generate remaining observations
  for(i=1;i<n;i++)
  { //printf("i=%d\n", i);
    simnestfactmvt(&m,grsize,cpar,&copcode,uvec);
    for(j=0;j<d;j++) 
    { //uvec[j]=pt_(&uvec[j],&df); 
      pchev_(nipol,qq,pp,pder, &nnew, &uvec[j], &ptem, &dertem, &ierr);
      uvec[j]=ptem;
    }
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

#ifdef AGAIN
/* inputs
     m0 = mgrp = #groups
     grsize = vector of group sizes, length mgrp; sum(mgrp)=d
       assume grsize[0] in group 1, next grsize[1] in group 2 etc.
     cpar = parameter vector with partial correlations with latent variables 
            and df for mvt
         cpar order grp[1],...,grp[m], then var[1],...,var[d]
   output
     zvec = random normal or t d-vector with nested-factor structure
*/
void simnestfactmvt(int *m0, int *grsize, double *cpar, int *copcode, double *zvec)
{ int cop,j,d,m,jg,ind1,ind2;
  double df,w1,*wgr,z,rho,rhe,alp2,denom,*rhpar,rhgr,rhgr1;
  double snorm(),rgamma(double,double); // in file rt.c
  cop=*copcode; m=*m0;
  for(jg=0,d=0;jg<m;jg++) d+=grsize[jg];
  wgr=(double *) malloc(m * sizeof(double));
  df=cpar[m+d];
  // convert to rhpar for bifactor model
  rhpar=(double *) malloc((2*d+1) * sizeof(double));
  ind1=0; 
  for(jg=0;jg<m;jg++)
  { ind2=ind1+grsize[jg];
    rhgr=cpar[jg]; rhgr1=sqrt(1.-rhgr*rhgr);
    for(j=ind1;j<ind2;j++)
    { rhpar[j]=cpar[m+j]*rhgr;
      rhpar[j+d]=cpar[m+j]*rhgr1;
    }
    ind1=ind2; 
  }
  //for(j=0;j<m+d; j++) printf("%f ", cpar[j]);
  //for(j=0;j<2*d; j++) printf("%f ", rhpar[j]);

  w1=snorm();  // latent variable
  for(jg=0;jg<m;jg++) wgr[jg]=snorm();  // subgroup latent variable
  ind1=0; 
  for(jg=0;jg<m;jg++)
  { ind2=ind1+grsize[jg];
    //printf("jg=%d, ind1=%d, ind2=%d\n", jg,ind1,ind2);
    for(j=ind1;j<ind2;j++)
    { z=snorm();
      rho=rhpar[j]; alp2=rhpar[d+j]*sqrt(1.-rho*rho); 
      rhe=sqrt(1.-rho*rho-alp2*alp2);
      zvec[j]=rho*w1+alp2*wgr[jg]+rhe*z;
    }
    ind1=ind2; 
  }
  //printf("\n");
  if(cop==BVT && df<300. && df>0.)  // bivariate t
  { denom=rgamma(df/2.,2.) / df; denom=sqrt(denom);
    for(j=0;j<d;j++) zvec[j]/=denom;
  }
  free(wgr); free(rhpar);
}
#endif
