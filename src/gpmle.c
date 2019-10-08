#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef MAIN
/* sample main program 
   gcc -DMAIN -o gpmle gpmle.c -lm
   gpmle < gpmle.in
   gpmle.in has sample data. 
   Linking to R: gcc -fpic -c gpmle.c ; gcc -shared -o gpmle.so gpmle.o
*/
#define N 100
main(int argc, char *argv[])
{ int i,n,iconv,maxitn;
  double y[N],av[7],xi,mu,sig;
  void gpmle(double [], int *, double *, double *, int *, int *, double []);
  if(argc==1) scanf("%d", &n);
  else { n=atoi(argv[1]); }
  for(i=0;i<n;i++) scanf("%lf", &y[i]);
  maxitn=30;
  gpmle(y,&n,&xi,&sig,&iconv,&maxitn,av);
  printf("iconv=%d, log-likelihood=%f\n", iconv, av[0]);
  printf("MLEs: xi,sigma:\n");
  printf("  %9.5f %9.5f\n", xi,sig);
  printf("Asymptotic covariance matrix:\n");
  printf("  %9.5f %9.5f\n", av[1],av[3]);
  printf("  %9.5f %9.5f\n", av[3],av[2]);
  exit(0);
}
#endif

/* maximum likelihood for 2-parameter generalized Pareto distribution */
/* input 
     y = vector that forms a decreasing sequence 
     m = sample size = dimension of y[]
     maxitn = max number of iterations
   output
     c = xi =tail index, 
     d = sigma = scale parameter, 
     iconv = convergence code (=1 for converged, =0 for not)
     av = asymptotic cov matrix stored as a vector (also log-likelihood)
     av[0]=log-lik, av[1]=var(xi), av[2]=var(sigma), av[3]=cov(xi,sigma),
*/
void gpmle(double y[], int *m,  double *cml, double *dml, int *iconv,
  int *maxitn, double av[])
{ double lcc,ldd,lcd,c,d,c1,c2,d2,xllk,xllkmax;
  double t1,t2,r,cmm,dmm,eta,eta1,dif,det,temp,tp2,s1,s2,s3;
  double xx,le,lee;
  int iter,i,mm,maxi;
  
  mm= *m; maxi= *maxitn;
  t1=0.; t2=0.;
  for(i=0;i<mm;i++)
  { xx=y[i];
    t1+=xx; t2+=xx*xx;
  }
  t1/=mm; t2/=mm;
  r=t2/(t1*t1);
  cmm=(r-2.)/(2.*(r-1.));
  /* if(*cml<5.) cmm= *cml;*/
  dmm=t1*(1.-cmm);
  eta1=cmm/dmm;
  if(1.+eta1*y[0]<=0.) eta1=.05-1./y[0];
  
  iter=0;
  xllkmax= -1.e20; dif=1.;
  // added line because of warnings
  le=0.; eta=eta1; s1=0.; s2=0.; s3=0.;
  while(iter<=maxi)
  { s1=0.; s2=0.; s3=0.;
    eta=eta1;
    for(i=0;i<mm;i++)
    { temp=1.+eta*y[i]; tp2=y[i]/temp;
      s1+=log(temp); s2+=tp2; s3+=tp2*tp2;
    }
    iter++;
    xllk= mm*log(eta*mm/s1)-s1- mm;
    if(xllk<xllkmax)
    { eta+=dif; dif/= -2.; eta1=eta-dif;
      while(1.+eta1*y[0]<=0) { eta1=(eta+eta1)/2.;}
    }
    else
    { xllkmax=xllk;
      le= mm/eta-(1.+ mm/s1)*s2;
      lee= - mm/(eta*eta) + mm*(s2*s2)/(s1*s1) +(1.+ mm/s1)*s3;
      dif=le/lee;
      eta1=eta-dif;
      /* printf("%3d %8.4f\n", iter,eta1);*/
      while(iter<=maxi && 1.+eta1*y[0]<=0)
      { dif/=2; iter++; eta1=eta-dif;}
      if(fabs(eta1-eta)<1.e-4 && iter<=maxi) { *iconv=1; goto A; }
    }
  }
  *iconv=0;
  if(fabs(le)<.01) *iconv=1;
  if(*iconv==0) return;
  
  A: c=s1/mm; *cml=c;
  d=c/eta; *dml=d;
  c2=c*c; d2=d*d; c1=1.+c;
  lcc= -2.*s1/(c2*c)+2.*s2/(d*c2)+c1*s3/(c*d2);
  ldd=(mm-2.*c1*s2/d+c*c1*s3/d2)/d2;
  lcd=(s2/d-c1*s3/d2)/d;
  det=lcc*ldd-lcd*lcd;
  av[1]= -ldd/det; av[2]= -lcc/det; av[3]= lcd/det;
  av[0]= xllkmax;
}
