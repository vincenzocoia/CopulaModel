#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef MAIN
/* sample main program 
   gcc -DMAIN -o gevmle gevmle.c -lm
   gevmle < gevmle.in
   gevmle.in is sample data, 
   first entry is number of data values,
   (or sample size specified on command line)
   remainder are data values
   Linking to R: gcc -fpic -c gevmle.c ; gcc -shared -o gevmle.so gevmle.o
*/
#define N 100
main(int argc, char *argv[])
{ int i,n,iconv,maxitn;
  double x[N],av[7],xi,mu,sig;
  void gevmle(double [], int *, double *, double *, double *, int *, int * , double []);
  if(argc==1) scanf("%d", &n);
  else { n=atoi(argv[1]); }
  
  for(i=0;i<n;i++) scanf("%lf", &x[i]);
  maxitn=30;
  gevmle(x,&n,&xi,&sig,&mu,&iconv,&maxitn,av);
  printf("iconv=%d, log-likelihood=%f\n", iconv, av[0]);
  printf("MLEs: xi,sigma,mu:\n");
  printf("  %9.5f %9.5f %9.5f\n", xi,sig,mu);
  printf("Asymptotic covariance matrix:\n");
  printf("  %9.5f %9.5f %9.5f\n", av[1],av[4],av[5]);
  printf("  %9.5f %9.5f %9.5f\n", av[4],av[2],av[6]);
  printf("  %9.5f %9.5f %9.5f\n", av[5],av[6],av[3]);
  exit(0);
}
#endif

/* maximum likelihood for 3-parameter generalized extreme value distribution */
/* inputs
     x[] = data vector 
     mptr = m = sample size = dimension of x
     maxitn = max number of iterations
   outputs
     c = xi= tail index, 
     a = sigma = scale parameter, 
     b = mu = location parameter 
     iconv = convergence code (=1 for converged, =0 for not)
     av = asymptotic cov matrix stored as a vector (also log-likelihood)
     av[0]=log-lik, av[1]=var(xi), av[2]=var(sigma), av[3]=var(mu),
     av[4]=cov(xi,sigma), av[5]=cov(xi,mu), av[6]=cov(sigma,mu)
*/
void gevmle(double *x, int *mptr, double *c, double *a, double *b, int *iconv,
  int *maxitn, double av[])
{ double ly,la,lb,lc,laa,lbb,lcc1,lcc2,lcc,lab,lac,lbc;
  double tol,y,yc,c02,c1,a02;
  double z,yi,yc1,yyc1,yyi,cdif,adif,bdif,xlk,xlkmax;
  double det,wab,wac,wbc,cn,an,bn,c0,a0,b0;
  double t1,t2;
  int i,iter,m;
  double gevllk(double *, int, double, double, double);
  
  m= *mptr;
  tol=5.e-4;
  iter=0;
  for(i=0,t1=0.,t2=0.;i<m;i++) { t1+=x[i]; t2+=x[i]*x[i]; }
  t1/=m; t2/=m;
  a0=sqrt(t2-t1*t1)/1.28255;
  b0=t1-a0*.57722;
  xlkmax=-1.e10; c0=-.03;
  for(cn=-.52; cn<=.53; cn+=.05)
  { xlk=gevllk(x,m,cn,a0,b0);
    if(xlk>xlkmax) { xlkmax=xlk; c0=cn; }
  }
  lb=0.; la=-m; lc=0.; laa=0.; lbb=0.; lcc1=0.; lcc2=0.;
  lab=0.; lac=0.; lbc=0.;
  // added line because of warnings
  lcc=0.; det=1.; wab=0.; wac=0.; wbc=0.;
  /*printf("%10.5f %10.5f %10.5f %10.5f\n", c0,a0,b0,xlkmax);*/
  while(iter< *maxitn)
  { lb=0.; la=-m; lc=0.; laa=0.; lbb=0.; lcc1=0.; lcc2=0.;
    lab=0.; lac=0.; lbc=0.;
    c02=c0*c0; c1=1.+c0;
    iter++;
    for(i=0;i<m;i++)
    { y=(x[i]-b0)/a0; z=1.+c0*y; ly=log(z);
      yi=1./z; yc=exp(-ly/c0); yc1=yc*yi; yyi=y*yi; yyc1=y*yc1;
      lc+=ly*(1.-yc)+c0*(yyc1-c1*yyi);
      la+=c1*yyi-yyc1;
      lb+=c1*yi-yc1;
      laa+=c1*yyi*(-1.+c0*yyi-yyc1)+yyc1;
      lbb+=c1*yi*(c0*yi-yc1);
      lab+=c1*yi*(c0*yyi-yyc1);
      lac+=yyi+c1*yyi*(yyc1/c0-yyi)-yyc1*ly/c02;
      lbc+=yi+c1*yi*(yyc1/c0-yyi)-yc1*ly/c02;
      lcc1+=2.*(yyi-yyc1)+c1*yyi*(c0*yyi-yyc1);
      lcc2+=2.*ly*(yyc1+yc-1.)-ly*ly*yc/c0;
    }
    a02=a0*a0; lc/=c02; la/=a0; lb/=a0;
    laa=laa/a02-la/a0; lbb/=a02; lab=lab/a02-lb/a0;
    lac/=a0; lbc/=a0; lcc=(lcc1+lcc2/c0)/c02;
    det=lcc*laa*lbb+2.*lac*lab*lbc;
    det=det-lbc*lbc*laa-lac*lac*lbb-lab*lab*lcc;
    wac=lab*lbc-lbb*lac;
    wbc=lac*lab-laa*lbc;
    wab=lac*lbc-lcc*lab;
    cdif=((laa*lbb-lab*lab)*lc+wac*la+wbc*lb)/det;
    adif=((lcc*lbb-lbc*lbc)*la+wac*lc+wab*lb)/det;
    bdif=((lcc*laa-lac*lac)*lb+wbc*lc+wab*la)/det;
    cn=c0-cdif; an=a0-adif; bn=b0-bdif;
    xlk= gevllk(x,m,cn,an,bn);
    /*printf("%5d %10.5f %10.5f %10.5f %10.5f %10.5f\n",
         iter,c0,a0,b0,xlkmax,xlk);*/
    if(xlk<xlkmax)  /* move in direction of gradient */
    { c0=cn+cdif; a0=an+adif; b0=bn+bdif;
      cdif=fabs(cdif)*(2.*(lc>0.)-1.);
      adif=fabs(adif)*(2.*(la>0.)-1.);
      bdif=fabs(bdif)*(2.*(lb>0.)-1.);
      do
      { cdif/=2.; adif/=2.; bdif/=2.;
        iter++;
        if(iter>= *maxitn) { *iconv= -1; goto A; }
        cn=c0+cdif; an=a0+adif; bn=b0+bdif;
        xlk=gevllk(x,m,cn,an,bn);
      }
      while(xlk<xlkmax);
    }
    xlkmax=xlk;
    if(fabs(cn-c0)<tol && fabs(an-a0)<tol && fabs(bn-b0)<tol)
    { *iconv=1; goto A; }
    c0=cn; a0=an; b0=bn;
  }
  *iconv=0;
  if(fabs(la)<1.e-2 && fabs(lb)<1.e-2 && fabs(lc)<1.e-2) *iconv=1;
  A: *c=c0; *a=a0; *b=b0;
  det*=-1.;
  av[1]=(laa*lbb-lab*lab)/det;
  av[2]=(lcc*lbb-lbc*lbc)/det;
  av[3]=(lcc*laa-lac*lac)/det;
  av[4]=wac/det;
  av[5]=wbc/det;
  av[6]=wab/det;
  av[0]=xlkmax;
}

/* inputs
     x[] = data vector 
     m = sample size = dimension of x 
     c = xi= tail index, 
     a = sigma = scale parameter, 
     b = mu = location parameter 
   output
     GEV log-likelihood
*/
double gevllk(double *x, int m, double c, double a, double b)
{ double xlik,z,ly,yc,c1;
  int i;
  xlik=-1.e10;
  if(a<=0.) return(xlik);
  c1=1.+c;
  xlik=-m*log(a);
  for(i=0;i<m;i++)
  { z=1.+c*(x[i]-b)/a;
    if(z<=0.) { xlik= -1.e10; return(xlik); }
    ly=log(z);
    yc=exp(-ly/c);
    xlik-=c1*ly/c+yc;
  }
  return(xlik);
}
