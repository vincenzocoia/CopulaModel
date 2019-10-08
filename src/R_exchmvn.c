#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define UB 6. // default upper limit for integration of common factor
//#define EPS 1.e-7
#define S 12  // max number of splits for Romberg extrapolation

int mm,kk,ksign;
double *ww,*xx,rs,r1,r32;

/* gcc -DMAIN2 -o R_exchmvn R_exchmvn.c qt.c -lm */
/* mvn rectangle probability and derivatives for positive exch case, 
   with Romberg integration */
/* version with pointers for link to R, zero indexes are used */
#ifdef MAIN2
main()
{ int m,i,k,ks;
  double rh,*a,*b,pr,eps;
  void r_exchmvn(int *, double *, double *, double *, double *, double *);
  void r_emvnd(int *, double *,double *,double *, int *,int *, double *,double *);
  void r_emvndrh(int *, double *, double *, double *, double *, double *);
  double heps,pr2,dera,derb,drh;
  
  eps=1.e-6;
  heps=1.e-4;
  scanf("%d", &m);
  while(m>0)
  { scanf("%lf", &rh);
    a=(double *) malloc(m * sizeof(double));
    b=(double *) malloc(m * sizeof(double));
    for(i=0;i<m;i++)  scanf("%lf", &a[i]);
    for(i=0;i<m;i++)  scanf("%lf", &b[i]);
    printf("m=%3d, rh=%6.2f\n", m,rh);
    for(i=0;i<m;i++) printf("%8.4f", a[i]);  printf("\n");
    for(i=0;i<m;i++) printf("%8.4f", b[i]);  printf("\n");
    r_exchmvn(&m,a,b,&rh,&eps,&pr);
    printf("exch.  : %9.5f\n", pr);
    a[0]+=heps;
    r_exchmvn(&m,a,b,&rh,&eps,&pr2);
    dera=(pr2-pr)/heps;
    a[0]-=heps;
    b[0]+=heps;
    r_exchmvn(&m,a,b,&rh,&eps,&pr2);
    derb=(pr2-pr)/heps;
    b[0]-=heps;
    printf("num deriv: dera1=%f, derb1=%f\n", dera,derb);
    k=1; ks=-1;
    r_emvnd(&m,a,b,&rh,&k,&ks,&eps,&dera);
    k=1; ks=1;
    r_emvnd(&m,a,b,&rh,&k,&ks,&eps,&derb);
    printf("integ:     dera1=%f, derb1=%f\n", dera,derb);
    
    a[m-1]+=heps;
    r_exchmvn(&m,a,b,&rh,&eps,&pr2);
    dera=(pr2-pr)/heps;
    a[m-1]-=heps;
    b[m-1]+=heps;
    r_exchmvn(&m,a,b,&rh,&eps,&pr2);
    derb=(pr2-pr)/heps;
    b[m-1]-=heps;
    printf("num deriv: dera1=%f, derb1=%f\n", dera,derb);
    k=m; ks=-1;
    r_emvnd(&m,a,b,&rh,&k,&ks,&eps,&dera);
    k=m; ks=1;
    r_emvnd(&m,a,b,&rh,&k,&ks,&eps,&derb);
    printf("integ:     dera1=%f, derb1=%f\n", dera,derb);

    rh+=heps;
    r_exchmvn(&m,a,b,&rh,&eps,&pr2);
    drh=(pr2-pr)/heps;
    rh-=heps;
    printf("num deriv: derrh=%f\n", drh);
    r_emvndrh(&m,a,b,&rh,&eps,&drh);
    printf("integ.   : derrh=%f\n", drh);
    free(a); free(b);
    scanf("%d", &m);
  }
}
#endif

/* gcc -DMAIN1 -o R_exchmvn R_exchmvn.c qt.c -lm */
/* mvn rectangle probability for positive exch case, */
/* version with pointers for link to R, zero indexes are used */
#ifdef MAIN1
main()
{ int m,i;
  double rh,*x,*w,pr,eps;
  void r_exchmvn(int *, double *, double *, double *, double *, double *);
  
  eps=1.e-6;
  scanf("%d", &m);
  while(m>0)
  { scanf("%lf", &rh);
    x=(double *) malloc(m * sizeof(double));
    w=(double *) malloc(m * sizeof(double));
    for(i=0;i<m;i++)  scanf("%lf", &w[i]);
    for(i=0;i<m;i++)  scanf("%lf", &x[i]);
    printf("m=%3d, rh=%f\n", m,rh);
    for(i=0;i<m;i++) printf("%8.4f", w[i]);  printf("\n");
    for(i=0;i<m;i++) printf("%8.4f", x[i]);  printf("\n");
    r_exchmvn(&m,w,x,&rh,&eps, &pr);
    printf("exch.  : %.10f\n", pr);
    free(x); free(w);
    scanf("%d", &m);
  }
}
#endif

/* version with all pointers and zero indices for interface to R directly */
/* inputs
     m = dimension
     w = m-vector of lower limits
     x = m-vector of upper limits
     rh = positive common correlation parameter
     eps = tolerance for convergence for numerical integration
   output
     pr = rectangle probability
*/
void r_exchmvn(int *m, double *w, double *x, double *rh, double *eps, double *pr)
{ double r_g(double);
  double romberg(double (*)(double), double, double, double);
  int i;
  extern int mm;
  extern double *ww,*xx,rs,r1;
  mm=*m; rs=sqrt(*rh); r1=sqrt(1.-(*rh));
  xx=(double *) malloc(mm * sizeof(double));
  ww=(double *) malloc(mm * sizeof(double));
  for(i=0;i<mm;i++) { ww[i]=w[i]; xx[i]=x[i];}
  *pr=romberg(r_g,-UB,UB,*eps);
  free(xx); free(ww);
}

/* integrand for rectangle probability */
/* z = real value */
double r_g(double z)
{ double pnorms(double),dnorms(double),a,b;
  extern int mm;
  extern double *ww,*xx,rs,r1;
  int i;
  double tem;
  for(i=0,tem=1.;i<mm;i++)
  { a=(ww[i]-rs*z)/r1; b=(xx[i]-rs*z)/r1;
    tem*=pnorms(b)-pnorms(a);
  }
  tem*=dnorms(z);
  return(tem);
}

/* rectangle probability = P(Z_j\in (a_j,b_j)): 
   derivative with respect to a_k or b_k,
   ks=-1 for a_k, ks=1 for b_k
*/
/* inputs
     m = dimension
     w = m-vector of lower limits
     x = m-vector of upper limits
     rh = positive common correlation parameter
     k = integer between 1 and m
     ks = -1 or 1, -1 for derivative wrt lower limit, 1 for upper limit
     eps = tolerance for convergence for numerical integration
   output
     deriv = derivative of rectangle probability wrt a_k or b_k
*/
void r_emvnd(int *m, double *w, double *x, double *rh, int *k, int *ks, 
  double *eps, double *deriv)
{ double r_gd(double),der;
  double romberg(double (*)(double), double, double, double);
  int i;
  extern int mm,kk,ksign;
  extern double *ww,*xx,rs,r1;
  mm=*m; kk=(*k-1); rs=sqrt(*rh); r1=sqrt(1.-(*rh)); ksign=*ks;
  xx=(double *) malloc(mm * sizeof(double));
  ww=(double *) malloc(mm * sizeof(double));
  for(i=0;i<mm;i++) { ww[i]=w[i]; xx[i]=x[i]; }
  der=romberg(r_gd,-UB,UB,*eps);
  free(xx); free(ww);
  *deriv= ksign*der;
}

/* integrand for emvnd */
double r_gd(double z)
{ double pnorms(double),dnorms(double),a,b;
  extern int mm,kk;
  extern double *ww,*xx,rs,r1;
  int i;
  double tem;
  for(i=0,tem=1.;i<mm;i++)
  { if(i!=kk)
    { a=(ww[i]-rs*z)/r1; b=(xx[i]-rs*z)/r1;
      tem*=pnorms(b)-pnorms(a);
    }
    else if(ksign==-1)
    { a=(ww[i]-rs*z)/r1; tem*=dnorms(a)/r1; }
    else
    { b=(xx[i]-rs*z)/r1; tem*=dnorms(b)/r1; }
  }
  tem*=dnorms(z);
  return(tem);
}


/* P(Z_j\in (a_j,b_j)): derivative with respect to rho */
/* inputs
     m = dimension
     w = m-vector of lower limits
     x = m-vector of upper limits
     rh = positive common correlation parameter
     eps = tolerance for convergence for numerical integration
   output
     deriv = derivative of rectangle probability wrt rho
*/
void r_emvndrh(int *m, double *w, double *x, double *rh, double *eps, 
  double *deriv)
{ double r_grh(double),der,tem,sum;
  double romberg(double (*)(double), double, double, double);
  double pnorms(double),dnorms(double);
  int i,k;
  double *t;
  extern int mm;
  extern double *ww,*xx,rs,r1,r32;
  mm=*m; rs=sqrt(*rh); r1=sqrt(1.-(*rh)); r32=r1*(1.-(*rh));
  xx=(double *) malloc(mm * sizeof(double));
  ww=(double *) malloc(mm * sizeof(double));
  t=(double *) malloc(mm * sizeof(double));
  for(i=0;i<mm;i++) { ww[i]=w[i]; xx[i]=x[i]; }
  if((*rh)>=0.) der=romberg(r_grh,-UB,UB,*eps);
  else /* rho=0 */
  { for(i=0;i<mm;i++) t[i]=pnorms(x[i])-pnorms(w[i]);
    for(k=0,sum=0.;k<mm;k++)
    { for(i=0,tem=1.;i<mm;i++)
      { if(i!=k) { tem*=t[i]; }
        else tem*=(x[i]*dnorms(x[i])-w[i]*dnorms(w[i]));
        // maybe check if x>10 or w<-10?
      }
      sum+=tem;
    }
    der=.5*sum;
  }
  free(xx); free(ww); free(t);
  *deriv=der;
}

/* integrand for emvndrh */
double r_grh(double z)
{ double pnorms(double),dnorms(double),a,b;
  extern int mm;
  extern double *ww,*xx,rs,r1,r32;
  int i,k;
  double tem,sum,tem2,*t;

  t=(double *) malloc(mm * sizeof(double));
  for(i=0;i<mm;i++) 
  { a=(ww[i]-rs*z)/r1; b=(xx[i]-rs*z)/r1;
    t[i]=pnorms(b)-pnorms(a);
  }
  for(k=0,sum=0.;k<mm;k++)
  { for(i=0,tem=1.;i<mm;i++)
    { if(i!=k) { tem*=t[i]; }
      else
      { a=(ww[i]-rs*z)/r1; b=(xx[i]-rs*z)/r1;
        tem2=dnorms(b)*(xx[i]-z/rs)-dnorms(a)*(ww[i]-z/rs);
        tem2*=.5/r32;
        tem*=tem2;
      }
    }
    sum+=tem;
  }
  tem=sum*dnorms(z);
  free(t);
  return(tem);
}

/*
double phi2(double z1, double z2, double rh)
{ double r1,tem;
  r1=1.-rh*rh;
  tem=(z1*z1+z2*z2-2.*rh*z1*z2)/r1;
  return( 0.1591549430918953*exp(-.5*tem)/sqrt(r1));
}*/

/* Romberg integration (one-dimensional integral) with relative convergence
   criterion, this integration method should be mentioned in any
   text on numerical methods; see also Numerical Recipes.
   Maybe can use the version which return a code for convergence.
*/
/* inputs
     g = integrand
     a = lower limit (finite)
     b = upper limit (finite)
     eps = tolerance for convergence for numerical integration
   output
     integral of g over (a,b)
*/
double romberg(double (*g)(double), double a, double b, double eps)
{ double t[S+1][S+1],h,fourj,sum,integ;
/* the 0th row and column of t[][] are not used */
  int m,k,i,j;
  h=b-a; m=1;
  t[1][1]=h*((*g)(a)+(*g)(b))/2.;
  for(k=2;k<=S;k++)
  { h/=2.; m*=2;
    /* for(i=1,sum=0.;i<=m;i+=2) sum+= g(a+i*h);*/
    for(i=1,sum=0.;i<=m;i+=2) sum+= (*g)(a+i*h);
    t[k][1]=t[k-1][1]*.5+sum*h;
    for(j=2,fourj=1.;j<=k;j++)
    { fourj*=4.;
      t[k][j]=t[k][j-1] + (t[k][j-1]-t[k-1][j-1])/(fourj-1.);
    }
    if(fabs((t[k][k]-t[k-1][k-1])/t[k][k]) <= eps)
    { integ=t[k][k]; return(integ);}
  }
  integ=t[S][S];
#ifdef DIAG 
  printf("*** convergence not reached ***\n");
#endif
  return(integ);
}
#undef UB 
//#undef EPS 
#undef S
