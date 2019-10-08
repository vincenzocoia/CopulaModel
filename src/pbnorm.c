#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* sample main program */
// gcc -DMAINB -o pbnorm pbnorm.c R_exchmvn.c qt.c -lm
// gcc -fpic -c pbnorm.c R_exchmvn.c qt.c
// gcc -shared -o pbnorm.so pbnorm.o R_exchmvn.o qt.o
#ifdef MAINB
main(int argc, char *argv[])
{ int n,i;
  double *z1,*z2,*r,*bcdf;
  void pbnorm(int *, double *, double *, double *, double *);
  scanf("%d", &n);
  z1=(double *) malloc(n * sizeof(double));
  z2=(double *) malloc(n * sizeof(double));
  r=(double *) malloc(n * sizeof(double));
  bcdf=(double *) malloc(n * sizeof(double));
  for(i=0;i<n;i++) scanf("%lf %lf %lf", &z1[i],&z2[i],&r[i]);
  pbnorm(&n,z1,z2,r,bcdf);
  for(i=0;i<n;i++) printf("%f %f %f %e\n", z1[i],z2[i],r[i],bcdf[i]);
  free(z1); free(z2); free(r); free(bcdf); 
}
#endif

/* front end of Donnelly's code to interface with R */
/* call to exchmvn when z1 or z2 is less than -5, as
   the extreme tail is not handle properly */
/* alternative is the code in library pbivnorm/mvtnorm */
/* cdf value is -1 for any illegal inputs */
/* inputs
     n0 = dimension of z1, z2, r, bcdf
     z1 = vector of first co-ordinates
     z2 = vector of second co-ordinates
     r = vector of correlation values
   output
     bcdf = vector of bivariate cdf values
*/
void pbnorm(int *n0, double *z1, double *z2, double *r, double *bcdf)
{ int n,i,m;
  double cdf1,cdf2,x1,x2,xmin,rh;
  double tem,lb[2],ub[2],eps;
  double alnorm(double, int);
  double min(double , double );
  double max(double , double );
  double bivnor(double , double , double );
  void r_exchmvn(int *, double *, double *, double *, double *, double *);
  n=*n0;
  eps=1.e-6;
  for(i=0;i<n;i++)
  { x1=z1[i]; x2=z2[i]; rh=r[i];
    // handle boundary cases first
    if(rh< -1 || rh>1) bcdf[i]=-1;
    // lower tail isn't handled in Donnelly's algorithm
    //else if(x1<=-6. || x2<=-6.) bcdf[i]=0.; /* temporary: isn't always good */
    else if(x1>=6.) bcdf[i]=alnorm(x2,0); // approx
    else if(x2>=6.) bcdf[i]=alnorm(x1,0);
    else if(rh==0.) 
    { cdf2=alnorm(x2,0); cdf1=alnorm(x1,0); bcdf[i]=cdf1*cdf2; }
    else if(rh==1.) { xmin=min(x1,x2); bcdf[i]=alnorm(xmin,0); }
    else if(rh==-1.) 
    { if(x1+x2<=0.) bcdf[i]=0.; 
      else { cdf2=alnorm(x2,0); cdf1=alnorm(x1,0); bcdf[i]=cdf1+cdf2-1.; }
    }
    else 
    { tem=bivnor(-x1,-x2,rh);
      // use alternative exchmvn if this returns<=2.e-9
      if(tem>2.e-9) bcdf[i]=tem;
      else
      { m=2; 
        if(rh>=0.) 
        { lb[0]=min(x1-5.,-6.); ub[0]=x1;
          lb[1]=min(x2-5.,-6.); ub[1]=x2;
          r_exchmvn(&m,lb,ub,&rh,&eps,&tem);
          bcdf[i]=tem;
        }
        else // Phi2(x1,x2,rho)=Phi(x1)-Phi2(x1,-x2,-rho)
        { rh=-rh;
          lb[0]=min(x1-5.,-6.); ub[0]=x1;
          lb[1]=min(-x2-5.,-6.); ub[1]=-x2;
          r_exchmvn(&m,lb,ub,&rh,&eps,&tem);
          bcdf[i]=alnorm(x1,0)-tem;
        }
      }
    }
  }
  return;
}

double min(double a, double b)
{ return (a<b)? a:b; }

double max(double a, double b)
{ return (a<b)? b:a; }

/* Bivariate normal survival function, translated from Fortran to C: 
   Donnelly, T.G. [1973],  Algorithm 462: bivariate normal distribution, 
   Communications of the association for computing machinery, 16, 638. 
*/
/* inputs
     ah = first co-ordinate
     ak = second co-ordinate
     r = correlation in the interval (-1,1)
   output
     bivariate normal survival function
*/
double bivnor(double ah, double ak, double r)
{  double twopi,b,gh,gk,rr,con,sqr,wh,wk,h2,a2,ex,h4,w2,sn,sp,ap,cn,t;
   double gw,g2,s2,s1,conex,sgn;
   double alnorm(double, int);
   int idig,i,is;
   twopi=6.283185307179587;
   b=0.;
   idig=9;
   gh=alnorm(ah,1)/2.; gk=alnorm(ak,1)/2.;
   if(r!=0.0)  rr=1.-r*r;
   else 
   { b=4.*gh*gk; if(b<0.) b=0.; if(b>1.) b=1.; return(b); }
   if(rr<0.) { return (-1.);}
   if(rr==0.0)  
   { if(r<0.) /* r=-1 */
     { if(ah+ak<0.) b=2.*(gh+gk)-1.;
       if(b<0.) b=0.; if(b>1.) b=1.; return(b); 
     }
     else /* r=1 */
     {  if(ah-ak<0.) b=2.*gk;
        else b=2.*gh;
        if(b<0.) b=0.; if(b>1.) b=1.; return(b);
     } 
   }
   /* rr != 0 , r!=0 */
   sqr=sqrt(rr);
   con=twopi*.5;
   for(i=1;i<=idig;i++) con/=10.;
   if(ah==0.0) 
   { if(ak==0.) 
     { b=atan(r/sqr)/twopi+.25; if(b<0.) b=0.; if(b>1.) b=1.; return(b); }
     else
     {  b+=gk;
        if(ah!=0.) { wh=-ah; wk=(ak/ah-r)/sqr; gw=2.*gh; is=-1;}
        else { wh=-ak; wk=(ah/ak-r)/sqr; gw=2.*gk; is=1;}
        goto L210;
     }
   }
   
      b=gh;
      if(ah*ak<0) { b-=.5;} 
      if(ah*ak!=0.) 
      { b+=gk;
        if(ah!=0.) { wh=-ah; wk=(ak/ah-r)/sqr; gw=2.*gh; is=-1;}
        else { wh=-ak; wk=(ah/ak-r)/sqr; gw=2.*gk; is=1;}
      }
      else { wh=-ah; wk=(ak/ah-r)/sqr; gw=2.*gh; is=-1;}
 L210: sgn=-1.; t=0.;
      if(wk!=0.) 
      { if(fabs(wk)==1.) { t=wk*gw*(1.-gw)*.5; goto L310;}
        else 
        { if(fabs(wk)>1.)
	  { sgn=-sgn; wh=wh*wk; g2=alnorm(wh,0); wk=1./wk;
	    if (wk<0) { b=b+.5;}
	    b=b-(gw+g2)*.5+ gw*g2;
	  }
	  h2=wh*wh; a2=wk*wk; h4=h2*.5; ex=0.0; if(h4<87.0) ex=exp(-h4);
	  w2=h4*ex; ap=1.; s2=ap-ex; sp=ap; s1=0.; sn=s1; conex=fabs(con/wk);
	  goto L290;
	}
      }
      else 
      {  goto L320; }

 L290: cn=ap*s2/(sn+sp); s1+=cn;
      while(fabs(cn)-conex>0) 
      {   sn=sp; sp+=1.; s2-=w2;
          if(fabs(w2)<=1.0e-15 || fabs(h4)<=1.0e-15) w2=0.0;
	  else { w2*=h4/sp;}
	  /*    underflow prevention   */
	  if(fabs(ap)<=1.0e-15 || fabs(a2)<=1.0e-15) { ap=0.0;}
	  else { ap=-ap*a2; }
          //goto L290;
	  cn=ap*s2/(sn+sp); s1+=cn;
      }
      //else
      t=(atan(wk)-wk*s1)/twopi; //goto L310;}
 L310: b+=sgn*t;
 L320: if(is<0) 
       {  if(ak!=0.) 
          {  wh=-ak; wk=(ah/ak-r)/sqr; gw=2.*gk; is=1; goto L210; }
          else { if(b<0.) b=0.; if(b>1.) b=1.; return(b); }
       }
       else { if(b<0.) b=0.; if(b>1.) b=1.; return(b); }
    
}

/* algorithm Applied Statistics 66 by I.D. Hill */
// survival function if upper=1, cdf if upper=0
/* inputs
     x = real value
     upper = flag for upper tail (1 for upper, 0 for lower)
   output
     lower or upper tail probability of standard normal
*/
double alnorm(double x, int upper)
{  int up;
   double ltone,utzero,con,prob,z,y;
   con=1.28; ltone=5.0; utzero=12.5;
   up=upper; z=x;
   if(z<0.0) { up=1-up; z=-z; }
   if(z<=ltone || (up==1 && z<=utzero)) 
   { y=0.5*z*z;
     if(z<=con) 
     { prob=0.5-z*(0.398942280444-0.399903438504*y/ 
         (y+5.75885480458-29.8213557808/ 
         (y+2.62433121679+48.6959930692/ 
         (y+5.92885724438))));
     }
     else
     { prob=0.398942280385*exp(-y)/ 
         (z-3.8052e-8+1.00000615302/ 
         (z+3.98064794e-4+1.98615381364/ 
         (z-0.151679116635+5.29330324926/ 
         (z+4.8385912808-15.1508972451/ 
         (z+0.742380924027+30.789933034/(z+3.99019417011))))));
     }
   }
   else { prob = 0.0; }
   if(up==0) prob=1.0-prob;
   return(prob);
}
