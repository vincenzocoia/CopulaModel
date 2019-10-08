#include <stdio.h>
#include <math.h>
/* Schervish, M.J. [1984]. Multivariate normal probabilities 
   with error bound. Appl. Statist. v 33, 81-94. */

/* f2() changed to ff(), in order to compile as a shared object
   for linking to R */

static double bcs[6]= {0.0, 1.0, 4.0, 6.0, 4.0, 1.0};
static double bcn[8]= {0.0,1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0};
static double ddo[26]= {0.0, -3.75043972, -3.3242574, -2.85697,
   -2.36675941, -2.3344142, -1.8891759, -1.7320508,
   -1.3556262,-1.15440539,-1.0,-0.74196378,-0.61670659,
   0.0, 0.61670659, 0.74196378, 1.0, 1.15440539,
   1.3556262, 1.7320508, 1.8891759, 2.3344142,
   2.36675941, 2.85697, 3.3242574, 3.75043972};
static double eo[14]= {0.0, - 2.85697, -2.3344142, -1.7320508, -1.3556262, 
   -1.0, -0.74196378, 0.0, 0.74196378, 1.0, 1.3556262, 1.7320508, 
   2.3344142, 2.85697};
static double coef[6][4]={0.0, 0.0, 0.0, 0.0,
        0.0, 0.311111111111111, 0.333333333333333, 0.5,
        0.0, 1.422222222222222, 0.0, 0.0,
        0.0, 0.533333333333333, 1.333333333333333, 0.0, 
        0.0, 1.422222222222222, 0.0, 0.0, 
        0.0, 0.311111111111111, 0.333333333333333, 0.5 };

/* inputs
     ub = vector of upper limit
     lb = vector of lower limit
     sig0 = vector from lower triangle of correlation matrix, by rows
     eps = tolerance for accuracy
     n = dimension of ub and lb
     inf0(i) = 0 if ith range is (lb(i),oo)
     inf0(i) = 1 if ith range is (-oo,ub(i))
     inf0(i) = 2 if ith range is (lb(i),ub(i))
   outputs
     prob = rectangle probability
     bound = error bound
     ifault = error code
              4 is sig0 is not positive definite
              1,2,3,5 : see Schervish's paper
              0 OK
*/
void mulnor(double ub[], double lb[], double sig0[], double eps, int n, 
   int inf0[], double *prob, double *bound, int *ifault)
{   double a[8], b[8], sig[22];
    double c[8], d[8], co[26], sd[3];
    double binc[8][6], bl[8][6], br[8][36][6], r[37][6], s[6][28],
        xss[7][7], pr2[7], prep[7], preb[7], cv[6][16], sinv[29],
        condl[6], xm[6], cond[6], beta[6][6], bh[6], fe[6], sigma[29],
        ep[6], del[4][6], bou4[6], bou5[6], ans[7], fact[6], prod[6], bint[7]; 
    /*bcs[6], bcn[8], eo[14], do[26];*/
    int inf[8];
    int intvl[6], ind[6], ksa[6], num[6], itype[6];
    int simps,ipr;
    double bou1, bou2, bou3, cup, det,
     eplos, epsi, ept, fac, fsa, rho, sigc, t, tem, temb, temp, wb,
     wl, wt, wu, x1, x2, xs, y1, y2, z;
    double alnorm(double, int); 
    double bivnor(double, double, double); 
    double ff(double, double); 
    double f3(double, double); 
    double f4(double, double);
    double f5(double, double); 
    double f6(double, double); 
    double gauinv(double, int *); 
    double phistar(double, double); 
    double max(double, double); 
    double min(double, double); 
    static double cons2=6.879833, cons3=4.517004, cons4=76.371214,
        epsmin=1.0e-8, epssim=6.0e-5;
    void invert(double [], double [], double [][16], int , double *, int *);
    int i,j,k,ij,ik,jk,l,jn,nm1,nm2,inj,mjs,ll,iend,igo,indl,ien,ijk,ni;
    int ity,njs,nl,numl,intl;

    /* initialization to avoid warnings */
    simps=0.; rho=0.; sd[1]=0.; sd[2]=0.; mjs=0;
    /* redefine new arrays for consistency with previous programs
       calling Fortran routine
    */
    for(i=1;i<=n;i++)
    {  inf[i]=inf0[i-1]; a[i]=ub[i-1]; b[i]=lb[i-1];}
    for(j=1;j<=(n*(n-1))/2;j++) { sig[j]=sig0[j-1];}

    *ifault = 0;
    if(eps <= epsmin) *ifault = 2;
    if(n <= 0 || n > 7) *ifault = 3;
    if(*ifault != 0) return;
    for(i=1;i<= n;i++)
    { if(inf[i] == 2 && a[i] < b[i]) *ifault = 100 + i;}
    if(*ifault != 0) return;
    *bound = eps;
    ept = eps * 0.15 / ((double) n);
    z = -gauinv(ept, ifault) + epsmin;
    if(*ifault != 0) return;
    cup = alnorm(z, 1);
    ik = 0; ij = 0;
    for(i=1;i<=n;i++)
    {  for(j=1;j<=i;j++)
       { ik++;
         if(i == j) sigma[ik] = 1.;
         else { ij++; sigma[ik] = sig[ij];}
       }
    }

    if(n > 2) 
    {  invert(sigma, sinv, cv, n, &det, ifault);
       if(*ifault != 0) return;
       simps = 1;
       if(det < 0.05 || eps <= epssim) simps = 0;
       *prob = 0.;
       det = sinv[1] * sinv[3] - sinv[2] * sinv[2];
       sd[1] = sqrt(sinv[3] / det);
       sd[2] = sqrt(sinv[1] / det);
       rho = -sinv[2] / (sd[1] * sd[2] * det);
       if(fabs(rho) > 1.) { *ifault = 4; return;}
    }
    nm2 = n - 2;
    nm1 = n - 1;
    eplos = 0.;
    for(l=1;l<= n;l++)
    {  c[l] = max(b[l], -z); d[l] = min(a[l], z);
       if(inf[l] == 0) d[l] = z;
       if(inf[l] == 1) c[l] = -z;
       if(a[l] > z || inf[l] == 0) eplos += cup;
       if(b[l] < -z || inf[l] == 1) eplos += cup;
       if(c[l] >= d[l]) return;
    }
    if(n == 1) { *prob = alnorm(d[1], 0) - alnorm(c[1], 0); return;}
    fac = 1.;
    ipr = 0;
    if(inf[1] == 1 && inf[2] == 1) 
    { ipr = 1; eplos -= 2. * cup; }
    else if(inf[1] == 0 && inf[2] == 0) 
    {  fac= -1.; ipr=1; d[1]=c[1]; d[2]=c[2]; eplos -= 2. * cup; }
    if(n == 2) 
    {  rho = sigma[2];
       if(fabs(rho) > 1.) { *ifault = 4; return;}
       y1 = -d[1] * fac; y2 = -d[2] * fac;
       *prob = bivnor(y1, y2, rho);
       if(ipr==1) return;
       x1 = -c[1]; x2 = -c[2];
       wl = bivnor(x1, x2, rho); wt = bivnor(x1, y2, rho);
       wb = bivnor(y1, x2, rho); 
       *prob = wl - wt - wb + *prob;
       return;
    }

    /* n>=3 */
    *ifault = 5;
    epsi = (eps - eplos) / ((double) nm2);

    for(l=1;l<= nm2;l++)
    {  cond[l] = 1. / sqrt(cv[l][1]);
       condl[l] = log(cond[l]);
       for(i=1;i<=l;i++)
       { beta[l][i] = 0.;
         for(j=1;j<=l;j++)
         {  jk = (l - i + 1) * (l - i) / 2 + j;
            if(j > l - i + 1) jk = j * (j - 1) / 2 + l - i + 1;
            jn = (n - l + j) * (n - l + j - 1) / 2 + n - l;
            beta[l][i] += sigma[jn] * cv[l][jk];
         }
       }

       k = n - l - 1;
       bh[k + 1] = beta[l][l];
       for(i=1;i<=k;i++)
       {  bh[i] = 0.;
          for(j=1;j<=l;j++)
          {  jn = (j + n - l) * (j + n - l - 1) / 2 + n - l;
             ijk = 1 + j * (j - 1) / 2;
             bh[i] += sigma[jn] * cv[l][ijk];
          }
       }
       k = 0; sigc = 0.;
       for(j=1;j<=l;j++)
       {  for(i=1;i<=j;i++)
          {  k++; sigc += bh[i] * bh[j] * sinv[k]; }
       }
       binc[1][l] = 1.;
       binc[2][l] = sqrt(sigc);
       binc[3][l] = 2. * sigc;
       binc[4][l] = cons2 * sigc * binc[2][l];
       binc[5][l] = 12. * sigc * sigc;
       if(simps==0) 
       {  binc[6][l] = sigc * sigc * binc[1][l] * cons3;
          binc[7][l] = pow(sigc, 3.) * cons4;
       }
       if(l >= nm2) 
       {  for(i=1;i<=nm2;i++)
          {  bh[i] = 0.;
             for(j=1;j<=nm2;j++)
             {  jk = (l - i + 1) * (l - i) / 2 + j;
                if(j > l - i + 1) jk = j * (j - 1) / 2 + l - i + 1;
                jn = (2 + j) * (1 + j) / 2 + 1;
                bh[i] += sigma[jn] * cv[l][jk];
             }
          }
       }
    }

    l = 1;
    if(simps==0) 
    {  for(i=1;i<= 25;i++) co[i] = ddo[i];
       iend = 25; ien = 7;
    }
    else 
    {  for(i=1;i<= 13;i++) co[i] = eo[i];
       iend = 13; ien = 5;
    }
    for(i=1;i<= nm1;i++)  xss[i][1] = 0.;
    xm[1] = 0.; prod[1] = 1.; pr2[1] = 1.;
    for(i=1;i<= nm2;i++)
    {  ni = n - i + 1; pr2[i + 1] = pr2[i] * (d[ni] - c[ni]); }

    bint[nm1] = 0.;

L90: intvl[l] = 2;
    ans[l] = 0.;
    bou4[l] = 0.; bint[l] = 0.; prep[l] = 0.; preb[l] = 0.;
    k = 1; nl = n - l + 1;
    s[l][1] = c[nl] - xm[l];
    s[l][iend+2] = d[nl] - xm[l];
    num[l] = iend + 2;
    for(i=1;i<=iend;i++)
    {  njs = i;
       if(s[l][1] < co[i] * cond[l]) break;
       num[l]--;
    }
    if(num[l] != 2) 
    {  for(i=njs;i<=iend;i++)
       {  mjs = iend - i + njs;
          if(s[l][iend + 2] >= co[mjs] * cond[l]) break;
          num[l]--;
       }
       if(num[l] != 2) 
       {  for(i=njs;i<=mjs;i++)
          { inj = i - njs + 2;   s[l][inj] = co[i] * cond[l]; }
       }
    }
    numl = num[l];
    s[l][numl] = s[l][iend + 2];
    ep[l] = epsi / pr2[l + 1];
    r[1][l] = s[l][2];
    ind[l] = 6;
    fe[l] = s[l][1];
    t = fe[l] / cond[l];
    bl[1][l] = phistar(t, condl[l]);
    bl[2][l] = bl[1][l] * fabs(t / cond[l]);
    bl[3][l] = bl[1][l] * ff(t, cond[l]);
    bl[4][l] = bl[1][l] * f3(t, cond[l]);
    bl[5][l] = bl[1][l] * f4(t, cond[l]);
    if(simps==0) 
    {  bl[6][l] = bl[1][l] * f5(t, cond[l]);
       bl[7][l] = bl[1][l] * f6(t, cond[l]);
    }

L100: t = r[k][l] / cond[l];
    br[1][k][l] = phistar(t, condl[l]);
    br[2][k][l] = br[1][k][l] * fabs(t / cond[l]);
    br[3][k][l] = br[1][k][l] * ff(t, cond[l]);
    br[4][k][l] = br[1][k][l] * f3(t, cond[l]);
    br[5][k][l] = br[1][k][l] * f4(t, cond[l]);
    if(simps==0) 
    {  br[6][k][l] = br[1][k][l] * f5(t, cond[l]);
       br[7][k][l] = br[1][k][l] * f6(t, cond[l]);
    }

L104: r[k + 1][l] = (fe[l] + r[k][l]) * 0.5;
    bou5[l] = ep[l] * (r[k][l] - s[l][1]);
    del[2][l] = r[k + 1][l] - fe[l];
    del[3][l] = 2. * del[2][l];
    bou1 = max(br[1][k][l], bl[1][l]) * binc[3][l] +
       2. * max(br[2][k][l], bl[2][l]) * binc[2][l] +
       max(br[3][k][l], bl[3][l]) * binc[1][l];
    bou3 = bou4[l] + bou1 * pow(del[3][l], 3.) * prod[l] / 12.;
    itype[l] = 3;
    if(bou3 <= bou5[l]) goto L200;
    bou1 = 0.;
    for(ij=1;ij<= 5;ij++)
    {  jk = 6 - ij;
       bou2 = max(br[ij][k][l], bl[ij][l]);
       bou1 += bou2 * binc[jk][l] * bcs[ij];
    }
    bou3 = bou4[l] + bou1 * pow(del[2][l], 5.) * prod[l] / 90.;
    itype[l] = 2;
    if(bou3 <= bou5[l]) goto L200;
    if(simps==0) 
    {  del[1][l] = 0.5 * del[2][l];
       bou1 = 0.; itype[l] = 1;
       for(ij=1;ij<= 7;ij++)
       {  jk = 8 - ij;
          bou2 = max(br[ij][k][l], bl[ij][l]);
          bou1 += bou2 * binc[jk][l] * bcn[ij];
       }
       bou3 = bou4[l] + bou1 * pow(del[1][l], 7.) * prod[l] * 8. / 945.;
       if(bou3 > bou5[l]) { k++; if(k > 35) return; else goto L100;}
    }
    else
    { k++; 
      if(k > 35) return; else goto L100;
    }

L200: bint[l] += bou3 - bou4[l];
    bou4[l] = bou3; ksa[l] = k;

    if(ind[l] == 6) 
    {  ind[l] = 5; xs = fe[l]; fact[l] = bl[1][l]; }
    else if(itype[l]>2) 
    {  ind[l] = 1; k = ksa[l]; xs = r[k][l]; fact[l] = br[1][k][l]; }
    else 
    {  if(itype[l]==2) 
       {  ind[l] = 3; k = ksa[l]; xs = r[k + 1][l];}
       else
       {  ind[l] = 4; k = ksa[l]; xs = 0.5 * (fe[l] + r[k + 1][l]);}
       t = xs / cond[l]; fact[l] = phistar(t, condl[l]);
    }

L203: xss[nm1][l+1] = xss[nm1][l] + bh[l] * (xs + xm[l]);
    for(ll=l;ll<= nm2;ll++)
    {  xss[ll][l+1] = xss[ll][l] + beta[ll][l] * (xs + xm[l]); }

    if(l != nm2) 
    {  xm[l+1] = xss[l][l+1];
       prod[l+1] = prod[l] * fact[l];
       l++; goto L90;
    }

    x1 = fac * (xss[nm1][nm1] - d[1]) / sd[1];
    x2 = fac * (xss[nm2][nm1] - d[2]) / sd[2];
    l = nm1;
    ans[l] = bivnor(x1, x2, rho);
    if(ipr!=1) 
    {  y1 = (xss[nm1][nm1] - c[1]) / sd[1];
       y2 = (xss[nm2][nm1] - c[2]) / sd[2];
       wu = bivnor(y1, y2, rho); wt = bivnor(x1, y2, rho);
       wb = bivnor(y1, x2, rho); ans[l] += wu - wt - wb;
    }

L310: if(l == 1) 
    {  *ifault = 0; *prob = ans[1]; *bound = bint[1] + eplos; return;}
    l--; indl = ind[l]; numl = num[l]; ity = itype[l];
    temp = fact[l] * ans[l + 1]; temb = bint[l + 1];
    fsa = 1.;
    if(indl == 1) 
    {  tem = temp; temp += prep[l]; temb += preb[l];
       prep[l] = tem; preb[l] = bint[l + 1];
       fsa = 2.;
    }
    ans[l] += coef[indl][ity] * temp * del[ity][l];
    bint[l] += coef[indl][ity] * temb * del[ity][l];
    ep[l] += del[ity][l] * coef[indl][ity] * (fsa*((double) nm2- l) 
     * epsi / pr2[l + 1] - temb) / (s[l][numl] - s[l][1]);

    if(indl != 1) 
    {  igo = indl - (1 + (itype[l] * (itype[l] - 1)) / 2);
       if(igo>=1 && igo <=4)
                         //goto (210, 208, 206, 205), igo;
       {  if(igo==1)
          {  ind[l] = 1; k = ksa[l]; xs = r[k][l]; fact[l] = br[1][k][l]; }
          else 
          {  if(igo==2)
             {  ind[l] = 2; k = ksa[l]; xs = 0.5 * (r[k][l] + r[k + 1][l]);}
             else if(igo==3)
             {  ind[l] = 3; k = ksa[l]; xs = r[k + 1][l];}
             else if(igo==4)
             {  ind[l] = 4; k = ksa[l]; xs = 0.5 * (fe[l] + r[k + 1][l]);}
             t = xs / cond[l]; fact[l] = phistar(t, condl[l]);
          }
          goto L203;
       }
    }

    k = ksa[l];
    for(i=1;i<=ien;i++) bl[i][l] = br[i][k][l];
    ind[l] = 5;
    fe[l] = r[k][l];
    if(k != 1) 
    {  k--; goto L104;}
    if(intvl[l] == num[l]) goto L310;
    intvl[l]++; intl = intvl[l]; r[1][l] = s[l][intl];
    goto L100;
}

/* inputs
     x = co-ordinate
   output
     univariate standard normal density
   Algorithm  AS 195.1, Appl Statist (1984), v 33 (1)
*/
double phistar(double x, double y)
{  double arg, tem;
   static double xlow= -87.0, sq2p= 0.91893853320467274;
   tem = 0.;
   arg = -0.5 * x * x - sq2p - y;
   if(arg > xlow) tem = exp(arg);
   return tem;
}

/* inverting a[ ] (vectorization of lower triangular of correlation matrix)
   and the n-2 lower right-hand principal minors */
/* inputs
     a = lower triangular of correlation matrix
     ai = inverse of a
     c = matrix to store coefficients
     n = dimension
   output
     d = determinant
     ifault = error code
*/
void invert(double a[], double ai[], double c[][16], int n, 
       double *d, int *ifault)
{  
/*    dimension a[1],ai[1],c[5][1],sum[28],b[28],bi[28]; */
   double sum[29],b[29],bi[29];
   double tol,su;
   int jt,ij,ik,jj,ke,ijt,jk,i,ijl,j,i1,j1,ij1,j2,j3,j4,j5,ii,l,kl,n2,k;
   tol=5.0e-08;
   *d=1.0;
   *ifault=4;
   jt=n*(n+1)/2+1;
   ij=jt-1;
   c[1][1]=1.0/a[ij];
   for(i=1;i<=ij;i++) sum[i]=0.;
   for(ij=1;ij<=n;ij++)
   {  ik=n-ij; jt--;
      if(a[jt]<=sum[jt]+tol) return;
      b[jt]=sqrt(a[jt]-sum[jt]);
      (*d)*=b[jt]*b[jt];
      if(ik==0) break;
      ijt=jt;
      for(jj=1;jj<=ik;jj++)
      {  jt--;  b[jt]=(a[jt]-sum[jt])/b[ijt]; }
      i=0;
      ke=ijt-1;
      for(jj=jt;jj<=ke;jj++)
      {  for(jk=jt;jk<=jj;jk++)
         {  i++;  sum[i]+=b[jj]*b[jk];}
      }
   }
   bi[1]=1./b[1];
   for(i=2;i<=n;i++)
   {  ijt=i*(i+1)/2; bi[ijt]=1./b[ijt];
      ijl=ijt-i; ij=i-1;
      for(j=1;j<=ij;j++)
      {  jt=ijt-j; kl=i-j+1; su=0.;
         for(k=kl;k<=i;k++)  su+=bi[ijl+k]*b[k*(k-1)/2+kl-1];
          bi[jt]=-su/b[kl*(kl-1)/2];
      }
   }
   if(n!=3) 
   {  n2=n-2;
      for(i1=2;i1<=n2;i1++)
      {  i=n-i1+1;
         j=0;
         for(j1=1;j1<=i1;j1++)
         {  ij1=((i+j1-1)*(i+j1-2))/2;
            for(j2=1;j2<=j1;j2++)
            {  j++; c[i1][j]=0.; j4=i+j2-1; j5=(j4*(j4-1))/2;
               for(j3=i;j3<=j4;j3++) c[i1][j]+=bi[ij1+j3]*bi[j5+j3];
            }
         }
      }
   }
   jt=0;
   for(j=1;j<=n;j++)
   {  jj=j*(j-1)/2;
      for(i=1;i<=j;i++)
      {  ii=i*(i-1)/2; jt=jt+1; ai[jt]=0.;
         for(l=1;l<=i;l++) ai[jt]+=bi[ii+l]*bi[jj+l];
      }
   }
   *ifault=0;
   return;
}

/* input
     p = value in (0,1)
   outputs
     ifault = error code
     standard normal quantile function at p
*/
double gauinv(double p, int *ifault)
{  static double p0=-0.322232431088, p1=-1.0, p2=-0.342242088547,
      p3=-0.204231210245e-01, p4=-0.453642210148e-04;
   static double  q0=0.99348462606e-01, q1=0.588581560495,
       q2=0.531103462366, q3=0.103537752850, q4=0.38560700634e-02;
   double tem,ps,yi;
   tem=0.0;
   ps=p;  if(ps>0.5) ps=1.0-ps;
   *ifault=1;
   if(ps<1.0e-20) return tem;
   *ifault=0;
   if(ps==0.5) return tem;
   yi=sqrt(log(1.0/(ps*ps)));
   tem=yi+((((yi*p4+p3)*yi+p2)*yi+p1)*yi+p0)
      /((((yi*q4+q3)*yi+q2)*yi+q1)*yi+q0);
   if(p<0.5) tem=-tem;
   return tem;
}

#ifdef ALONE
double min(double a, double b)
{  return (a<b)? a:b; }

double max(double a, double b)
{  return (a<b)? b:a; }
#endif

double ff(double x, double y) 
{ return fabs((x*x - 1.) / (y*y)); }

double f3(double x, double y) 
{ return fabs(-x * (x*x - 3.) / (y*y*y)); }

double f4(double x, double y) 
{ return fabs((3. + x*x * (-6. + x*x)) / pow(y,4.)); }

double f5(double x, double y) 
{ return fabs(x * (-15. + x*x * (10. - x*x))) / pow(y, 5.); }

double f6(double x, double y) 
{ return fabs(-15. + x*x* (45. - x*x* (-15. + x*x))) / pow(y, 6.); }

