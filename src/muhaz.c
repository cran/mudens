/*
Produced by PROMULA.FORTRAN to C V7.16 on 12/28/06 at 15:07:19
*/
#define LPROTOTYPE
#define TRUE 1
#define FALSE 0
#include "muhaz.h"
#include <math.h>
#include "fortran.h"
/*
Global variables
*/
double Xhpilot[1000];
int Xalpha,Xbeta,Xnu,Xdens;

/* Function atpos
*/
int atpos(double *v, int *n, double *x)
{
/*
     Determines the position of the scalar x, within the vector v
*/
static int atpos,i;
   if(*x < v[0]) {
      atpos = 0;
      return atpos;
   }
   if(*x > v[*n-1]) {
      atpos = *n;
      return atpos;
   }
   for(i=1; i<=*n; i++) {
      if(*x-v[i-1] >= 0.0e0) atpos = i;
   }
   return atpos;
}

/* Function bsmoth
*/
void bsmoth(int *n, 
			double *z, 
			double *bopt, 
			int *m, 
			double *zz, 
			double *bsm, 
			double *b, 
			int *kflag, 
			double *endl, 
			double *endr
			)
{
/*
     Computes the smoothed bandwidths at zz, using optimal bandwidths
     at z, and smoothing bandwidth b
*/
#define zero 0.0e0
#define one 1.0e0
static double T1;
static int ilo,ihi,i,j;
static double z0,sum1,sum2,ker,q,u;
   for(j=0; j<*m; j++) {
      z0 = zz[j];
/*
 --- compute for "active" points only
*/
      ibnds(z,n,&z0,b,&ilo,&ihi);
      sum1 = sum2 = zero;
      for(i=ilo-1; i<ihi; i++) {
/*
 --- u in [-1, 1], for any i in [ilo, ihi] (this is what ibnds does)
*/
         u = (z0-z[i])/ *b;
         if((*kflag == 0 || *endl+*b <= z0) && (z0 <= *endr-*b)) {
/*
 --- INTERIOR point; no boundary correction
*/
            T1 = one;
            ker = kernel(&T1,&u);
         }
         else if(*endl <= z0 && z0 < *endl+*b) {
/*
     LEFT boundary
*/
            q = (z0-*endl)/ *b;
            ker = kernel(&q,&u);
         }
         else {
/*
     RIGHT boundary
*/
            q = (*endr-z0)/ *b;
            T1 = -u;
            ker = kernel(&q,&T1);
         }
         sum1 += (ker*bopt[i]);
         sum2 += ker;
      }
      bsm[j] = sum1/sum2;
   }
   return;
}
#undef zero
#undef one

/* Function func
*/
void func(int *n,
		  double *x,
		  double *delta,
		  double *z,
		  double *b,
		  double *endl,
		  double *endr,
		  double *q,
		  double *y,
		  double *bb,
		  double *vv,
		  double *bpilot,
		  int *kflag
		  )
{
static double negy,zz,newz,fz,k;
   newz = *z-*b**y;
   fz = hazden(n,x,delta,&newz,bpilot,endl,endr,kflag);
   negy = *y;
   if(*endr-*b < *z && *z <= *endr) negy = -*y;
   k = kernel(q,&negy);
   *bb = fz*k;
   zz = surfct(x,delta,n,&newz);
   *vv = k*k*fz/zz;
   return;
}

/* Function gets
 Computes the value of the survival function at x0
 Uses BINARY search
     x --> table (matrix) computed using Kaplan-Meier
           1st column: survival times
           2nd column: survival function values for the corresp. time
     n --> number of unique survival times
     x0 --> survival time for which the survival value is desired
*/

double gets(double **x, int *n, double *x0)
{
static double gets;
static int ilo,ihi,ihf;
   if(*x0 < x[0][0]) gets = 1.0e0;
   else if(*x0 >= x[0][*n-1]) gets = x[1][*n-1];
   else {
/*
     binary search (halving at each step)
*/
      ilo = 1;
      ihi = *n;
S100:
      if(ihi-ilo == 1) {
         gets = x[1][ilo-1];
         return gets;
      }
      ihf = (ilo+ihi)/2;
      if(x[0][ihf-1] < *x0) {
         ilo = ihf;
         goto S100;
      }
      else if(x[0][ihf-1] > *x0) {
         ihi = ihf;
         goto S100;
      }
      else gets = x[1][ihf-1];
   }
   return gets;
}

/* Function glmin
 Computes optimal GLOBAL bandwidth, by minimizing the
 Integrated Mean Squared Error (IMSE).
                         Algorithm
 For each bandwidth, compute IMSE.  Take as OPTIMAL bw the one yielding
 the smallest IMSE.  If optimal bw is the first, i.e. bw(1), it means
 that no minimum was reached.  In this case, set optimal global bw to
 the largest, i.e. bw(gridb)
*/
void glmin(int *n,
		   double *x,
		   double *delta,
		   double *z,
		   int *gridz,
		   double *bw,
		   int *gridb,
		   double *endl,
		   double *endr,
		   double *bpilot,
		   double *imsemn,
		   double *globlb,
		   double *glmse,
		   int *kflag
		   )
{
#define zero 0.0e0
#define huge 1.0e30
static int i,j;
static double imse,mse,bias,var;
   *imsemn = huge;
   *globlb = bw[*gridb-1];
   for(j=0; j<*gridb; j++) {
      imse = zero;
      for(i=0; i<*gridz; i++) {
         msemse(n,&z[i],endl,endr,x,delta,&bw[j],&mse,&bias,&var,bpilot,&
           Xhpilot[i],kflag);
         imse += mse;
      }
      if(imse > zero && imse < *imsemn) {
         *imsemn = imse;
         *globlb = bw[j];
      }
      glmse[j] = imse;
   }
   return;
}
#undef zero
#undef huge

/* Function h
*/
double h(int *j)
{
static double fact[8] = {
   1.0e0,1.0e0,2.0e0,6.0e0,24.0e0,120.0e0,720.0e0,5040.0e0
};
static double h;
   h = (double)((2.0F*(float)*j+Xalpha+Xbeta+1.0F)/pow(2.0F,(Xalpha+
     Xbeta+1)))*fact[*j]*fact[*j+1+Xalpha+Xbeta-1]/(fact[*j+1+Xalpha-1]*fact[*j
     +1+Xbeta-1]);
   return h;
}

/* Function hazden
*/
double hazden(int *n,
			  double *x,
			  double *delta,
			  double *z,
			  double *b,
			  double *endl,
			  double *endr,
			  int *kflag
			  )
{
#define zero 0.0e0
#define one 1.0e0
static double T1;
static double hazden;
static int i,ilo,ihi;
static double q,u,ker,fz,fz1;
/*
 compute indices of "active" points
*/
   ibnds(x,n,z,b,&ilo,&ihi);
   hazden = fz = fz1 = zero;
   for(i=ilo; i<=ihi; i++) {
      if(delta[i-1] == one) {
         u = (*z-x[i-1])/ *b;
/*
 --- u in [-1, 1], for any i in [ilo, ihi] (this is what ibnds does)
*/
         if((*kflag == 0 || *endl+*b <= *z) && (*z <= *endr-*b)) {
/*
 --- INTERIOR point; no correction
*/
            T1 = one;
            ker = kernel(&T1,&u);
         }
         else if(*endl <= *z && *z < *endl+*b) {
/*
 --- LEFT boundary correction
*/
            q = (*z-*endl)/ *b;
            ker = kernel(&q,&u);
         }
         else {
/*
 --- RIGHT boundary correction
*/
            if(*kflag == 1) {
/*
 --- LEFT only boundary correction; like INTERIOR
*/
               T1 = one;
               ker = kernel(&T1,&u);
            }
            else {
               q = (*endr-*z)/ *b;
               if(-q <= u) {
                  T1 = -u;
                  ker = kernel(&q,&T1);
               }
               else goto S100;
            }
         }
         if(Xdens > 0) fz += (ker/(double)*n);
         else fz += (ker/(double)(*n-i+1));
      }
S100:;
   }
   hazden = fz/pow(*b,(double)(Xnu+1));
/*
 --- do not allow negative hazards
*/
   if(hazden < zero) hazden = zero;
   return hazden;
}
#undef zero
#undef one

/* Function ibnds
*/
void ibnds(double *x,
		   int *n,
		   double *z,
		   double *b,
		   int *ilo,
		   int *ihi
		   )
{
/*
 Ensures that all components of x, with indices ilo...ihi,
 are within (z-b, z+b) ("active points", i.e. yielding weights
 in the kernel smoothing)
     x --> vector of observation points
     n --> length of x
     z --> grid point for which "active" observation points are compute
     b --> bandwidth; "active" points are those in [z-b, z+b]
     ilo <-- index of the first "active" observation point
     ihi <-- index of the last "active" observation point
*/
static int i;
static double what;
   what = *z-*b;
   for(i=1; i<=*n; i++) {
      if(what < x[i-1]) {
         *ilo = i;
         goto S150;
      }
   }
   *ilo = *n+1;
S150:
   what = *z+*b;
   if(what >= x[*n-1]) {
      *ihi = *n;
      return;
   }
   for(i= *n; i>=*ilo; i--) {
      if(what > x[i-1]) {
         *ihi = i;
         return;
      }
   }
   *ihi = 0;
   return;
}

/* Function intgrl
*/
void intgrl(int *n,
			double *x,
			double *delta,
			double *z,
			double *b,
			double *endl,
			double *endr,
			double *q,
			double *r,
			double *s,
			double *valueb,
			double *valuev,
			double *bpilot,
			int *kflag
			)
{
#define jmax 6
#define epsi (1.0e-3)
#define huge 1.0e30
static int D2;
static int j;
static double oldb,oldv;
   oldb = oldv = -huge;
   for(j=1,D2=(jmax-j+1); D2>0; D2--,j+=1) {
      try(n,x,delta,z,b,endl,endr,q,r,s,valueb,valuev,&j,bpilot,kflag);
      if(fabs(*valueb-oldb) <= epsi*fabs(oldb) && fabs(*valuev-oldv) <= epsi*
        fabs(oldv)) return;
      oldb = *valueb;
      oldv = *valuev;
   }
   return;
}

/* Function kapmei
*/
#undef jmax
#undef epsi
#undef huge
void kapmei(double *times, double *delta, int *n, double **x, int *count)
{
/*
 Computes Kaplan-Meier estimates
*/
#define one 1.0e0
static int i,j,equals,lsteql,atrisk,events;
static double prob;
   atrisk = *n;
   *count = lsteql = 0;
   prob = one;
/*
 Loop over all observations
*/
   i = 1;
S10:
   if(i >= *n) goto S200;
/*
     find how many equal times
*/
   equals = 0;
   events = (int)delta[i-1];
   for(j=i+0; j<*n; j++) {
      if(times[j] == times[i-1]) {
         equals += 1;
         events = events+(int)delta[j];
      }
      else goto S150;
   }
S150:
   *count += 1;
   atrisk -= lsteql;
   lsteql = equals+1;
   x[0][*count-1] = times[i-1];
   x[1][*count-1] = prob*(one-(double)events/(double)atrisk);
   prob = x[1][*count-1];
   i += lsteql;
   goto S10;
S200:
   return;
}

/* Function kernel
*/
#undef one
double kernel(double *q, double *x)
{
#define one 1.0e0
#define two 2.0e0
static int K1 = 0;
static int K2 = 1;
static int K3 = 2;
static int K4 = 3;
static double kernel;
static int k;
static double c1,c2,ker,temp1,temp2;
   k = Xnu+2;
   ker = 0.0e0;
   c1 = two*(one+*x)/(one+*q)-one;
   c2 = (one-*q)/(one+*q);
   if(Xnu == 0) ker = h(&K1)+h(&K2)*pzero(&K2,&c1)*pzero(&K2,&c2);
   else if(Xnu == 1) ker = h(&K2)*pzero(&K2,&c1)*pone(&K2,&c2)+h(&K3)*pzero(&K3,
     &c1)*pone(&K3,&c2);
   else if(Xnu == 2) ker = h(&K3)*pzero(&K3,&c1)*ptwo(&K3,&c2)+h(&K4)*pzero(&K4,
     &c1)*ptwo(&K4,&c2);
   temp1 = temp2 = one;
   if(Xalpha > 0) temp1 = pow(*q-*x,(double)Xalpha);
   if(Xbeta >= 0) temp2 = pow(one+*x,(double)Xbeta);
   ker = temp1*temp2*pow(two/(one+*q),(double)(Xalpha+Xbeta+1+Xnu))*ker;
   if(fifmod(Xnu,2) == 1) ker = -ker;
   kernel = ker;
   return kernel;
}

/* Function knncen
*/
#undef one
#undef two
void knncen(double *times,
			double *status,
			int *n,
			double *z,
			int *nz,
			int *k,
			double *bw
			)
{
/*
     k-th nearest neighbor for survival data with censorship
     z  --> grid points where bandwidths are calculated
     nz --> number of grid points
     For each grid point z0 in z, computes the bandwidth bw,
     so that in each [z0-bw, z0+bw] there are k neighbors
     (uncensored observations, i.e. status=1)
*/
static int i,j,iv,ilo,ihi,ipos,count;
static double z0,tcopy[20000],td[20000];
/*
     "Clean" the survival times vector,
     i.e. eliminate all censored observations
*/
   count = 0;
   for(i=0; i<*n; i++) {
      if(status[i] != 0.0e0) {
/*
 --- UNCENSORED observation
*/
         count += 1;
         tcopy[count-1] = times[i];
      }
   }
/*
     Compute bandwidth for each grid point z, so that in each vicinity
     of z, there are k survival times (uncensored)
*/
   for(i=0; i<*nz; i++) {
      z0 = z[i];
      ipos = atpos(tcopy,&count,&z0);
      ilo = fifmax0(ipos-*k,1);
      ihi = fifmin0(ipos+*k,count);
      iv = 0;
      for(j=ilo-1; j<ihi; j++) {
         iv += 1;
         td[iv-1] = fabs(tcopy[j]-z0);
      }
      sorter(td,&iv);
      bw[i] = td[*k-1];
   }
   return;
}

/* Function knnhad
     This procedure provides hazard function or density estimates,
     with boundary modifications and local/global bandwidth choice,
     for censored and uncensored case.
 Nearest methods algorithms to compute the bandwidth to be used with
 modified kernel polynomials

      PARAMETERS:

      n --> number of observations
      x --> vector of survival times
      delta --> censoring indicator for hazard estimation
                   0 - censored observation
                   1 - uncensored observation

      bwchoi --> method to determine the bandwidth
                   1 - nearest neighbors, eliminating censored observat
                   2 - nearest neighbors, with censored observations
      gridz --> number of points in the minimization grid z(gridz)
      z --> the minimization grid to determine optimal bw.
            gridz equidistant points between startz and endz;
            Number of points gridz influences computing time strongly.
            Usually: gridz < m
      m --> number of points in output (estimation) grid
      zz <-- output (estimation) grid, computed as m equidistant points
             between startz and endz at which curve is estimated
      bpilot --> initial bandwidth which is used to  estimate bias and
      endl --> assumed left endpoint of function to be estimated.
               Boundary kernels are used in [endl,endl+b).
      endr --> assumed right endpoint of function to be estimated.
               Boundary kernels are used in (endr-b,endr],
               b being the bandwidth.
      bsmo --> bandwidth for smoothing of local bandwidths,
               Small value recommended for local bandwidth choice.
               Specification not necessary for global bandwidth choice.
      kflag --> 0 - unmodified kernels (Epanechnikov)
                1 - boundary corrected kernels
      fzz <--  function estimate at zz(i), i=1, m
      kmin <-> minimum number of neighbors.
               OUTPUT: optimum number of neighbors
      kmax --> maximum number of neighbors.
      bopt <-- chosen bandwidth at each z(j), j=1, gridz, for local cho
               Note that final bandwidth is used on grid zz(m)
               and is obtained from bopt(gridz) by smoothing
               with bandwidth bsmo.
      bopt1 <-- smoothed bandwidth used at output point zz(i)
      kimse <-- IMSE at each k
*/
void knnhad(int *n,
			double *x,
			double *delta,
			int *niu,
			int *alph,
			int *bet,
			int *dns,
			int *bwchoi,
			int *gridz,
			double *z,
			int *m,
			double *zz,
			double *bpilot,
			double *endl,
			double *endr,
            double *bsmo,
			int *kflag,
			double *fzz,
			int *kmin,
			int *kmax,
			double *bopt,
			double *bopt1,
			double *kimse
			)
{
static int i;
/*
 --- Save parameters to the common storage
*/
   Xnu = *niu;
   Xalpha = *alph;
   Xbeta = *bet;
   Xdens = *dns;
/*
 --- Compute the hazard estimates for bw=bpilot.
     They will be used in msemse
*/
   for(i=0; i<*gridz; i++) {
      Xhpilot[i] = hazden(n,x,delta,&z[i],bpilot,endl,endr,kflag);
   }
   if(*bwchoi == 1) {
/*
     --- simple nearest neighbor approach
*/
      knnmin(x,delta,n,z,gridz,endl,endr,bpilot,bopt,kmin,kmax,kimse,kflag);
   }
   else if(*bwchoi == 2) {
/*
     --- modified nearest neighbor approach
*/
      olafmn(x,delta,n,z,gridz,endl,endr,bpilot,bopt,kmin,kmax,kimse,kflag);
   }
   else {
/*
     for future algorithms ...
*/
   }
/*
       --- smooth the local bandwidths: BOPT ==> BOPT1(i)
*/
   bsmoth(gridz,z,bopt,m,zz,bopt1,bsmo,kflag,endl,endr);
/*
     compute the estimates, using optimal bandwidths
*/
   for(i=0; i<*m; i++) {
      fzz[i] = hazden(n,x,delta,&zz[i],&bopt1[i],endl,endr,kflag);
   }
   return;
}

/* Function knnmin
     Computes the bandwidth at each grid point z
     First it finds optimum number of neighbors, minimizing the IMSE
*/
void knnmin(double *x,
			double *delta,
			int *n,
			double *z,
			int *gridz,
			double *endl,
			double *endr,
			double *bpilot,
			double *bopt,
			int *kmin,
			int *kmax,
			double *kimse,
			int *kflag
			)
{
#define zero 0.0e0
#define huge 1.0e5
static int D2;
static int k,i,kopt;
static double imse,imsemn,bias,var,mse,bwi,zi;
   if(*kmin == *kmax) {
      knncen(x,delta,n,z,gridz,kmin,bopt);
      return;
   }
   imsemn = huge;
/*
     For the moment, the maximum number of neighbors to be considered
     is half of the non-censored observations
*/
   for(k= *kmin,D2=(*kmax-k+1); D2>0; D2--,k+=1) {
/*
     compute the bandwidths bopt for k neighbors
*/
      knncen(x,delta,n,z,gridz,&k,bopt);
      imse = zero;
/*
     compute MSE at each gridpoint
*/
      for(i=0; i<*gridz; i++) {
         zi = z[i];
         bwi = bopt[i];
         msemse(n,&zi,endl,endr,x,delta,&bwi,&mse,&bias,&var,bpilot,&Xhpilot[i],
           kflag);
         imse += mse;
      }
      if(imse < imsemn) {
         kopt = k;
         imsemn = imse;
      }
      kimse[k-*kmin] = imse;
   }
/*
     kopt is returned in kmin
*/
   *kmin = kopt;
/*
     compute the bandwidths for kopt
*/
   knncen(x,delta,n,z,gridz,&kopt,bopt);
   return;
}

/* Function ksi
*/
#undef zero
#undef huge
double ksi(int *j)
{
static double ksi,a,b,c,d;
   a = -2.0e0;
   b = *j+Xalpha-1;
   c = *j+Xbeta-1;
   d = 2.0e0*(double)*j+(double)Xalpha+(double)Xbeta;
   ksi = a*b*c*d;
   return ksi;
}

/* Function loclmn
 Computes optimal LOCAL bandwidths at each gridpoint of z.
                         Algorithm
 For each gridpoint in z, computes MSE for each bandwidth in the grid.
 The OPTIMAL bandwidth is the one yielding the smallest MSE.
 NOTE: If no minimum was found for MSE, then set bopt to largest bw
*/
void loclmn(int *n,
			double *x,
			double *delta,
			double *z,
			int *gridz,
			double *bw,
			int *gridb,
			double *bopt,
			double *endl,
			double *endr,
			double *bpilot,
			double *msemin,
			double *biasmn,
            double *varmin,
			int *kflag
			)
{
#define zero 0.0e0
#define huge 1.0e30
static double mse,bias,var,amin;
static int i,j;
   for(i=0; i<*gridz; i++) {
      amin = huge;
      bopt[i] = bw[*gridb-1];
      for(j=0; j<*gridb; j++) {
         msemse(n,&z[i],endl,endr,x,delta,&bw[j],&mse,&bias,&var,bpilot,&
           Xhpilot[i],kflag);
         if(mse > zero && mse < amin) {
            amin = mse;
            bopt[i] = bw[j];
            biasmn[i] = bias;
            varmin[i] = var;
         }
      }
      if(amin == huge) {
/*
     --- all mse were 0.  Set bopt to the largest bw
*/
         msemin[i] = zero;
      }
/*
         call intpr('i=', 2, i, 1)
         call dblepr('bw=', 3, bopt(i), 1)
*/
      msemin[i] = amin;
   }
   return;
}
#undef zero
#undef huge

/* Function locolf
     Computes array of bandwidths at each grid point in xgrid
     x --> matrix generated by a Kaplan-Meier survival estimation
           1st column: unique survival times
           2nd column: corresponding survival function values
     nobs --> initial number of observations (censored and uncensored)
     xgrid --> vector of survival times, where the bandwidths are
               to be computed
     ngrid --> number of grid points
     n --> number of unique survival times
     k --> number of neighbors
     bw <-- vector at bandwidths, computed at each grid point
*/
void locolf(double **x,
			int *nobs,
			double *xgrid,
			int *ngrid,
			int *n,
			int *k,
			double *bw
			)
{
static int i;
   for(i=0; i<*ngrid; i++) {
      bw[i] = oneolf(x,n,&xgrid[i],nobs,k);
   }
   return;
}
/* Function luo
*/
double luo(int *j)
{
static double luo;
   luo = (double)(2.0F*(float)*j*(float)(*j+Xalpha+Xbeta)*(float)(2**j+Xalpha+
     Xbeta-2));
   return luo;
}

/* Function msemse
*/
void msemse(int *n,
			double *z,
			double *endl,
			double *endr,
			double *x,
			double *delta,
			double *b,
			double *mse,
			double *bias,
			double *var,
			double *bpilot,
			double *fz,
			int *kflag
			)
{
/*
 Computes MSE at z, for bw=b, using (2.4) from Mueller, pg. 64
*/
#define one 1.0e0
static double r,s,q,valueb,valuev;
   if((*kflag == 0 || *endl+*b <= *z) && (*z <= *endr-*b)) {
/*
     --- INTERIOR point; no correction
*/
      q = one;
      r = -one;
      s = one;
   }
   else if(*endl <= *z && *z < *endl+*b) {
/*
     --- LEFT boundary correction
*/
      q = (*z-*endl)/ *b;
      r = -one;
      s = q;
   }
   else {
/*
     --- RIGHT boundary correction
*/
      if(*kflag == 1) {
/*
--- actually, only LEFT boundary correction; this is like INTERIOR
*/
         q = one;
         r = -one;
         s = one;
      }
      else {
/*
 --- indeed, RIGHT boundary correction
*/
         q = (*endr-*z)/ *b;
         r = -q;
         s = one;
      }
   }
   intgrl(n,x,delta,z,b,endl,endr,&q,&r,&s,&valueb,&valuev,bpilot,kflag);
   *bias = valueb-*fz;
   *var = valuev/(double)*n/pow(*b,(double)(Xnu+1));
   *mse = *bias**bias+*var;
   return;
}

/*
Function newhad
     This procedure provides hazard function or density estimates
     for survival data, with boundary modifications and local/global
     bandwidth choice, for censored and uncensored case.
        INPUT DATA NEED NOT BE ORDERED
        ** VERSION JULY 92
        ** MODIFIED DEC 93
        ** NO RESPONSIBILITY IS ASSUMED FOR CORRECTNESS OF CODE
        ** COPYRIGHT H.G.MUELLER & J.L.WANG, DIVISION OF
        ** STATISTICS, UNIVERSITY OF CALIFORNIA, DAVIS, CA 95616 USA
        ** PROCEDURE DESCRIBED IN:
        **  MUELLER,H.G.,WANG,J.L.: HAZARD RATE ESTIMATION UNDER RANDOM
        **  CENSORING WITH VARYING KERNELS AND BANDWIDTHS,
        **  BIOMETRICS 50, 61-76, 1994.

 ***** Modified: Dan M. Serachitopol, Nov 1997
      PARAMETERS:

      n --> number of observations
      x --> vector of observations (survival times)
      delta --> censoring vector for the observations
                   0 -  censored observation
                   1 -  uncensored observation
      local -->  local or global bandwidth choice by minimizing
                 direct convolution type estimates of mse/imse.
                   0 - global minimizer
                   1 - local minimizer
      z --> the minimization grid
      gridz --> number of points for the minimization grid.
                NOTE: gridz influences computing time strongly.
                      Usually gridz < m
      niu   --> order of derivative to be estimated (=0,1 or 2)
                niu=1 or niu=2 only allowed for dens=0 or dens=1
      alph, bet --> determine choice of kernel and boundary kernels.
                    Kernels for left boundary interval are given by:
                          (x+1)**bet*(q-x)**alph*P(x),
                    where P(x) is a polynomial, with support on [-1,q]
                    The situation is symmetric for the right boundary.
                    For alph=bet=mu, smooth optimum boundary kernels
                    are obtained.
                    For alpha=mu-1, beta=mu, boundary kernels with
                    improved variance behavior are computed.
                    Typical choice: alpha=0, beta=1; alpha=1, beta=2.
          NOTE: Order of kernels is nu+2.
                Epanechnikov kernel can be continued by choices:
                   alpha=beta=mu=1 or alpha=0, beta=1.
      dns   --> type of curve to be estimated
                   0 - hazard estimation
                   1 - density estimation
      zz --> estimation grid grid, where the hazards are estimated
      m --> number of points in the estimation grid
      bpilot --> initial bandwidth which is used to  estimate bias/var
      bw    --> the bandwidths grid, scanned in the MSE minimization
      gridb --> number of bandwidths in the grid
                NOTE: If gridb=1, startb is used as a global optimal
                      bandwidth to compute the hazard estimates
      endl, endr --> bounds of the  function to be estimated.
                     Corrected boundary kernels are used in:
                     LEFT  boundary: [endl, endl+b]
                     RIGHT boundary: [endr-b, endr]
      bsmo --> bandwidth for smoothing the local optimal bandwidths.
               Small value recommended for local bandwidth choice.
               Specification not necessary for global bandwidth choice.
      kflag --> 0 - no boundary correction
                1 - LEFT boundary correction, ONLY
                2 - LEFT and RIGHT boundary correction
              OUTPUT
      fzz <-- hazard estimates at zz(1:m)
      bopt <-- optimal local bandwidth at z(1:gridz)
      bopt1 <-- bandwidth used to compute estimated hazards.
                It is obtained from bopt, smoothing with the bandwidth
      msemin <-- minimum MSE at each z, for local choice
      biasmn <-- minimum bias at each z, for local choice
      varmin  <-- minimum variance at each z, for local choice
      imsemn <-- minimum IMSE,  for global and local choice
      globlb  <-- optimal global bandwidth, resulting from MSE minimiza
      glmse <-- IMSE at bandwidth grid, for global choice
      b <-- bandwidth grid for MSE minimization, for global choice

*/
#undef one
void newhad(int* n, 
			double* x,
			double* delta,
			int* local,
			double* z, 
			int* gridz, 
			int* niu,
			int* alph,
			int* bet,
			int* dns,
			double* zz,
			int* m, 
			double* bpilot, 
			double* bw, 
			int* gridb, 
			double* endl,
			double* endr, 
			double* bsmo, 
			int* kflag, 
			double* fzz,
			double* bopt,
			double* bopt1,
			double* msemin,
			double* biasmn,
			double* varmin,
			double* imsemn,
			double* globlb,
			double* glmse
			)
{
static int i;
/*
 --- Save parameters to the common storage
*/
   Xnu = *niu;
   Xalpha = *alph;
   Xbeta = *bet;
   Xdens = *dns;
   if(*gridb == 1) {
      *globlb = bw[0];
      goto S90;
   }
/*
 --- Compute the hazard estimates for bw=bpilot.
     They will be used in msemse
*/
   for(i=0; i<*gridz; i++) {
      Xhpilot[i] = hazden(n,x,delta,&z[i],bpilot,endl,endr,kflag);
   }
   if(*local == 1) {
      loclmn(n,x,delta,z,gridz,bw,gridb,bopt,endl,endr,bpilot,msemin,biasmn,
        varmin,kflag);
/*
 --- compute IMSE
*/
      *imsemn = 0.0e0;
      for(i=0; i<*gridz; i++) {
         *imsemn += msemin[i];
      }
/*
 --- smooth the optimal local bandwidths
*/
      bsmoth(gridz,z,bopt,m,zz,bopt1,bsmo,kflag,endl,endr);
   }
   else glmin(n,x,delta,z,gridz,bw,gridb,endl,endr,bpilot,imsemn,globlb,glmse,
     kflag);
S90:
   for(i=0; i<*m; i++) {
      if(*gridb == 1 || *local == 0) fzz[i] = hazden(n,x,delta,&zz[i],globlb,
        endl,endr,kflag);
      else fzz[i] = hazden(n,x,delta,&zz[i],&bopt1[i],endl,endr,kflag);
   }
   return;
}

/* Function olafbw
*/
void olafbw(double *times, double *delta, int *n, double *z, int *gridz, int *k, double *bopt)
{
/*
 Computes local bandwidths at each grid point in z
*/
static int count;
static double a[2][20000];
double *x[2];
x[0] = a[0];
x[1] = a[1];
/*
     Call Kaplan-Meier to compute the survival function
*/
   kapmei(times,delta,n,x,&count);
/*
     now compute the bandwidths
*/
   locolf(x,n,z,gridz,&count,k,bopt);
   return;
}

/* Function olafmn
*/
void olafmn(double *x,
			double *delta,
			int *n,
			double *z,
			int *gridz,
			double *endl,
			double *endr,
			double *bpilot,
			double *bopt,
			int *kmin,
			int *kmax,
			double *kimse,
			int *kflag
			)
{
/*
     Computes the bandwidth at each grid point z
     First it finds optimum number of neighbors, minimizing the IMSE
*/
#define zero 0.0e0
#define huge 1.0e5
static int D2;
static int k,i,kopt;
static double imse,imsemn,bias,var,mse,bwi,zi;
   if(*kmin == *kmax) {
      olafbw(x,delta,n,z,gridz,kmin,bopt);
      return;
   }
   imsemn = huge;
/*
     For the moment, the maximum number of neighbors to be considered
     is half of the non-censored observations
*/
   for(k= *kmin,D2=(*kmax-k+1); D2>0; D2--,k+=1) {
/*
     compute the bandwidths bopt for k neighbors
*/
      olafbw(x,delta,n,z,gridz,&k,bopt);
      imse = zero;
/*
     compute MSE at each gridpoint
*/
      for(i=0; i<*gridz; i++) {
         zi = z[i];
         bwi = bopt[i];
         msemse(n,&zi,endl,endr,x,delta,&bwi,&mse,&bias,&var,bpilot,&Xhpilot[i],
           kflag);
         imse += mse;
      }
      if(imse < imsemn) {
         kopt = k;
         imsemn = imse;
      }
      kimse[k-*kmin] = imse;
   }
/*
     kopt is returned in kmin
*/
   *kmin = kopt;
/*
     compute the bandwidths for kopt
*/
   olafbw(x,delta,n,z,gridz,&kopt,bopt);
   return;
}

/* Function oneolf
	Parameters:
     n --> unique survival times
     x0 --> grid point for which a nearest neighbor bandwidth is to be
            calculated
    nobs --> number of observations (censored and uncensored)
     k --> number of neighbors within the bandwidth
 Computes bandwidth at grid point x0, using the survival matrix x
     x --> survival matrix computed using Kaplan-Meier
           first column: distinct survival times
           second column: survival values
*/

#undef zero
#undef huge
double oneolf(double **x, int *n, double *x0, int *nobs, int *k)
{
#define epsi (1.0e-5)
#define one 1.0e0
#define onemor (one+epsi)
#define oneles (one-epsi)
static double T1,T2;
static double oneolf;
static int i,ix0,ilo,ihi,count;
static double dx[200000],Const,bw,bw0,r,r0,ds;
   ix0 = atpos(&x[0][0],n,x0);
   ilo = fifmax0(1,ix0-*k);
   ihi = fifmin0(*n,ix0+*k);
   count = 0;
/*
     compute the distances from x0 to nearest neighbors
*/
   for(i=ilo-1; i<ihi; i++) {
      count += 1;
      dx[count-1] = fabs(x[0][i]-*x0);
   }
   sorter(dx,&count);
/*
     compute the largest distance so that: S(x0-r)-S(x0+r-0) <= (k-1)/n
*/
   Const = onemor*(double)(*k-1)/(double)*nobs;
   bw = -99.99e0;
   for(i=0; i<count; i++) {
      r = dx[i];
      T1 = *x0-r;
      T2 = *x0+r;
      ds = gets(x,n,&T1)-gets(x,n,&T2);
      if(ds > Const) goto S300;
      else bw = r;
   }
S300:
/*
     Decide from bw, bw+, r-
*/
   bw0 = onemor*bw;
   T1 = *x0-bw0;
   T2 = *x0+bw0;
   ds = gets(x,n,&T1)-gets(x,n,&T2);
   if(ds > Const) {
      oneolf = bw;
      return oneolf;
   }
   r0 = oneles*r;
   T1 = *x0-r0;
   T2 = *x0+r0;
   ds = gets(x,n,&T1)-gets(x,n,&T2);
   if(ds > Const) oneolf = bw0;
   else oneolf = r0;
   return oneolf;
}

/* Function pone
*/
#undef epsi
#undef one
#undef onemor
#undef oneles
double pone(int *j, double *x)
{
static int K1 = 2;
static int K2 = 1;
static double pone,pj1;
   pj1 = (double)(0.5F*(float)(Xalpha+Xbeta+2));
   if(*j == 1) pone = pj1;
   else {
      pone = sigma(&K1)*pzero(&K2,x)+(sigma(&K1)**x+tao(&K1))*pj1+ksi(&K1);
      pone /= luo(&K1);
   }
   return pone;
}

/* Function ptwo
*/
double ptwo(int *j, double *x)
{
/*
 Called with j=2 or j=3
*/
static int K1 = 2;
static int K2 = 1;
static int K3 = 3;
static double ptwo,pj2;
   pj2 = 2.0e0*sigma(&K1)*pone(&K2,x)/luo(&K1);
   if(*j == 2) {
      ptwo = pj2;
      return ptwo;
   }
/*
 --- j = 3
*/
   ptwo = (2.0*sigma(&K3)*pone(&K1,x)+(sigma(&K3)**x+tao(&K3))*pj2)/luo(&K3);
   return ptwo;
}

/* Function pzero
*/
double pzero(int *j, double *x)
{
#define half (1.0e0/2.0e0)
static int K1 = 2;
static int K2 = 3;
static double pzero,pj1,pj2;
   pj1 = half*((double)(Xalpha+Xbeta+2)**x+(double)Xalpha-(double)Xbeta);
   if(*j == 1) {
      pzero = pj1;
      return pzero;
   }
   pj2 = ((sigma(&K1)**x+tao(&K1))*pj1+ksi(&K1))/luo(&K1);
   if(*j == 2) {
      pzero = pj2;
      return pzero;
   }
/*
     --- j = 3
*/
   pzero = ((sigma(&K2)**x+tao(&K2))*pj2+ksi(&K2)*pj1)/luo(&K2);
   return pzero;
}

/* Function sigma
*/
#undef half
double sigma(int *j)
{
extern int Xalpha,Xbeta;
static double sigma,temp;
   temp = (double)(2.0F*(float)*j+Xalpha+Xbeta);
   sigma = temp*(temp-1.0)*(temp-2.0);
   return sigma;
}

/* Function sorter
     Sort the vector v, increasingly
     v - vector to be sorted
	 n - length of v
*/
void sorter(double *v, int *n)
{
static unsigned int qdone;
static int i;
static double temp;
   if(*n == 1) return;
S100:
   qdone = TRUE;
   for(i=1; i<=*n-1; i++) {
      if(v[i-1] > v[i]) {
         qdone = FALSE;
         temp = v[i-1];
         v[i-1] = v[i];
         v[i] = temp;
      }
   }
   if(!qdone) goto S100;
   return;
}

/* Function surfct
 Computes the empirical survival function of the UNCENSORED observation
 In the original code from Mueller, he uses ALL observations
*/
double surfct(double *x, double *delta, int *n, double *xx)
{
static double surfct;
static int index,i;
   index = 0;
   for(i=0; i<*n; i++) {
      if(x[i] <= *xx && delta[i] == 1) index += 1;
   }
   surfct = 1.0e0-(double)index/(double)(*n+1);
   return surfct;
}

/* Function tao
*/ 
double tao(int *j)
{
extern int Xalpha,Xbeta;
static double tao;
   tao = (double)((2.0F*(float)*j+Xalpha+Xbeta-1)*(float)(Xalpha*Xalpha-Xbeta*
     Xbeta));
   return tao;
}

/* Function try
*/ 
void try(int *n,
		 double *x,
		 double *delta,
		 double *z,
		 double *b,
		 double *endl,
		 double *endr,
		 double *q,
		 double *r,
		 double *s,
		 double *valueb,
		 double *valuev,
		 int *iterat,
		 double *bpilot,
		 int *kflag
		 )
{
#define zero 0.0e0
#define half 0.5e0
extern void func();
static int i,it;
static double sumb,sumv,del,xx,br,bs,bxx,vr,vs,vxx,tnm;
   if(*iterat == 1) {
      func(n,x,delta,z,b,endl,endr,q,r,&br,&vr,bpilot,kflag);
      func(n,x,delta,z,b,endl,endr,q,s,&bs,&vs,bpilot,kflag);
      *valueb = half*(*s-*r)*(br+bs);
      *valuev = half*(*s-*r)*(vr+vs);
   }
   else {
      it = fifipow(2,*iterat-2);
      tnm = (double)it;
      del = (*s-*r)/tnm;
      xx = *r+half*del;
      sumb = sumv = zero;
      for(i=1; i<=it; i++) {
         func(n,x,delta,z,b,endl,endr,q,&xx,&bxx,&vxx,bpilot,kflag);
         sumb += bxx;
         sumv += vxx;
         xx += del;
      }
      *valueb = half*(*valueb+(*s-*r)*sumb/tnm);
      *valuev = half*(*valuev+(*s-*r)*sumv/tnm);
   }
   return;
}
#undef zero
#undef half
