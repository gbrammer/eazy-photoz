/*
 * This is from SM
 * 
 */
#include <stdio.h>
#include <stdlib.h>

#define REAL double

  

/*********************************************/
/*
 * Spline package for interpolation.
 *
 * Written by Jim Gunn, modified by Robert Lupton
 */
#include <stdio.h>
#define d12 .0833333333
static double ddb, dde;		/* second derivatives at each end of spline */

/*****************************************************/
/*
 * returns index i of first element of t >= x
 */
int
circe(double x,double *t, long n)
{
   register long lo, hi, mid ;
   double tm ;
   lo = 0 ;
   hi = n-1 ;
   if ( x >= t[0] && x <= t[n-1])  {
      while ( hi - lo > 1 ) {
	 mid = ( lo + hi )/2. ;
	 tm = t[mid] ;
	 if ( x < tm ) hi = mid ;
	 else  lo = mid ;
      }
      return(hi);
   } else {
      return(-1);
   }
}

/**********************************************************/
/*
 * solves tridiagonal linear system . Diag elts mii=ai, subdiag mii-1=bi,
 * superdiag mii+1=ci.
 * it is understood that b0 and cn-1 are zero, but are not referenced.
 * f is rhs, x soln, may be same array; all arrays are trashed
 *
 * Re-written; JEG version incorrect.
 */
void
tridi(double *a,double *b,double *c,double *f,double *x, long n)
{
   long i;

   c[0] /= a[0];
   for(i=1;i < n;i++) {
      c[i] /= (a[i] - b[i]*c[i-1]);
   }
   f[0] /= a[0];
   for(i=1;i < n;i++) {
      f[i] = (f[i] - b[i]*f[i-1])/(a[i] - b[i]*c[i-1]);
   }
   x[n - 1] = f[n - 1];
   for(i = n - 2;i >= 0;i--) {
      x[i] = f[i] - c[i]*x[i+1];
   }
}

/*********************************************************/
/*
 * sets up spline derivative array k for a given x and y array of length
 * n points, n-1 intervals, for given estimates for the second derivatives
 * at the endpoints, q2b and q2e; "natural" boundary conditions for q2b=q2e=0
 */
int
initialise_spline(double *x,double *y,double *k,long n,double q2b,double q2e)
{
   REAL hio, hip, dio, dip ;
   register REAL *a, *b, *c, *f;
   long i, ip ;

   ddb = q2b; dde = q2e;		/* save end second derivatives */
   if((a = (double *)malloc(4*(unsigned)n*sizeof(double))) == NULL) {
      fprintf(stderr,"Can't malloc space in initialise_spline\n");
      exit(1);
   }
   b = a + n ;
   c = b + n ;
   f = c + n ;

   hio = 0. ;
   dio = 0. ;

/*	BIG BUG in borland compiler: ip is used before its computation !!! */
#ifdef	__BORLANDC__
	for( i=0, ip = 1 ; i < n ; i++, ip++ ) {
      hip = ( ip < n ? x[ip] - x[i] : 0. ) ;
      dip = ( ip < n ? (y[ip] - y[i])/hip : 0. ) ;
      b[i] = ( ip < n ? hip : hio ) ;
      a[i] = 2.*( hip + hio ) ;
      c[i] = ( i > 0 ? hio : hip ) ;
      f[i] = 3.*(hip*dio + hio*dip ) ;
      if( i == 0 ) f[0] = 3.* hip * dip  - hip * hip * q2b * 0.5 ;
      else if ( i == n-1 )
	 f[n-1] = 3.* hio* dio + hio* hio* q2e* 0.5  ;
      dio = dip ;
      hio = hip ;
   }
#else
	for( i=0 ; i < n ; i++ ) {
      hip = ((ip = i+1) < n ? x[ip] - x[i] : 0. ) ;
      dip = ( ip < n ? (y[ip] - y[i])/hip : 0. ) ;
      b[i] = ( ip < n ? hip : hio ) ;
      a[i] = 2.*( hip + hio ) ;
      c[i] = ( i > 0 ? hio : hip ) ;
      f[i] = 3.*(hip*dio + hio*dip ) ;
      if( i == 0 ) f[0] = 3.* hip * dip  - hip * hip * q2b * 0.5 ;
      else if ( i == n-1 )
	 f[n-1] = 3.* hio* dio + hio* hio* q2e* 0.5  ;
      dio = dip ;
      hio = hip ;
   }
#endif
   tridi(a,b,c,f,k,n) ;
   free((char *)a);
   return(0);
}

/*********************************************************/
/*
 * general spline evaluator; xyk are the x,y, and derivative arrays, of
 * dimension n, dx the argument = (x-xi-1), m the index of the next GREATER
 * abscissa.
 */
int
interpolate_spline(double z,double *val,double *x,double *y,double *k,long n)
{
   REAL h, t, d, a, b,
	 dx ;
   long m ;

   if((m = circe( z, x, n)) < 0) {	/* out of bounds */
      if(z < x[0]) {
	 dx = z - x[0];
	 *val = y[0] + dx*(k[0] + 0.5*dx*ddb);
      } else {
	 dx = z - x[n-1];
	 *val = y[n-1] + dx*(k[n-1] + 0.5*dx*dde);
      }
      return(-1);
   }
   dx = z - x[m-1];
   
   h = x[m] - x[m-1];
   d = (y[m] - y[m-1])/h;
   t = dx/h;
   a = (k[m-1] - d)*(1 - t);
   b = (k[m]-d) *t ;
   *val = t*y[m] + (1 - t)*y[m-1] + h*t*(1 - t)*(a - b);
   return(0);
}

/****
* Simple linear interpolation
*
* 'istart' speeds up the routine significantly when interpolating
* over a full grid of 'xval's.  Set istart=0 for arbitrary interpolations 
* where you're not sure that x[istart] < xval.
* 
* The routine as it is assumes that xval < max(x).  Change the 'while' statements
* for more general behavior.
*
*****/
int interpol(double xval,double *out, double *x, double *y, long n, long *istart)
{
    double x1,x2,y1,y2;
    long i;
    
    //// Start and index, istart.  Increment i until you
    //// find the x[i] position that is greater than 
    //// the interpolated position.  
    i=*istart;
    while (x[i]<xval) ++i;
    
    //// x[i] is now greater than xval because the 
    //// expression (x[i]<xval) is false, assuming
    //// that xval < max(x).
    //// -or-
    //// i=istart if x[istart] > xval, which shouldn't happen.
      
    /*
    while (x[i]<xval) {
          ++i;
        if ( i==n ) {
            fprintf(stderr,"interpol() out of range: xmax=%lf, target=%lf\n",x[n-1],xval);
            return(1);
        }
    }
    */
      
    if (i > *istart) --i;
    x1 = x[i];
    x2 = x[i+1];
    y1 = y[i];
    y2 = y[i+1];
    
    *out = ((y2-y1)/(x2-x1))*(xval-x1)+y1;
    *istart = i;

    return(0);
}

