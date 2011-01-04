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

