#include "defs.h"

/// Search through array 'chi2' and find points that are less than 
/// their N=2x'Ntest' nearest neighbors.
void prior_findmin(double *chi2, long NZ, long *mins, long *Nmin)
{
      long i,j;
      int Ntest,count,flag;
      double chimin,chimin2;
      Ntest=5;
            
      count=0;
      i=Ntest;
      /////// Find local minima
      while (i<NZ-Ntest) {
            flag=0;
            for (j=i-Ntest;j<i;++j) 
                  if (chi2[j]<=chi2[i]) ++flag;
            
            for (j=i+1;j<=i+Ntest;++j)
                  if (chi2[j]<=chi2[i]) ++flag;
            
            if (flag==0) {
                  mins[count]=i;
                  //printf("test %d %lf %ld\n",count,chi2[i],mins[count]);
                  ++count;
                  chimin=chi2[i];
            }
            ++i;
      }
      
      //// check for minima at the edges
      flag=-1;
      chimin2=chi2[Ntest];
      for (i=0;i<Ntest;++i) 
            if (chi2[i]<chimin2) {
                  flag=i;
                  chimin2=chi2[i];
            }
      if (flag>-1) {
            mins[count]=flag;
            ++count;
      }

      flag=-1;                  
      chimin2=chi2[NZ-Ntest-1];
      for (i=NZ-Ntest;i<NZ;++i) 
            if (chi2[i]<chimin2) {
                  flag=i;
                  chimin2=chi2[i];
            }
      if (flag>-1) {
            mins[count]=flag;
            ++count;
      }

      if (count==0) {
            chimin=1.e30;
            mins[0]=0;
            for (i=0;i<Ntest;++i) 
                  if (chi2[i]<chimin) {
                    mins[0]=i;
                    chimin=chi2[i];
                  }
            for (i=NZ-Ntest;i<NZ;++i) 
                  if (chi2[i]<chimin) {
                    mins[0]=i;
                    chimin=chi2[i];
                  }
                  
            *Nmin=1;
      } else *Nmin=count;      
}

void get_prior_file_size(long *NZ_prior, long *NK_prior)
{
      char *arg,buff[1024];
      char delim[]=" \n\t,\r";
      long nzp,nkp;
      int lineflag;
      FILE *fp;
      
      //// PRIOR_FILE = "kprior_zmax7.dat";
      
      if (!(fp = fopen(PRIOR_FILE,"r"))) {
            fprintf(stderr,"Error opening prior definition file, %s!\n",
               PRIOR_FILE);
            exit(1);
      }
      
      //// set lineflag=1 when you find a column where the first column is called 'z'
      lineflag=0;
      while (lineflag == 0 && fgets(buff,1024,fp)) {
          arg = strtok(buff,delim);
          while (arg[0] == '#') ++arg;
          if (arg[0] != '\0') {
            if (!(strcmp(arg,"z"))) lineflag=1;
          } else {
            arg = strtok(NULL,delim);
            if (!(strcmp(arg,"z"))) lineflag=1;
          }         
      }
      if (lineflag == 0) {
          fprintf(stderr,"No header line found in PRIOR_FILE.\n");
          fprintf(stderr,"The first non-# column in the header line needs to be z.\n");
          exit(1);
      }
      
      /// Count K bins
      nkp=0;
      while ( (arg=strtok(NULL,delim)) ) {
        ++nkp;
//        printf("arg:   %s\n",arg);
      }
      
      /// Count z lines
      nzp = 0;
      while(fgets(buff,1024,fp)) {
            ++nzp;
      }
      fclose(fp);
      
      *NZ_prior=nzp;
      *NK_prior=nkp;
            
}

////// Format of prior file: header line showing mag bins
////// additional lines each with z, p(z|Mag)
//////
////// e.g.
//////  z   18   20    22
////// 0.5  0.9  0.6  0.3
////// 2.5  0.0  0.4  0.5
////// ...
//////
////// each column is therefore the redshift distribution 
////// at the magnitude defined in the first line
void read_prior_file(double **priorkz, double *priorzval, double *klim, 
                     long nzp, long nkp)
{
      char *arg,buff[1024];
      char delim[]=" \n\t,\r";
      int i,j,lineflag;
      double dtemp;
      FILE *fp;
            
      if (!(fp = fopen(PRIOR_FILE,"r"))) {
            fprintf(stderr,"Error opening prior definition file, %s!\n",
               PRIOR_FILE);
            exit(1);
      }
      
      ///// Skip # comment lines and read Kmag header line                       
      lineflag=0;
      while (lineflag == 0 && fgets(buff,1024,fp)) {
          arg = strtok(buff,delim);
          while (arg[0] == '#') ++arg;
          if (arg[0] != '\0') {
            if (!(strcmp(arg,"z"))) lineflag=1;
          } else {
            arg = strtok(NULL,delim);
            if (!(strcmp(arg,"z"))) lineflag=1;
          }         
      }
      
      ////// Read K mags from header line
      j=0;
      while ( (arg=strtok(NULL,delim)) ) {
            sscanf(arg,"%lf",&dtemp);
            klim[j] = dtemp;
            //printf("klim: %lf %ld %d\n",klim[j],nkp,j);
            ++j;
      }
      
      ///// Read in prior grid
      j=0;
      while(fgets(buff,1024,fp)) {
            arg = strtok(buff,delim);
            sscanf(arg,"%lf",&dtemp);
            priorzval[j]=dtemp;
            //printf("z: %lf\n",priorzval[j]);
            i=0;
            while ( (arg=strtok(NULL,delim)) ) {
                  sscanf(arg,"%lf",&dtemp);
                  priorkz[i][j]=dtemp;
                  //printf("dtemp: %e\n", priorkz[i][j]);
                  ++i;
            }            
            ++j;
      }                  

}

///////// Interpolate the prior grid to the redshift grid
///////// defined in zphot.param.
void interpolate_prior(double **priorkz_first, double **priorkz_out, double *priorzval,
                       long nzp, long nkp)
{
    long mag,iz,istart;

    for (mag=0;mag<nkp;++mag) {
        
        
        istart=0;
        for (iz=0;iz<NZ;++iz) {
            //////// No extrapolation
            if (ztry[iz] > priorzval[nzp-1] || ztry[iz] < priorzval[0]) {
            
                priorkz_out[mag][iz] = 1.;
     //           printf("extrapolate: %lf %lf %lf\n",ztry[iz],priorzval[nzp-1],priorzval[0]);
            } else { 
                interpol(ztry[iz],&priorkz_out[mag][iz],
                               priorzval,priorkz_first[mag],nzp,&istart);  
      //          printf("prior %ld %lf %lf\n",mag,ztry[iz],priorkz_out[mag][iz]);
            }
            //// Enforce minimum prior
            //if (priorkz_out[mag][iz] < 1.e-20) priorkz_out[mag][iz] = 1.e-20;
        }  
    }
    
}

//// Multiply p(z) by a wiggle centered at z=1.61 to get rid of z=1.6 NMBS bump
void tweak_z16_gauss (double *pz) {
    long i,j; 
    double tweak,c236,ctot,cfact;
        
    for (i=0;i<NZ-1;++i) {

        ctot=0;
        for (j=0;j<6;++j) ctot+=coeffs[i][j];
        c236 = coeffs[i][1]+coeffs[i][2]+coeffs[i][5];
        cfact = c236/ctot;
        cfact = 1;

        tweak = 1-0.99*exp(-1.*pow(ztry[i]-1.61,2)/pow(0.1,2))*cfact;
        tweak += 0.05*exp(-1.*pow(ztry[i]-1.44,2)/pow(0.07,2))*cfact;
        tweak += 0.05*exp(-1.*pow(ztry[i]-1.78,2)/pow(0.07,2))*cfact;

        pz[i]*=tweak;
    }
    
}

//// Multiply X-band prior to observed p(z) [from chi2] and find where p(z) maximized
void apply_prior(double **priorkz, long nzp, long nkp, long first_best,
                 double *chi2, double kflux, long *bestz, double *pzout)
{

      long i,Kbinlo,Kbinhi;
      double maxprob,kmag,norm,factor_lo,factor_hi;
      
      //// Handle negative fluxes in the prior band
      //// by setting the magnitude to very low value.  Shouldn't 
      //// be necessary because neg fluxes here are handled in main.c.
      if (kflux>0) 
        kmag = PRIOR_ABZP-2.5*log10(kflux);
      else
        kmag = 40;
      // printf("kmag, kflux: %lf %lf\n",kmag,kflux);

      //// get K column
      Kcolumn=0;
      while (klim[Kcolumn]<kmag && Kcolumn<(nkp-1)) ++Kcolumn;
      // printf("Kmag_col: %lf %lf\n",kmag,klim[Kcolumn]); exit(1);

      /// Interpolate prior as function of mag
      if (Kcolumn == 0) {
                Kbinlo = 0;
                Kbinhi = 1;
                factor_lo = 1;
                factor_hi = 0;
      } else {
                Kbinlo = Kcolumn-1;
                Kbinhi = Kcolumn;
                factor_lo = klim[Kbinhi]-kmag;
                factor_lo/=(klim[Kbinhi]-klim[Kbinlo]);
                factor_hi = kmag-klim[Kbinlo];
                factor_hi/=(klim[Kbinhi]-klim[Kbinlo]);
      }
      //printf("interp prior: %lf %lf %lf\n",klim[Kbinlo],kmag,klim[Kbinhi]);

      ////// find minimum chisq
      // chimin=1.e30;
      // chimin_i=0;
      // for (i=0;i<NZ;++i) if (chi2[i]<chimin) {
      //       chimin=chi2[i];
      //       chimin_i = i;
      // }
      //// Have first_best = izbest from main.c
                        
	  //// Subtract the minimum chisq from whole array so maximum probability is 
      //// one and won't have problems with the probabilities being too small.
      //// 'pzout' is a probability, 'chi2' is chisq
      norm = 0;
      pzout[0] = exp(-0.5*(chi2[0]-chi2[first_best])/CHI2_SCALE)*(priorkz[Kbinlo][0]*factor_lo+priorkz[Kbinhi][0]*factor_hi);
      
      ///////// Multiply by probabilities by (1+z)
      //pzout[0] *= (1+ztry[0]);
      
      maxprob = pzout[0];
      *bestz = 0;
      
      for (i=1;i<NZ;++i) {
        //pzout[i] = exp(-0.5*(chi2[i]-chi2[*bestz])/CHI2_SCALE)*priorkz[Kcolumn][i];
        pzout[i] = exp(-0.5*(chi2[i]-chi2[first_best])/CHI2_SCALE)*(priorkz[Kbinlo][i]*factor_lo+priorkz[Kbinhi][i]*factor_hi);
        //pzout[i] *= (1+ztry[i]);
      }
           
      // tweak_z16_gauss(pzout);
      
      norm = 0;
      for (i=1;i<NZ;++i) {
          if (pzout[i]>maxprob) {
      		  maxprob=pzout[i];
      		  *bestz=i;
      	  }       
          norm += (ztry[i]-ztry[i-1])*(pzout[i-1]+pzout[i]);        
      }
      for (i=0;i<NZ;++i) pzout[i] /= (norm/2.);  //// Normalize total probability to unity.
      
}