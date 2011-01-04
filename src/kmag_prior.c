#include "defs.h"

/// Search through array 'probz' and find points that are less than 
/// their N=2x'Ntest' nearest neighbors.
void prior_findmin(double *probz, long NZ, long *mins, long *Nmin)
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
                  if (probz[j]<=probz[i]) ++flag;
            
            for (j=i+1;j<=i+Ntest;++j)
                  if (probz[j]<=probz[i]) ++flag;
            
            if (flag==0) {
                  mins[count]=i;
                  //printf("test %d %lf %ld\n",count,probz[i],mins[count]);
                  ++count;
                  chimin=probz[i];
            }
            ++i;
      }
      
      //// check for minima at the edges
      flag=-1;
      chimin2=probz[Ntest];
      for (i=0;i<Ntest;++i) 
            if (probz[i]<chimin2) {
                  flag=i;
                  chimin2=probz[i];
            }
      if (flag>-1) {
            mins[count]=flag;
            ++count;
      }

      flag=-1;                  
      chimin2=probz[NZ-Ntest-1];
      for (i=NZ-Ntest;i<NZ;++i) 
            if (probz[i]<chimin2) {
                  flag=i;
                  chimin2=probz[i];
            }
      if (flag>-1) {
            mins[count]=flag;
            ++count;
      }

      if (count==0) {
            chimin=1.e30;
            mins[0]=0;
            for (i=0;i<Ntest;++i) 
                  if (probz[i]<chimin) {
                    mins[0]=i;
                    chimin=probz[i];
                  }
            for (i=NZ-Ntest;i<NZ;++i) 
                  if (probz[i]<chimin) {
                    mins[0]=i;
                    chimin=probz[i];
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
        }  
    }
    
}
          
void apply_prior(double **priorkz, double *priorzval, double *klim, 
                     long nzp, long nkp,
                 double *probz, double kflux, long *bestz, double *pzout)
{

  //    long *mins,Nmin,j,chimin_i;
      long i;
      double maxprob,prob_i,kmag,norm;
      
      //// Handle negative fluxes in the prior band
      //// by setting the magnitude to very low value
      if (kflux>0) 
        kmag = PRIOR_ABZP-2.5*log10(kflux);
      else
        kmag = 40;

      // printf("kmag, kflux: %lf %lf\n",kmag,kflux);

//      mins = malloc(sizeof(long)*nzp);
      
//      prior_findmin(probz,NZ,mins,&Nmin);

      ////// find minimum chisq
/*      chimin=1.e30;
      chimin_i=0;
      for (i=0;i<Nmin;++i) if (probz[mins[i]]<chimin) {
            chimin=probz[mins[i]];
            chimin_i = mins[i];
      }
*/      
      ///// force minima to be within 3-sigma contour
/*      j=-1;
      for (i=0;i<Nmin;++i) if (probz[mins[i]]<=chimin+9) mins[++j]=mins[i];
      Nmin=j+1;
*/

      ///// Only adjust bestz if original zphot>1.5 (TURNED OFF)
      //if (ztry[chimin_i]>1.5) *bestz=chimin_i; else {
      
            //// get K column
            Kcolumn=0;
            while (klim[Kcolumn]<kmag && Kcolumn<(nkp-1)) ++Kcolumn;
            // printf("Kmag_col: %lf %lf\n",kmag,klim[Kcolumn]); exit(1);
            
            // for (i=0;i<nzp;++i) printf("prior: %lf %lf\n",priorzval[i],priorkz[Kcolumn][i]); exit(1);
            
/*            maxprob=-1;
            for (i=0;i<Nmin;++i) {
                  ///// scale chisq to min(chisq)=1
                  // prob_i=exp(-1*probz[mins[i]]/chimin)*priorkz[Kcolumn][mins[i]];
                  ///// Likelihood = exp(-1/2 * chisq)
                  prob_i=exp(-0.5*probz[mins[i]]/CHI2_SCALE)*priorkz[Kcolumn][mins[i]];
                  if (prob_i>maxprob) {
                        maxprob=prob_i;
                        *bestz=mins[i];
                  }
            }
*/            
      //}
      
      //// Don't force 3-sigma contour or that best is one of the original peaks
      maxprob = -1;
      for (i=0;i<NZ;++i) {
      	prob_i = exp(-0.5*probz[i]/CHI2_SCALE)*priorkz[Kcolumn][i];
      	if (prob_i>maxprob) {
      		maxprob=prob_i;
      		*bestz=i;
      	}
      }

	  //// Subtract the minimum chisq from whole array so maximum probability is 
      //// one and won't have problems with the probabilities being too small.
      //// Then normalize integrated probability to equal 1.
      //// 'pzout' is a probability, 'probz' is chisq
      norm = 0;
      pzout[0] = exp(-0.5*(probz[0]-probz[*bestz])/CHI2_SCALE)*priorkz[Kcolumn][0];
      for (i=1;i<NZ;++i) {
        pzout[i] = exp(-0.5*(probz[i]-probz[*bestz])/CHI2_SCALE)*priorkz[Kcolumn][i];
        norm+=(ztry[i]-ztry[i-1])*(pzout[i-1]+pzout[i]);
      }
      for (i=0;i<NZ;++i) pzout[i]/=norm/2.;
      
//      free(mins);
}
                  
      
      

      
      
            
            
