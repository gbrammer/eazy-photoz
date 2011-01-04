#include "defs.h"

double CHI2_SCALE;

// ** control parameters from zphot.params

int NTEMP;  // to be defined from template file
char TEMPLATES_FILE[1024];
int NTEMPL;

//redshift grid
int NZ;
double Z_STEP;
int Z_STEP_TYPE; 
double Z_MIN;
double Z_MAX; 


// filters
char FILTERS_RES[1024];
int SMOOTH_FILTERS;
double SMOOTH_SIGMA;
int FILTER_FORMAT;

// input data
char CATALOG_FILE[1024];
double NOT_OBS_THRESHOLD;
int N_MIN_COLORS;
char TEMP_ERR_FILE[1024];
char WAVELENGTH_FILE[1024];

// output
char OUTPUT_DIRECTORY[1024];
char MAIN_OUTPUT_FILE[1024];
int OBS_SED_FILE;
int TEMP_SED_FILE;
int POFZ_FILE;
int VERBOSE_LOG;
int BINARY_OUTPUT;

int PRINT_ERRORS;

// cosmology
double H0;
double OMEGA_M;
double OMEGA_L;

// template caching
int DUMP_TEMPLATE_CACHE;
int USE_TEMPLATE_CACHE;
char CACHE_FILE[1024];

// template fitting
int TEMPLATE_COMBOS;
double NMF_TOLERANCE;
double TEMP_ERR_A2;
double SYS_ERR;
int APPLY_IGM;
int FIX_ZSPEC; 

char PRIOR_FILE[1024];
int APPLY_PRIOR;
double PRIOR_ABZP;
int PRIOR_FILTER;

// start of filter throughput data linked list
filt_data filt_thru;

long nfilter,*usefilt,nusefilt,totfiltid,totcolumn;

double *totflux;

long nobj;
char **objid;
int idsize;

double **fnu, **efnu, *zspec, *ra, *dec;
double *lambdac,*ztry, *zspec;
double *templ,**tempf,*tempage;
double temp_err[NTEMPLMAX],**temp_errf, temp_scale[NTEMPLMAX];
double ***tempfilt, *temp_err_a, *tnorm;

int **temp_combine;

filt_data **pusefilt;

int CHECK_ZP_OFFSETS;
double ZP_OFFSET_TOL;

int *zpfilterid, zpcontinue, *filt_to_col;
double *zpfactor;

FILE *fplog;

double *dasum,*dbsum;

int iz,*idtemp1,*idtemp2a,*idtemp2b,**idtempall,ntemp_all,ngoodfilters,ncols;
double *pz1,*pz2,*pzall,*pzuse,*pzout,*pzsum,*zinverse;
double *atemp1,*atemp2a,*atemp2b,**coeffs;

double **priorkz,*klim;
long Kcolumn;

int NZ_resid, NF_resid;
double *zresid,*dzresid,**resid_table,**resid_int;

int32_t nusefilt32, NTEMP32, NZ32, nobj32, NZ32, *izsave32, NTEMPL32;

int main() {

	char coeff_file[1024],pz_file[1024];
	FILE *cf, *pzf;
		
    long i,j,iobj,izbest,iz_zm1,iz_zm2, *izsave;
    double zm1,zm2;
    FILE *fp;
    char ofile[1024],oparfile[1024],outfilename[1024];
    char strfmt[64];
    
    int getzStatus;

    double weighted_z (double *pofz);
    double weighted_z2 (double *pofz);
    double odds (double *pofz);
    // double calculate_mass (long zidx, double **coeffs);

    double l68,u68,l95,u95,l99,u99;
    
    double **priorkz_raw, *priorzval;
    long NZ_prior,NK_prior,izprior;

    
    time_t then,now;

    time(&then);
  
    printf("Reading in program parameters...\n\n");
    getparams();

    if (VERBOSE_LOG) {
        strcpy(oparfile,OUTPUT_DIRECTORY);
        strcat(oparfile,"/");
        strcat(oparfile,MAIN_OUTPUT_FILE);
        strcat(oparfile,".param");
        //printf("%s\n",oparfile);
        if (!(fplog = fopen(oparfile,"w"))) {
            fprintf(stderr,"\n\nuh,oh. Couldn't open %s.  Does %s exist?\n",oparfile,OUTPUT_DIRECTORY);
            exit(1);
        }
        printparams_logfile();
    }
        
    printf("Initializing...");
    init();
    
    printf("  Phew!! Ready to go.\n");

    ////// Some options don't apply IF redshift fixed to z_spec --> turn them off
    if (FIX_ZSPEC) {
        if (POFZ_FILE) fprintf(stderr,"FIX_ZSPEC set, ignoring POFZ_FILE ...\n");
        if (APPLY_PRIOR) fprintf(stderr,"FIX_ZSPEC set, ignoring APPLY_PRIOR ...\n");
        if (PRINT_ERRORS) fprintf(stderr,"FIX_ZSPEC set, ignoring PRINT_ERRORS ...\n");  
        APPLY_PRIOR=0;
        PRINT_ERRORS=0;
    }
    
    /////// Initialize Prior
    if (APPLY_PRIOR && FIX_ZSPEC==0) {
        
        Kcolumn=0;
        //// Get size of grid
        get_prior_file_size(&NZ_prior, &NK_prior);
        printf("Read PRIOR_FILE NZ: %ld NK: %ld \n",NZ_prior,NK_prior);
        
        klim = malloc(sizeof(double)*NK_prior);
        priorzval = malloc(sizeof(double)*NZ_prior);
      
        priorkz_raw = malloc(sizeof(double *)*NK_prior);
        for (i=0;i<NK_prior;++i) priorkz_raw[i] = malloc(sizeof(double)*NZ_prior);

        priorkz = malloc(sizeof(double *)*NK_prior);
        for (i=0;i<NK_prior;++i) priorkz[i] = malloc(sizeof(double)*NZ);
      
        //// Read the grid
        read_prior_file(priorkz_raw,priorzval,klim,NZ_prior,NK_prior);
        
        //// Interpolate Prior grid to redshift grid
        interpolate_prior(priorkz_raw,priorkz,priorzval,NZ_prior,NK_prior);
        
        free(priorkz_raw);        
    }
    
    pzuse = malloc(sizeof(double)*NZ);
    if (pzuse==NULL) {
        fprintf(stderr,"pzuse: memory allocation error!\n");
        exit(1);
    }
    pzout = malloc(sizeof(double)*NZ);
    if (pzout==NULL) {
        fprintf(stderr,"pzout: memory allocation error!\n");
        exit(1);
    }

    pzsum = malloc(sizeof(double)*NZ);
    if (pzsum==NULL) {
        fprintf(stderr,"pzsum: memory allocation error!\n");
        exit(1);
    }

    zinverse = malloc(sizeof(double)*NZ);
    if (zinverse==NULL) {
        fprintf(stderr,"zinverse: memory allocation error!\n");
        exit(1);
    }
        
    //////// Initialize p(z), template(z), normalization(z) vectors 
    //////// depending on linear combination options
    if (TEMPLATE_COMBOS <= 2) {
        pz1 = malloc(sizeof(double)*NZ);
        if (pz1==NULL) {
            fprintf(stderr,"pz1: memory allocation error!\n");
            exit(1);
        }
        for (i=0;i<NZ;++i) pz1[i] = 1.e30;
        idtemp1 = malloc(sizeof(int)*NZ);
        if (idtemp1==NULL) {
            fprintf(stderr,"idtemp1: memory allocation error!\n");
            exit(1);
        }
        atemp1 = malloc(sizeof(double)*NZ);
        if (atemp1==NULL) {
            fprintf(stderr,"atemp1: memory allocation error!\n");
            exit(1);
        }       
    }
    
    if (TEMPLATE_COMBOS == 2) {
        pz2 = malloc(sizeof(double)*NZ);
        if (pz2==NULL) {
            fprintf(stderr,"pz2: memory allocation error!\n");
            exit(1);
        }
        for (i=0;i<NZ;++i) pz2[i] = 1.e30;
        idtemp2a = malloc(sizeof(int)*NZ);
        if (idtemp2a==NULL) {
            fprintf(stderr,"idtemp2a: memory allocation error!\n");
            exit(1);
        }
        idtemp2b = malloc(sizeof(int)*NZ);
        if (idtemp2b==NULL) {
            fprintf(stderr,"idtemp2b: memory allocation error!\n");
            exit(1);
        }
        atemp2a = malloc(sizeof(double)*NZ);
        if (atemp2a==NULL) {
            fprintf(stderr,"atemp2a: memory allocation error!\n");
            exit(1);
        }       
        atemp2b = malloc(sizeof(double)*NZ);
        if (atemp2b==NULL) {
            fprintf(stderr,"atemp2b: memory allocation error!\n");
            exit(1);
        }       
    }

    coeffs = malloc(sizeof(double *)*NZ);
    for (i=0;i<NZ;++i) coeffs[i] = malloc(sizeof(double)*NTEMP);
    if (coeffs==NULL) {
    	fprintf(stderr,"coeffs: memory allocation error!\n");
    	exit(1);
    }
    
    if (TEMPLATE_COMBOS == 99) {
        pzall = malloc(sizeof(double)*NZ);
        if (pzall==NULL) {
            fprintf(stderr,"pzall: memory allocation error!\n");
            exit(1);
        }
        for (i=0;i<NZ;++i) pzall[i] = 1.e30;
        ///// Get template normalizations out of getphotz
        idtempall = malloc(sizeof(int *)*NZ);
        for (i=0;i<NZ;++i) idtempall[i] = malloc(sizeof(int)*NTEMP);
        if (idtempall==NULL) {
            fprintf(stderr,"idtempall: memory allocation error!\n");
            exit(1);
        }
    }
    
    ////// Initialise output file
    strcpy(ofile,OUTPUT_DIRECTORY);
    strcat(ofile,"/");
    strcat(ofile,MAIN_OUTPUT_FILE);
    strcat(ofile,".zout");
    if (!(fp = fopen(ofile,"w"))) {
        fprintf(stderr,"Oops...can't open the output file, %s\n",
            ofile);
        exit(1);
    }
 
    printf("Getting photometric redshifts...\n");

    ////// Header line for MAIN_OUTPUT_FILE, add columns for PRIOR output if necessary
    printf("# id z_spec");
    fprintf(fp,"# id z_spec");
    
    switch (TEMPLATE_COMBOS) {
      case 1:
        printf(" z_1 z_m1 chi_1 temp_1");
        fprintf(fp," z_1 z_m1 chi_1 temp_1");
        if (APPLY_PRIOR && FIX_ZSPEC==0) {
            printf(" z_p chi_p temp_p z_m2 odds");
            fprintf(fp," z_p chi_p temp_p z_m2 odds");
        }
        ncols=4+5*(APPLY_PRIOR)+6*(PRINT_ERRORS);
        break;
      case 2:
        printf(" z_2 z_m1 chi_2 temp_2a temp_2b");
        fprintf(fp," z_2 z_m1 chi_2 temp2a temp_2b");
        if (APPLY_PRIOR && FIX_ZSPEC==0) {
            printf(" z_p chi_p temp_pa temp_pb z_m2 odds");
            fprintf(fp," z_p chi_p temp_pa temp_pb z_m2 odds");
        }
        ncols=5+6*(APPLY_PRIOR)+6*(PRINT_ERRORS);
        break;

      case 99:  
        printf(" z_a z_m1 chi_a");
        fprintf(fp," z_a z_m1 chi_a");
        if (APPLY_PRIOR && FIX_ZSPEC==0) {
            printf(" z_p chi_p z_m2 odds");
            fprintf(fp," z_p chi_p z_m2 odds");
        }
        ncols=3+4*(APPLY_PRIOR)+6*(PRINT_ERRORS);
        break;
    }
    if (PRINT_ERRORS && FIX_ZSPEC==0) {
        fprintf(fp," l68 u68  l95 u95  l99 u99 ");
    }
    
    fprintf(fp," nfilt");
    ncols+=1;
    
    printf("\n");
    fprintf(fp,"\n");

    ////// VERSION INFORMATION
    fprintf(fp,"# EAZY v1.00 (July 9, 2008)\n");
    
    //////// Apply zero point corrections
    if (GET_ZP_OFFSETS) zeropoints();
        
    //// Dump coefficients to binary file
    izsave = malloc(sizeof(long)*nobj);
    if (BINARY_OUTPUT) {
        nusefilt32 = (int32_t) nusefilt;
        NTEMP32 = (int32_t) NTEMP;
        NZ32 = (int32_t) NZ;
        nobj32 = (int32_t) nobj;
        
    	sprintf(coeff_file,"%s/%s.coeff",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE);
    	cf = fopen(coeff_file,"w");
        fwrite(&nusefilt32,sizeof(int32_t),1,cf);
        fwrite(&NTEMP32,sizeof(int32_t),1,cf);
        fwrite(&NZ32,sizeof(int32_t),1,cf);
        fwrite(&nobj32,sizeof(int32_t),1,cf);
    
    
    //// Dump p(z) to binary
    //if (POFZ_FILE) {
    	sprintf(pz_file,"%s/%s.pz",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE);
    	pzf = fopen(pz_file,"w");
        fwrite(&NZ32,sizeof(int32_t),1,pzf);
        fwrite(&nobj32,sizeof(int32_t),1,pzf);
 	}    	
        
    ////////
    //// Main loop to get photz for each object in the catalog
    ////////    
    for (iobj=0; iobj<nobj; ++iobj) { 
      izsave[iobj] = -1;
      ngoodfilters = 0;
      for(j=0;j<nusefilt;++j)
        if (fnu[iobj][j] > NOT_OBS_THRESHOLD && efnu[iobj][j] > 0) {
          ++ngoodfilters;
      }
      
      //////// If too few filters, skip  
      if (ngoodfilters < N_MIN_COLORS) {
        printf("%s: only %d filters available, skipping...\n",objid[iobj],ngoodfilters);
        fprintf(fp,"%s %lf",objid[iobj],zspec[iobj]);
        for (j=0;j<ncols;++j) fprintf(fp," -99");
        fprintf(fp,"\n");
		for (j=0;j<NZ;++j) pzout[i] = 0;
		for (j=0;j<NTEMP;++j) coeffs[0][j] = 0;
		if (BINARY_OUTPUT) {
			fwrite(pzout,sizeof(double),NZ,pzf);
			fwrite(coeffs[0],sizeof(double),NTEMP,cf);
		}
		izsave[iobj] = -1;
      } else { 

        /////////////////  Call out to getphotz to get Chisq as a function of z //////////////////
        ////////// Unused arrays are not malloc-ed, just those appropriate for TEMPLATE_COMBOS
        getzStatus = getphotz(iobj, pz1, idtemp1, atemp1, 
                       pz2, idtemp2a, idtemp2b, atemp2a, atemp2b,
                       pzall, idtempall, coeffs, &ntemp_all);
                
        if (getzStatus == 1) {
            printf("%s: oops, got singular matrix for template normalization!  skipping...\n",objid[iobj]);
            fprintf(fp,"%s %lf",objid[iobj],zspec[iobj]);
            for (j=0;j<ncols;++j) fprintf(fp," -99.9"); 
            fprintf(fp,"\n");
            continue; //// next object
        }     
        //////////////////////////////////////////////////////////////////////////////////////////

        /////// Find chisq minima and print results to the terminal and MAIN_OUTPUT_FILE
        /////// Also reset chisq vectors for when FIX_ZSPEC==true and not all Pofz values are replaced

        //// add spaces to id column so all the other columns line up in output file
        strcpy(strfmt,objid[iobj]);
        i=0; while (strfmt[i] != '\0') ++i;
        --i; while (i<=idsize && i < 62) strfmt[++i] = ' '; // max size of objid: IDBUFF=64
        strfmt[i] = '\0';
            
        switch (TEMPLATE_COMBOS) {
          case 1:
            izbest=0;
            for (i=0;i<NZ;++i) {
                pzuse[i] = pz1[i];
                pzout[i] = pzuse[i];
                if (pz1[i]<pzuse[izbest]) izbest=i;
                pz1[i]=1.e30;
            }
            fprintf(fp,"%s %7.4lf %6.3lf %6.3lf %e %5d ",
                strfmt,zspec[iobj],ztry[izbest],weighted_z(pzuse),pzuse[izbest],idtemp1[izbest]+1);
            printf("%s %7.4lf %6.3lf %6.3lf %e %5d ",
                strfmt,zspec[iobj],ztry[izbest],weighted_z(pzuse),pzuse[izbest],idtemp1[izbest]+1);
            break;
          case 2:
            izbest=0;
            for (i=0;i<NZ;++i) {
                pzuse[i] = pz2[i];
                pzout[i] = pzuse[i];
                if (pz2[i]<pzuse[izbest]) izbest=i;
                pz2[i]=1.e30;
             }            
            fprintf(fp,"%s %7.4lf %6.3lf %6.3lf %e %5d %5d ",
                strfmt,zspec[iobj],ztry[izbest],weighted_z(pzuse),pzuse[izbest],idtemp2a[izbest]+1,idtemp2b[izbest]+1);
            printf("%s %7.4lf %6.3lf %6.3lf %e %5d %5d ",
                strfmt,zspec[iobj],ztry[izbest],weighted_z(pzuse),pzuse[izbest],idtemp2a[izbest]+1,idtemp2b[izbest]+1);
            break;
          case 99:
            izbest=0;
            for (i=0;i<NZ;++i) {
                pzuse[i] = pzall[i];
                pzout[i] = pzuse[i];
                if (pzall[i]<pzuse[izbest]) izbest=i;
                pzall[i]=1.e30;
            }
            fprintf(fp,"%s %7.4lf %6.3lf %6.3lf %e ",
                strfmt,zspec[iobj],ztry[izbest],weighted_z(pzuse),pzuse[izbest]);
            printf("%s %7.4lf %6.3lf %6.3lf %e ",
                strfmt,zspec[iobj],ztry[izbest],weighted_z(pzuse),pzuse[izbest]);
            break;
            
          izprior = izbest;  //// For when APPLY_PRIOR=0
        }
        
        izprior=izbest;                               
        ////////  Apply prior, redshift of maximum likelihood at ztry[izprior]
        ////////  (could be the same as the minimum chisq *before* applying the prior: ztry[izbest])
        if (APPLY_PRIOR && FIX_ZSPEC==0) {

            if (fnu[iobj][PRIOR_FILTER] < NOT_OBS_THRESHOLD || efnu[iobj][PRIOR_FILTER] < NOT_OBS_THRESHOLD) {
                izprior = izbest;     ///// Do nothing and prior result will be the same as first result
                ///// but we still need to convert chi-sq ("pzuse") to likelihoods ("pzout")
                for (i=0;i<NZ;++i) 
                    pzout[i] = exp(-0.5*(pzuse[i]-pzuse[izbest])/CHI2_SCALE);                
            } else {
                apply_prior(priorkz,priorzval,klim,NZ,NK_prior,pzuse,fnu[iobj][PRIOR_FILTER],&izprior,pzout);
            }
            
            switch (TEMPLATE_COMBOS) {
              case 1:
                fprintf(fp,"%7.3f %e %5d %7.3f %7.3f",ztry[izprior],pzuse[izprior],idtemp1[izprior]+1,weighted_z2(pzout),odds(pzout));
                printf("%7.3f %e %5d %7.3f %7.3f",ztry[izprior],pzuse[izprior],idtemp1[izprior]+1,weighted_z2(pzout),odds(pzout));
                break;
              case 2:
                fprintf(fp,"%7.3f %e %5d %5d %7.3f %7.3f",ztry[izprior],pzuse[izprior],idtemp2a[izprior]+1,idtemp2b[izprior]+1,weighted_z2(pzout),odds(pzout));
                printf("%7.3f %e %5d %5d %7.3f %7.3f",ztry[izprior],pzuse[izprior],idtemp2a[izprior]+1,idtemp2b[izprior]+1,weighted_z2(pzout),odds(pzout));
                break;
              case 99:
                fprintf(fp,"%7.3f %e %7.3f %7.3f",ztry[izprior],pzuse[izprior],weighted_z2(pzout),odds(pzout));
                printf("%7.3f %e %7.3f %7.3f",ztry[izprior],pzuse[izprior],weighted_z2(pzout),odds(pzout));
                break;
            }
        } else {        	
        	for (i=0;i<NZ;++i) pzout[i] = exp(-0.5*(pzuse[i]-pzuse[izbest])/CHI2_SCALE);
        }
        
        //////// Print confidence intervals
        if (PRINT_ERRORS && FIX_ZSPEC==0) {        
            //get_errors(&l68, &u68, &l95, &u95, &l99, &u99);
            //l95=l68;
            //u95=u68;
            get_errors(&l68, &u68, &l95, &u95, &l99, &u99);
            fprintf(fp,"   %6.3lf %6.3lf  %6.3lf %6.3lf  %6.3lf %6.3lf ",
                   l68,u68,l95,u95,l99,u99);            
        }
        
        fprintf(fp," %4d ",ngoodfilters);
        
        printf("\n");
        fprintf(fp,"\n");
        
        if (POFZ_FILE && FIX_ZSPEC==0 && BINARY_OUTPUT==0) {
            sprintf(outfilename,"%s/%s.pz",OUTPUT_DIRECTORY,objid[iobj]);
            pz_output(outfilename,iobj);    
        }

		if (BINARY_OUTPUT) {  
		   //// POFZ
           if (APPLY_PRIOR) 
            	fwrite(pzout,sizeof(double)*NZ,1,pzf);
        	else
            	fwrite(pzuse,sizeof(double)*NZ,1,pzf);
		   //// Coeffs, OBS_SED
		   fwrite(coeffs[izprior],sizeof(double)*NTEMP,1,cf);

        }
        

//////  Uncomment these lines and the commented lines in the IF statements below to print template filter fluxes and the
//////  full best-fit templates at the 'marginalized' redshifts, rather than at the probability peaks.

//        zm1 = weighted_z(pzuse);
//        iz_zm1=0;
//        while (ztry[iz_zm1] < zm1 && iz_zm1 < (NZ-1)) ++iz_zm1;
//
//        zm2 = weighted_z2(pzout);
//        iz_zm2=0;
//        while (ztry[iz_zm2] < zm2 && iz_zm2 < (NZ-1)) ++iz_zm2;


		if (TEMP_SED_FILE && BINARY_OUTPUT==0) {
            sprintf(outfilename,"%s/%s.temp_sed",OUTPUT_DIRECTORY,objid[iobj]);
            temp_sed_output(outfilename,iobj,izbest,izprior);
////            temp_sed_output(outfilename,iobj,iz_zm1,iz_zm2);
        }
        
        izsave[iobj] = izprior;
        
	    if (OBS_SED_FILE && BINARY_OUTPUT==0) {
			sprintf(outfilename,"%s/%s.obs_sed",OUTPUT_DIRECTORY,objid[iobj]);
            obs_sed_output(outfilename,iobj,izbest,izprior);
////            obs_sed_output(outfilename,iobj,iz_zm1,iz_zm2);
   	    }

      }
      	  
    }
    fclose(fp);
  
  	if (BINARY_OUTPUT) {  //// Coeffs, OBS_SED
        fclose(pzf);
        izsave32 = malloc(sizeof(int32_t)*nobj);
        for (i=0;i<nobj;++i) {
            izsave32[i] = izsave[i];
        }
        fwrite(izsave32,sizeof(int32_t)*nobj,1,cf);
  		fwrite(tnorm,sizeof(double)*NTEMP,1,cf);
  		fclose(cf);
  	}
  	
    time(&now);
        
    if (VERBOSE_LOG) {
    	fprintf(fplog,"# EAZY v1.00 (July 9, 2008)\n");
        fprintf(fplog,"# Took %ld seconds.\n",now-then);
        fclose(fplog);
    }
    
    printf("All done.  Took %ld seconds.\n",now-then);
    
    return(0);
}

/////// Input is chi-squared
double weighted_z (double *pofz) {
    
    long i;
    double this,next,num,denom,dz,cmin;
    
    cmin = 1.e30;
    for (i=0; i<NZ; ++i) if (pofz[i] < cmin) cmin = pofz[i];

    num=0.;
    denom=0.;
    this = exp(-0.5*(pofz[0]-cmin)/CHI2_SCALE);
    for (i=1; i<NZ; ++i) {
        dz = ztry[i]-ztry[i-1];
        next = exp(-0.5*(pofz[i]-cmin)/CHI2_SCALE);
        
        num += (this*ztry[i-1]+next*ztry[i])*dz;
        denom += (this+next)*dz;
        
        this = next;
        // printf(" w(z) %lf %le %lf\n",dz,exp(-0.5*pofz[i]),ztry[i]);
    }
    
    return (num/denom);  //// factor of 2 from trapezoid rule cancels out
}

/////// Input is a probability, not chi-squared
double weighted_z2 (double *pofz) {
    
    long i;
    double this,next,num,denom,dz,cmin;
    
    cmin = 1.e30;
    for (i=0; i<NZ; ++i) if (pofz[i] < cmin) cmin = pofz[i];

    num=0.;
    denom=0.;
    this = pofz[0];
    for (i=1; i<NZ; ++i) {
        dz = ztry[i]-ztry[i-1];
        next = pofz[i];
        
        num += (this*ztry[i-1]+next*ztry[i])*dz;
        denom += (this+next)*dz;
        
        this = next;
    }
    
    return (num/denom);  //// factor of 2 from trapezoid rule cancels out
}

double odds (double *pofz) {

    long i;
    double this,next,num,denom,dz,cmin;
    double zbest,outodds,limit;
    
    limit = 0.2;
    
    cmin = 1.e30;
    for (i=0; i<NZ; ++i) if (pofz[i] < cmin) cmin = pofz[i];

    num=0.;
    denom=0.;
    this = pofz[0];
    for (i=1; i<NZ; ++i) {
        dz = ztry[i]-ztry[i-1];
        next = pofz[i];
        
        num += (this*ztry[i-1]+next*ztry[i])*dz;
        denom += (this+next)*dz;
        
        this = next;
    }
    
    zbest =num/denom;  //// factor of 2 from trapezoid rule cancels out

    outodds=0.;
    i=0;
    while (ztry[i] < (zbest-limit*(1+zbest)) && i < NZ) ++i;
    if (i == 0) ++i;
    this = pofz[i-1];
    while (ztry[i] < (zbest+limit*(1+zbest)) && i < NZ) {
        dz = ztry[i]-ztry[i-1];
        next = pofz[i];
        outodds += (this+next)*dz;        
        this = next;
        ++i;
    }
    if (i < NZ-1) {
        dz = ztry[i]-ztry[i-1];
        next = pofz[i];
        outodds += (this+next)*dz;
    }
            
    return (outodds/denom);

}

//////////  Get confidence intervals by (linear) interpolating the redshift grid
//////////  at particular cumulative probabilities
void get_errors(double *l68, double *u68, double *l95, double *u95, double *l99, double *u99)
{

    double pztotal;
    long i,istart;
        
    /////// Get total normalization
    pztotal = 0.;
    zinverse[0] = ztry[0]; /// not inverse here but the middle of the integration trapezoid
    for (i=1;i<NZ;++i) {
    	zinverse[i] = (ztry[i]+ztry[i-1])/2.;
    	pztotal += (ztry[i]-ztry[i-1])*(pzout[i]+pzout[i-1]);                   
    }
    
    ///// get lower limits, factors of 1/2 from the trapezoid rule cancel out in pzsum/pztotal
    pzsum[0] = 0.;
    for (i=1;i<NZ;++i) {
    	pzsum[i] = pzsum[i-1]+(ztry[i]-ztry[i-1])*(pzout[i]+pzout[i-1])/pztotal;
//    	printf("\nps %13.6f",pzsum[i]);
    }
    istart=0;
    interpol(0.00134990,l99,pzsum,zinverse,NZ,&istart);
    interpol(0.022750,l95,pzsum,zinverse,NZ,&istart);
    interpol(0.158655,l68,pzsum,zinverse,NZ,&istart);

    ///// get upper limits, invert ztry because interpol assumes monotonically *increasing* abscissa,
    ///// so you need to reverse 'pzsum'.
    zinverse[0] = ztry[NZ-1];
    for (i=1;i<NZ;++i) {
    	zinverse[i] = (ztry[NZ-i]+ztry[NZ-i-1])/2.;
    }
    pzsum[0] = 0.;
    for (i=1;i<NZ;++i) pzsum[i] = pzsum[i-1]+(ztry[NZ-i]-ztry[NZ-i-1])*(pzout[NZ-i]+pzout[NZ-i-1])/pztotal;

    istart=0;
    interpol(0.00134990,u99,pzsum,zinverse,NZ,&istart);
    interpol(0.022750,u95,pzsum,zinverse,NZ,&istart);
    interpol(0.158655,u68,pzsum,zinverse,NZ,&istart);
         
}   

	
