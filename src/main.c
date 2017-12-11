#include "defs.h"

char EAZY_VERSION[] = "2015-05-08";

double CHI2_SCALE;
long RandSeed=-1;

char zphot_param_file[64], zphot_translate_file[64], zphot_zeropoint_file[64];

// ** control parameters from zphot.params

int NTEMP;
int NTEMP_REST;  // to be defined from template file
char TEMPLATES_FILE[1024];
char RF_TEMPLATES_FILE[1024];
char ZBIN_FILE[1024];
int NTEMPL;

//redshift grid
int NZ,NZCOMOVE;
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
int MAGNITUDES;
double NOT_OBS_THRESHOLD;
int N_MIN_COLORS;
char TEMP_ERR_FILE[1024];
char WAVELENGTH_FILE[1024];
char LAF_FILE[1024];
char DLA_FILE[1024];

// output
char OUTPUT_DIRECTORY[1024];
char MAIN_OUTPUT_FILE[1024];
int OBS_SED_FILE;
int TEMP_SED_FILE;
int POFZ_FILE;
int VERBOSE_LOG;
int BINARY_OUTPUT;

int PRINT_ERRORS;

//// IGM parametrization from Inoue et al. 2014
double *lam1, *ALAF1, *ALAF2, *ALAF3, *ADLA1, *ADLA2;
int NA;

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
double SCALE_2175_BUMP;

char PRIOR_FILE[1024];
int APPLY_PRIOR;
double PRIOR_ABZP;
int PRIOR_FILTER, PRIOR_FILTER_IDX;
char REST_FILTERS[1024],Z_COLUMN[512];
int USE_ZSPEC_FOR_REST;
double RF_PADDING;
int RF_ERRORS;

int READ_ZBIN, ZBIN_OPENED;

// start of filter throughput data linked list
filt_data filt_thru;

long nfilter,*okfilt,nusefilt,nrestfilt,totfiltid,totcolumn;

double *totflux;

long nobj;
char **objid;
int idsize;

double **fnu, **efnu, *zspec, *ra, *dec;
double *lambdac,*lambdac_sort,*ztry, *comovingdist,*zcomove,*zspec;
double *templ,**tempf,**tempf_rest,*tempage;
double temp_err[NTEMPLMAX],**temp_errf, temp_scale[NTEMPLMAX];
double ***tempfilt, ***tempfilt_rest, *temp_err_a, *tnorm;

long *lc_sort_idx;

int **temp_combine;

filt_data **pusefilt;

int GET_ZP_OFFSETS;
double ZP_OFFSET_TOL;

int *zpfilterid, zpcontinue, *filt_to_col;
double *zpfactor;

FILE *fplog, *fprf;

double *dasum,*dbsum;

int iz,*idtemp1,*idtemp2a,*idtemp2b,**idtempall,ntemp_all,ngoodfilters,ncols;
double *pz1,*pz2,*pzall,*chi2fit,*pzout,*pzsum,*zinverse;
double *atemp1,*atemp2a,*atemp2b,**coeffs;

double **priorkz,*klim;
long Kcolumn;
int32_t *klim_idx;

double *z_a, *z_p, *z_m1, *z_m2, *z_peak;

int32_t nusefilt32, NTEMP32, nobj32, NZ32, NK32,*izsave32, NTEMPL32;

//// Get command line options for zphot.param, zphot.translate, and zphot.zeropoint files
//// http://www.gnu.org/software/hello/manual/libc/Example-of-Getopt.html
int processCommandLine(int argc, char **argv) {

    int c;
    
    opterr = 0;

    //// Default files
    strcpy(zphot_param_file,"zphot.param");
    strcpy(zphot_translate_file,"zphot.translate");
    strcpy(zphot_zeropoint_file,"zphot.zeropoint");
  
    while ((c = getopt (argc, argv, "z:p:t:")) != -1)
      switch (c) {
        case 'p':
          strcpy(zphot_param_file,optarg);
          break;
        case 't':
          strcpy(zphot_translate_file,optarg);
          break;
        case 'z':
          strcpy(zphot_zeropoint_file,optarg);
          break;
        case '?':
          if (optopt == 'z' || optopt == 'p' || optopt == 't')
            fprintf (stderr, "EAZY: Option -%c requires a filename.\n", optopt);
          else if (isprint (optopt)) {
            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            fprintf (stderr, "EAZY accepts the following command line options (default filename):\n");
            fprintf (stderr,"   -p (zphot.param)\n   -t (zphot.translate)\n   -z (zphot.zeropoint)\n");
          } else
            fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
          exit(1);
        default:
          abort ();
    }    
    return 0;
}

int main(int argc, char **argv) {

	char coeff_file[1024],pz_file[1024],zbin_file[1024];
	FILE *cf, *pzf;
		
    long i,j,iobj,izbest,iz_zm1,iz_zm2, *izsave,ipeakmax;
    double zm1,zm2,zpeak_best,zpeak_prob,pztot;
    FILE *fp, *fpzbin; //, *fpasciicoeff;
    char ofile[1024],outfilename[1024];
    char strfmt[64];
    
    int getzStatus,oneFiltGTzero;

    double weighted_z (double *pofz);
    double weighted_z2 (double *pofz);
    double odds (double *pofz);
    void best_peak (double *pofz, double *zpeak_best, double *zpeak_prob, long *ipeakmax);
    // void find_zpvd (double *pofz, long ipeakmax, double *zpvd);
    void rest_frame_colors();
    float ran1(long *idum);
    
    // double calculate_mass (long zidx, double **coeffs);

    double l68,u68,l95,u95,l99,u99;
    
    double **priorkz_raw, *priorzval,q_z,ztest;
    long NZ_prior,NK_prior,izprior,izpeak;
    
    double zmc;
    long izmc;

    time_t then,now,middle;

    time(&then);
    
    processCommandLine(argc,argv);
    printf("EAZY, version %s\n\n", EAZY_VERSION);
    printf("Control files:\n    param = %s\ntranslate = %s\nzeropoint = %s\n\n",zphot_param_file,zphot_translate_file,zphot_zeropoint_file);
    // exit(1);
                
    printf("Reading in program parameters...\n\n");
    getparams();

    printf("Initializing...");

    init();
             
    // printf("  Phew!! Ready to go.\n");

    ////// Some options don't apply IF redshift fixed to z_spec --> turn them off
    if (FIX_ZSPEC) {
        if (POFZ_FILE) fprintf(stderr,"FIX_ZSPEC set, ignoring POFZ_FILE ...\n");
        if (APPLY_PRIOR) fprintf(stderr,"FIX_ZSPEC set, ignoring APPLY_PRIOR ...\n");
        if (PRINT_ERRORS) fprintf(stderr,"FIX_ZSPEC set, ignoring PRINT_ERRORS ...\n");  
        APPLY_PRIOR=0;
        PRINT_ERRORS=0;
    }

    /////// Make sure BINARY_OUTPUT is set if computing RF color errors
    if (RF_ERRORS && isdigit(REST_FILTERS[0]) && BINARY_OUTPUT==0 && READ_ZBIN==0)  {
        fprintf(stderr,"RF_ERRORS set, forcing BINARY_OUTPUT=1 ...\n");
        BINARY_OUTPUT=1;
    }
    
    /////// Initialize Prior
    if (APPLY_PRIOR && FIX_ZSPEC==0) {
        
        Kcolumn=0;
        //// Get size of grid
        get_prior_file_size(&NZ_prior, &NK_prior);
        printf("Read PRIOR_FILE NZ: %ld NK: %ld \n",NZ_prior,NK_prior);
        
        klim = malloc(sizeof(double)*NK_prior);
        klim_idx = malloc(sizeof(int32_t)*nobj);
        
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

    okfilt = malloc(sizeof(long)*nusefilt);
    if(okfilt == NULL) {
        fprintf(stderr, "okfilt failed, out of memory\n");
        exit(1);
    }

    chi2fit = malloc(sizeof(double)*NZ);
    if (chi2fit==NULL) {
        fprintf(stderr,"chi2fit: memory allocation error!\n");
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

    //////// Apply zero point corrections
    if (GET_ZP_OFFSETS) zeropoints();

  if (!READ_ZBIN) { 

    ////// Initialise output file
    strcpy(ofile,OUTPUT_DIRECTORY);
    strcat(ofile,"/");
    strcat(ofile,MAIN_OUTPUT_FILE);
    strcat(ofile,".zout");
    if (!(fp = fopen(ofile,"w"))) {
        fprintf(stderr,"Oops...can't open the output file, %s\n",ofile);
        exit(1);
    }

    // strcpy(ofile,OUTPUT_DIRECTORY);
    // strcat(ofile,"/");
    // strcat(ofile,MAIN_OUTPUT_FILE);
    // strcat(ofile,".coeff.txt");
    // if (BINARY_OUTPUT==0) if (!(fpasciicoeff = fopen(ofile,"w"))) {
    //     fprintf(stderr,"Oops...can't open the output file, %s\n",ofile);
    //     exit(1);
    // }
    
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
    if (PRINT_ERRORS && FIX_ZSPEC==0) fprintf(fp," l68 u68  l95 u95  l99 u99 ");
    
    
    fprintf(fp," nfilt");
    ncols+=1;
    if (FIX_ZSPEC==0) {
        fprintf(fp," q_z ");
        ncols+=1;
        fprintf(fp,"z_peak peak_prob z_mc");
        ncols+=3;
    }
    
    printf("\n");
    fprintf(fp,"\n");

    ////// VERSION INFORMATION
    fprintf(fp,"# EAZY %s\n",EAZY_VERSION);
            
    //// Dump coefficients to binary file
    izsave = malloc(sizeof(long)*nobj);
    if (BINARY_OUTPUT && !READ_ZBIN) {
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
       	sprintf(pz_file,"%s/%s.pz",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE);
    	pzf = fopen(pz_file,"w");
        fwrite(&NZ32,sizeof(int32_t),1,pzf);
        fwrite(&nobj32,sizeof(int32_t),1,pzf);

        //// zbin file
        sprintf(zbin_file,"%s/%s.zbin",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE);
    	fpzbin = fopen(zbin_file,"w");
        fwrite(&nobj32,sizeof(int32_t),1,fpzbin);
        
 	}    	
            
    ////////
    //// Main loop to get photz for each object in the catalog
    ////////    
    for (iobj=0; iobj<nobj; ++iobj) { 
        
      izsave[iobj] = -1;
      ngoodfilters = 0;
      oneFiltGTzero=0;
      for(j=0;j<nusefilt;++j) {
          if (fnu[iobj][j] > NOT_OBS_THRESHOLD && efnu[iobj][j] > 0) {
              if ((fnu[iobj][j]/efnu[iobj][j]) > 1) ++oneFiltGTzero;              
              ++ngoodfilters;
              okfilt[j] = 1;
              //printf("%lf %lf %lf\n", fnu[iobj][j],efnu[iobj][j],fnu[iobj][j]/efnu[iobj][j]);
          } else okfilt[j] = 0;
      }
      
      //////// If too few filters, skip  
      if (oneFiltGTzero <= 1 || ngoodfilters < N_MIN_COLORS) {
          if (oneFiltGTzero <= 1) 
              printf("%s: Only one or no filters with flux/err > 1, skipping...\n",objid[iobj]);
          else 
              printf("%s: only %d filters available, skipping...\n",objid[iobj],ngoodfilters);
          fprintf(fp,"%s %lf",objid[iobj],zspec[iobj]);
          for (j=0;j<ncols;++j) fprintf(fp," -99");
          fprintf(fp,"\n");
          
          // if (BINARY_OUTPUT==0) {
          //     for (j=0;j<NTEMP;++j) fprintf(fpasciicoeff," -99");
          //     fprintf(fpasciicoeff,"\n");
          // }
          
          for (j=0;j<NZ;++j) pzout[j] = 0;
          for (j=0;j<NTEMP;++j) coeffs[0][j] = 0;
          if (BINARY_OUTPUT && !READ_ZBIN) {
              fwrite(pzout,sizeof(double),NZ,pzf);
              fwrite(coeffs[0],sizeof(double),NTEMP,cf);
          }
          izsave[iobj] = -1;

      } else { 

        /////////////////  Call out to getphotz to get Chisq as a function of z //////////////////
        ////////// Unused arrays are not malloc-ed, just those appropriate for TEMPLATE_COMBOS
        getzStatus = getphotz(iobj, pz1, idtemp1, atemp1, 
                       pz2, idtemp2a, idtemp2b, atemp2a, atemp2b,
                       pzall, idtempall, coeffs, &ntemp_all, -1);
                
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
                chi2fit[i] = pz1[i];
                pzout[i] = chi2fit[i];
                if (pz1[i]<chi2fit[izbest]) izbest=i;
                pz1[i]=1.e30;
                z_a[iobj] = ztry[izbest];
            }
            fprintf(fp,"%s %7.4lf %6.3lf %6.3lf %e %5d ",
                strfmt,zspec[iobj],ztry[izbest],weighted_z(chi2fit),chi2fit[izbest],idtemp1[izbest]+1);
            printf("%s %7.4lf %6.3lf %6.3lf %e %5d ",
                strfmt,zspec[iobj],ztry[izbest],weighted_z(chi2fit),chi2fit[izbest],idtemp1[izbest]+1);
            break;
          case 2:
            izbest=0;
            for (i=0;i<NZ;++i) {
                chi2fit[i] = pz2[i];
                pzout[i] = chi2fit[i];
                if (pz2[i]<chi2fit[izbest]) izbest=i;
                pz2[i]=1.e30;
                z_a[iobj] = ztry[izbest];
             }            
            fprintf(fp,"%s %7.4lf %6.3lf %6.3lf %e %5d %5d ",
                strfmt,zspec[iobj],ztry[izbest],weighted_z(chi2fit),chi2fit[izbest],idtemp2a[izbest]+1,idtemp2b[izbest]+1);
            printf("%s %7.4lf %6.3lf %6.3lf %e %5d %5d ",
                strfmt,zspec[iobj],ztry[izbest],weighted_z(chi2fit),chi2fit[izbest],idtemp2a[izbest]+1,idtemp2b[izbest]+1);
            break;
          case 99:
            izbest=0;
            for (i=0;i<NZ;++i) {
                chi2fit[i] = pzall[i];
                pzout[i] = chi2fit[i];
                if (pzall[i]<chi2fit[izbest]) izbest=i;
                pzall[i]=1.e30;
                z_a[iobj] = ztry[izbest];
            }
            z_m1[iobj] = weighted_z(chi2fit);
            
            fprintf(fp,"%s %7.4lf %6.3lf %6.3lf %e ",
                strfmt,zspec[iobj],ztry[izbest],z_m1[iobj],chi2fit[izbest]);
            printf("%s %7.4lf %6.3lf %6.3lf %e ",
                strfmt,zspec[iobj],ztry[izbest],z_m1[iobj],chi2fit[izbest]);
            break;
            
            z_p[iobj] = z_p[iobj];
            z_m2[iobj] = z_m1[iobj];
            
            izprior = izbest;  //// For when APPLY_PRIOR=0
        }
        
        izprior=izbest;                               
        ////////  Apply prior, redshift of maximum likelihood at ztry[izprior]
        ////////  (could be the same as the minimum chisq *before* applying the prior: ztry[izbest])
        if (APPLY_PRIOR && FIX_ZSPEC==0) {

            if (fnu[iobj][PRIOR_FILTER_IDX] < NOT_OBS_THRESHOLD || efnu[iobj][PRIOR_FILTER_IDX] < NOT_OBS_THRESHOLD) {
                izprior = izbest;     ///// Do nothing and prior result will be the same as first result
                ///// but we still need to convert chi-sq ("chi2fit") to likelihoods ("pzout")
                pzout[0] = exp(-0.5*(chi2fit[0]-chi2fit[izbest])/CHI2_SCALE); //*(1+ztry[0]);   
                pztot = 0;
                for (i=1;i<NZ;++i) {
                    pzout[i] = exp(-0.5*(chi2fit[i]-chi2fit[izbest])/CHI2_SCALE); //*(1+ztry[i]);
                    pztot+=(ztry[i]-ztry[i-1])*(pzout[i]+pzout[i-1]);  
                }
                for (i=1;i<NZ;++i) pzout[i] /= pztot/2.;
                klim_idx[iobj] = -1;        
                //fprintf(stderr,"%ld Bad K\n",iobj);
                
            } else {
                //printf("\n\nAPPLY PRIOR!\n\n");
                apply_prior(priorkz,NZ,NK_prior,izbest,chi2fit,fnu[iobj][PRIOR_FILTER_IDX],&izprior,pzout);
                klim_idx[iobj] = (int32_t) Kcolumn;
                pztot=0.; 
                for (i=0;i<NZ;++i) pztot+=pzout[i];
                if ((pztot == 0) | isnan(pztot)) {
                    printf("Sum p(z) = 0, check this object and the prior: ");
                    pzout[0] = 1.;
                    izprior = 0;
                    chi2fit[0] = 1.e8;
                    for (i=1;i<NZ;++i) pzout[i] = 0;
                }
            }
            //printf("\nPZOUT: %lf\n", pztot);
            if ((pztot == 0) | isnan(pztot)) {
                printf("Sum p(z) = 0, check this object and the prior: ");
                pzout[0] = 1.;
                izprior = 0;
                chi2fit[0] = 1.e8;
                for (i=1;i<NZ;++i) pzout[i] = 0;
            }
            
            z_p[iobj] = ztry[izprior];
            z_m2[iobj] = weighted_z2(pzout);
            switch (TEMPLATE_COMBOS) {
              case 1:
                fprintf(fp,"%7.3f %e %5d %7.3f %7.3f",ztry[izprior],chi2fit[izprior],idtemp1[izprior]+1,z_m2[iobj],odds(pzout));
                printf("%7.3f %e %5d %7.3f %7.3f",ztry[izprior],chi2fit[izprior],idtemp1[izprior]+1,z_m2[iobj],odds(pzout));
                break;
              case 2:
                fprintf(fp,"%7.3f %e %5d %5d %7.3f %7.3f",ztry[izprior],chi2fit[izprior],idtemp2a[izprior]+1,idtemp2b[izprior]+1,z_m2[iobj],odds(pzout));
                printf("%7.3f %e %5d %5d %7.3f %7.3f",ztry[izprior],chi2fit[izprior],idtemp2a[izprior]+1,idtemp2b[izprior]+1,z_m2[iobj],odds(pzout));
                break;
              case 99:
                fprintf(fp,"%7.3f %e %7.3f %7.3f",ztry[izprior],chi2fit[izprior],z_m2[iobj],odds(pzout));
                printf("%7.3f %e %7.3f %7.3f",ztry[izprior],chi2fit[izprior],z_m2[iobj],odds(pzout));
                break;
            }
        } else {        	
        	//for (i=0;i<NZ;++i) pzout[i] = exp(-0.5*(chi2fit[i]-chi2fit[izbest])/CHI2_SCALE);
            pzout[0] = exp(-0.5*(chi2fit[0]-chi2fit[izbest])/CHI2_SCALE); //*(1+ztry[0]);   
            pztot = 0;
            for (i=1;i<NZ;++i) {
                pzout[i] = exp(-0.5*(chi2fit[i]-chi2fit[izbest])/CHI2_SCALE); //*(1+ztry[i]);
                pztot+=(ztry[i]-ztry[i-1])*(pzout[i]+pzout[i-1]);  
            }
            for (i=1;i<NZ;++i) pzout[i] /=pztot/2.;
        }
        
        //////// Print confidence intervals
        get_errors(&l68, &u68, &l95, &u95, &l99, &u99);
        if (PRINT_ERRORS && FIX_ZSPEC==0) {        
            //get_errors(&l68, &u68, &l95, &u95, &l99, &u99);
            //l95=l68;
            //u95=u68;
            fprintf(fp,"   %6.3lf %6.3lf  %6.3lf %6.3lf  %6.3lf %6.3lf ",
                   l68,u68,l95,u95,l99,u99);            
        }
        
        fprintf(fp," %4d ",ngoodfilters);
        
        zpeak_best = 0.;
        izpeak=izprior;
        if (FIX_ZSPEC==0) {
            
            //// Find peak with largest integrated probability
            best_peak (pzout, &zpeak_best, &zpeak_prob, &ipeakmax);
            ztest=100;
            if (zpeak_best > 0) {
                for (i=0;i<NZ;++i) if (fabs(ztry[i]-zpeak_best) < ztest) {
                    ztest = fabs(ztry[i]-zpeak_best);
                    izpeak=i;
                }
            }       
            
            izprior = izpeak;
            if (izpeak <= 0) {
                izpeak = 0;
                zpeak_best = ztry[0];
            }
            if (izpeak >= (NZ-1)){
                izpeak = ztry[NZ-1];
                zpeak_best = ztry[NZ-1];
            }
            
            //// testing
            //fprintf(stderr,"%ld Kcolumn main: %ld\n",iobj,klim_idx[iobj]);
            //for (i=0;i<NZ;++i) printf("\n zzp %ld %ld %8.3f %13.5e %13.5e %13.5e",iobj,klim_idx[iobj],ztry[i],chi2fit[i],pzout[i],priorkz[Kcolumn][i]);
            
            //// Q_z quality parameter
            q_z = chi2fit[izprior]/(ngoodfilters-3)*(u99-l99)/odds(pzout);
            fprintf(fp," %13.5e",q_z);
            
            z_peak[iobj] = zpeak_best;
            fprintf(fp," %7.4f %7.3f",zpeak_best,zpeak_prob);
            
            //// z_mc from Wittman (astro-ph:0905.0892).  Use cumulative probability, "pzsum", from get_errors()
            izmc=0;
            interpol(ran1(&RandSeed),&zmc,pzsum,ztry,NZ,&izmc);
            fprintf(fp," %7.4f",zmc);
            
            // find_zpvd (pzout, ipeakmax, &zpvd);
            // fprintf(fp," %7.3f",zpvd);

        }
        
        printf("\n");
        fprintf(fp,"\n");
        
        if (POFZ_FILE && FIX_ZSPEC==0 && BINARY_OUTPUT==0 && !READ_ZBIN) {
            sprintf(outfilename,"%s/%s_%s.pz", OUTPUT_DIRECTORY, MAIN_OUTPUT_FILE, objid[iobj]);
            pz_output(outfilename,iobj);    
        }

		if (BINARY_OUTPUT && !READ_ZBIN) {  
		   //// POFZ
            //            if (APPLY_PRIOR) 
            //              fwrite(pzout,sizeof(double)*NZ,1,pzf);
            // else
            //              fwrite(chi2fit,sizeof(double)*NZ,1,pzf);
           //// Always print pz without the prior, and also print prior(K,z) to apply yourself
           //for (i=0;i<NZ;++i) chi2fit[i] = exp(-0.5*chi2fit[i]);
           fwrite(chi2fit,sizeof(double)*NZ,1,pzf);
		   //// Coeffs, OBS_SED
		   fwrite(coeffs[izprior],sizeof(double)*NTEMP,1,cf);
        }
        
        // if (BINARY_OUTPUT==0) {
        //     for (i=0;i<NTEMP;++i) fprintf(fpasciicoeff," %11.4le",coeffs[izprior][i]);
        //     fprintf(fpasciicoeff,"\n");
        // }
        
//////  Uncomment these lines and the commented lines in the IF statements below to print template filter fluxes and the
//////  full best-fit templates at the 'marginalized' redshifts, rather than at the probability peaks.
//        zm1 = weighted_z(chi2fit);
//        iz_zm1=0;
//        while (ztry[iz_zm1] < zm1 && iz_zm1 < (NZ-1)) ++iz_zm1;
//
//        zm2 = weighted_z2(pzout);
//        iz_zm2=0;
//        while (ztry[iz_zm2] < zm2 && iz_zm2 < (NZ-1)) ++iz_zm2;


		if (TEMP_SED_FILE && BINARY_OUTPUT==0 && !READ_ZBIN) {
            sprintf(outfilename, "%s/%s_%s.temp_sed", OUTPUT_DIRECTORY, MAIN_OUTPUT_FILE, objid[iobj]);
            temp_sed_output(outfilename,iobj,izbest,izprior);
////            temp_sed_output(outfilename,iobj,iz_zm1,iz_zm2);
        }
        
        izsave[iobj] = izprior;
        
	    if (OBS_SED_FILE && BINARY_OUTPUT==0 && !READ_ZBIN) {
			sprintf(outfilename, "%s/%s_%s.obs_sed", OUTPUT_DIRECTORY, MAIN_OUTPUT_FILE, objid[iobj]);
            obs_sed_output(outfilename,iobj,izbest,izprior);
////            obs_sed_output(outfilename,iobj,iz_zm1,iz_zm2);
   	    }

      }
      	  
    }
    fclose(fp);
    //if (BINARY_OUTPUT==0) fclose(fpasciicoeff);
    
  	if (BINARY_OUTPUT && !READ_ZBIN) {  //// Coeffs, OBS_SED
        
        if (FIX_ZSPEC==0 && APPLY_PRIOR==1) {
            NK32 = (int32_t) NK_prior;
            fwrite(&NK32,sizeof(int32_t),1,pzf);
            fwrite(klim,sizeof(double)*NK_prior,1,pzf);
            for (i=0;i<NK_prior;++i) fwrite(priorkz[i],sizeof(double)*NZ,1,pzf);
            fwrite(klim_idx,sizeof(int32_t)*nobj,1,pzf);
            // for (i=0;i<nobj;++i) {
            //     NK32 = (int32_t) klim_idx[i];
            //     fwrite(&NK32,sizeof(int32_t),1,pzf);
            // }            
        }
        fclose(pzf);
        
        izsave32 = malloc(sizeof(int32_t)*nobj);
        for (i=0;i<nobj;++i) {
            izsave32[i] = izsave[i];
        }
        fwrite(izsave32,sizeof(int32_t)*nobj,1,cf);
  		fwrite(tnorm,sizeof(double)*NTEMP,1,cf);
  		fclose(cf);
        fwrite(z_a,sizeof(double)*nobj,1,fpzbin);
        fwrite(z_p,sizeof(double)*nobj,1,fpzbin);
        fwrite(z_m1,sizeof(double)*nobj,1,fpzbin);
        fwrite(z_m2,sizeof(double)*nobj,1,fpzbin);
        fwrite(z_peak,sizeof(double)*nobj,1,fpzbin);
        fclose(fpzbin);      
  	}
  } else {
      
      printf("Using pre-computed redshifts from %s for rest-frame colors.\n",ZBIN_FILE);
           
  }
    
    time(&middle);
    
    if (VERBOSE_LOG && !READ_ZBIN) {
  	    fprintf(fplog,"# EAZY %s\n",EAZY_VERSION);
        fprintf(fplog,"# Took %ld seconds.\n",middle-then);
        fclose(fplog);
        printf("Done fitting redshifts.  Took %ld seconds.\n",middle-then);
    }
        
    if (nrestfilt > 0) {
        fprintf(stderr,"\nFitting rest-frame color ");
        rest_frame_colors();
        fprintf(stderr," Done\n\n");
        time(&now);     
        printf("EAZY %s\nAll done.  Took %ld seconds.\n",EAZY_VERSION,now-then);
    }
    
    return(0);
}

/////// Input is chi-squared, convert in place
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

//// Find multiple peaks, distinguished as being above some threshold value, 'thresh'
void best_peak (double *pofz, double *zpeak_best, double *zpeak_prob, long *ipeakmax)
{
    long step,i;
    double zpeak_i, prob_i; //, zpeak_best, zpeak_prob;
    double thresh=1.e-2,pmax;
    long ipeak;
    
    double z_best, z_prob;
    long imax;
    
    //double psum, psum_i;
    
    // psum=0.;
    // for (i=1;i<NZ;++i) {
    //     psum_i = (ztry[i]-ztry[i-1])*(pofz[i]+pofz[i-1]); 
    //     psum+=psum_i; 
    //     printf("psum: %.4lf %.4lf\n",ztry[i], psum/2.);
    // }
    
    pmax=0.;
    for (i=0;i<NZ;++i) if (pofz[i] > pmax) pmax = pofz[i];
    if (thresh > pmax) thresh = 0.1 * pmax;
    
    z_prob = -1;
    z_best = -1;
    step=0;
    ipeak=0;
    while (step < NZ) {
        //// Start below threshold.  When threshold is reached, start integrating until the threshold is reached again.
        while (pofz[step] < thresh && step < NZ-1) ++step;
        prob_i = 0;
        zpeak_i = 0;
        pmax=0.;
        if (step==0) ++step;
        while (pofz[step] > thresh && step < NZ-1) {
            prob_i += (ztry[step]-ztry[step-1])*(pofz[step]+pofz[step-1]);
            zpeak_i += (ztry[step]-ztry[step-1])*(pofz[step]*ztry[step]+pofz[step-1]*ztry[step-1]);
            if (pofz[step] > pmax) {
                ipeak = step;
                pmax = pofz[ipeak];
            }
            ++step;
        }
        prob_i /= 2.; //// Trapezoid rule
        //printf("\n\n z_peak: %ld %.4lf %.4lf    %.4lf %.4lf  %.3lf\n",step, z_best, z_prob, zpeak_i/2./prob_i, prob_i, psum/2.);
        if (prob_i > z_prob) {
            z_best = zpeak_i / 2. / prob_i; 
            z_prob = prob_i*1.;
            imax = ipeak*1;
        }
        if (step == NZ-1) ++step;
    }   
    
    *zpeak_best = z_best;
    *zpeak_prob = z_prob;
    *ipeakmax = imax;
    
    // step=0;test=0;
    // for (step=1;step<NZ;++step) test+=(ztry[step]-ztry[step-1])*(pofz[step]+pofz[step-1]);
    // fprintf(stderr,"TEST: %lf\n",test);
           
}

//// Calculate a "z_m2" like redshift around the peak found in z_peak, but only where p > thresh*pmax
void find_zpvd (double *pofz, long ipeakmax, double *zpvd) {

    double thresh=0.8;
    double pmax, num, denom,pzi;
    long i, ipstart, ipstop;
    
    pmax = pofz[ipeakmax];
    thresh *= pmax;
    
    thresh = 0.02;
    
    i = ipeakmax; while (i > 1 && pofz[i] > thresh) --i; ipstart = i;
    i = ipeakmax; while (i < (NZ-1) && pofz[i] > thresh) ++i; ipstop = i;
    
    if (ipstart == ipstop) {
        *zpvd = ztry[ipeakmax];
    } else {
        num = 0.; denom=0;
        i=ipstart+1;
        num = (ztry[i]-ztry[i-1])*((pofz[i]-thresh)*ztry[i]+(0.)*ztry[i-1]);
        denom = (ztry[i]-ztry[i-1])*(pofz[i]+0-thresh);
        
        for (i=ipstart+2; i<ipstop;++i) {
            pzi = ((pofz[i]-thresh) > 0) ? (pofz[i]-thresh) : 0;
            num += (ztry[i]-ztry[i-1])*((pzi-thresh)*ztry[i]+(pofz[i-1]-thresh)*ztry[i-1]);
            denom += (ztry[i]-ztry[i-1])*(pzi+pofz[i-1]-2*thresh);
        }
        num += (ztry[i]-ztry[i-1])*(0.*ztry[i]+(pofz[i-1]-thresh)*ztry[i-1]);
        denom += (ztry[i]-ztry[i-1])*(0+pofz[i-1]-thresh);
        
        *zpvd = num / denom;
        //printf("zpvd: %ld %ld %lf %7.3f\n",ipstart,ipstop,pmax,*zpvd);
        
    }
}

double odds (double *pofz) 
{

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
    
    zbest = num/denom;  //// factor of 2 from trapezoid rule cancels out

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
    for (i=1;i<NZ;++i) pztotal += (ztry[i]-ztry[i-1])*(pzout[i]+pzout[i-1]);                   
	
    ///// get upper limits, invert ztry because interpol assumes monotonically *increasing* abscissa,
    ///// so you need to reverse 'pzsum'.
    zinverse[0] = ztry[NZ-1];
    for (i=1;i<NZ;++i) zinverse[i] = (ztry[NZ-i]+ztry[NZ-i-1])/2.;
    
    pzsum[0] = 0.;
    for (i=1;i<NZ;++i) pzsum[i] = pzsum[i-1]+(ztry[NZ-i]-ztry[NZ-i-1])*(pzout[NZ-i]+pzout[NZ-i-1])/pztotal;

    istart=0;
    interpol(0.00134990,u99,pzsum,zinverse,NZ,&istart);
    interpol(0.022750,u95,pzsum,zinverse,NZ,&istart);
    interpol(0.158655,u68,pzsum,zinverse,NZ,&istart);
    
    ///// get lower limits, factors of 1/2 from the trapezoid rule cancel out in pzsum/pztotal
    zinverse[0] = ztry[0]; /// not inverse here but the middle of the integration trapezoid
    for (i=1;i<NZ;++i) zinverse[i] = (ztry[i]+ztry[i-1])/2.;

    pzsum[0] = 0.;
    for (i=1;i<NZ;++i) pzsum[i] = pzsum[i-1]+(ztry[i]-ztry[i-1])*(pzout[i]+pzout[i-1])/pztotal;

    istart=0;
    interpol(0.00134990,l99,pzsum,zinverse,NZ,&istart);
    interpol(0.022750,l95,pzsum,zinverse,NZ,&istart);
    interpol(0.158655,l68,pzsum,zinverse,NZ,&istart);
         
}   

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
        int j;
        long k;
        static long iy=0;
        static long iv[NTAB];
        float temp;

        if (*idum <= 0 || !iy) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ;
                        *idum=IA*(*idum-k*IQ)-IR*k;
                        if (*idum < 0) *idum += IM;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ;
        *idum=IA*(*idum-k*IQ)-IR*k;
        if (*idum < 0) *idum += IM;
        j=iy/NDIV;
        iy=iv[j];
        iv[j] = *idum;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

	
