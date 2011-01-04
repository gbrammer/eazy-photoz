#include "defs.h"
#define BUFFSIZE 2048

/*************************************
*
* Read in catalog, filters, templates
*
**************************************/

char **col_data;
int ncol,  id_col, spec_col, *filt_defnum, *err_to_col; // *filt_to_col

//////// Read in zeropoint file
int read_zeropoint_file() {

    char buff[256],filt_col[256], factor[256],zptype[8];
    int i,MAG;    
    FILE *fp;
    double factor_i;
    
    //// No translation file
    if (!(fp = fopen("zphot.zeropoint","r"))) {
        zpcontinue=1;
        for (i=0;i<ncol;++i) zpfactor[i] = 1.;
        return(0);
    }
    
    zpcontinue = 0;
    
    while(fgets(buff,256,fp)) {     
        sscanf(buff,"%s %s",filt_col,factor);
        if (!strcmp(filt_col,"CONTINUE")) zpcontinue=1; 
        if (filt_col[0] !='\0' && factor[0]!='\0' && filt_col[0] != '#') {
            
            factor_i = atof(factor);
            
            if (filt_col[0] == 'M') {
                if (!isdigit(filt_col[1])) {
                    fprintf(stderr,"Bad line in zphot.zeropoint: %s.\n",buff);
                    exit(1);
                }
                MAG=1;
                filt_col[0]='F';
                strcpy(zptype,"mag");
                factor_i = pow(10,-0.4*factor_i);
            } else strcpy(zptype,"flux");
                        
            for (i=0;i<ncol;++i) {
                if (!strcmp(filt_col,col_data[i])) {
                    if (factor_i < 0) 
                        printf("zphot.zeropoint: (flux) scaling %lf < 0 for column %s, ignoring.\n",factor_i,filt_col);
                    else {
                        printf("zphot.zeropoint: will apply a %s correction of %s to column %s\n",
                                zptype,factor,filt_col);
                        zpfactor[i] = factor_i;
                    }
                }
            }
        }
    }
    
    fclose(fp);
    return(0);    
}

//////// Read in translation file, zphot.translate, and swap columns as defined
int translate_columns(char **col_data, int ncol) {

    char buff[256],from[256], to[256];
    int i;    
    FILE *fp;
    
    //// No translation file
    if (!(fp = fopen("zphot.translate","r"))) {
        return(0);
    }
    
    while(fgets(buff,256,fp)) {     
        sscanf(buff,"%s %s",from,to);
        if (from[0] !='\0' && to[0]!='\0') for (i=0;i<ncol;++i) {
            if (!strcmp(from,col_data[i])) {
                printf("zphot.translate: translating column name... %s -> %s\n",col_data[i],to);
                strcpy(col_data[i],to);
            }
        }
    }
    fclose(fp);
    return(0);    
}
 
void get_column_defs(char *buff) {
    char temp[BUFFSIZE];
    char *arg,*cp;
    char  delim[]=" \n\t,\r";
    int i,j,k, nerror,ifilt;
    double avlambda, totfilt,dl;
    double *filt_smooth,g1,g2;
     filt_data *pfilt_thru;

    strcpy(temp,buff);
    arg=strtok(temp,delim);
    if (!arg) {
        fprintf(stderr,"First line of catalog file is blank? Quitting...\n");
        exit(1);
    }
    if (*arg != '#') {
        fprintf(stderr,"Invalid format: first line of catalog file doesn't start with a '#'.\n");
        exit(1);
    }

    ncol=0;
    spec_col=-1;
    id_col=-1;
    
    ////// Handle case where there's no space between the # and the first 
    ////// column header or where there are multiple ##### and no space 
    while (arg[0]=='#') ++arg;
    if (arg[0]!='\0') ++ncol;
  
    while  ( (arg=strtok(NULL,delim)) ) ++ncol;
    printf("Found %d fields in catalog data file header...\n",ncol);
  
    filt_smooth = malloc(sizeof(double)*NTEMPLMAX);
    if(filt_smooth == NULL) {
        fprintf(stderr, "filt_smooth: out of memory!\n");
        exit(1);
    }

    zpfactor = malloc(sizeof(double)*ncol);
    if(zpfactor == NULL) {
        fprintf(stderr, "zpfactor: out of memory!\n");
        exit(1);
    }
    for (i=0;i<ncol;++i) zpfactor[i] = -1;

    col_data = malloc(sizeof(char *)*ncol);
    if(col_data == NULL) {
        fprintf(stderr, "col_data i: out of memory!\n");
        exit(1);
    }
    for (i=0;i<ncol;++i) {
        col_data[i] = malloc(sizeof(char)*BUFFSIZE);
        if (col_data[i]==NULL){
            fprintf(stderr, "col_data j: out of memory!\n");
            exit(1);
        }
    }

    arg=strtok(buff,delim); // throw away leading #
    while (arg[0]=='#') ++arg;
    if (arg[0]!='\0') {
        j=1;
        strcpy(col_data[0],arg);
    } else j=0;

    for (i=j;i<ncol;++i) {
        arg=strtok(NULL,delim);
        strcpy(col_data[i],arg);
    }
    
    translate_columns(col_data,ncol);
    if (GET_ZP_OFFSETS) read_zeropoint_file();
    
    nusefilt = 0;
    nerror = 0;
    for (i=0;i<ncol;++i) {
        
        ///// test that F and E columns have numbers following F or E
        if (col_data[i][0]=='F' && isdigit(col_data[i][1])) ++nusefilt;
        if (col_data[i][0]=='E' && isdigit(col_data[i][1])) ++nerror;
        
        if (!(strcmp("z_spec",col_data[i]))) {
            if (spec_col > -1) {
                fprintf(stderr,"Multiple z_spec columns found!\n");
                fprintf(stderr,"Check columns %d and %d.\n",spec_col,i);
                exit(1);
            } else spec_col = i;
        }
    }
  
    if (FIX_ZSPEC==1 && spec_col==-1) {
        fprintf(stderr,"\nError: FIX_ZSPEC is 'y' and no z_spec column found.\n\n");
        exit(1);
    }
  
    printf("%ld filter data columns found.\n",nusefilt);

    if (nerror != nusefilt) {
        fprintf(stderr, "Warning: the number of data columns [%ld] is different\n",
                         nusefilt);
        fprintf(stderr, " from the number of error columns [%d]. Please fix!!\n",
                        nerror);
        exit(1);
    }


    filt_to_col = malloc(sizeof(int )*nusefilt);
    if(filt_to_col == NULL) {
        fprintf(stderr, "filt_to_col i: out of memory!\n");
        exit(1);
    }

    err_to_col= malloc(sizeof(int )*nusefilt);
    if(err_to_col == NULL) {
        fprintf(stderr, "efilt_to_col i: out of memory!\n");
        exit(1);
    }
  
    filt_defnum = malloc(sizeof(int )*nusefilt);
    if(filt_defnum == NULL) {
        fprintf(stderr, "filter_defnum: out of memory!\n");
        exit(1);
    }
    
    zpfilterid = malloc(sizeof(int)*nusefilt);
    if(zpfilterid == NULL) {
        fprintf(stderr, "zpfilterid: out of memory!\n");
        exit(1);
    }
    
    totfiltid=0;
    j = 0;
    for(i=0; i<ncol;++i) {
        if (!(strcmp(col_data[i],"id"))) {
            if (id_col > -1) {
                fprintf(stderr,"Multiple id column definitions??\n");
                fprintf(stderr,"Check  col. %d vs. %d!  Quitting...\n",
                id_col,i);
                exit(1);
            } else id_col = i;
        } else if (col_data[i][0] == 'T') {
            /////// 'TOTxxx' column (awkward?)
            cp = col_data[i];
            ++cp;
            if (cp[0] == 'O') {
                ++cp;
                if (cp[0] == 'T' && isdigit(cp[1])) {
                    ++cp;
                    totfiltid = atoi(cp);
                    if ( (totfiltid < 1) || (totfiltid > nfilter) ) {
                        fprintf(stderr, "Bad filter number, %ld, for TOTAL column %d!\n",
                                totfiltid,i+1);
                        fprintf(stderr,"N.B.  There are only %ld filters defined in the filter file.\n",
                                nfilter);
                        exit(1);
                    }                         
                    totcolumn = i;
                }
            }
        } else if (col_data[i][0] == 'F' && isdigit(col_data[i][1])) {
            cp = col_data[i];
            ++cp;
            ifilt = atoi(cp);
            if ( (ifilt < 1) || (ifilt> nfilter) ) {
                fprintf(stderr, "Bad filter number, %d, for DATA column %d!\n",
                                ifilt,i+1);
                fprintf(stderr,"N.B.  There are only %ld filters defined in the filter file.\n",
                                nfilter);
                exit(1);
            } 
            for(k=0;k<j;++k) {
                if (ifilt == filt_defnum[k]) {
                    fprintf(stderr,"Duplicate filter definition for filter %d!?\n",
                                   ifilt);
                    fprintf(stderr,"Check column %d vs. column %d.\n",i,filt_to_col[k]);
                    exit(1);
                }
            }
            filt_defnum[j] = ifilt;
            filt_to_col[j] = i;
            
            if (zpfactor[i] > 0) zpfilterid[j] = ifilt; else zpfilterid[j]=-1;
            
            ++j;
        }
    }

    // now get error columns
    for (i=0;i<nusefilt;++i) err_to_col[i]=-1;

    for(i=0; i<ncol;++i) {
        if (col_data[i][0] == 'E' && isdigit(col_data[i][1])) {
            cp = col_data[i];
            ++cp;
            ifilt = atoi(cp);
            if ((ifilt < 1)|| (ifilt>nfilter)) {
                fprintf(stderr, "Bad filter number, %d, for error column %d!\n",
                                ifilt,i);
                fprintf(stderr,"N.B.  There are only %ld filters defined in the filter file.",
                               nfilter);
                exit(1);
            } 
            j = 0;
            while ((j<nusefilt)&&(ifilt!=filt_defnum[j])) {
                ++j;
            }
            if (j==nusefilt) {
                fprintf(stderr,"The filter number, %d, in error column %d doesn't match\n",
                                ifilt,i);
                fprintf(stderr,"any of the filters numbers defined by the data columns!\n");
                exit(1);
            }
            if (err_to_col[j]> -1) {
                fprintf(stderr,"Duplicate error column definition!?\n");
                fprintf(stderr,"Check column %d vs. column %d.\n",i,err_to_col[j]);
	            exit(1);
            }
            err_to_col[j] = i;
        }
    }

    // for convenience, now get pointers to used filters
    pusefilt = malloc(sizeof(filt_data *)*nusefilt);
    if (!(pusefilt)) {
        fprintf(stderr,"pusefilt: memory allocation error!\n");
        exit(1);
     }

    // lambdac is mean wavelength of filter
    lambdac = malloc(sizeof(double)*nusefilt);
    if (!(lambdac)) {
        fprintf(stderr,"malloc prob\n");
        exit(1);
    }
    for(j=0;j<nusefilt;++j) {
        pfilt_thru = &filt_thru;
        // walk the list
        for(i=1;i<filt_defnum[j];++i) pfilt_thru = pfilt_thru->pnext;
        pusefilt[j] = pfilt_thru;
    }
  
    for(j=0;j<nusefilt;++j) {
        avlambda = 0;
        totfilt = 0;
        for (i=1;i<pusefilt[j]->ndat;++i) {
            dl = pusefilt[j]->plambda[i]-pusefilt[j]->plambda[i-1];
            avlambda += dl*(pusefilt[j]->pthrough[i-1]*pusefilt[j]->plambda[i-1]+
                            pusefilt[j]->pthrough[i]*pusefilt[j]->plambda[i]);
            totfilt  += dl*(pusefilt[j]->pthrough[i-1]+pusefilt[j]->pthrough[i]);
        }  
        lambdac[j] = avlambda/totfilt;
    }
    
    /////// Smooth filters with a gaussian of width [SMOOTH_SIGMA] \AA
    if (SMOOTH_FILTERS==1 && ((!USE_TEMPLATE_CACHE && !DUMP_TEMPLATE_CACHE) 
           || DUMP_TEMPLATE_CACHE)) {

        printf("Smoothing filters with gaussian sig=%lf A.\n",SMOOTH_SIGMA);
        if (VERBOSE_LOG) fprintf(fplog,"# Smoothing filters with gaussian sig=%lf A.\n#\n",SMOOTH_SIGMA);

        for(j=0;j<nusefilt;++j) {
    
            for (i=0;i<pusefilt[j]->ndat;++i) {
                filt_smooth[i]=0.;
                g1 = 1./sqrt(2.*3.14159)/SMOOTH_SIGMA*
                     exp(-0.5*pow(pusefilt[j]->plambda[0]-pusefilt[j]->plambda[i],2)/SMOOTH_SIGMA/SMOOTH_SIGMA);
                for (k=1;k<pusefilt[j]->ndat;++k) {
                    dl = pusefilt[j]->plambda[k]-pusefilt[j]->plambda[k-1];
                    g2 = 1./sqrt(2.*3.14159)/SMOOTH_SIGMA*
                            exp(-0.5*pow(pusefilt[j]->plambda[k]-pusefilt[j]->plambda[i],2)/SMOOTH_SIGMA/SMOOTH_SIGMA);
                    filt_smooth[i]+=dl*(g1*pusefilt[j]->pthrough[k-1]+g2*pusefilt[j]->pthrough[k]);
                    g1=g2;
                }
            }
            for (i=0;i<pusefilt[j]->ndat;++i) {
                pusefilt[j]->pthrough[i]=filt_smooth[i]/2.;
            }
        }
    }
  

    if (id_col > -1) printf("\n\nID column found: %d\n",id_col);
    if (spec_col > -1) printf("z_spec column found: %d\n",spec_col);

    printf("\n\nSummary of filter columns found: \n");

    for (i=0;i<nusefilt;++i) {

        printf("Filter #%d, RES#%d: %s - lambda_c=%lf\n",i+1,filt_defnum[i],pusefilt[i]->filtname,lambdac[i]);
        printf("     [flux col: %d, error col.: %d]\n\n",filt_to_col[i],err_to_col[i]);
    
        if (VERBOSE_LOG)  {
            fprintf(fplog,"#  Filter #%d, RES#%d: %s - lambda_c=%lf\n",
                          i+1,filt_defnum[i],pusefilt[i]->filtname,lambdac[i]);
            fprintf(fplog,"#       [flux col: %d, error col.: %d]\n",
	        filt_to_col[i],err_to_col[i]);
        }    
    }

    ////// Check PRIOR_FILTER and change to zero index column of catalog data
    if (APPLY_PRIOR) {
        j=0;
        while (j<nusefilt) {
            if (filt_defnum[j] == PRIOR_FILTER) {
                printf("PRIOR_FILTER, %d, found: %s\n\n",PRIOR_FILTER,pusefilt[j]->filtname);
                PRIOR_FILTER = j;
                break;
            }
            ++j;
        } 
        if (j==nusefilt) {
            fprintf(stderr,"\nPRIOR_FILTER, %d, is not in the catalog.\n",PRIOR_FILTER);
            exit(1);
        }
    } 
    
    free(filt_smooth);    
}

void readfilters() { 
    FILE *fp;
    long fnum, i, itrash;
    char buff[180];
    double lambda, filterf;
    filt_data *pfilt_thru;

    if (!(fp = fopen(FILTERS_RES,"r"))) {
        fprintf(stderr,"Error opening filter definition file, %s!\n",
                       FILTERS_RES);
        exit(1);
    }
 
    pfilt_thru = & filt_thru;
    pfilt_thru->ndat = 0;
    pfilt_thru->pnext = 0;

    printf("Filters defined in file: (# -- ID string) \n");
    while(fgets(buff,180,fp)) {     
        ++nfilter;
        if (nfilter > 1) {
            pfilt_thru->pnext = malloc(sizeof(filt_thru));
            pfilt_thru=pfilt_thru->pnext;
            pfilt_thru->pnext = 0;
        } 
        
        sscanf(buff,"%ld %s",&fnum,pfilt_thru->filtname);
        printf(" %ld -- %s\n",nfilter,pfilt_thru->filtname);
        //      printf("Number of filter data lines: %ld\n",fnum);
      
        pfilt_thru->plambda = malloc(sizeof(lambda)*fnum);
        pfilt_thru->pthrough = malloc(sizeof(filterf)*fnum);
        pfilt_thru->ndat = fnum; 
        if (fnum > NTEMPLMAX) {
            fprintf(stderr,"Error: number of filter points, %ld, exceeds maximum allowed [%ld]!\n",
                    fnum,NTEMPLMAX);
            fprintf(stderr,"If necessary, increase size of NTEMPLMAX variable.\n");
            exit(1);
        }

      
        for (i=0;i<fnum;++i) {
            if (!(fgets(buff,180,fp))) {
                printf("ran out of lines! exiting...");
                exit(1);
            }
        
            sscanf(buff,"%ld %lf %lf",&itrash,&lambda,&filterf);
        
            ///// Convert energy/wavelength to #photons/wavelength
            // if (FILTER_FORMAT == 0) filterf=filterf*lambda;  
            ///// DON'T DO THIS HERE, DO IT IN INTEGRATION STEP
            
            pfilt_thru->plambda[i]= lambda;
            if (filterf < 0) filterf=0.0;
            pfilt_thru->pthrough[i]= filterf;
      
        }
    } 
    fclose(fp);
    printf("%ld filters read in!\n",nfilter);
}
	  
//// Read catalog data file
void readdata() {
    FILE *fp;
    long i,j,k;
    char buff[BUFFSIZE], *arg;
    double f,ef,total_scale;
    int firstline = 1;
    char  delim[]=" \n\t,\r";
    long IDBUFF=64,idsize_max;
    
    printf("\nProcessing data catalog [%s]...\n",CATALOG_FILE);

    if (!(fp = fopen(CATALOG_FILE,"r"))) {
        fprintf(stderr,"Hey, I can't read from %s!\n", CATALOG_FILE);
        exit(1);
    }
    j = 0;
    while (fgets(buff,BUFFSIZE,fp)) {
        // printf("buff %s %ld\n",buff,j);
        if (firstline)  {
            get_column_defs(buff);
            firstline = 0;
        }
        if (buff[0] !='#') ++j;
    }
    fclose(fp);
 
    // ** Emulate sm -- assume lines starting with "#" are comments to be ignored!
    nobj = j;
  
    if (nobj<1) {
        fprintf(stderr,
            "Hm... a null catalog file is not very interesting. I'm quitting."); 
        exit(1);
    }
      
    printf("Catalog file contains %ld objects.\n",nobj);
    fnu = malloc(sizeof(double *)*nobj);
    if(fnu == NULL) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }

    efnu = malloc(sizeof(double *)*nobj);
    if(efnu == NULL) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
	
    if (totfiltid > 0) {
        totflux = malloc(sizeof(double)*nobj);
        if (totflux == NULL) {
            fprintf(stderr,"totflux: out of memory\n");
            exit(1);
        }
    }
     	
    objid = malloc(sizeof(char *)*nobj);
    if(objid == NULL) {
        fprintf(stderr, "objid i: out of memory!\n");
        exit(1);
    }
    for (i=0;i<nobj;++i) {
        objid[i] = malloc(sizeof(char)*IDBUFF);
        if (objid[i]==NULL){
            fprintf(stderr, "objid j: out of memory!\n");
            exit(1);
        }
    }
    idsize_max=0;

    zspec = malloc(sizeof(f)*nobj);
    if(zspec == NULL) {
        fprintf(stderr, "zspec: out of memory\n");
        exit(1);
    }
    for (i=0;i<nobj;++i) zspec[i] = -1.0;
		
    printf("Again, I have %ld filters per object\n",nusefilt);
    
    for (i=0;i<nobj; ++i) {
        fnu[i] = malloc(sizeof(f)*nusefilt);
        if(fnu[i] == NULL) {
            fprintf(stderr, "out of memory\n");
            exit(1);
		}
        efnu[i] = malloc(sizeof(ef)*nusefilt);
        if(efnu[i] == NULL) {
            fprintf(stderr, "out of memory\n");
            exit(1);
        }
    }   
  
    printf("Catalog data space successfully allocated!\n");
   
    /////
    ///// Now read the suckers in!
    /////
    if (!(fp = fopen(CATALOG_FILE,"r"))) {
        fprintf(stderr, "Hey, I can't read from %s\n",CATALOG_FILE);
        exit(1);
    }  
     
    j = 0;
    while (fgets(buff,BUFFSIZE,fp))
        if (buff[0]!='#') {   // ignore comment "#" lines
        for (i=0;i<ncol;++i) {
            if (i)    
                arg=strtok(NULL,delim);
            else
                arg=strtok(buff,delim);   
            if (!arg) {
                fprintf(stderr,
                  "Something is wrong with the data in line %ld: \n",
                  j+1);
                fprintf(stderr, "%s\n", buff);
                fprintf(stderr, "Missing or corrupted columns?\n");
                exit(1);
            }
            strcpy(col_data[i],arg);
        }
        if (id_col > -1) strcpy(objid[j],col_data[id_col]); else strcpy(objid[j],"-1");
        idsize=0;
        while (objid[j][idsize] != '\0' && idsize < IDBUFF) ++idsize;
        if (idsize > idsize_max) idsize_max=idsize;
        
        if (spec_col >-1) zspec[j] = atof(col_data[spec_col]);
        for (i=0;i<nusefilt;++i) {
            fnu[j][i] = atof(col_data[filt_to_col[i]]);
            efnu[j][i] = atof(col_data[err_to_col[i]]);
        }
        ////// Read in absolute value of total fluxes
        if (totfiltid > 0) totflux[j] = fabs(atof(col_data[totcolumn]));
        
        // printf("%ld objid %ld zspec %lf\n",j,objid[j],zspec[j]);
        ++j;
    }
//    fprintf(stderr,"IDSIZE_MAX= %ld\n",idsize_max);
    idsize=idsize_max;
    
    fclose(fp);
    if (j!=nobj) {
        printf("Hoi! The number of objects I read doesn't match what I expected. Try again...\n");
        exit(1);
    }
    
    ///////
    /////// Scale to total fluxes
    ///////
    if (totfiltid > 0) {
        k=0;
        while (filt_defnum[k] != totfiltid && k < nusefilt) ++k;
        if (k == nusefilt) {
            fprintf(stderr,"No color flux column column found for total flux filter #%ld.  Skipping correction...\n",totfiltid);
        } else {
            printf("Applying total flux correction = column#%ld / column#%d...",totcolumn+1,filt_to_col[k]+1);    
            for (j=0;j<nobj;++j) if (fnu[j][k] > NOT_OBS_THRESHOLD && totflux[j] > 0) {
                total_scale = totflux[j] / fnu[j][k];
                //printf("HERE total: %lf = %lf / %lf\n",total_scale,totflux[j],fnu[j][k]);
                for (i=0;i<nusefilt;++i) {
                    fnu[j][i] *= total_scale;
                    efnu[j][i] *= total_scale;
                }
            }
            printf("Done.\n");
        }
    }
    printf("Data read in successfully!\n\n");
}

//// Read templates and interpolate to wavelength grid
void readtemplates() {
    FILE *fp, *tp;
    long i,j,k,istart;
    char buff[BUFFSIZE], fname[512], *arg;
    char  delim[]=" \n\t,\r";
    double *lambda_conv;
    double tlam[NTEMPLMAX],tf[NTEMPLMAX],tlamk,tfk,ts[NTEMPLMAX];
    double *templmid,*tempfmid,numsum,h;
    long ntlam;
    int itemp;
  
    ///////////////
    ////  Wavelength grid
    ///////////////
    printf("Opening template wavelength grid file, %s\n", WAVELENGTH_FILE);
    if (!(fp = fopen(WAVELENGTH_FILE,"r"))) {
        fprintf(stderr,"Oops... open: %s failed...\n",WAVELENGTH_FILE);
        exit(1);
    }
    i=0;
    while (fgets(buff,BUFFSIZE,fp)) if (buff[0]!='#') ++i;
    fclose(fp);
    NTEMPL = i;
    printf("Found %d wavelength grid points...\n",NTEMPL);

    templ =  malloc(sizeof(double)*NTEMPL);
    if (!(templ)) {
        fprintf(stderr,"templ: memory allocation error!\n");
        exit(1);
    }
    //// Wavelength midpoints
    templmid = malloc(sizeof(double)*NTEMPL);
    if (!(templmid)) {
        fprintf(stderr,"templmid: memory allocation error!\n");
        exit(1);
    }
    tempfmid = malloc(sizeof(double)*NTEMPL);        
    if (!(tempfmid)) {
        fprintf(stderr,"tempfmid: memory allocation error!\n");
        exit(1);
    }


    //// now read them in!
    if (!(fp = fopen(WAVELENGTH_FILE,"r"))) {
        fprintf(stderr,"Oops... re-open %s failed...\n",WAVELENGTH_FILE);
        exit(1);
    }
    i=0;
    while (fgets(buff,BUFFSIZE,fp)) if (buff[0]!='#') {
        arg=strtok(buff,delim);
        templ[i] =atof(arg);
        ++i;
    }
    fclose(fp);

    //// now do a sanity check on wavelength grid.
    //// Make sure it is monotonic increasing,with no duplicate
    //// values caused by precision problems!! 

    for(i=0;i<NTEMPL;++i) {
        if (templ[i] <0) {
            fprintf(stderr,"Hey! Negative wavelength value \n");
            fprintf(stderr,"for wavelength  # %ld: %lf?\n", i+1,templ[i]);
            exit(1);
        }
    }
    for(i=1;i<NTEMPL;++i) {
        if (templ[i] <= templ[i-1]) {
            fprintf(stderr,
                "Wavelength values not monotonically increasing in file!\n");
            fprintf(stderr,
                "Wavelength  # %ld [%lf} is not greater than preceeding value [%lf].\n",
                i+1,templ[i],templ[i-1]);
            fprintf(stderr,
                "This might be caused by insufficient precision for the wavelength values.\n");
            exit(1);
        }
    }
    
    ///// Check that wavelength grid contains
    ///// all filters over desired redshift range
    j=0;
    k=0;
    for (i=0;i<nusefilt;++i) {
        if (lambdac[i]<lambdac[j]) j=i; // min
        if (lambdac[i]>lambdac[k]) k=i; // max
    } 
      
    if (lambdac[j] < templ[0]*(1+Z_MAX)) {
        fprintf(stderr,"WARNING: Shortest wavelength filter falls off the wavelength grid at z >%6.3f.\n",
              lambdac[j]/templ[0]-1.);
    }
    if (lambdac[k] > templ[NTEMPL-1]*(1+Z_MIN)) {
        fprintf(stderr,"WARNING:  Longest wavelength filter falls off the wavelength grid at z <%6.3f.\n",
              lambdac[k]/templ[NTEMPL-1]-1.);
    }

    ///////////////
    ////  Template error file
    ///////////////    
    printf("\nOpening template error profile file, %s: ", TEMP_ERR_FILE);
    if (!(tp = fopen(TEMP_ERR_FILE,"r"))) {
        fprintf(stderr,"Oops! Opening %s failed.\n",TEMP_ERR_FILE);
        exit(1);
    }

    k=0;
    while (fgets(buff,BUFFSIZE,tp) && k < NTEMPLMAX) 
        if (buff[0]!='#') {
        	
        	arg=strtok(buff,delim);
        	tlam[k] = atof(arg);
 			arg=strtok(NULL,delim);
 			if (!(arg)) {
 				fprintf(stderr,"Template error file has only one column at line %ld.\n",k);
 				exit(1);
 			} else tf[k] = atof(arg);
 			
 			arg = strtok(NULL,delim);
 			if (arg) ts[k] = atof(arg); else ts[k] = 1.;
 			
 			//tf[k] = tfk;
            //sscanf(buff,"%lf %lf\n",&tlamk,&tfk);
            //tlam[k] = tlamk;
            //tf[k] = tfk;
            ++k;
    }
    
    if (k == NTEMPLMAX) {
        fprintf(stderr,"Too many lines [>%ld] in template error file!\n",NTEMPLMAX);
        exit(1);
    }
    printf(" %ld lines read in...\n",k);
    fclose(tp);
    ntlam = k;

    //// Interpolate template error file to standard wavelength grid
    //// Template error function is set to zero where the template error
    //// is not defined on the wavelength grid.
    k=0;
    while (templ[k] < tlam[0] && k < NTEMPL) {
        temp_err[k] = 0.0;
        temp_scale[k] = 0.0;
        ++k;
    }
    if (k > 0) fprintf(stderr,"Warning: template error function will be set to zero BELOW %le AA\n",tlam[0]);
    
    istart = 0;    
    while (templ[k] < tlam[ntlam-1] && k < NTEMPL) {
        interpol(templ[k],&(temp_err[k]),tlam,tf,ntlam,&istart);
        if (temp_err[k]<0) temp_err[k] = 0.0;
        interpol(templ[k],&(temp_scale[k]),tlam,ts,ntlam,&istart);
        if (temp_scale[k]<0) temp_err[k] = 1.0;
        ++k;
    }
    if (k < NTEMPL) {
        fprintf(stderr,"Warning: template error function will be set to zero BEYOND %le AA\n",templ[k]);
        while (k < NTEMPL) {
            temp_err[k] = 0.0;
            temp_scale[k] = 1.0;
            ++k;
        }
    }

    ///////////////
    ////  Templates file
    /////////////// 
    printf("Opening template list file, %s ...\n",TEMPLATES_FILE);
    if (!(fp = fopen(TEMPLATES_FILE,"r"))) {
        fprintf(stderr,"Oops! Opening %s failed.\n",TEMPLATES_FILE);
        exit(1);
    }

    j = 0;
    while (fgets(buff,BUFFSIZE,fp)) if (buff[0]!='#') ++j;
    NTEMP = j;
    fclose(fp);
    printf("File contains %d templates...\n\n",NTEMP);
    if (NTEMP == 0) {
        fprintf(stderr,"Empty template file, %s?\n",TEMPLATES_FILE);
        exit(1);
    }
    
    tnorm = malloc(sizeof(double)*NTEMP);
    if (!(tnorm)) {
        fprintf(stderr,"tnorm: memory allocation error!\n");
        exit(1);
    }
    
    lambda_conv = malloc(sizeof(double)*NTEMP);
    if (!(lambda_conv)) {
        fprintf(stderr,"lambda_conv: memory allocation error!\n");
        exit(1);
    }

    temp_err_a = malloc(sizeof(double)*NTEMP);
    if (temp_err_a==NULL) {
        fprintf(stderr,"temp_erra: memory allocation error!\n");
        exit(1);
    }

    tempage = malloc(sizeof(double)*NTEMP);
    if (tempage == NULL) {
        fprintf(stderr,"tempage: memory allocation error!\n");
        exit(1);
    }

    temp_combine = malloc(sizeof(int *)*NTEMP);
    if (temp_combine == NULL) {
        fprintf(stderr,"temp_combine i: out of memory!\n");
        exit(1);
    }
    for(i=0;i<NTEMP;++i) {
        temp_combine[i] = malloc(sizeof(int)*NTEMP);
        if (temp_combine[i] == NULL) {
            fprintf(stderr,"temp_combine j: out of memory!\n");
            exit(1);
        }
    }
    for (i=0;i<NTEMP;++i)
        for(j=0;j<NTEMP;++j) temp_combine[i][j] = 0;

    tempf = malloc(sizeof(double *)*NTEMP);
    if (tempf == NULL) {
        fprintf(stderr,"tempf i: out of memory!\n");
        exit(1);
    }
    for(i=0;i<NTEMP;++i) {
        tempf[i] = malloc(sizeof(double)*NTEMPL);
        if (tempf[i] == NULL) {
            fprintf(stderr,"tempf j: out of memory!\n");
            exit(1);
        }
    }
    
    if (!(fp = fopen(TEMPLATES_FILE,"r"))) {
        fprintf(stderr,"Oops... re-open of %s failed...\n", TEMPLATES_FILE);
        exit(1);
    }

    j=0;
    while (fgets(buff,BUFFSIZE,fp))
      if (buff[0]!='#') {   // ignore comment "#" lines
        
        arg = strtok(buff,delim);
        if (!arg) {
            fprintf(stderr,
              "Something is wrong with the line for template  %ld: \n",
              j+1);
            fprintf(stderr, "%s\n", buff);
            fprintf(stderr, "Blank line?\n");
            exit(1);
        } 
      
        k = atoi(arg)-1;
        if (j!=k) {
            fprintf(stderr,"Bad template number in template line %ld,\n", j+1);
            fprintf(stderr,"    ?? %s\n",buff);
            exit(1);
        }

        arg = strtok(NULL,delim);
        if (!(arg)) {
            fprintf(stderr,"No template file name in template line %ld,\n",j+1);
            fprintf(stderr,"    ?? %s\n",buff);
            exit(1);
        }
        strcpy(fname,arg);
        if (NTEMP < 40) printf("Working on template file %s\n",fname);
            
        arg=strtok(NULL,delim);
        if (!arg) {  
            fprintf(stderr,"%s: Prematurely truncated line for template # %ld!?: \n",
                TEMPLATES_FILE,j+1);
            fprintf(stderr, " -> %s?\n", buff);
            exit(1);
        }
        lambda_conv[j] = atof(arg);
        if (NTEMP < 40) printf("   Wavelength conversion factor to be used: %lf\n",
            lambda_conv[j]);
      
        arg=strtok(NULL,delim);
        if (!arg) {  
            fprintf(stderr,"%s: Prematurely truncated line for template # %ld!?: \n",
                TEMPLATES_FILE,j+1);
            fprintf(stderr, " no age -> %s?\n", buff);
            exit(1);
        }
        tempage[j] = atof(arg);
        if (NTEMP < 40) printf("   Template model age: %lf\n",
            tempage[j]);
      
        arg=strtok(NULL,delim);
        if (!arg) {  
	       fprintf(stderr,"%s: Prematurely truncated line for template # %ld!?: \n",
		      TEMPLATES_FILE,j+1);
	       fprintf(stderr, " no temp_err_a -> %s?\n", buff);
	       exit(1);
        }
        temp_err_a[j] = atof(arg);
        if (NTEMP < 40) printf("   Individual template error amplitude:  %lf\n",temp_err_a[j]);
      
        if (TEMPLATE_COMBOS == 1 || TEMPLATE_COMBOS == 99) {
            if (NTEMP < 40) printf(" Template combination mode: %d\n",TEMPLATE_COMBOS);       
            if (VERBOSE_LOG > 0 && NTEMP < 40) {
                fprintf(fplog,"#  Template %ld: %s\n",j+1,fname);
                fprintf(fplog,"#     %lf   %lf   %lf\n",
                    lambda_conv[j], tempage[j], temp_err_a[j]);
            }        
        } else if (TEMPLATE_COMBOS == 2) {
          i=0;
          arg=strtok(NULL,delim);
          while (arg) {
            itemp = atoi(arg);
            if ((itemp<1)||(itemp>NTEMP)) {
                fprintf(stderr,"Bad template number to combine with? %d\n",itemp);
                exit(1);
            }
            itemp = itemp-1;
            temp_combine[j][itemp] = 1;
            ++i;
            if (i>NTEMP-1) {
                fprintf(stderr,
                  "Too many templates listed to combine with template %ld!\n",j);
                exit(1);
            }
            arg=strtok(NULL,delim);
          }
      
          if (NTEMP < 40) {
            printf("   Combined fit with template(s) # "); 
            for (i=0;i<NTEMP;++i) {
              if (temp_combine[j][i]) { 
                printf(" %ld ",i+1);
              }
            }
            printf("\n\n");
          }
    
          if (VERBOSE_LOG > 0 && NTEMP < 40) {
            fprintf(fplog,"#  Template %ld: %s\n",j+1,fname);
            fprintf(fplog,"#     %lf   %lf   %lf\n",
                lambda_conv[j], tempage[j], temp_err_a[j]);
            fprintf(fplog,"#     combines with:");       
            for (i=0;i<NTEMP;++i) if (temp_combine[j][i])  
                fprintf(fplog," %ld ",i+1);
            fprintf(fplog,"\n");
          }
        } else {
            for (i=0;i<NTEMP;++i) 
                temp_combine[j][i] = 1;
            if (NTEMP < 40) printf(" Combined with _ALL_ other templates \n\n");

            if (VERBOSE_LOG > 0 && NTEMP < 40) {
                fprintf(fplog,"#  Template %ld: %s\n",j+1,fname);
                fprintf(fplog,"#     %lf   %lf   %lf\n",
                    lambda_conv[j], tempage[j], temp_err_a[j]);
                fprintf(fplog,"#     combines with _ALL_ other templates\n");       
            }

        }

        ///////////////
        ////  Individual templates
        ///////////////
        if (!(tp = fopen(fname,"r"))) {
            fprintf(stderr,"  Can't open the template file!  Bad name?\n");
            exit(1);
        }
                
        //// load template into temporary array
        k=0;
        while (fgets(buff,BUFFSIZE,tp)) 
            if (buff[0]!='#'){
                sscanf(buff," %lf %lf\n",&tlamk,&tfk);
                tlam[k] = tlamk*lambda_conv[j];
                tf[k] = tfk;
                ++k;                
                if (k>NTEMPLMAX-1) {
                    fprintf(stderr,
                    "Too many lines [>%ld] in template SED file!\n",NTEMPLMAX);
                    exit(1);
                }
        }
        ntlam = k;
        if (NTEMP < 40) printf("%ld SED lines read in...\n",ntlam);
        fclose(tp);


        //// now do a sanity check on wavelength grid.
        //// Make sure it is monotonic increasing,with no duplicate
        //// values caused by precision problems!! Otherwise spline
        //// and everything fails!!

        for(i=0;i<ntlam;++i) {
            if (tlam[i] <0) {
                fprintf(stderr,"Hey! Negative wavelength value for wavelength  # %ld: %lf?\n", i+1,tlam[i]);
                exit(1);
            }
        }
        for(i=1;i<ntlam;++i) {
            if (tlam[i] <= tlam[i-1]) {
                fprintf(stderr,
                    "Wavelength values not monotonically increasing in file!\n");
                fprintf(stderr,
                    "Wavelength  # %ld [%lf} is not greater than preceeding value [%lf].\n",
                    i+1,tlam[i],tlam[i-1]);
                fprintf(stderr, 
                    "This might be caused by insufficient precision for the wavelength values.\n");
                exit(1);
            }
        }

        ///// Check that max wavelength in wavelength grid is in template
        if (tlam[ntlam-1] < templ[NTEMPL-1]) {
            fprintf(stderr,"!! WARNING !! template %s stops at %lf, while the grid extends to %lf \n",
                                 fname,tlam[ntlam-1],templ[NTEMPL-1]);   
        }

        ///// Check that min wavelength in wavelength grid is in template
        if (tlam[0] > templ[0]) {
            fprintf(stderr,"!! WARNING !! template %s starts at %lf, while the grid starts at %lf \n",
                                 fname,tlam[0],templ[0]);   
        }
        
        /////// Now interpolate to standard wavelength grid, CONSERVING FLUX.  templmid
        /////// is the linear-interpolated midpoint between points on the master grid, tlam
        /////// is the grid read from the template file.
        ///////
        /////// lambda -->
        /////// master grid:       |       |       |       |       |       |       |       |
        /////// templmid:          x   |       |       |       |       |       |       |
        /////// tlam:          ^ ^  ^    ^ ^ ^   ^  ^ ^     ^     ^^^^   ^^      ^^^^    ^^
        ///////
        ///////                         _       _       ----     
        ///////                           _       __         --  
        ///////                             _        _             etc.
        ///////                               _       _
        ///////
        /////// Points on the master grid are integrations of the tlam grid at the steps indicated
        /////// above.
        
        i=0;
        /////// Get wavelength points at the midpoints of the master grid
        tempfmid[0] = 0.;
        templmid[0] = templ[0];
        for (k=1;k<NTEMPL;++k) {
            templmid[k] = (templ[k]+templ[k-1])/2.;

            //// Set master template flux to zero where tlam falls of the master grid
            if (templ[k] < tlam[0] || templ[k] > tlam[ntlam-1]) {
                tempfmid[k] = 0.;
                tempf[j][k] = 0.;
            } else {
                interpol(templmid[k],&tempfmid[k],tlam,tf,ntlam,&i);            
            }
        }
        
        ////// Rebin template grid to master wavelength grid, conserving template flux
        i=0;
        for (k=1;k<NTEMPL-1;++k) {
            
            if (templ[k] < tlam[0] || templ[k] > tlam[ntlam-1]) {
                tempf[j][k] = -100.0;   //// How to best handle when template doesn't cover full master grid?
                continue; //// next k
            }
                        
            numsum=0.;
            //// Go to where tlam is greater than the first midpoint
            while (tlam[i] < templmid[k] && i < ntlam) ++i;
            istart=i;
            
            /////// First point
            if (tlam[i] < templmid[k+1]) {
                h = tlam[i]-templmid[k];
                numsum+=h*(tf[i]+tempfmid[k]);
                //printf("  first t %lf %le\n",tlam[i],tf[i]);
                ++i;
            }
            
            if (i==0) ++i;
                    
            /////// Template points between master grid points
            while (tlam[i] < templmid[k+1] && i < ntlam) {
                h = tlam[i]-tlam[i-1];
                numsum+=h*(tf[i]+tf[i-1]);
                //printf("  mid   t %lf %le\n",tlam[i],tf[i]);
                ++i;
            }
            
            //// If no template points between master grid points, then just use interpolated midpoints
            if ( i == istart ) {
                h = templmid[k+1]-templmid[k];
                numsum=h*(tempfmid[k+1]+tempfmid[k]);
            } else {  
                ///// Last point              
                --i;
                h = templmid[k+1]-tlam[i];
                numsum+=h*(tempfmid[k+1]+tf[i]);
            }
            
            tempf[j][k] = numsum*0.5/(templmid[k+1]-templmid[k]); //*temp_scale[k];
        }
      
        ++j; 
    }

    fclose(fp);

    if (TEMPLATE_COMBOS == -2) TEMPLATE_COMBOS = 2;  /// FOR CASE statements in MAIN and GETPHOTZ
    
    //// Now normalize all templates to 5500 A flux

    k=0;
    while (templ[k]<5500.0 && k<NTEMPL) ++k;
    if (k==NTEMPL) {
        fprintf(stderr,"Wavelength file problem! Can't find 5500 A point??\n");
        exit(1);
    }

    for(i=0;i<NTEMP;++i) {
        tnorm[i] = tempf[i][k];
        
        for (j=0;j<NTEMPL;++j) {
            tempf[i][j]=tempf[i][j]/tnorm[i];
            /////// watch out for negative spline values at large lambda
            //if (tempf[i][j] < 0.0) tempf[i][j] = 0.0;
        }
    }

    free(templmid);
    free(tempfmid);
    free(lambda_conv);
    printf("Templates loaded.\n");
}
 


void getigmfactors (double ztarg, double *daz, double *dbz)
{ 
    double l2=1216.0,
         l3=1026.0,
         l4=973.0,
         l5=950.0,
         a2=3.6e-3,
         a3=1.7e-3,
         a4=1.2e-3,
         a5=9.3e-4,
         a,b,c,d;

    ///////// CHANGE LONG TO DOUBLE //////
    double lam1, lam2;
    double dl;

    if (ztarg <= 0.0) {
        *daz = 1;
        *dbz = 1;
        return;
    }
    
    lam1 =  (1050.0*(1+ztarg)+0.5);
    lam2 =  (1170.0*(1+ztarg)+0.5);
    
    *daz = 0.0;
    for (dl=lam1;dl<=lam2; dl+=1/4.) {
        a = exp(-a2*pow((dl/l2),3.46));
        *daz += a;
    }
    *daz = 1.0-*daz/((lam2-lam1)*4.);

    lam1  = (920.0*(1+ztarg)+0.5);
    lam2 =  (1015.0*(1+ztarg)+0.5);
  
    *dbz = 0;
    for (dl=lam1;dl<=lam2; dl+=1/4.)  {
        a=a3*pow(dl/l3,3.46);
        b=a4*pow(dl/l4, 3.46);
        c=a5*pow(dl/l5, 3.46);
        d=exp(-(a+b+c));
        *dbz += d;
    }

    *dbz = 1.0-*dbz/((lam2-lam1)*4.);

}

//////////////    
///// Generate template/redshift grid integrated through filters
//////////////
void makegrid() {

    char tempcache[1024], tempcache_info[1024],temp_file[1024];
    long i,j,k,ktemp, ifilt,NFILT;
    int **falls_off, **falls_off_lowz;
    double  filtsum, tempsum;
    double flamtofnu;
    double dw;
    
    double *igm_corr,*templz,*convert_detector;
    double *sedz, *filt_int; 

    FILE *fp, *tf;
  
    printf("\n\nConverting template SEDs to filter system...\n");

    //// get memory for filter/template/redshift grid
    tempfilt = malloc(sizeof(double *)*NZ);
    if(tempfilt == NULL) {
        fprintf(stderr, "tempfilt i: out of memory\n");
        exit(1);
    }
  
    for(i=0;i<NZ;++i) {
        tempfilt[i] = malloc(sizeof(double *)*NTEMP);
        if(tempfilt[i] == NULL) {
            fprintf(stderr, "tempfilt j: out of memory\n");
            exit(1);
        }
    }

    falls_off = malloc(sizeof(int *)*NTEMP);
    for (i=0;i<NTEMP;++i) {
        falls_off[i] = malloc(sizeof(int)*nusefilt);
        if(falls_off[i] == NULL) {
            fprintf(stderr, "falls_off: out of memory\n");
            exit(1);
        }
        for (j=0;j<nusefilt;++j) falls_off[i][j] = 0.;
        
    }
    falls_off_lowz = malloc(sizeof(int *)*NTEMP);
    for (i=0;i<NTEMP;++i) {
        falls_off_lowz[i] = malloc(sizeof(int)*nusefilt);
        if(falls_off_lowz[i] == NULL) {
            fprintf(stderr, "falls_off_lowz: out of memory\n");
            exit(1);
        }
        for (j=0;j<nusefilt;++j) falls_off_lowz[i][j] = 0.;
        
    }
        
    for(i=0;i<NZ;++i)
      for(j=0;j<NTEMP;++j)  {
        tempfilt[i][j] = malloc(sizeof(double)*nusefilt);
        if(tempfilt[i][j] == NULL) {
            fprintf(stderr, "tempfilt k: out of memory\n");
            exit(1);
        }
    }

    //// get memory for redshifted template errors  
    temp_errf = malloc(sizeof(double *)*NZ);
    if(temp_errf == NULL) {
        fprintf(stderr, "temp_errf i: out of memory\n");
        exit(1);
    }
  
    for(i=0;i<NZ;++i) {
        temp_errf[i] = malloc(sizeof(double)*nusefilt);
        if(temp_errf[i] == NULL) {
            fprintf(stderr, "temp_errf j: out of memory\n");
            exit(1);
        }
    }

    //// now compute template error at >(redshifted)< filter wavelengths
    for (i=0;i<NZ;++i) 
      for (j=0;j<nusefilt;++j) {
        for (k=0;k<NTEMPL;++k) if (templ[k]>=lambdac[j]/(1.+ztry[i])) break;
        temp_errf[i][j] = temp_err[k];
    }
  
    strcpy(tempcache,OUTPUT_DIRECTORY);
    strcat(tempcache,"/");
    strcat(tempcache,CACHE_FILE);
    //strcpy(tempcache,"tempfilt.dat");
    
    if (USE_TEMPLATE_CACHE && !DUMP_TEMPLATE_CACHE) {
        printf("Loading cache...");        
        if (!(fp = fopen(tempcache,"r"))) {
            fprintf(stderr,"uh,oh. Something went wrong opening cache file!\n");
            perror("help");
            exit(1);
        }

		for (i=0;i<4;++i) fread(&j,sizeof(long),1,fp); /// dummy NFILT NTEMP NZ NOBJ 
        for (i=0;i<NZ;++i)
          for (j=0;j<NTEMP;++j) {
            ktemp=fread(tempfilt[i][j],sizeof(filtsum)*nusefilt,1,fp);
            if (ktemp!=1) {
                fprintf(stderr, 
                    "Cache file corrupted! Couldn't read filters for %ld, %ld\n",
                    i,j);
            	exit(1);
            }
        }   
        fclose(fp);
        printf("Template cache loaded in.\n");
    } else {   //////// Do it the hard way and generate template grid
    
      igm_corr = malloc(sizeof(double)*NTEMPLMAX);
      if(igm_corr == NULL) {
        fprintf(stderr, "igm_corr: out of memory\n");
        exit(1);
      }
      
      templz = malloc(sizeof(double)*NTEMPLMAX);
      if(templz == NULL) {
        fprintf(stderr, "templz: out of memory\n");
        exit(1);
      }

      convert_detector = malloc(sizeof(double)*NTEMPLMAX);
      if(convert_detector == NULL) {
        fprintf(stderr, "convert_detector: out of memory\n");
        exit(1);
      }

      sedz = malloc(sizeof(double)*NTEMPLMAX);
      if(sedz == NULL) {
        fprintf(stderr, "sedz: out of memory\n");
        exit(1);
      }

      filt_int = malloc(sizeof(double)*NTEMPLMAX);
      if(filt_int == NULL) {
        fprintf(stderr, "filt_int: out of memory\n");
        exit(1);
      }
      
      fprintf(stderr,"Generating template grid ");
      
      //// Step through the redshift grid
      for (i=0;i<NZ;++i) {

        //// Progress indicator
        if (NZ >= 10) 
          if ((i % (NZ/10))==0) fprintf(stderr,">");

        for (j=0;j<NTEMPL;++j) 
            templz[j] = templ[j]*(1.0+ztry[i]);

        /////// These factors for integrating f_nu*dlambda.  Multiply by templz**2 for integrating f_nu*dnu
        if  (FILTER_FORMAT == 1) {
            ///// FILTER_FORMAT == 1 for a photon-counting detector, and integrating in F_nu
            for (j=0;j<NTEMPL;++j) 
                convert_detector[j] = 1./templz[j];
        } else {
            ///// FILTER_FORMAT == 0 for an energy-counting detector, and integrating in F_nu
            for (j=0;j<NTEMPL;++j) 
                convert_detector[j] = 1./(templz[j]*templz[j]);
        }
            
        ////// Get IGM factors if APPLY_IGM is set AND convert templates to F_nu
        if (APPLY_IGM) {
      
            //getigmfactors(ztry[i],&dasum,&dbsum);
           
            for (j=0;j<NTEMPL;++j) {

                flamtofnu = templz[j]/(5500.*(1.0+ztry[i]));
                flamtofnu = flamtofnu*flamtofnu;
      
                if ((templ[j]>=912.0)&&(templ[j]<1026.0)) 
                    igm_corr[j] = 1.0-dbsum[i];
                else 
                    igm_corr[j] = 1.0;
                
                //// lyman limit
                if ((templ[j]<912.0)) igm_corr[j] = 0. ;
      
                if ((templ[j]>=1026.0)&&(templ[j]<1216.0))
                    igm_corr[j] = 1.0 - dasum[i];
      
                igm_corr[j] = igm_corr[j]*flamtofnu;
      
            } /// NTEMPL
        } else {  ////// No IGM Correction
            for (j=0;j<NTEMPL;++j) {
                flamtofnu = templz[j]/(5500.*(1.0+ztry[i]));
                flamtofnu = flamtofnu*flamtofnu;
                igm_corr[j] = flamtofnu;
            }
        } ///// IF APPLY_IGM                

        for (ktemp = 0;ktemp<NTEMP;++ktemp) {  
            
            //////  apply IGM correction and change to f_nu
            for (j=0;j<NTEMPL;++j) 
                sedz[j] = tempf[ktemp][j]*igm_corr[j];
               
            ////// Now convolve template with filter responses
            for (ifilt=0;ifilt<nusefilt;++ifilt) {
                
                NFILT = pusefilt[ifilt]->ndat;
                
                ///// Interpolate filter to redshifted wavelength grid without extrapolation
                /////    AND
                ///// integrate template through filter responses
                /////
                ///// temp = rebinned f_nu template to master grid from readtemplates()
                ///// filt_int = interpolated filter curve at master_grid*(1+z)
                /////
                /////                              Integrate_trapezoid_rule(master_grid*(1+z), temp x filt_int)
                ///// Flux[z][template][filter] =                    -----------------------
                /////                                 Integrate_trapezoid_rule(master_grid*(1+z), filt_int)
              
                j=0;
                while (templz[j] <= pusefilt[ifilt]->plambda[0] && j < NTEMPL) ++j;  //// skip steps where filter not defined
                                
                ///// Check to make sure template defined at all filter points
                if (tempf[ktemp][j] < -99) {
                    if (!falls_off[ktemp][ifilt]) {
                        fprintf(stderr,"WARNING: Filter %d falls off template %ld at z>%6.3f, will ignore this template at those redshifts.\n",filt_defnum[ifilt],ktemp+1,ztry[i]);
                    }
                    falls_off[ktemp][ifilt]=1;
                    tempfilt[i][ktemp][ifilt] = -100;
                    continue;
                }
                k=j;
                while (templz[k] <= pusefilt[ifilt]->plambda[NFILT-1] && k < NTEMPL) ++k;
                if (tempf[ktemp][k] < -99 || k == NTEMPL) {
                    falls_off_lowz[ktemp][ifilt] = i;
                    tempfilt[i][ktemp][ifilt] = -100;
                    continue;
                }
                
                ////// Everything OK, start interpolation                
                k=0;
                interpol(templz[j],&filt_int[j],pusefilt[ifilt]->plambda,pusefilt[ifilt]->pthrough,NFILT,&k);
                
                //if (i==0 && ifilt==0) printf("HERE start %lf %lf %lf %lf\n",templz[j],pusefilt[ifilt]->plambda[k],pusefilt[ifilt]->pthrough[k],filt_int[j]);
                
                filtsum = 0.0;
                tempsum = 0.0;
                ++j;
                while (templz[j] <= pusefilt[ifilt]->plambda[NFILT-1] && j < NTEMPL) {
                    interpol(templz[j],&filt_int[j],pusefilt[ifilt]->plambda,pusefilt[ifilt]->pthrough,NFILT,&k);
                    //if (i==0 && ifilt==0) printf("HERE       %lf %lf %lf %lf\n",templz[j],pusefilt[ifilt]->plambda[k],pusefilt[ifilt]->pthrough[k],filt_int[j]);

                    dw = 1./templz[j-1] - 1./templz[j];    ///// Integrate in dnu
                    dw = templz[j]-templz[j-1];            ///// Integrate in dlambda
                    
                    tempsum += dw*(filt_int[j]*sedz[j]*convert_detector[j]+filt_int[j-1]*sedz[j-1]*convert_detector[j-1]);
                    filtsum += dw*(filt_int[j]*convert_detector[j]+filt_int[j-1]*convert_detector[j-1]);
               
                    //////////// Integrate_trapezoid_rule(master_grid*(1+z), 1/[master_grid*(1+z)] x temp x filt_int)
               //     dw = templz[j]-templz[j-1];
               //     filtsum += dw*(filt_int[j]/templz[j]+filt_int[j-1]/templz[j-1]);
               //     tempsum += dw*(filt_int[j]/templz[j]*sedz[j]+filt_int[j-1]/templz[j-1]*sedz[j-1]);
                
                    ++j;
                    
                } 
                tempfilt[i][ktemp][ifilt] = tempsum / filtsum; //// Full template cache

                // printf("z temp filt flux %lf %ld %ld %13.5e\n",ztry[i],ktemp,ifilt,tempfilt[i][ktemp][ifilt]);
        
            } //// End Filter For 
        } //// End template For
    } //// End redshfit For

    for (ktemp = 0;ktemp<NTEMP;++ktemp) {  
        for (ifilt=0;ifilt<nusefilt;++ifilt) {
            if (falls_off_lowz[ktemp][ifilt] > 0)
                fprintf(stderr,"WARNING: Filter %d falls off template %ld AND/OR master grid at z<%6.3f, will ignore this template at those redshifts.\n",filt_defnum[ifilt],ktemp+1,ztry[falls_off_lowz[ktemp][ifilt]]);
        }
    }
    
    fprintf(stderr," Done.\n");  //// Template cache done
    
    free(igm_corr);
    free(templz);
    free(sedz);
    free(filt_int);
    free(falls_off);
    free(falls_off_lowz);
    
	if (BINARY_OUTPUT) {
    	sprintf(temp_file,"%s/%s.temp_sed",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE);
    	tf = fopen(temp_file,"w");

        NTEMP32 = (int32_t) NTEMP;
        NZ32 = (int32_t) NZ;
        NTEMPL32 = (int32_t) NTEMPL;
        
    	fwrite(&NTEMP32,sizeof(int32_t),1,tf);
    	fwrite(&NTEMPL32,sizeof(int32_t),1,tf);
    	fwrite(&NZ32,sizeof(int32_t),1,tf);
    	fwrite(templ,sizeof(double)*NTEMPL,1,tf);
    	for (j=0;j<NTEMP;++j)
    		fwrite(tempf[j],sizeof(double)*NTEMPL,1,tf);
    	fwrite(dasum,sizeof(double)*NZ,1,tf);
    	fwrite(dbsum,sizeof(double)*NZ,1,tf);
    	fclose(tf);
    }

    if (DUMP_TEMPLATE_CACHE || BINARY_OUTPUT) {
    
        fprintf(stdout,"Dumping filter templates ...");

        if (!(fp = fopen(tempcache,"w"))) {
            fprintf(stderr,"\n\nuh,oh. Something went wrong with writing cache! Does ./%s exist?\n",OUTPUT_DIRECTORY);
            perror("help");
            exit(1);
        }
        
        nusefilt32 = (int32_t) nusefilt;
        NTEMP32 = (int32_t) NTEMP;
        NZ32 = (int32_t) NZ;
        nobj32 = (int32_t) nobj;
        
        fwrite(&nusefilt32,sizeof(int32_t),1,fp);
        fwrite(&NTEMP32,sizeof(int32_t),1,fp);
        fwrite(&NZ32,sizeof(int32_t),1,fp);
        fwrite(&nobj32,sizeof(int32_t),1,fp);
        for (i=0;i<NZ;++i) 
          for (j=0;j<NTEMP;++j) 
            fwrite(tempfilt[i][j],sizeof(filtsum)*nusefilt,1,fp);
				
        fwrite(lambdac,sizeof(double)*nusefilt,1,fp);

		fwrite(ztry,sizeof(double)*NZ,1,fp);
                
        for (i=0;i<nobj;++i)
        	fwrite(fnu[i],sizeof(double)*nusefilt,1,fp);

        for (i=0;i<nobj;++i)
        	fwrite(efnu[i],sizeof(double)*nusefilt,1,fp);
        
        fclose(fp);
        
        fprintf(stdout," Done.\n\n");
    
    	////////////  Dump info of binary cache grid to [tempcache].info
    	///// Read in cache with IDL:
    	///// tempfilt = dblarr(NFILT,NTEMP,NZ)
    	///// openr,lun,'tempfilt.dat',/get_lun ; [you might need /swap_endian on Intel Macs with Universal (pre-Intel) IDL versions
    	///// readu,lun,tempfilt
    	///// close,/all
    	sprintf(tempcache_info,"%s/%s.readbin.pro",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE);
//    	strcpy(tempcache_info,tempcache);
//    	strcat(tempcache_info,".info");
 		fp = fopen(tempcache_info,"w");
		fprintf(fp,"; s[0:3] = NFILT NTEMP NZ NOBJ \n");
		fprintf(fp,"; %ld %d %d %ld\n; Read with IDL:\n",nusefilt,NTEMP,NZ,nobj);
		fprintf(fp,"; Template & catalog fluxes\n openr,lun,'%s',/swap_if_big_endian,/get_lun ; might not need 'swap_endian'\n",tempcache);
		fprintf(fp," s = lonarr(4) & readu,lun,s\n");
		fprintf(fp," NFILT=s[0] & NTEMP = s[1] & NZ = s[2] & NOBJ = s[3]\n");
		fprintf(fp," tempfilt = dblarr(NFILT,NTEMP,NZ)\n");
		fprintf(fp," lc = dblarr(NFILT) ; central wavelengths\n");
		fprintf(fp," zgrid = dblarr(NZ)\n");
		fprintf(fp," fnu = dblarr(NFILT,NOBJ)\n");
		fprintf(fp," efnu = dblarr(NFILT,NOBJ)\n");
		fprintf(fp," readu,lun,tempfilt,lc,zgrid,fnu,efnu\n close,/all\n");
		
		if (BINARY_OUTPUT) {
		    //// OBS_SED -> .coeff
			fprintf(fp,"; \n; Coeffs (OBS_SED)\n openr,lun,'%s/%s.coeff',/get_lun,/swap_if_big_endian\n",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE);
			fprintf(fp," s = lonarr(4) & readu,lun,s\n");
			fprintf(fp," NFILT=s[0] & NTEMP = s[1] & NZ = s[2] & NOBJ = s[3]\n");
			fprintf(fp," coeffs = dblarr(NTEMP,NOBJ)\n");
			fprintf(fp," izbest = lonarr(NOBJ)\n");
			fprintf(fp," tnorm = dblarr(NTEMP)\n");
			fprintf(fp," readu,lun,coeffs,izbest,tnorm \n close,/all\n");

			//// TEMP_SED
			fprintf(fp,"; \n; Full templates (TEMP_SED)\n openr,lun,'%s/%s.temp_sed',/get_lun,/swap_if_big_endian\n",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE);
			fprintf(fp," s = lonarr(3) & readu,lun,s\n");
			fprintf(fp," NTEMP=s[0] & NTEMPL = s[1] & NZ = s[2]\n");
			fprintf(fp," templam = dblarr(NTEMPL)\n");
			fprintf(fp," temp_seds = dblarr(NTEMPL,NTEMP)\n");
			fprintf(fp," da = dblarr(NZ)\n");
			fprintf(fp," db = dblarr(NZ)\n");
			fprintf(fp," readu,lun,templam,temp_seds,da,db \n close,/all\n");

			//// POFZ_FILE 
			fprintf(fp,"; \n; POFZ\n openr,lun,'%s/%s.pz',/get_lun,/swap_if_big_endian\n",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE);
			fprintf(fp," s = lonarr(2) & readu,lun,s\n");
			fprintf(fp," NZ = s[0] & NOBJ = s[1]\n");
			fprintf(fp," pz = dblarr(NZ,NOBJ)\n");
			fprintf(fp," readu,lun,pz \n close,/all\n");
			
			fprintf(fp," end\n");
		}
		
//		fprintf(fp,"#"); j=0; i=0; for (k=0;k<nusefilt;++k) fprintf(fp," %lf",tempfilt[i][j][k]); fprintf(fp,"\n");
//		fprintf(fp,"#"); j=1; i=0; for (k=0;k<nusefilt;++k) fprintf(fp," %lf",tempfilt[i][j][k]); fprintf(fp,"\n");
				
 		//fprintf(fp,"# z \n");
 		//for (i=0;i<NZ;++i) fprintf(fp," %lf\n",ztry[i]);
 		
 		fclose(fp);
 		
        if (USE_TEMPLATE_CACHE == 0 && BINARY_OUTPUT==0) exit(0);

    }
  }
  printf("\nTemplate grid computed/loaded!\n");

}

//// Initialize z grid, read in filters/data/templates      
void init() {
    long i;
    double zval; 
  
    readfilters();   
    readdata();   
    readtemplates();

    //// define redshift grid
    printf("\n\nMaking redshift grid...\n");
    NZ = 0;
    zval = Z_MIN;
    //// figure out by brute force how z values there are...
    while (zval <= Z_MAX) {
        if (Z_STEP_TYPE)
            zval += Z_STEP*(1.+zval);
        else
            zval += Z_STEP;
        ++NZ;
    }
    printf("Number of redshift points: %d\n",NZ);
    
    //// now allocate and fill redshift grid   
    ztry = malloc(sizeof(double)*NZ);
    if (ztry == NULL) {
        fprintf(stderr, "ztry: out of memory\n");
        exit(1);
    }
    zval = Z_MIN;
    ztry[0] = zval;
    for (i=1;i<NZ;++i) {
        if (Z_STEP_TYPE)
            zval += Z_STEP*(1.+zval);
        else
        zval += Z_STEP;
        ztry[i] = zval;
    }
    
    ///// Compute IGM factors
    dasum = malloc(sizeof(double)*NZ);
    if (dasum == NULL) {
        fprintf(stderr, "dasum: out of memory\n");
        exit(1);
    }
    dbsum = malloc(sizeof(double)*NZ);
    if (dbsum == NULL) {
        fprintf(stderr, "dbsum: out of memory\n");
        exit(1);
    }
    for (i=0;i<NZ;++i) {
        getigmfactors(ztry[i],&dasum[i],&dbsum[i]);
//        printf("igm: %lf %le %le\n",ztry[i],dasum[i],dbsum[i]);
    }

    makegrid();
        
}
