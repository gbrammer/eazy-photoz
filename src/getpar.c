#include "defs.h"

struct param_data {
     char name[256];
     int defined;
     char default_value[1024];
     char comment[1024];
     };

#define NKEY 51

void assign_par(int i, char *arg) {
  switch (i) {
  case 0: 
    strcpy(FILTERS_RES,arg);
    printf("FILTERS_RES: %s\n",FILTERS_RES);
    break;
  case 1:
    strcpy(TEMPLATES_FILE,arg);
    printf("TEMPLATES_FILE: %s\n",TEMPLATES_FILE);
    break;
  case 2: 
    strcpy(WAVELENGTH_FILE,arg);
    printf("WAVELENGTH_FILE: %s\n",WAVELENGTH_FILE);
    break;
  case 3:
    strcpy(TEMP_ERR_FILE,arg);
    printf("TEMP_ERR_FILE: %s\n",TEMP_ERR_FILE);
    break;
  case 4:
    strcpy(CATALOG_FILE,arg);
    printf("CATALOG_FILE: %s\n",CATALOG_FILE);
    break;
  case 5:
    NOT_OBS_THRESHOLD=atof(arg);
    printf("NOT_OBS_THRESHOLD: %lf\n",NOT_OBS_THRESHOLD);
    break;
  case 6:
    N_MIN_COLORS=atoi(arg);
    printf("N_MIN_COLORS: %d\n",N_MIN_COLORS);
    break;
  case 7:
    strcpy(OUTPUT_DIRECTORY,arg);
    printf("OUTPUT_DIRECTORY: %s\n",OUTPUT_DIRECTORY);
    break;
  case 8:
    strcpy(MAIN_OUTPUT_FILE,arg);
    printf("MAIN_OUTPUT_FILE: %s\n",MAIN_OUTPUT_FILE);
    break;
  case 9:
    if (*arg=='n' || *arg=='N' || *arg=='0') OBS_SED_FILE = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') OBS_SED_FILE = 1;
    else {
      fprintf(stderr,"\n OBS_SED_FILE:  bad value,  %s\n",arg);
      exit(1);
    }
    printf("OBS_SED_FILE: %c\n",*arg);
    break;
  case 10:
    if (*arg=='n' || *arg=='N' || *arg=='0') TEMP_SED_FILE = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') TEMP_SED_FILE = 1;
    else {
      fprintf(stderr,"\n TEMP_SED_FILE: bad value, %s\n",arg);
      exit(1);
    }
    printf("TEMP_SED_FILE: %c\n",*arg);
    break;
  case 11:
    if (*arg=='n' || *arg=='N' || *arg=='0') POFZ_FILE = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') POFZ_FILE = 1;
    else {
      fprintf(stderr,"\n POFZ_FILE:  bad value, %s\n",arg);
      exit(1);
    }
    printf("POFZ_FILE: %c\n",*arg);
    break;
  case 12:
    H0 = atof(arg);
    printf("H0: %lf\n",H0);
    break;
  case 13:
    OMEGA_M = atof(arg);
    printf("OMEGA_M: %lf\n",OMEGA_M);
    break;
  case 14:
    OMEGA_L = atof(arg);
    printf("OMEGA_L: %lf\n",OMEGA_L);
    break;
  case 15:
    if (*arg=='n' || *arg=='N' || *arg=='0') DUMP_TEMPLATE_CACHE = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') DUMP_TEMPLATE_CACHE = 1;
    else {
      fprintf(stderr,"\n DUMP_TEMPLATE_CACHE:  bad value, %s\n",arg);
      exit(1);
    }
    printf("DUMP_TEMPLATE_CACHE: %c\n",*arg);
    break;
  case 16:
    if (*arg=='n' || *arg=='N' || *arg=='0') USE_TEMPLATE_CACHE = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') USE_TEMPLATE_CACHE = 1;
    else {
      fprintf(stderr,"\n USE_TEMPLATE_CACHE: bad value, %s\n",arg);
      exit(1);
    }
    printf("USE_TEMPLATE_CACHE: %c\n",*arg);
    break;
  case 17:
    Z_MIN = atof(arg);
    printf("Z_MIN: %lf\n",Z_MIN);
    break;
  case 18:
    Z_MAX = atof(arg);
    printf("Z_MAX: %lf\n",Z_MAX);
    break;
  case 19:
    Z_STEP = atof(arg);
    printf("Z_STEP: %lf\n",Z_STEP);
    break;
  case 20:
    Z_STEP_TYPE = atoi(arg);
    printf("Z_STEP_TYPE: %d\n",Z_STEP_TYPE);
    break;
  case 21:
    if (*arg == 'a') TEMPLATE_COMBOS=99; else TEMPLATE_COMBOS = atoi(arg);
    /*if (*arg=='n' || *arg=='N' || *arg=='0') TEMPLATE_COMBOS = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') TEMPLATE_COMBOS = 1;
    else if (*arg=='a' || *arg=='A' || *arg=='2') TEMPLATE_COMBOS = 2;
    else*/
    if (TEMPLATE_COMBOS !=1 && TEMPLATE_COMBOS !=2 && TEMPLATE_COMBOS != -2 && TEMPLATE_COMBOS !=99) {
      fprintf(stderr,"\n TEMPLATE_COMBOS:  bad value, %s\n",arg);
      exit(1);
    }
    printf("TEMPLATE_COMBOS: %c\n",*arg);
    break;
  case 22:
    TEMP_ERR_A2 = atof(arg);
    printf("TEMP_ERR_A2: %lf\n",TEMP_ERR_A2);
    break;
  case 23:
    SYS_ERR = atof(arg);
    printf("SYS_ERR: %lf\n",SYS_ERR);
    break;
  case 24:
    if (*arg=='n' || *arg=='N' || *arg=='0') APPLY_IGM = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') APPLY_IGM = 1;
    else {
      fprintf(stderr,"\n APPLY_IGM:  bad value, %s\n",arg);
      exit(1);
    }
    printf("APPLY_IGM: %c\n",*arg);
    break;
  case 25:
    if (*arg=='n' || *arg=='N' || *arg=='0') FIX_ZSPEC = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') FIX_ZSPEC = 1;
    else {
      fprintf(stderr,"\n FIX_ZSPEC:  bad value, %s\n",arg);
      exit(1);
    }
    printf("FIX_ZSPEC: %c\n",*arg);
    break;
  case 26:
    if (*arg=='n' || *arg=='N' || *arg=='0') APPLY_PRIOR = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') APPLY_PRIOR = 1;
    else {
      fprintf(stderr,"\n APPLY_PRIOR:  bad value, %s\n",arg);
      exit(1);
    }
    printf("APPLY_PRIOR: %c\n",*arg);
    break;
  case 27:
    strcpy(PRIOR_FILE,arg);
    printf("PRIOR_FILE: %s\n",PRIOR_FILE);
    break;
  case 28:
    PRIOR_ABZP = atof(arg);
    printf("PRIOR_ABZP: %lf\n",PRIOR_ABZP);
    break;
  case 29:
    PRIOR_FILTER = atoi(arg);
    printf("PRIOR_FILTER: %d\n",PRIOR_FILTER);
    break;
  case 30:
    NMF_TOLERANCE = atof(arg);
    printf("NMF_TOLERANCE: %lf\n",NMF_TOLERANCE);
    break;
  case 31:
    if (*arg=='n' || *arg=='N' || *arg=='0') PRINT_ERRORS = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') PRINT_ERRORS = 1;
    else {
      fprintf(stderr,"\n PRINT_ERRORS:  bad value, %s\n",arg);
      exit(1);
    }
    printf("PRINT_ERRORS: %c\n",*arg);
    break;
  case 32:
    if (*arg=='n' || *arg=='N' || *arg=='0') VERBOSE_LOG = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') VERBOSE_LOG = 1;
    else {
      fprintf(stderr,"\n VERBOSE_LOG:  bad value, %s\n",arg);
      exit(1);
    }
    printf("VERBOSE_LOG: %c\n",*arg);
    break;
  case 33:
    if (*arg=='n' || *arg=='N' || *arg=='0') SMOOTH_FILTERS = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') SMOOTH_FILTERS = 1;
    else {
      fprintf(stderr,"\n SMOOTH_FILTERS:  bad value, %s\n",arg);
      exit(1);
    }
    printf("SMOOTH_FILTERS: %c\n",*arg);
    break;
  case 34:
    SMOOTH_SIGMA = atof(arg);
    printf("SMOOTH_SIGMA: %lf\n",SMOOTH_SIGMA);
    break;
  case 35:
    strcpy(CACHE_FILE,arg);
    printf("CACHE_FILE: %s\n",CACHE_FILE);
    break;
  case 36:
    if (*arg=='0') FILTER_FORMAT=0;
    else if (*arg=='1') FILTER_FORMAT=1;
    else {
        fprintf(stderr,"\n FILTER_FORMAT: bad value, %s\n",arg);
        exit(1);
    }
    printf("FILTER_FORMAT: %d\n",FILTER_FORMAT);
    break;
  case 37:
    if (*arg=='n' || *arg=='N' || *arg=='0') GET_ZP_OFFSETS = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') GET_ZP_OFFSETS = 1;
    else {
      fprintf(stderr,"\n GET_ZP_OFFSETS:  bad value, %s\n",arg);
      exit(1);
    }
    printf("GET_ZP_OFFSETS: %c\n",*arg);
    break;
  case 38:
    ZP_OFFSET_TOL = atof(arg);
    printf("ZP_OFFSET_TOL: %lf\n",ZP_OFFSET_TOL);
    break;
  case 39:
    CHI2_SCALE = atof(arg);
    printf("CHI2_SCALE: %lf\n",CHI2_SCALE);
    break;
  case 40:
    if (*arg=='n' || *arg=='N' || *arg=='0') BINARY_OUTPUT = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') BINARY_OUTPUT = 1;
    else {
      fprintf(stderr,"\n BINARY_OUTPUT:  bad value, %s\n",arg);
      exit(1);
    }
    printf("BINARY_OUTPUT: %c\n",*arg);
    break;
  case 41: 
    strcpy(REST_FILTERS,arg);
    printf("REST_FILTERS: %s\n",REST_FILTERS);
    break;
  case 42: 
    strcpy(Z_COLUMN,arg);
    printf("Z_COLUMN: %s\n",Z_COLUMN);
    break;      
  case 43:
    if (*arg=='n' || *arg=='N' || *arg=='0') USE_ZSPEC_FOR_REST = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') USE_ZSPEC_FOR_REST = 1;
    else {
        fprintf(stderr,"\n USE_ZSPEC_FOR_REST:  bad value, %s\n",arg);
        exit(1);
    }
    printf("USE_ZSPEC_FOR_REST: %c\n",*arg);
    break;
  case 44:
    if (*arg=='n' || *arg=='N' || *arg=='0') MAGNITUDES = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') MAGNITUDES = 1;
      else {
          fprintf(stderr,"\n MAGNITUDES:  bad value, %s\n",arg);
          exit(1);
      }
    printf("MAGNITUDES: %c\n",*arg);
    break;
  case 45:
    SCALE_2175_BUMP = atof(arg);
    printf("SCALE_2175_BUMP: %lf\n",SCALE_2175_BUMP);
    break;
  case 46:
    if (*arg=='n' || *arg=='N' || *arg=='0') {
        READ_ZBIN = 0;
        strcpy(ZBIN_FILE,"no");
    } else if (*arg=='y' || *arg=='Y' || *arg=='1') {
        READ_ZBIN = 1;
        strcpy(ZBIN_FILE,"yes");
    } else {
        strcpy(ZBIN_FILE,arg);
        READ_ZBIN = 1;
        //fprintf(stderr,"\n READ_ZBIN:  bad value, %s\n",arg);
        //exit(1);
      }
      printf("READ_ZBIN: %s %d\n",arg,READ_ZBIN);
      break;
  case 47:
    RF_PADDING = atof(arg);
    printf("RF_PADDING: %lf\n",RF_PADDING);
    break;
  case 48:
    if (*arg=='n' || *arg=='N' || *arg=='0') RF_ERRORS = 0;
    else if (*arg=='y' || *arg=='Y' || *arg=='1') RF_ERRORS = 1;
    else {
        fprintf(stderr,"\n RF_ERRORS:  bad value, %s\n",arg);
        exit(1);
    }
    printf("RF_ERRORS: %c\n",*arg);
    break;
  case 49:
    strcpy(LAF_FILE,arg);
    printf("LAF_FILE: %s\n",LAF_FILE);
    break;
  case 50:
    strcpy(DLA_FILE,arg);
    printf("DLA_FILE: %s\n",DLA_FILE);
    break;
  default:
    fprintf(stderr,"I shouldn't be here!!\n");
    fprintf(stderr,"%d -- %s\n",i,arg);
    exit(1);
    break;
  }
}
	 

void getparams()
{  char buff[1024];
   FILE *fp;
   char *arg;

   int ndef = 0;
   int nline = 0;
   int argcount,i=0;

   // Parameter definitions
   struct param_data params[NKEY] = {
     // USE hyperz format as base
     //   {"AOVSED", 0, "NIL"},  // Vega SED
     {"FILTERS_RES", 0 , "FILTER.RES.latest","Filter transmission data"},            // 0 
     {"TEMPLATES_FILE", 0, "templates/eazy_v1.2_dusty.spectra.param","Template definition file"},       // 1 
     {"WAVELENGTH_FILE",0,"templates/EAZY_v1.1_lines/lambda_v1.1.def","Wavelength grid definition file"},    // 2
     {"TEMP_ERR_FILE",0,"templates/TEMPLATE_ERROR.eazy_v1.0","Template error definition file"},     // 3 
     {"CATALOG_FILE", 0, "hdfn_fs99_eazy.cat","Catalog data file"},                      // 4
     {"NOT_OBS_THRESHOLD",0, "-90","Ignore flux point if <NOT_OBS_THRESH"},  // 5
     {"N_MIN_COLORS",0,"5","Require N_MIN_COLORS to fit"},                    // 6
     {"OUTPUT_DIRECTORY", 0, "OUTPUT","Directory to put output files in"},      // 7
     {"MAIN_OUTPUT_FILE", 0, "photz","Main output file, .zout"},                 // 8
     {"OBS_SED_FILE", 0, "n","Write out observed SED/object, .obs_sed"},      // 9
     {"TEMP_SED_FILE",0,"n","Write out best template fit/object, .temp_sed"}, // 10
     {"POFZ_FILE",0,"n","Write out Pofz/object, .pz"},                        // 11
     {"H0",0,"70.0","Hubble constant (km/s/Mpc)"},                                      // 12
     {"OMEGA_M",0,"0.3","Omega_matter"},                                      // 13
     {"OMEGA_L",0,"0.7","Omega_lambda"},                                      //14
     {"DUMP_TEMPLATE_CACHE",0,"n","Write binary template cache"},             // 15
     {"USE_TEMPLATE_CACHE",0,"n","Load in template cache"},                   // 16
     {"Z_MIN",0,"0.01","Minimum redshift"},                                   //17
     {"Z_MAX",0,"6.0","Maximum redshift"},                                    //  18  
     {"Z_STEP",0,"0.01","Redshift step size"},                                //19
     {"Z_STEP_TYPE",0,"1"," 0 = ZSTEP, 1 = Z_STEP*(1+z)"},                    //20
     {"TEMPLATE_COMBOS",0,"a","Template combination options: "},// 21
     {"TEMP_ERR_A2",0,"0.50","Template error amplitude"}, // 22
     {"SYS_ERR",0,"0.00","Systematic flux error (% of flux)"},                // 23
     {"APPLY_IGM",0,"y","Apply Madau 1995 IGM absorption"},                   // 24
     {"FIX_ZSPEC",0,"n","Fix redshift to catalog zspec"},                     // 25
     {"APPLY_PRIOR",0,"y","Apply apparent magnitude prior"},                          // 26
     {"PRIOR_FILE",0,"templates/prior_K_extend.dat","File containing prior grid"},       // 27
     {"PRIOR_ABZP",0,"25.0","AB zeropoint of fluxes in catalog.  Needed for calculating apparent mags!"},             // 28
     {"PRIOR_FILTER",0,"28","Filter from FILTER_RES corresponding to the columns in PRIOR_FILE"}, //29
     {"NMF_TOLERANCE",0,"1.e-4","Tolerance for non-negative combinations (TEMPLATE_COMBOS=a)"}, //30
     {"PRINT_ERRORS",0,"y","Print 68, 95 and 99% confidence intervals"}, //31
     {"VERBOSE_LOG",0,"y","Dump information from the run into [MAIN_OUTPUT_FILE].param"}, //32
     {"SMOOTH_FILTERS",0,"n","Smooth filter curves with Gaussian"}, //33
     {"SMOOTH_SIGMA",0,"100.","Gaussian sigma (in Angstroms) to smooth filters"}, //34
     {"CACHE_FILE",0,"photz.tempfilt","Template cache file (in OUTPUT_DIRECTORY)"}, //35
     {"FILTER_FORMAT",0,"1","Format of FILTERS_RES file -- 0: energy-  1: photon-counting detector"}, //36
     {"GET_ZP_OFFSETS",0,"n","Look for zphot.zeropoint file and compute zeropoint offsets"}, //37
     {"ZP_OFFSET_TOL",0,"1.e-4","Tolerance for iterative fit for zeropoint offsets [not implemented]"}, //38
     {"CHI2_SCALE",0,"1.0","Scale ML Chi-squared values to improve confidence intervals"}, //39
     {"BINARY_OUTPUT",0,"y","Save OBS_SED, TEMP_SED, PZ in binary format to read with e.g IDL"}, //40
     {"REST_FILTERS",0,"---","Comma-separated list of rest frame filters to compute"}, //41
     {"Z_COLUMN",0,"z_peak","Redshift to use for rest-frame color calculation (z_a, z_p, z_m1, z_m2, z_peak)"}, //42
     {"USE_ZSPEC_FOR_REST",0,"y","Use z_spec when available for rest-frame colors"}, //43
     {"MAGNITUDES",0,"n","Catalog photometry in magnitudes rather than f_nu fluxes"}, //44
     {"SCALE_2175_BUMP",0,"0.00","Scaling of 2175A bump.  Values 0.13 (0.27) absorb ~10 (20) % at peak."}, //45
     {"READ_ZBIN",0,"n","Get redshifts from OUTPUT_DIRECTORY/MAIN_OUTPUT_FILE.zbin rather than fitting them."}, //46
     {"RF_PADDING",0,"1000.","Padding (Ang) for choosing observed filters around specified rest-frame pair."}, //47
     {"RF_ERRORS",0,"n","Compute RF color errors from p(z)"}, //48
     {"LAF_FILE",0,"templates/LAFcoeff.txt","File containing the Lyman alpha forest data from Inoue"}, //49
     {"DLA_FILE",0,"templates/DLAcoeff.txt","File containing the damped Lyman absorber data from Inoue"} //50
   };


    // delimeter/whitespace characters for token search
    char delim[]=" \n\t";

    fp = fopen(zphot_param_file,"r");
    ///////// Generate zphot.param.default if zphot.param not found
    if (!(fp)) {
        
        if (strcmp(zphot_param_file,"zphot.param")) {
            fprintf(stderr,
                "<%s> not found.\nRun EAZY with no options and no zphot.param file to generate <zphot.param.default>.\n",zphot_param_file);            
            exit(1);
        } else 
            fprintf(stderr,"<%s> not found.  Generating <zphot.param.default>\n",zphot_param_file);
        
        fp = fopen("zphot.param.default","w");
        fprintf(fp,"#### EAZY Default parameters\n\n");
        fprintf(fp,"## Filters\n");
        i=0; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // FILTERS_RES
        i=36; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // FILTER_FORMAT
        i=33; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // SMOOTH_FILTERS
        i=34; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // SMOOTH_SIGMA
     
        fprintf(fp,"\n## Templates\n");
        i=1; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // TEMPLATES_FILE
        i=21; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // TEMPLATE_COMBOS
              fprintf(fp,"%-39s #         1 : one template at a time\n"," "); // TEMPLATE_COMBOS
              fprintf(fp,"%-39s #         2 : two templates, read allowed combinations from TEMPLATES_FILE\n"," "); // TEMPLATE_COMBOS
              fprintf(fp,"%-39s #        -2 : two templates, all permutations\n"," "); // TEMPLATE_COMBOS
              fprintf(fp,"%-39s # a <or> 99 : all templates simultaneously\n"," "); // TEMPLATE_COMBOS
        i=30; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // NMF_TOLERANCE
        i=2; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // WAVELENGTH_FILE
        i=3; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // TEMP_ERR_FILE
        i=22; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // TEMP_ERR_A2
        i=23; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // SYS_ERR
        i=24; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // APPLY_IGM
        i=49; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // LAF_FILE
        i=50; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // DLA_FILE
        i=45; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // SCALE_2175_BUMP

        i=15; fprintf(fp,"\n%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // DUMP_TEMP 
        i=16; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // USE_TEMP
        i=35; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // CACHE_FILE
        
        fprintf(fp,"\n## Input Files\n");
        i=4; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment);  // CATALOG_FILE
        i=44; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment);  // MAGNITUDES
        i=5; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // NOT_OBS_THRESHOLD
        i=6; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // N_MIN_COLORS

        fprintf(fp,"\n## Output Files\n");
        i=7; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // OUTPUT_DIRECTORY
        i=8; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // MAIN_OUTPUT_FILE
        i=31; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // PRINT_ERRORS
        i=39; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // CHI2_SCALE
        i=32; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // VERBOSE_LOG
        i=9; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // OBS_SED_FILE
        i=10; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // TEMP_SED_FILE
        i=11; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // POFZ_FILE
        i=40; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // BINARY_OUTPUT

        fprintf(fp,"\n## Redshift / Mag prior\n");
        i=26; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // APPLY_PRIOR
        i=27; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // PRIOR_FILE
        i=29; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // PRIOR_FILTER
        i=28; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // PRIOR_ABZP
        
        fprintf(fp,"\n## Redshift Grid\n");
        i=25; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // FIX_ZSPEC
        i=17; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // Z_MIN
        i=18; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // Z_MAX
        i=19; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // Z_STEP
        i=20; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // Z_STEP_TYPE

        fprintf(fp,"\n## Zeropoint Offsets\n");
        i=37; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // GET_ZP_OFFSETS
        i=38; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // ZP_OFFSET_TOL

        fprintf(fp,"\n## Rest-frame colors\n");
        i=41; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // REST_FILTERS
        i=47; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // RF_PADDING
        i=48; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // RF_ERRORS
        i=42; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // Z_COLUMN
        i=43; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // USE_ZSPEC_FOR_REST
        i=46; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // READ_ZBIN
        
        fprintf(fp,"\n## Cosmology\n");
        i=12; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // H0
        i=13; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // OMEGA_M
        i=14; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,params[i].default_value,params[i].comment); // OMEGA_L
      
        fclose(fp);
        exit(1);
    }

    while(fgets(buff,1024,fp)) {
        ++nline;
        argcount = 0;
        arg = strtok(buff,delim);
        i = NKEY;
        while (arg) {
            if (*arg=='#') break;
            if (argcount == 0) { // look for key match
                i = 0;
                while (strcmp(params[i].name,arg)&&i<NKEY) ++i;
                if (i<NKEY) {
                //	   printf("    %s:  ", params[i].name);
                    ++argcount;
                } else {
                    fprintf(stderr,
                        "getphotz: unrecognized text on line %d!\n",nline);
                    fprintf(stderr," ?? %s ??\n",arg);
                    // exit(1); //// Don't die if unrecognized keywords found
                }
            } else if (argcount==1) {
                ++argcount;
                if (params[i].defined) 
                    fprintf(stderr,"getphotz:  Warning!! %s redefined on line %d\n",
                                    params[i].name,nline);
                assign_par(i,arg);
                params[i].defined = 1;
                ++ndef;
            } else {
                fprintf(stderr,
                    "getphotz:  unrecognized text on line %d!\n",nline);
                fprintf(stderr," ?? %s ??\n",arg);
                // exit(1); //// Don't die if unrecognized keywords found
            }
            //printf("argcount %d\n", argcount);
            arg = strtok(NULL,delim);
        } 
    }
   
    printf("\n %d parameters defined in zphot.param.\n",ndef);
 //   if (ndef<NKEY) {
        printf(
          "\nDefault values will be assumed for the following parameters: \n");
        for (i=0;i<NKEY;++i) {
            if(!(params[i].defined)) {
                assign_par(i,params[i].default_value);
            }
        }
 //   }
    
}

//////// put a pointer to the file in the arguments to 
//////// the subroutine and dump verbose information from main/getphotz ?
void printparams_logfile(FILE *fp) 
{
    int i;
    //FILE *fp;
    //char outparfile[1024];
    
    //Parameter definitions
   struct param_data params[NKEY] = {
     // USE hyperz format as base
     //   {"AOVSED", 0, "NIL"},  // Vega SED
     {"FILTERS_RES", 0 , "FILTER.RES.latest","Filter transmission data"},            // 0 
     {"TEMPLATES_FILE", 0, "templates/eazy_v1.2_dusty.spectra.param","Template definition file"},       // 1 
     {"WAVELENGTH_FILE",0,"template/EAZY_v1.1_lines/lambda_v1.1.def","Wavelength grid definition file"},    // 2
     {"TEMP_ERR_FILE",0,"templates/TEMPLATE_ERROR.eazy_v1.0","Template error definition file"},     // 3 
     {"CATALOG_FILE", 0, "hdfn_fs99_eazy.cat","Catalog data file"},                      // 4
     {"NOT_OBS_THRESHOLD",0, "-90","Ignore flux point if <NOT_OBS_THRESH"},  // 5
     {"N_MIN_COLORS",0,"5","Require N_MIN_COLORS to fit"},                    // 6
     {"OUTPUT_DIRECTORY", 0, "OUTPUT","Directory to put output files in"},      // 7
     {"MAIN_OUTPUT_FILE", 0, "photz","Main output file, .zout"},                 // 8
     {"OBS_SED_FILE", 0, "n","Write out observed SED/object, .obs_sed"},      // 9
     {"TEMP_SED_FILE",0,"n","Write out best template fit/object, .temp_sed"}, // 10
     {"POFZ_FILE",0,"n","Write out Pofz/object, .pz"},                        // 11
     {"H0",0,"70.0","Hubble constant (km/s/Mpc)"},                                      // 12
     {"OMEGA_M",0,"0.3","Omega_matter"},                                      // 13
     {"OMEGA_L",0,"0.7","Omega_lambda"},                                      //14
     {"DUMP_TEMPLATE_CACHE",0,"n","Write binary template cache"},             // 15
     {"USE_TEMPLATE_CACHE",0,"n","Load in template cache"},                   // 16
     {"Z_MIN",0,"0.01","Minimum redshift"},                                   //17
     {"Z_MAX",0,"6.0","Maximum redshift"},                                    //  18  
     {"Z_STEP",0,"0.01","Redshift step size"},                                //19
     {"Z_STEP_TYPE",0,"1"," 0 = ZSTEP, 1 = Z_STEP*(1+z)"},                    //20
     {"TEMPLATE_COMBOS",0,"a","Template combination options: "},// 21
     {"TEMP_ERR_A2",0,"0.50","Template error amplitude"}, // 22
     {"SYS_ERR",0,"0.00","Systematic flux error (% of flux)"},                // 23
     {"APPLY_IGM",0,"y","Apply Madau 1995 IGM absorption"},                   // 24
     {"FIX_ZSPEC",0,"n","Fix redshift to catalog zspec"},                     // 25
     {"APPLY_PRIOR",0,"y","Apply apparent magnitude prior"},                          // 26
     {"PRIOR_FILE",0,"templates/prior_K_extend.dat","File containing prior grid"},       // 27
     {"PRIOR_ABZP",0,"25.0","AB zeropoint of fluxes in catalog.  Needed for calculating apparent mags!"},             // 28
     {"PRIOR_FILTER",0,"28","Filter from FILTER_RES corresponding to the columns in PRIOR_FILE"}, //29
     {"NMF_TOLERANCE",0,"1.e-4","Tolerance for non-negative combinations (TEMPLATE_COMBOS=a)"}, //30
     {"PRINT_ERRORS",0,"y","Print 68, 95 and 99% confidence intervals"}, //31
     {"VERBOSE_LOG",0,"y","Dump information from the run into [MAIN_OUTPUT_FILE].param"}, //32
     {"SMOOTH_FILTERS",0,"n","Smooth filter curves with Gaussian"}, //33
     {"SMOOTH_SIGMA",0,"100.","Gaussian sigma (in Angstroms) to smooth filters"}, //34
     {"CACHE_FILE",0,"photz.tempfilt","Template cache file (in OUTPUT_DIRECTORY)"}, //35
     {"FILTER_FORMAT",0,"1","Format of FILTERS_RES file -- 0: energy-  1: photon-counting detector"}, //36
     {"GET_ZP_OFFSETS",0,"n","Look for zphot.zeropoint file and compute zeropoint offsets"}, //37
     {"ZP_OFFSET_TOL",0,"1.e-4","Tolerance for iterative fit for zeropoint offsets [not implemented]"}, //38
     {"CHI2_SCALE",0,"1.0","Scale ML Chi-squared values to improve confidence intervals"}, //39
     {"BINARY_OUTPUT",0,"y","Save OBS_SED, TEMP_SED, PZ in binary format to read with e.g IDL"}, //40
     {"REST_FILTERS",0,"---","Comma-separated list of rest frame filters to compute"}, //41
     {"Z_COLUMN",0,"z_peak","Redshift to use for rest-frame color calculation (z_a, z_p, z_m1, z_m2, z_peak)"}, //42
     {"USE_ZSPEC_FOR_REST",0,"y","Use z_spec when available for rest-frame colors"}, //43
     {"MAGNITUDES",0,"n","Catalog photometry in magnitudes rather than f_nu fluxes"}, //44
     {"SCALE_2175_BUMP",0,"0.00","Scaling of 2175A bump.  Values 0.13 (0.27) absorb ~10 (20) % at peak."}, //45
     {"READ_ZBIN",0,"n","Get redshifts from OUTPUT_DIRECTORY/MAIN_OUTPUT_FILE.zbin rather than fitting them."}, //46
     {"RF_PADDING",0,"1000.","Padding (Ang) for choosing observed filters around specified rest-frame pair."}, //47
     {"RF_ERRORS",0,"n","Compute RF color errors from p(z)"}, //48
     {"LAF_FILE",0,"templates/LAFcoeff.txt","File containing the Lyman alpha forest data from Inoue"}, //49
     {"DLA_FILE",0,"templates/DLAcoeff.txt","File containing the damped Lyman absorber data from Inoue"} //50
   };


    ///// Print used parameters to MAIN_OUTPUT_FILE.param  
    /*strcpy(outparfile,OUTPUT_DIRECTORY);
    strcat(outparfile,"/");
    strcat(outparfile,MAIN_OUTPUT_FILE);
    strcat(outparfile,".param");
    fp = fopen(outparfile,"w");
    */
    
    fprintf(fp,"################   Run parameters (can feed this file back to EAZY)  ####################\n");
    fprintf(fp,"## Filters\n");
      i=0; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,FILTERS_RES,params[i].comment); // FILTERS_RES
      i=36; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,FILTER_FORMAT,params[i].comment); // FILTER_FORMAT
      i=33; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,SMOOTH_FILTERS,params[i].comment); // SMOOTH_FILTERS
      i=34; fprintf(fp,"%-20s %-18.2f # %s\n",params[i].name,SMOOTH_SIGMA,params[i].comment); // SMOOTH_SIGMA
     
    fprintf(fp,"\n## Templates\n");
      i=1; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,TEMPLATES_FILE,params[i].comment); // TEMPLATES_FILE
      i=21; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,TEMPLATE_COMBOS,params[i].comment); // TEMPLATE_COMBOS
      i=30; fprintf(fp,"%-20s %-18.2e # %s\n",params[i].name,NMF_TOLERANCE,params[i].comment); // NMF_TOLERANCE
      i=2; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,WAVELENGTH_FILE,params[i].comment); // WAVELENGTH_FILE
      i=3; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,TEMP_ERR_FILE,params[i].comment); // TEMP_ERR_FILE
      i=22; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,TEMP_ERR_A2,params[i].comment); // TEMP_ERR_A2
      i=23; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,SYS_ERR,params[i].comment); // SYS_ERR
      i=24; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,APPLY_IGM,params[i].comment); // APPLY_IGM
      i=49; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,LAF_FILE,params[i].comment); // LAF_FILE
      i=50; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,DLA_FILE,params[i].comment); // DLA_FILE
      i=45; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,SCALE_2175_BUMP,params[i].comment); // SCALE_2175_BUMP

      i=15; fprintf(fp,"\n%-20s %-18i # %s\n",params[i].name,DUMP_TEMPLATE_CACHE,params[i].comment); // DUMP_TEMP 
      i=16; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,USE_TEMPLATE_CACHE,params[i].comment); // USE_TEMP
      i=35; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,CACHE_FILE,params[i].comment); // CACHE_FILE
    
    fprintf(fp,"\n## Input Files\n");
      i=4; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,CATALOG_FILE,params[i].comment);  // CATALOG_FILE
      i=44; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,MAGNITUDES,params[i].comment);  // MAGNITUDES
      i=5; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,NOT_OBS_THRESHOLD,params[i].comment); // NOT_OBS_THRESHOLD
      i=6; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,N_MIN_COLORS,params[i].comment); // N_MIN_COLORS

    fprintf(fp,"\n## Output Files\n");
      i=7; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,OUTPUT_DIRECTORY,params[i].comment); // OUTPUT_DIRECTORY
      i=8; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,MAIN_OUTPUT_FILE,params[i].comment); // MAIN_OUTPUT_FILE
      i=31; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,PRINT_ERRORS,params[i].comment); // PRINT_ERRORS
      i=39; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,CHI2_SCALE,params[i].comment); // CHI2_SCALE
      i=32; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,VERBOSE_LOG,params[i].comment); // VERBOSE_LOG
      i=9; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,OBS_SED_FILE,params[i].comment); // OBS_SED_FILE
      i=10; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,TEMP_SED_FILE,params[i].comment); // TEMP_SED_FILE
      i=11; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,POFZ_FILE,params[i].comment); // POFZ_FILE
      i=40; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,BINARY_OUTPUT,params[i].comment); // BINARY_OUTPUT

    fprintf(fp,"\n## Redshift / Mag prior\n");
      i=26; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,APPLY_PRIOR,params[i].comment); // APPLY_PRIOR
      i=27; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,PRIOR_FILE,params[i].comment); // PRIOR_FILE
      i=29; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,PRIOR_FILTER,params[i].comment); // PRIOR_FILTER
      i=28; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,PRIOR_ABZP,params[i].comment); // PRIOR_ABZP
    
    fprintf(fp,"\n## Redshift Grid\n");
      i=25; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,FIX_ZSPEC,params[i].comment); // FIX_ZSPEC
      i=17; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,Z_MIN,params[i].comment); // Z_MIN
      i=18; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,Z_MAX,params[i].comment); // Z_MAX
      i=19; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,Z_STEP,params[i].comment); // Z_STEP
      i=20; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,Z_STEP_TYPE,params[i].comment); // Z_STEP_TYPE

    fprintf(fp,"\n## Zeropoint Offsets\n");
      i=37; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,GET_ZP_OFFSETS,params[i].comment); // GET_ZP_OFFSETS
      i=38; fprintf(fp,"%-20s %-18.3e # %s\n",params[i].name,ZP_OFFSET_TOL,params[i].comment); // ZP_OFFSET_TOL

    fprintf(fp,"\n## Rest-frame colors\n");
      i=41; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,REST_FILTERS,params[i].comment); // REST_FILTERS
      i=47; fprintf(fp,"%-20s %-18.0f # %s\n",params[i].name,RF_PADDING,params[i].comment); // RF_PADDING
      i=48; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,RF_ERRORS,params[i].comment); // RF_ERRORS
      i=42; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,Z_COLUMN,params[i].comment); // Z_COLUMN
      i=43; fprintf(fp,"%-20s %-18i # %s\n",params[i].name,USE_ZSPEC_FOR_REST,params[i].comment); // USE_ZSPEC_FOR_REST
      i=46; fprintf(fp,"%-20s %-18s # %s\n",params[i].name,ZBIN_FILE,params[i].comment); // READ_ZBIN

    fprintf(fp,"\n## Cosmology\n");
      i=12; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,H0,params[i].comment); // H0
      i=13; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,OMEGA_M,params[i].comment); // OMEGA_M
      i=14; fprintf(fp,"%-20s %-18.3f # %s\n",params[i].name,OMEGA_L,params[i].comment); // OMEGA_L
      
    fprintf(fp,"#\n####################################\n#\n");
    
    //fclose(fp);
   
}
