#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>

extern double chiscale;

extern char EAZY_VERSION[];

extern char zphot_param_file[64], zphot_translate_file[64], zphot_zeropoint_file[64];

//// Templates 
extern int NTEMP;
extern int NTEMP_REST;  // to be defined from template file
extern char TEMPLATES_FILE[1024];
extern char RF_TEMPLATES_FILE[1024];
extern char ZBIN_FILE[1024];

//// Maximum number of lines in template file
#define NTEMPLMAX (long)1.e5 

extern int NTEMPL;

//// Redshift grid
extern int NZ, NZCOMOVE;
extern double Z_STEP;
extern int Z_STEP_TYPE; 
extern double Z_MIN;
extern double Z_MAX; 

//// For output, set to standard 32bit integer
int32_t nusefilt32, NTEMP32, nobj32, NZ32, NK32,*izsave32, NTEMPL32;

//// Filters
extern char FILTERS_RES[1024];
extern int SMOOTH_FILTERS;
extern double SMOOTH_SIGMA;
extern int FILTER_FORMAT;

typedef struct 
{ long ndat;
  char filtname[180];
  double *plambda;
  double *pthrough;
  void *pnext;
} filt_data;

//// Input data
extern char CATALOG_FILE[1024];
extern double NOT_OBS_THRESHOLD;
extern int N_MIN_COLORS;
extern int MAGNITUDES;
extern char WAVELENGTH_FILE[1024];
extern char TEMP_ERR_FILE[1024];
extern char LAF_FILE[1024];
extern char DLA_FILE[1024];

//// Output
extern char OUTPUT_DIRECTORY[1024];
extern char MAIN_OUTPUT_FILE[1024];
extern int OBS_SED_FILE;
extern int TEMP_SED_FILE;
extern int POFZ_FILE;
extern int BIGDUMP_FILE;
extern int PRINT_ERRORS;
extern double CHI2_SCALE;
extern int BINARY_OUTPUT;

extern int VERBOSE_LOG;
extern FILE *fplog, *fprf;

//// Cosmology
extern double H0;
extern double OMEGA_M;
extern double OMEGA_L;

//// Template caching
extern int DUMP_TEMPLATE_CACHE;
extern int USE_TEMPLATE_CACHE;
extern char CACHE_FILE[1024];

//// Template fitting
extern int TEMPLATE_COMBOS;
extern double NMF_TOLERANCE;
extern double SYS_ERR;
extern double TEMP_ERR_A2;
extern int APPLY_IGM;  
extern int FIX_ZSPEC;  
extern double SCALE_2175_BUMP;

extern filt_data filt_thru;
extern filt_data **pusefilt;

extern int iz,*idtemp1,*idtemp2a,*idtemp2b,**idtempall,ntemp_all,ngoodfilters,ncols;
extern double *pz1,*pz2,*pzall,*chi2fit,*pzout,*atemp1,*atemp2a,*atemp2b,**coeffs;

///// Compute zeropoint offsets
extern int GET_ZP_OFFSETS;
extern double ZP_OFFSET_TOL;
extern int *zpfilterid, zpcontinue, *filt_to_col;
extern double *zpfactor;

///// Rest-frame filters
extern char REST_FILTERS[1024],Z_COLUMN[512];
extern int USE_ZSPEC_FOR_REST;
extern int READ_ZBIN, ZBIN_OPENED;
extern double RF_PADDING;
extern int RF_ERRORS;

extern int *filt_defnum;
extern long nfilter;
extern long *okfilt;
extern long nusefilt,nrestfilt;
extern long totfiltid;
extern long totcolumn;
extern double *totflux;
extern long nobj, *nobjfilt;
extern double **fnu, **efnu, *zspec, *lambdac,*lambdac_sort, *ra, *dec;
extern char **objid;
extern int idsize;
extern long *lc_sort_idx;

//// Fluxes of redshifted templates integrated through filters
extern double ***tempfilt;
extern double ***tempfilt_rest;
extern double *templ,**tempf,**tempf_rest,*tempage;
extern double  *temp_err_a, **temp_errf;
extern double  temp_err[NTEMPLMAX], temp_scale[NTEMPLMAX];
extern double *tnorm;
extern int **temp_combine;
extern double *ztry,*zspec,*comovingdist,*zcomove;

extern double *z_a, *z_p, *z_m1, *z_m2, *z_peak;

//// IGM parametrization from Inoue et al. 2014
extern double *lam1, *ALAF1, *ALAF2, *ALAF3, *ADLA1, *ADLA2;
extern int NA;
void read_Inoue_coeffs();
double tLSLAF();
double tLCLAF();
double tLSDLA();
double tLCDLA();

void getparams();
void printparams_logfile(); 

//// Spline interpolation
int initialise_spline(double *x,double *y,double *k,long n, double q2b,double q2e);
int interpolate_spline(double z,double *val,double *x,double *y,double *k,long n);

//// Linear interpolation
int interpol(double xval,double *out, double *x, double *y, long n, long *istart);

//// makegrid.c
void init();

int getphotz(long iobj, double *pz1, int *idtemp1, double *atemp1,
                         double *pz2, int *idtemp2a, int *idtemp2b, double *atemp2a, double *atemp2b,
                         double *pzall, int **idtempall, double **coeffs, int *ntemp_all, long fixidx);

extern double *dasum,*dbsum;
void getigmfactors (double ztarg, double *daz, double *dbz);

int zeropoints();

void pz_output(char *filename, long iobj);
void temp_sed_output(char *filename, long iobj, long izbest, long izprior);
void obs_sed_output(char *filename, long iobj, long izbest, long izprior);

void get_errors(double *l68, double *u68, double *l95, double *u95, double *l99, double *u99);

double cosmotl(double z, double lambda);

int compare_for_sort (const void *X, const void *Y);

//// Prior
extern char PRIOR_FILE[1024];
extern int APPLY_PRIOR;
extern double PRIOR_ABZP;
extern int PRIOR_FILTER, PRIOR_FILTER_IDX;
extern long Kcolumn;
extern int32_t *klim_idx;
extern double **priorkz,*klim;

void prior_findmin(double *probz, long NZ, long *mins, long *Nmin);
void get_prior_file_size(long *NZ_prior, long *NK_prior);
void read_prior_file(double **priorkz, double *priorzval, double *klim,
        long nzp, long nkp);
void interpolate_prior(double **priorkz_first, double **priorkz_out, double *priorzval,
                       long nzp, long nkp);
void apply_prior(double **priorkz, long nzp, long nkp, long first_best,
        double *probz, double kflux, long *bestz, double *pzout);

