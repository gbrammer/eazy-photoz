#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>

extern double chiscale;

//// Templates 
extern int NTEMP;  // to be defined from template file
extern char TEMPLATES_FILE[1024];

//// Maximum number of lines in template file
#define NTEMPLMAX (long)1.e5 

extern int NTEMPL;

//// Redshift grid
extern int NZ;
extern double Z_STEP;
extern int Z_STEP_TYPE; 
extern double Z_MIN;
extern double Z_MAX; 

//// For output, set to standard 32bit integer
int32_t nusefilt32, NTEMP32, NZ32, nobj32, NZ32, *izsave32, NTEMPL32;

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
extern char WAVELENGTH_FILE[1024];
extern char TEMP_ERR_FILE[1024];

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
extern FILE *fplog;

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

extern filt_data filt_thru;
extern filt_data **pusefilt;

int iz,*idtemp1,*idtemp2a,*idtemp2b,**idtempall,ntemp_all,ngoodfilters,ncols;
double *pz1,*pz2,*pzall,*pzuse,*pzout,*atemp1,*atemp2a,*atemp2b,**coeffs;

///// Compute zeropoint offsets
int GET_ZP_OFFSETS;
double ZP_OFFSET_TOL;
int *zpfilterid, zpcontinue, *filt_to_col;
double *zpfactor;

extern long nfilter;
extern long *usefilt;
extern long nusefilt;
extern long totfiltid;
extern long totcolumn;
extern double *totflux;
extern long nobj, *nobjfilt;
extern double **fnu, **efnu, *zspec, *lambdac, *ra, *dec;
extern char **objid;
extern int idsize;

//// Fluxes of redshifted templates integrated through filters
extern double ***tempfilt;
extern double *templ,**tempf,*tempage;
extern double  *temp_err_a, **temp_errf;
extern double  temp_err[NTEMPLMAX], temp_scale[NTEMPLMAX];
extern double *tnorm;
extern int **temp_combine;
extern double *ztry,*zspec;

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
                         double *pzall, int **idtempall, double **coeffs, int *ntemp_all);

extern double *dasum,*dbsum;
void getigmfactors (double ztarg, double *daz, double *dbz);

int zeropoints();

void pz_output(char *filename, long iobj);
void temp_sed_output(char *filename, long iobj, long izbest, long izprior);
void obs_sed_output(char *filename, long iobj, long izbest, long izprior);

void get_errors(double *l68, double *u68, double *l95, double *u95, double *l99, double *u99);

double cosmotl(double z, double lambda);

//// Prior
extern char PRIOR_FILE[1024];
extern int APPLY_PRIOR;
extern double PRIOR_ABZP;
extern int PRIOR_FILTER;
extern long Kcolumn;
extern double **priorkz,*klim;

void prior_findmin(double *probz, long NZ, long *mins, long *Nmin);
void get_prior_file_size(long *NZ_prior, long *NK_prior);
void read_prior_file(double **priorkz, double *priorzval, double *klim,
        long nzp, long nkp);
void interpolate_prior(double **priorkz_first, double **priorkz_out, double *priorzval,
                       long nzp, long nkp);
void apply_prior(double **priorkz, double *priorzval, double *klim, 
        long nzp, long nkp,
        double *probz, double kflux, long *bestz, double *pzout);

