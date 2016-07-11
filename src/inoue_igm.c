#include "defs.h"

// #include <stdio.h>
// #include <math.h>
// #include <string.h>
// #include <stdlib.h>
// #include <time.h>
// #include <ctype.h>
// #include <unistd.h>

/* c translation of the fortran IGM code from Inoue et al. (2014):

http://adsabs.harvard.edu/abs/2014arXiv1402.0677I
http://www.las.osaka-sandai.ac.jp/~inoue/ANAIGM/ANAIGM.tar.gz

*/

void read_Inoue_coeffs() {
    /* Read the coefficients files */
    
    FILE *tlaf, *tdla;
    char *arg,buff[1024];
    char delim[]=" \n\t,\r";
    int j, ix;
    double value;
    
    NA=39;
    lam1 = malloc(sizeof(double)*NA);    
    ALAF1 = malloc(sizeof(double)*NA);    
    ALAF2 = malloc(sizeof(double)*NA);    
    ALAF3 = malloc(sizeof(double)*NA);    
    ADLA1 = malloc(sizeof(double)*NA);    
    ADLA2 = malloc(sizeof(double)*NA);    

    if (!(tlaf = fopen(LAF_FILE,"r"))) {
        fprintf(stderr,"  Can't open the LAF file!  Bad name?\n");
        exit(1);
    }

    if (!(tdla = fopen(DLA_FILE,"r"))) {
        fprintf(stderr,"  Can't open the DLA file!  Bad name?\n");
        exit(1);
    }
    
    j=0;
    while(fgets(buff,1024,tlaf)) {
        sscanf(buff,"%d %lf %lf %lf %lf\n", &ix, &lam1[j], &ALAF1[j], &ALAF2[j], &ALAF3[j]);
        ++j;
    }
    
    j=0;
    while(fgets(buff,1024,tdla)) {
        sscanf(buff,"%d %lf %lf %lf\n", &ix, &lam1[j], &ADLA1[j], &ADLA2[j]);
        ++j;
    }
    
}

double tLSLAF(double zS, double lobs) {
    /* Lyman series, LAF */
    
    double tLSLAF_value, z1LAF, z2LAF;
    int j;
    
    z1LAF = 1.2;
    z2LAF = 4.7;
    
    tLSLAF_value = 0.0; // ! Initialize
    for (j=0;j<NA;++j) {
       if ((lobs < lam1[j]*(1.0+zS)) & (lobs > lam1[j])) {
          if (lobs < lam1[j]*(1.0+z1LAF)) {
              tLSLAF_value += ALAF1[j] * pow(lobs / lam1[j], 1.2);
          } else if (lobs < lam1[j]*(1.0+z2LAF)) {
              tLSLAF_value += ALAF2[j] * pow(lobs / lam1[j], 3.7);
          } else 
              tLSLAF_value += ALAF3[j] * pow(lobs / lam1[j], 5.5);
       }
    }
    
    return tLSLAF_value;
}

double tLSDLA(double zS, double lobs) {
    /* Lyman series, DLA */
    
    double tLSDLA_value, z1DLA;
    int j;
    z1DLA = 2.0;
    
    tLSDLA_value = 0.0; // ! Initialize
    for (j=0;j<NA;++j) {
       if ((lobs < lam1[j]*(1.0+zS)) & (lobs > lam1[j])) {
          if (lobs < lam1[j]*(1.0+z1DLA)) {
              tLSDLA_value += ADLA1[j] * pow(lobs / lam1[j], 2.0);
          } else 
              tLSDLA_value += ADLA2[j] * pow(lobs / lam1[j], 3.0);
       }
    }
    
    return tLSDLA_value;
}

double tLCDLA(double zS, double lobs) {
    /* Lyman continuum, DLA */
    double tLCDLA_value, z1DLA, lamL;
    z1DLA = 2.0;
    lamL = 911.8;
    
    if (lobs > lamL*(1.0+zS)) {
        tLCDLA_value = 0.0;
    } else if (zS < z1DLA) {
       tLCDLA_value = 0.2113 * pow(1.0+zS, 2) - 0.07661 * pow(1.0+zS, 2.3) * pow(lobs/lamL, (-3e-1)) 
           - 0.1347 * pow(lobs/lamL, 2);
    } else {
       if (lobs > lamL*(1.0+z1DLA)) {
          tLCDLA_value = 0.04696 * pow(1.0+zS, 3) 
               - 0.01779 * pow(1.0+zS, 3.3) * pow(lobs/lamL, (-3e-1)) 
               - 0.02916 * pow(lobs/lamL, 3);
       } else {
          tLCDLA_value = 0.6340 + 0.04696 * pow(1.0+zS, 3) 
               - 0.01779 * pow(1.0+zS, 3.3) * pow(lobs/lamL, (-3e-1)) 
               - 0.1347 * pow(lobs/lamL, 2) - 0.2905 * pow(lobs/lamL, (-3e-1));
       }
    }
    
    return tLCDLA_value;
    
}

double tLCLAF(double zS, double lobs) {
    /* Lyman continuum, Ly-a forest */
    double tLCLAF_value, z1LAF, z2LAF, lamL;
    
    z1LAF = 1.2;
    z2LAF = 4.7;
    lamL = 911.8;
    
    if (lobs > lamL*(1.0+zS)) {
        tLCLAF_value = 0.0;
    } else if (zS < z1LAF) {
        tLCLAF_value = 0.3248 * (pow(lobs/lamL, 1.2) - pow(1.0+zS, -9e-1) * pow(lobs/lamL, 2.1));
    } else if (zS < z2LAF) {
       if (lobs > lamL*(1.0+z1LAF)) {
           tLCLAF_value = 2.545e-2 * (pow(1.0+zS, 1.6) * pow(lobs/lamL, 2.1) - pow(lobs/lamL, 3.7));
       } else {
          tLCLAF_value = 2.545e-2 * pow(1.0+zS, 1.6) * pow(lobs/lamL, 2.1) 
              + 0.3248 * pow(lobs/lamL, 1.2) - 0.2496 * pow(lobs/lamL, 2.1);
       }
    } else {
       if (lobs > lamL*(1.0+z2LAF)) {
           tLCLAF_value = 5.221e-4 * (pow(1.0+zS, 3.4) * pow(lobs/lamL, 2.1) - pow(lobs/lamL, 5.5));
       } else if (lobs > lamL*(1.0+z1LAF)) {
          tLCLAF_value = 5.221e-4 * pow(1.0+zS, 3.4) * pow(lobs/lamL, 2.1) 
              + 0.2182 * pow(lobs/lamL, 2.1) - 2.545e-2 * pow(lobs/lamL, 3.7);
       } else {
          tLCLAF_value = 5.221e-4 * pow(1.0+zS, 3.4) * pow(lobs/lamL, 2.1) 
              + 0.3248 * pow(lobs/lamL, 1.2) - 3.140e-2 * pow(lobs/lamL, 2.1);
       }
    }
    
    return tLCLAF_value;
    
}

// int main() {
//     double zS = 5.;
//     double lrest, lobs, tau;
//     int i;
//     
//     read_Inoue_coeffs();
//     // printf("%f  %e %e %e    %e %e\n", lam1[10], ALAF1[10], ALAF2[10], ALAF3[10], ADLA1[10], ADLA2[10]);
//     
//     for (i=0; i<=1300; ++i) {
//         lrest = 100. + i;
//         lobs = lrest*(1.+zS);
//         //tau = tLSLAF(zS, lobs);
//         //tau = tLCLAF(zS, lobs); 
//         //tau = tLSDLA(zS, lobs);
//         tau = tLSLAF(zS, lobs) + tLCLAF(zS, lobs) + tLSDLA(zS, lobs) + tLCDLA(zS, lobs); 
//         printf("%7.1f  %13.4le\n", lrest, exp(-tau));
//     }
// }