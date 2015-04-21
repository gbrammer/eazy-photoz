#include "defs.h"

int iz;

double cosmotl(double z, double lambda)
{ // calculates time as a function of z in flat lambda universes
  // (omega_l + omega_m = 1, omega_l>0)
  //  arguments: $1=z, $2=omega_lambda

    if (z <= 0) return (1);
    
    return (asinh(pow((1.0+z),-1.50)*sqrt(lambda/(1.0-lambda)))
          /asinh(sqrt(lambda/(1.0-lambda))));
}
  
int getphotz(long iobj, double *pz1, int *idtemp1, double *atemp1, 
                         double *pz2, int *idtemp2a, int *idtemp2b, double *atemp2a, double *atemp2b,
                         double *pzall, int **idtempall, double **coeffs, int *ntemp_all, long fixidx)
{
    void nonneg_fact(double **amatrix, double *bvector, double *coeffs_z, long NTEMP,
                 double toler);

    long i,j,k,NTEMP_I;
    double chisum, age_univ, te,ztest;
    static int ifirst=1, *oktemp;
    static double tempsumfbest,*sigi2;
    static double **amatrix,*bvector,*coeffs_z;
    int itemp1,itemp2;
    double a11,a12,a21,a22, b1,b2, ai11,ai12,ai21,ai22,det;
    double c1,c2,  tempchi, d1,s1, a1;
    double pzmin,pzmin2;
    int izspec,NZ_use;
  
//    char templogfile[1024];
//    FILE *tlog;
    
    c1=0.;
    c2=0.;  
    if (ifirst) {
    
        ////// Take out and make it a global variable
        // okfilt = malloc(sizeof(ifirst)*nusefilt);
        // if(okfilt == NULL) {
        //     fprintf(stderr, "okfilt failed, out of memory\n");
        //     exit(1);
        // }
        /* // Use in place
        tempsumfbest = malloc(sizeof(c1)*nusefilt);
        if(tempsumfbest == NULL) {
            fprintf(stderr, "tempsumfbest failed, out of memory\n");
            exit(1);
        }
        */
        oktemp = malloc(sizeof(int)*NTEMP);
        if(oktemp == NULL) {
            fprintf(stderr, "oktemp failed, out of memory\n");
            exit(1);
        }
        sigi2 = malloc(sizeof(double)*nusefilt);
        if(sigi2 == NULL) {
            fprintf(stderr, "sigi2 failed, out of memory\n");
            exit(1);
        }

        if (TEMPLATE_COMBOS == 99) {
            amatrix = malloc(sizeof(double *)*NTEMP);
            for (i=0;i<NTEMP;++i) 
                amatrix[i] = malloc(sizeof(double)*NTEMP);
            if(amatrix == NULL) {
                fprintf(stderr,"amatrix failed, out of memory\n");
                exit(1);
            }
    
            bvector = malloc(sizeof(double)*NTEMP);
            if(bvector == NULL) {
                fprintf(stderr,"bvector failed, out of memory\n");
                exit(1);
            }
    
            coeffs_z = malloc(sizeof(double)*NTEMP);
            if(coeffs_z == NULL) {
                fprintf(stderr,"coeffs_z failed, out of memory\n");
                exit(1);
            }
        }
        ifirst = 0;
    }  

    pzmin = 1e30;
    pzmin2 = 1e30;

    izspec=0;
    if (FIX_ZSPEC || fixidx >= 0) {
        //while ( (ztry[izspec] < zspec[iobj]) && (izspec < (NZ-1)) ) ++izspec;
 
        if (fixidx < 0) {
            ztest=100;
            for (iz=0;iz<NZ;++iz) if (fabs(ztry[iz]-zspec[iobj]) < ztest) {
                ztest = fabs(ztry[iz]-zspec[iobj]);
                izspec=iz;
            }
        } else izspec=fixidx; 
        NZ_use = izspec+1;
        
    } else NZ_use = NZ;
        
    if (TEMPLATE_COMBOS <= 2) {
      for (iz=izspec;iz<NZ_use;++iz) {

        age_univ = cosmotl(ztry[iz],OMEGA_L)*3e19/H0/3.14e7/1e9;
        //  printf("cosmo time: %lf\n",age_univ);

        //// choose which fitting templates are allowed    
        //// implement age constraint for now
        for (itemp1=0;itemp1<NTEMP;++itemp1) {
            oktemp[itemp1] = 1;
            ///////   Added line to allow fixing fitting to template
            // if (tempage[itemp1]<age_univ && (tuse[iobj]==itemp1 || tuse[iobj]==-1)) 
            if (tempage[itemp1]<age_univ)
                oktemp[itemp1]=1;
            else
                oktemp[itemp1]=0;
        }        
         
        pz1[iz] = 1e30;
        for (itemp1=0; itemp1<NTEMP; ++itemp1) {
            
            //// Check if filters fall off templates
            chisum=100;
            for (i=0;i<nusefilt;++i) {
                if (tempfilt[iz][itemp1][i] < -90.0) {
                    chisum = -100;
                    //fprintf(stderr,"HERE skipping template %d because of filt %ld at z=%lf\n",itemp1,i,ztry[iz]);
                }
            }
            
            if (oktemp[itemp1] && chisum > 0) {
                                    
                //  first compute min. chi1 for single template;
                //  get best template normalization
                
                d1 = 0;
                s1 = 0;
                // printf("first loop\n");
                for (i=0;i<nusefilt;++i) {
                    if (okfilt[i]) {
                        te = temp_err_a[itemp1]*temp_errf[iz][i]*TEMP_ERR_A2;
                        sigi2[i] = (SYS_ERR*SYS_ERR+te*te)
                            *fnu[iobj][i]*fnu[iobj][i];
                        sigi2[i] += efnu[iobj][i]*efnu[iobj][i];
                        d1 += tempfilt[iz][itemp1][i]*fnu[iobj][i]/sigi2[i];
                        s1 += tempfilt[iz][itemp1][i]*tempfilt[iz][itemp1][i]/sigi2[i];
                    }
                }
                a1 = d1/s1;
                chisum = 0.0;
                for (i=0;i<nusefilt;++i) {
                    if (okfilt[i]) {
                        tempsumfbest = a1*tempfilt[iz][itemp1][i]; 
                        tempchi = (tempsumfbest-fnu[iobj][i]);
                        tempchi = tempchi*tempchi/sigi2[i];
                        chisum += tempchi;
                    }
                }
                    
                /////// pass out results
                if (chisum < pz1[iz]) {
                    pz1[iz] = chisum;
                    idtemp1[iz] = itemp1;
                    atemp1[iz] = a1;
                }
            } 
        }
		for (itemp1=0; itemp1<NTEMP; ++itemp1) coeffs[iz][itemp1] = 0.;
		coeffs[iz][idtemp1[iz]] = atemp1[iz];

        ////// Two templates
        if (TEMPLATE_COMBOS == 2) {
            
            ///// Default two template fit is the best
            ///// *single* template fit.  Handles the case
            ///// where only negative combinations of two 
            ///// templates improves the fit.
            pz2[iz] = pz1[iz];
            idtemp2a[iz] = idtemp1[iz];
            idtemp2b[iz] = idtemp1[iz];
            atemp2a[iz] = atemp1[iz];
            atemp2b[iz] = 0.;
            
            for (itemp1=0; itemp1<NTEMP; ++itemp1){
                if (oktemp[itemp1]) {
                  for(itemp2=itemp1+1;itemp2<NTEMP; ++itemp2){

                    //// Check if filters fall off templates
                    chisum=100;
                    for (i=0;i<nusefilt;++i) {
                        if (tempfilt[iz][itemp1][i] < -90.0 || tempfilt[iz][itemp2][i] < -90.0) {
                            chisum = -100;
                        }
                    }

                    if (oktemp[itemp2]&&temp_combine[itemp1][itemp2] && chisum > 0) {
                        // printf("itemp 1 2 %d %d\n",itemp1,itemp2);
                        a11 =0;
                        a12 = 0;
                        a21 = 0;
                        a22 = 0;
                        b1 = 0;
                        b2 = 0;
                        chisum = 0.0;
                        for (i=0; i<nusefilt; ++i)
                          if (okfilt[i]) {
                            //printf("itemp1 itemp2 %d %d\n",itemp1,itemp2);
                            te = TEMP_ERR_A2*temp_errf[iz][i];
                            sigi2[i] = (SYS_ERR*SYS_ERR+te*te)
                                *fnu[iobj][i]*fnu[iobj][i];
                            sigi2[i] += efnu[iobj][i]*efnu[iobj][i];
	   
                            a11 += tempfilt[iz][itemp1][i]*tempfilt[iz][itemp1][i]/sigi2[i];
                            a12 += tempfilt[iz][itemp2][i]*tempfilt[iz][itemp1][i]/sigi2[i];
                            a21 += tempfilt[iz][itemp1][i]*tempfilt[iz][itemp2][i]/sigi2[i];
                            a22 += tempfilt[iz][itemp2][i]*tempfilt[iz][itemp2][i]/sigi2[i];
                            b1 += fnu[iobj][i]*tempfilt[iz][itemp1][i]/sigi2[i];
                            b2 += fnu[iobj][i]*tempfilt[iz][itemp2][i]/sigi2[i];
                        }
                        //// now invert "a" matrix;
                        det = a11*a22 - a21*a12;
                        ai11 = a22/det;
                        ai12 = -a21/det;
                        ai21 = -a12/det;
                        ai22 = a11/det;

                        // printf("determinant %lf\n",det);
                        if (det == 0 ) {
                            // printf("huh, how did I get a singular matrix!? \n");
                            return(1);
                        }
                        c1 = ai11*b1 + ai12*b2;
                        c2 = ai21*b1 + ai22*b2;
                        chisum = 0.0;
                        for (i=0;i<nusefilt;++i) 
                          if (okfilt[i]) {
                            tempsumfbest = c1*tempfilt[iz][itemp1][i] 
                                              + c2*tempfilt[iz][itemp2][i];
                            tempchi = (tempsumfbest-fnu[iobj][i]);
                            tempchi = tempchi*tempchi/sigi2[i];
                            chisum += tempchi;
                        }
                    }
                    //printf("and the chi^2 is? %lf\n",chisum);
                    ////// Only accept fit IF both coefficients are non-negative
                    if (chisum < pz2[iz] && c1>0 && c2>0) {
                        pz2[iz] = chisum;
                        idtemp2a[iz] = itemp1;
                        idtemp2b[iz] = itemp2;
                        atemp2a[iz] = c1;
                        atemp2b[iz] = c2;
                    }	    
                  }
	            }	            
            }
        	for (itemp1=0; itemp1<NTEMP; ++itemp1) coeffs[iz][itemp1] = 0.;   
        	coeffs[iz][idtemp2a[iz]] = atemp2a[iz];
        	coeffs[iz][idtemp2b[iz]] = atemp2b[iz];
        }
          
      }  //// For (iz...)    
    } else {
    ////// Fit linear combinations of *all* templates
    ////// in SPECTRA_FILE using Sha, Saul & Lee algorithm

      for (iz=izspec;iz<NZ_use;++iz) {

        //// Check if filters fall off templates
        chisum=100;
        for (itemp1=0;itemp1<NTEMP;++itemp1) {
            for (i=0;i<nusefilt;++i) {
                if (tempfilt[iz][itemp1][i] < -90.0) {
                    chisum = -100;
                }
            }
        }
        if (chisum < 0) {
            pzall[iz] = 9.99e30;
            continue;
        }
        
        age_univ = cosmotl(ztry[iz],OMEGA_L)*3e19/H0/3.14e7/1e9;
        //    printf("cosmo time: %lf\n",age_univ);

        // choose which fitting templates are allowed    
        // implement age constraint for now
        NTEMP_I = 0.;
        for (itemp1=0;itemp1<NTEMP;++itemp1) {
            ///////   Added line to allow fixing fitting to template
            // if (tempage[itemp1]<age_univ && (tuse[iobj]==itemp1 || tuse[iobj]==-1)) 
            if (tempage[itemp1]<age_univ) {
                  oktemp[NTEMP_I]=itemp1;
                  ++NTEMP_I;
            }
        }

        // fprintf(tlog,"%7.3f ",ztry[iz]);
        
        /////  Set up sigma^2, now only use filters with detections
        for (i=0;i<nusefilt;++i) if (okfilt[i]) {
            te = temp_errf[iz][i]*TEMP_ERR_A2;  ////// Force individual temp_err_a = 1
            sigi2[i] = (SYS_ERR*SYS_ERR+te*te)
                             *fnu[iobj][i]*fnu[iobj][i];
            sigi2[i] += efnu[iobj][i]*efnu[iobj][i];
            // fprintf(tlog,"%13.5e %16.5e ",fnu[iobj][i],sigi2[i]);
        }
         
        ///// "B" vector for template normalization
        for (i=0; i<NTEMP_I; ++i) {
            bvector[i] = 0.;

            for (j=0;j<nusefilt;++j) 
              // should be if (okfilt[j]) ??
              if (fnu[iobj][j] > 0 && efnu[iobj][j] > 0 && okfilt[j])
                bvector[i]-=fnu[iobj][j]*tempfilt[iz][oktemp[i]][j]/sigi2[j];
                    
             //printf("bvector %ld %e\n",itemp1,bvector[itemp1]);
        }  
      
        ///// "A" matrix for template normalization    
        for (i=0;i<NTEMP_I;++i) {
            for (j=0;j<NTEMP_I;++j) {
                amatrix[i][j] = 0.;
                for (k=0;k<nusefilt;++k) {
                    // if (fnu[iobj][k] > NOT_OBS_THRESHOLD && efnu[iobj][k] > 0) 
                    if (okfilt[k])
                        amatrix[i][j]+=tempfilt[iz][oktemp[i]][k]*
                                       tempfilt[iz][oktemp[j]][k]/sigi2[k];                    
                }    
                ///printf("am[%ld,%ld]=%e\n",i,j,amatrix[i][j]);
                  
            }
            //coeffs_z[i] = 1.;
            if (bvector[i] < 0) {
                coeffs_z[i] = 1.;
            } else {
                coeffs_z[i] = 0.;
            }
        }
      
        ///////////// Calculate template normlizations /////////////////    
        nonneg_fact(amatrix, bvector, coeffs_z, NTEMP_I, NMF_TOLERANCE);
        ////////////////////////////////////////////////////////////////     
        
        ////// Compute Chi squared                
        chisum = 0.0;
        for (i=0;i<nusefilt;++i) {
            if (okfilt[i]) {
                tempsumfbest = 0.;
                for (j=0;j<NTEMP_I;++j) tempsumfbest+=coeffs_z[j]*tempfilt[iz][oktemp[j]][i]; 

                // fprintf(tlog,"%13.5e ",tempsumfbest);
      		
                tempchi = (tempsumfbest-fnu[iobj][i]);
                tempchi = tempchi*tempchi/sigi2[i];
                chisum += tempchi;
            }
        }
        // fprintf(tlog,"\n");
        
        /////// Pass out results
        pzall[iz] = chisum;
        *ntemp_all = NTEMP_I;
        for (i=0;i<NTEMP_I;++i) {
            idtempall[iz][i] = oktemp[i];
            coeffs[iz][oktemp[i]] = coeffs_z[i];
        }
      }
//      fclose(tlog);
    }
    
    /*free(okfilt);
    free(oktemp);
    free(tempsumfbest);
    free(sigi2);
    free(bvector);
    
    free(amatrix);
    free(coeffs_z);
    */
    return(0);
}

///////  "Multiplicative updates for nonnegative quadratic programming ..." 
///////    - F. Sha, L. Saul, D. Lee 2006 (earlier version: 2002) [via Google Scholar]
///////
///////
///////   N = N_templates
///////   A = [N x N] matrix where
///////     A[i][j] = sum_over_filters_k(templateI[k]*templateJ[k]/sig2[k])
///////
///////   b = N element vector where
///////     b[i] = -1*sum_over_filters_k(templateI[k]*sourceFlux[k]/sig2[k])
///////
///////   sig2[k] is the squared error on the photometry in filter, k
///////
///////   Difference in Chi^2 between tol=1.e-3 and tol=1.e-7 is only few% and 1e-3 is much faster
///////   Best-fit redshifts are identical at these tolerances (but not at tol=1.e-2).
///////
void nonneg_fact(double **amatrix, double *bvector, double *coeffs, long NTEMP,
                 double toler)
{
    long i,j,itcount,MAXITER;
    double tolnum,toldenom,tol;
    double vold,av;
     
    MAXITER=100000; //// 1.e5

//    av = malloc(sizeof(double)*NTEMP);
//    vold = malloc(sizeof(double)*NTEMP);
      
    tol = 100;
      
    itcount=0;
    while (tol>toler && itcount<MAXITER) {
      
        tolnum=0.;
        toldenom = 0.;
        tol=0;
        for (i=0;i<NTEMP;++i) {
                  
            vold = coeffs[i];
            av = 0.;
                  
            for (j=0;j<NTEMP;++j) av+=amatrix[i][j]*coeffs[j];
            
            ////// Update coeffs in place      
            coeffs[i]*=-1.*bvector[i]/av;
                  
            tolnum+=fabs(coeffs[i]-vold);
            toldenom+=vold;
            // tol+=fabs(coeffs[i]-vold[i])/vold[i];
        }
        // tol/=NTEMP; 
        tol = tolnum/toldenom;
        // printf("      %ld   %le\n",itcount,tol);
        ++itcount;
    }
    //printf("    %ld\n",itcount);
    
    //free(av);
    //free(vold);                  
    //printf(" MAXITER: %ld\n",itcount);            
}                          
                 
                    

