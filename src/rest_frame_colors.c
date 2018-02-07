#include "defs.h"

int compare_doubles (const void *a, const void *b) {
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

void init_rf_file(FILE *fprf) {
    
    long i,j;
    
    fprintf(fprf,"# id z DM nfilt_fit chi2_fit ");
    for (i=0;i<nrestfilt;++i) fprintf(fprf,"  L%d",filt_defnum[nusefilt+i]);
    
    if (RF_ERRORS) {
        fprintf(fprf, "  p16 p84 p025 p975");
    }
    
    fprintf(fprf,"\n#\n");

    for (i=0;i<nrestfilt;++i) fprintf(fprf,"# %d: %s, %13.5le\n",filt_defnum[nusefilt+i],pusefilt[nusefilt+i]->filtname,lambdac[nusefilt+i]);

    fprintf(fprf,"#\n# L: eazy-interpolated color indices\n");
    
    fprintf(fprf,"#\n# Rest-frame colors computed using templates in %s\n#\n",RF_TEMPLATES_FILE);
    if (FIX_ZSPEC) fprintf(fprf,"# z = z_spec\n#\n"); else {
        fprintf(fprf,"# z = %s",Z_COLUMN);
        if (USE_ZSPEC_FOR_REST) fprintf(fprf,", will use z_spec when available");
        fprintf(fprf,"\n#\n");
    }
    
    //// will need sorted array later
    qsort (lambdac_sort, nusefilt, sizeof (double), compare_doubles);
    for (i=0;i<nusefilt;++i) for (j=0;j<nusefilt;++j) if (lambdac[i] == lambdac_sort[j]) lc_sort_idx[j] = i;
    //for (i=0;i<nusefilt;++i) printf("sort: %lf %lf %lf %ld\n",lambdac[i],lambdac[lc_sort_idx[i]],lambdac_sort[i],lc_sort_idx[i]); exit(1);
    
}

/// Rest-frame colors
void rest_frame_colors() {
    
    int getzStatus, nfilt_weight;
    long i,j,k,iobj,*izuse;
    char outparfile[1024],outrffile[1024],coeffile[1024],strfmt[64];
    double lcmin,lcmax,dlc,tau;
    double *zuse,ztest,zrest;
    double comove, DM;
    double chi2_fit,temp_flux_filt1, temp_flux_filt2;
    
    dlc = RF_PADDING;
    
    void interrest(long iobj, long izrest, int i_rf, double *interrest_flux, int *flag);
    // double irest_flux;
    // int irest_flag;
    
    int32_t NOBJ32, NRESTFILT32, NTEMP32;
    FILE *fpcoeff;
    
    //// Full rest-frame color grid
    char pz_file[1024];
    FILE *pzf;    
    int32_t pzdummy;
    double *chi2_fit_full, *colors_z, color_blue, color_red, pztot;
    double *colors_sort, *colors_cumprob, color_lo, color_hi;
    
    long NZ_prior, NK_prior, izbest, izprior, *colors_sort_idx, icolor;
        
    //// Open ".pz" file for getting p(z)
    if (RF_ERRORS) {
        sprintf(pz_file,"%s/%s.pz",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE);
        pzf = fopen(pz_file,"r");
        fread(&pzdummy,sizeof(int32_t),1,pzf); // NZ
        fread(&pzdummy,sizeof(int32_t),1,pzf); // nobj
        
        chi2_fit_full = malloc(sizeof(double)*NZ);
        colors_z = malloc(sizeof(double)*NZ);
        colors_sort = malloc(sizeof(double)*NZ);
        colors_cumprob = malloc(sizeof(double)*NZ);
        colors_sort_idx = malloc(sizeof(long)*NZ);
        
        if (APPLY_PRIOR) {
            get_prior_file_size(&NZ_prior, &NK_prior);
        }
    
    }

    /////////////////////////////////////////////////////////////////
    ////
    ////            New output param file for RF filters
    ////
    /////////////////////////////////////////////////////////////////
    if (nrestfilt == 1) {
        lcmin=lambdac[nusefilt];
        lcmax=lambdac[nusefilt];
        // dlc=1000;
        sprintf(outparfile,"%s/%s.%d.rf.param",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE,filt_defnum[nusefilt]);
        sprintf(outrffile,"%s/%s.%d.rf",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE,filt_defnum[nusefilt]);
        sprintf(coeffile,"%s/%s.%d.coeff",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE,filt_defnum[nusefilt]);
    } else {
        // dlc=1000;
        if (lambdac[nusefilt] > lambdac[nusefilt+1]) {
            lcmin=lambdac[nusefilt+1];
            lcmax=lambdac[nusefilt];
            if ((lcmin-dlc) < 0) dlc=(lcmin-1);

            sprintf(outparfile,"%s/%s.%d-%d.rf.param",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE,filt_defnum[nusefilt+1],filt_defnum[nusefilt]);
            sprintf(outrffile,"%s/%s.%d-%d.rf",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE,filt_defnum[nusefilt+1],filt_defnum[nusefilt]);
            sprintf(coeffile,"%s/%s.%d-%d.coeff",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE,filt_defnum[nusefilt+1],filt_defnum[nusefilt]);
        } else {
            lcmin=lambdac[nusefilt];
            lcmax=lambdac[nusefilt+1];
            if ((lcmin-dlc) < 0) dlc=(lcmin-1);
            
            sprintf(outparfile,"%s/%s.%d-%d.rf.param",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE,filt_defnum[nusefilt],filt_defnum[nusefilt+1]);
            sprintf(outrffile,"%s/%s.%d-%d.rf",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE,filt_defnum[nusefilt],filt_defnum[nusefilt+1]);
            sprintf(coeffile,"%s/%s.%d-%d.coeff",OUTPUT_DIRECTORY,MAIN_OUTPUT_FILE,filt_defnum[nusefilt],filt_defnum[nusefilt+1]);
        }
    }

    //printf("outparfile: %s\n",outparfile);
    if (!(fplog = fopen(outparfile,"w"))) {
        fprintf(stderr,"\n\nuh,oh. Couldn't open %s.  Does %s exist?\n",outparfile,OUTPUT_DIRECTORY);
        exit(1);
    }
    printparams_logfile(fplog);
    fclose(fplog);

    if (!(fplog = fopen(outrffile,"w"))) {
        fprintf(stderr,"\n\nuh,oh. Couldn't open %s.  Does %s exist?\n",outparfile,OUTPUT_DIRECTORY);
        exit(1);
    }
    init_rf_file(fplog);

    if (BINARY_OUTPUT) {
        if (!(fpcoeff = fopen(coeffile,"w"))) {
            fprintf(stderr,"\n\nuh,oh. Couldn't open %s.  Does %s exist?\n",coeffile,OUTPUT_DIRECTORY);
            exit(1);
        }
        NOBJ32 = (int32_t) nobj;
        NRESTFILT32 = (int32_t) nrestfilt;
        NTEMP32 = (int32_t) NTEMP;
        fwrite(&NOBJ32,sizeof(int32_t),1,fpcoeff);
        fwrite(&NRESTFILT32,sizeof(int32_t),1,fpcoeff);
        fwrite(&NTEMP32,sizeof(int32_t),1,fpcoeff);
        //// Write z=0 fluxes of templates used for rest-frame fluxes in one or both rest-frame filters
        for (i=0;i<NTEMP_REST;++i) fwrite(&tempfilt_rest[0][i][nusefilt],sizeof(double),1,fpcoeff); //coeffs[izuse[iobj]][i]*tempfilt_rest[0][i][nusefilt]
        if (nrestfilt == 2) for (i=0;i<NTEMP_REST;++i) fwrite(&tempfilt_rest[0][i][nusefilt+1],sizeof(double),1,fpcoeff);
    }
    
    /////////////////////////////////////////////////////////////////
    //////
    //////     Remake template error function to downweight fits outside of
    //////     range covered by desired rest-frame color filters
    //////
    /////////////////////////////////////////////////////////////////

    //dlc = 1.e5; //// Uncomment this line to use all filters in the fit without downweighting
    
    tau=0.05;
    i=0;
    while (templ[i] < (lcmin-dlc) && i < NTEMPL) {
        temp_err[i] = pow(10,pow(log10(templ[i])-log10(lcmin-dlc),2)/tau)-1;
        ++i;
    }
    while (templ[i] < (lcmax+dlc) && i < NTEMPL) {
        temp_err[i] = 0.;
        ++i;
    }
    while (i < NTEMPL) {
        temp_err[i] = pow(10,pow(log10(templ[i])-log10(lcmax+dlc),2)/tau)-1;
        ++i;      
    }
    for (i=0;i<NTEMPL;++i) {
        if (temp_err[i] > 1.e3) temp_err[i] = 1.e3;
        // printf("9999 %le %le\n",templ[i],temp_err[i]);
    }
    
    //// now compute template error at >(redshifted)< filter wavelengths
    for (i=0;i<NZ;++i) {
      for (j=0;j<nusefilt;++j) {
        for (k=0;k<NTEMPL;++k) if (templ[k]>=lambdac[j]/(1.+ztry[i])) break;
        // temp_errf[i][j] = temp_err[k];  //// New template error to downweight filters outside of desired RF color
        temp_errf[i][j] = sqrt(pow(temp_errf[i][j],2)+temp_err[k]*temp_err[k]); //// Quadrature combination of original + new
      }
    }

    /////////////////////////////////////////////////////////////////
    ////
    ////       Find which redshifts to use
    ////
    /////////////////////////////////////////////////////////////////
    zuse = malloc(sizeof(double)*nobj);
    izuse = malloc(sizeof(long)*nobj);
    
    for (i=0;i<nobj;++i) zuse[i] = z_a[i];
    
    if (!strcmp(Z_COLUMN,"z_p")) for (i=0;i<nobj;++i) zuse[i] = z_p[i];
    if (!strcmp(Z_COLUMN,"z_m1")) for (i=0;i<nobj;++i) zuse[i] = z_m1[i];
    if (!strcmp(Z_COLUMN,"z_m2")) for (i=0;i<nobj;++i) zuse[i] = z_m2[i];
    if (!strcmp(Z_COLUMN,"z_peak")) for (i=0;i<nobj;++i) zuse[i] = z_peak[i];
    
    if (USE_ZSPEC_FOR_REST) for (i=0;i<nobj;++i) if (zspec[i] > 0) zuse[i] = zspec[i];
    
    for (i=0;i<nobj;++i) {
        zrest=zuse[i];
        ztest=100;
        if (zrest < 0) izuse[i] = -1; else {
            for (j=0;j<NZ;++j) if (fabs(ztry[j]-zrest) < ztest) {
                ztest = fabs(ztry[j]-zrest);
                izuse[i]=j;
            }
        }       
    }
    
    /////////////////////////////////////////////////////////////////
    ////
    ////      Main loop over objects
    ////
    /////////////////////////////////////////////////////////////////
    for (iobj=0;iobj<nobj;++iobj) {
        
        //// Progress indicator
        if (nobj >= 10) 
          if ((iobj % (nobj/10))==0) fprintf(stderr,">");

        zrest = zuse[iobj];
        
        //// Read p(z) from the .pz file
        if (RF_ERRORS) fread(chi2_fit_full,sizeof(double)*NZ,1,pzf);
        
        strcpy(strfmt,objid[iobj]);
        i=0; while (strfmt[i] != '\0') ++i;
        --i; while (i<=idsize && i < 62) strfmt[++i] = ' '; // max size of objid: IDBUFF=64
        strfmt[i] = '\0';
        
        // printf("%s %lf\n",strfmt,zrest);
        
        //// Find filters above "not-observed" threshold
        ngoodfilters = 0;
        for(j=0;j<nusefilt;++j) {
            if (fnu[iobj][j] > NOT_OBS_THRESHOLD && efnu[iobj][j] > 0) {
                ++ngoodfilters;
                okfilt[j] = 1;
            } else okfilt[j] = 0;
        }
        
        if (ngoodfilters < N_MIN_COLORS || izuse[iobj] < 0) {
            fprintf(fplog,"%s %8.5f %6.2f %3d %13.5e",strfmt,zrest, -99.,-99,-99.);
            
            for (i=0;i<nrestfilt;++i) fprintf(fplog,"%14.5e",-99.);
            
            if (RF_ERRORS) {
                fprintf(fplog, "  -99 -99 -99 -99");
            }
            
            fprintf(fplog,"\n");
            if (BINARY_OUTPUT) fwrite(coeffs[0],sizeof(double)*NTEMP_REST,1,fpcoeff);            
            continue; //// next object
        } 
            
        //// Do fit for coefficients
        getzStatus = getphotz(iobj, pz1, idtemp1, atemp1, 
                       pz2, idtemp2a, idtemp2b, atemp2a, atemp2b,
                       pzall, idtempall, coeffs, &ntemp_all, izuse[iobj]);
            
        if (getzStatus == 1) {
            printf("%s: oops, got singular matrix for template normalization!  skipping...\n",objid[iobj]);
            fprintf(fplog,"%s %8.5f %6.2f %3d %13.5e",strfmt,zrest, -99.,-99,-99.);    
            for (i=0;i<nrestfilt;++i) fprintf(fplog,"%14.5e",-99.);
            fprintf(fplog,"\n");        
            if (BINARY_OUTPUT) fwrite(coeffs[0],sizeof(double)*NTEMP_REST,1,fpcoeff);            
            continue; //// next object
        }     

        chi2_fit = 0.;
        switch (TEMPLATE_COMBOS) {
            case 1:
                chi2_fit = pz1[izuse[iobj]];
                break;
            case 2:
                chi2_fit = pz2[izuse[iobj]];
                break;
            case 99:
                chi2_fit = pzall[izuse[iobj]];
                break;
        }
        
        //// dH = c / H0 = 3000 / h Mpc
        //// dC = dH * \int_0^z 1/eta dz = dH * "comove"
        //// dM = dC (for Omega_k = 0)
        //// dL = (1+z) * dM
        //// DM = 5 * log10( dL / 10pc) 
        //DM = 5. * log10( (1.+zrest) * comove * 1e5 * 3000. / (H0/100.) ) - 2.5 * log10( 1.+zrest ) ; // the last factor from InterRest, but wrong??
        i=0; interpol(zrest,&comove,zcomove,comovingdist,NZCOMOVE,&i);        
        DM = 5. * log10( (1.+zrest) * 3000. / (H0/100.) * 1.e6 * comove / 10. / sqrt(1.+zrest) ); // last factor of sqrt(1+z) from Rudnick 03
        // DM = 5. * log10( (1.+zrest) * 3000. / (H0/100.) * 1.e6 * comove / 10. ); // 
        if (zrest==0) DM = 0;
            
        ///// Number of filters with reasonable weights in the fit
        nfilt_weight = 0;
        for (i=0;i<nusefilt;++i) if (temp_errf[izuse[iobj]][i] < 0.1 && fnu[iobj][i] > NOT_OBS_THRESHOLD && efnu[iobj][i] > 0) ++nfilt_weight;
                        
        fprintf(fplog,"%s %8.5f %6.2f %3d %13.5e",strfmt,zrest, DM,nfilt_weight,chi2_fit);

        temp_flux_filt1=0.;
        for (i=0;i<NTEMP_REST;++i) {
            temp_flux_filt1+=coeffs[izuse[iobj]][i]*tempfilt_rest[0][i][nusefilt];
            //fprintf(fplog,"\ncoeff %le\n",coeffs[izuse[iobj]][i]);
        }
        fprintf(fplog,"%14.5e",temp_flux_filt1);
              
        //// Don't run InterRest routine here, run InterRest directly if you want it  
        // interrest(iobj,izuse[iobj],0,&irest_flux, &irest_flag);
        // fprintf(fplog,"%14.5e",irest_flux);
        
        if (nrestfilt == 2) {
            temp_flux_filt2=0.;
            for (i=0;i<NTEMP_REST;++i) 
                temp_flux_filt2+=coeffs[izuse[iobj]][i]*tempfilt_rest[0][i][nusefilt+1];
            fprintf(fplog,"%14.5e",temp_flux_filt2);        

            // interrest(iobj,izuse[iobj],1,&irest_flux, &irest_flag);
            // fprintf(fplog,"%14.5e",irest_flux);
        }
           
        if (BINARY_OUTPUT) fwrite(coeffs[izuse[iobj]],sizeof(double)*NTEMP_REST,1,fpcoeff);            
        
        //// Do fit for coefficients at *all* redshifts
        if (RF_ERRORS) {
            
            izbest = 0;
            if (APPLY_PRIOR) {
                apply_prior(priorkz, NZ, NK_prior, izbest, chi2_fit_full,
                            fnu[iobj][PRIOR_FILTER_IDX], &izprior, pzout);
            } else {
                for (i=0;i<NZ;++i) {
                    pzout[i] = exp(-0.5*(chi2_fit_full[i]-chi2_fit_full[izbest])/CHI2_SCALE);
                }
            }
            
            //// Normalization of p(z), make pzout = pzout * dz so can sum directly
            pztot = 0.; 
            pzout[0] *= 0;
            for (i=1; i<NZ; ++i) {
                pzout[i] *= ztry[i]-ztry[i-1];
                pztot += pzout[i];
            }
            //fprintf(stderr, "PZTOT: %.3f\n", pztot);
            
            //printf("ID: %s, zp=%.3f, %.3f\n", strfmt, ztry[izprior], z_p[iobj]);

            if ((pztot == 0) | isnan(pztot)) {
                fprintf(stderr, "\n[%s] Sum p(z) = 0, check this object and the prior\n", strfmt);
                fprintf(fplog, " %.3f  %.3f", 0., 0.);
                fprintf(fplog, " %.3f  %.3f\n", 0., 0.);
                continue;
            }
                                
            //// Get fit coefficients on the redshift grid
            getzStatus = getphotz(iobj, pz1, idtemp1, atemp1, 
                        pz2, idtemp2a, idtemp2b, atemp2a, atemp2b,
                        pzall, idtempall, coeffs, &ntemp_all, -1);
            
            //// Compute the RF color on the redshift grid
            for (i=0; i<NZ; ++i) {
                color_blue=0.;
                color_red=0.;
                for (j=0;j<NTEMP_REST;++j) {
                    color_blue += coeffs[i][j]*tempfilt_rest[0][j][nusefilt];
                    if (nrestfilt == 2) {
                        color_red += coeffs[i][j]*tempfilt_rest[0][j][nusefilt+1];
                    } else {
                        color_red = 1.;
                    }
                }                        
                //colors_z[i] = -2.5*log10(color_blue/color_red);
                //colors_sort[i] = -2.5*log10(color_blue/color_red);
                colors_z[i] = (color_blue/color_red);
                colors_sort[i] = (color_blue/color_red);
            }
            
            //// Sort by colors
            qsort (colors_sort, NZ, sizeof (double), compare_doubles);
            for (i=0;i<NZ;++i) for (j=0;j<NZ;++j) if (colors_z[i] == colors_sort[j]) colors_sort_idx[j] = i;
            
            //// Cumulative probability from p(z) on sorted colors
            colors_cumprob[0] = pzout[colors_sort_idx[0]]/pztot;
            for (i=1;i<NZ;++i) colors_cumprob[i] = colors_cumprob[i-1]+pzout[colors_sort_idx[i]]/pztot;
            
            //for (i=0;i<NZ;++i) printf("%d  %.4f  %.4e\n", iobj, colors_sort[i], colors_cumprob[i]);
            
            //// Interpolate cumulative probability to get confidence intervals
            icolor=0;
            interpol(0.16, &color_lo, colors_cumprob, colors_sort, NZ, &icolor);
            interpol(0.84, &color_hi, colors_cumprob, colors_sort, NZ, &icolor);
            
            fprintf(fplog, " %.3f  %.3f", color_lo, color_hi);
            
            icolor=0;
            interpol(0.025, &color_lo, colors_cumprob, colors_sort, NZ, &icolor);
            interpol(0.975, &color_hi, colors_cumprob, colors_sort, NZ, &icolor);
            
            fprintf(fplog, " %.3f  %.3f", color_lo, color_hi);
            
            // printf("\n\nColor: %.2f (%.2f  -- %.2f) / tot:%.2f\n\n", 
            //        -2.5*log10(temp_flux_filt1/temp_flux_filt2), color_lo, 
            //        color_hi, colors_cumprob[i-1]);
        }
        
        fprintf(fplog,"\n");
        
    }
    
    fclose(fplog);
    
    if (BINARY_OUTPUT) {
        fwrite(tnorm,sizeof(double)*NTEMP_REST,1,fpcoeff);
        fclose(fpcoeff);            
    }
}

/// Rest-frame colors with uncertainties, draw from p(z)
void draw_rest_frame_colors() {

}

void interrest(long iobj, long izrest, int i_rf, double *interrest_flux, int *flag) {
    
    /// Rest-frame colors
    double lam_rf_z, cobs, clza, clzb, *ctemp, *ctemp_sort,csum,mlz,mlze,cfactor,ecobs;
    int i_rf_lo, i_rf_hi, rf_flag,i_rf_temp,icta, ictb;
    long ict;
    
    ctemp = malloc(sizeof(double)*NTEMP_REST);
    ctemp_sort = malloc(sizeof(double)*NTEMP_REST);

    lam_rf_z = lambdac[nusefilt+i_rf]*(1+ztry[izrest]);
            
    //// rf_flag bits, 2^: 
    ////     0 = rf_filter bluer than bluest obs filter
    ////     1 =         redder       reddest 
    ////     2 = observed color bluer than bluest template color
    ////     3 =                redder     reddest
    rf_flag = 0; 
               
    //// Find red and blue filters on either side of the redshifted rest-frame filter
    ict=0; while (lam_rf_z > lambdac_sort[ict] && ict < nusefilt) ++ict;
    //// use sorted
    if (ict==0) {
        i_rf_lo = 0;
        i_rf_hi = 1;
        rf_flag+=pow(2,0);
    }
    if (ict==nusefilt) {
        i_rf_lo = nusefilt-2;
        i_rf_hi = nusefilt-1;
        rf_flag+=pow(2,1);
    }
    if (ict>0 && ict < nusefilt) {
        i_rf_lo = ict-1;
        i_rf_hi = ict;
    }
    //printf("ilo ihi %d %d\n",i_rf_lo,i_rf_hi);
    
    //// Expand around filters where object is undetected (flux < 0) 
    while (fnu[iobj][lc_sort_idx[i_rf_lo]] < 0 && i_rf_lo > 0) --i_rf_lo;
    while (fnu[iobj][lc_sort_idx[i_rf_hi]] < 0 && i_rf_hi < (nusefilt-1)) ++i_rf_hi;
    i_rf_lo = lc_sort_idx[i_rf_lo];           
    i_rf_hi = lc_sort_idx[i_rf_hi];  
             
    //printf("filters: %8.1lf %8.1lf %8.1lf %5.2lf %d %ld\n",lambdac[i_rf_lo], lam_rf_z, lambdac[i_rf_hi],ztry[izrest],i_rf_hi,ict);
    //printf("observed: %lf %lf\n",fnu[iobj][i_rf_lo],fnu[iobj][i_rf_hi]);
                
    if (rf_flag > -1) {
        // printf("lc: %lf %lf %lf\n",lambdac[i_rf_lo],lam_rf_z,lambdac[i_rf_hi]);
        if (fnu[iobj][i_rf_lo] <= 0 || fnu[iobj][i_rf_hi] <=0) {
            *interrest_flux = -99.;
            *flag = -1;
        } else {    
                    
            //// find the templates whose colors bracket the observed colors
            for (i_rf_temp = 0;i_rf_temp<NTEMP_REST;++i_rf_temp) {
                ctemp[i_rf_temp] = -2.5*log10(tempfilt_rest[izrest][i_rf_temp][i_rf_lo]/tempfilt_rest[izrest][i_rf_temp][i_rf_hi]);
                ctemp_sort[i_rf_temp] = ctemp[i_rf_temp];
                // printf("ct  %lf %ld \n",ctemp[i_rf_temp],NTEMP_REST);
            }
                     
            cobs = -2.5*log10(fnu[iobj][i_rf_lo]/fnu[iobj][i_rf_hi]);
            
            //// sort template colors
            qsort (ctemp_sort, NTEMP_REST, sizeof (double), compare_doubles);
            //for (ict = 0;ict<NTEMP;++ict) printf("ctemp %lf %lf %lf\n",ctemp[ict],ctemp_sort[ict],cobs);
            
            ict = 0;
            while (ctemp_sort[ict] < cobs && ict < NTEMP_REST) ++ict;
            ///////////////// ict is the middle point
            // printf("ICT: %d\n",ict);
                    
            if (ict == 0) {
                //// extrapolation off blue end of template color distribution
                icta = 0; while (ctemp[icta] != ctemp_sort[0] && icta < NTEMP_REST) ++icta;
                ictb = 0; while (ctemp[ictb] != ctemp_sort[1] && ictb < NTEMP_REST) ++ictb;
                rf_flag+=pow(2,2);
                //if (ctemp_sort[0]-cobs > ctemp_sort[NTEMP-1]-ctemp_sort[0]) rf_flag+=pow(2,3);
                rf_flag+=pow(2,2);
            }
                    
            if (ict == NTEMP_REST) {
                //// extrapolation off red end
                icta = 0; while (ctemp[icta] != ctemp_sort[NTEMP_REST-2] && icta < NTEMP_REST) ++icta;
                ictb = 0; while (ctemp[ictb] != ctemp_sort[NTEMP_REST-1] && ictb < NTEMP_REST) ++ictb;
                //if (cobs-ctemp_sort[NTEMP_REST-1] > ctemp_sort[NTEMP_REST-1]-ctemp_sort[0]) rf_flag+=pow(2,3);
                rf_flag+=pow(2,3);
            }
                    
            //// Find indices of original array in sorted array
            if (ict > 0 && ict < NTEMP_REST) {
                icta = 0; while (ctemp[icta] != ctemp_sort[ict-1] && icta < NTEMP_REST) ++icta;
                ictb = 0; while (ctemp[ictb] != ctemp_sort[ict] && ictb < NTEMP_REST) ++ictb;
            }
            
            clza = -2.5*log10(tempfilt_rest[izrest][icta][i_rf_lo]/tempfilt_rest[0][icta][nusefilt+i_rf]);
            clzb = -2.5*log10(tempfilt_rest[izrest][ictb][i_rf_lo]/tempfilt_rest[0][ictb][nusefilt+i_rf]);
            
            //for (ict=0;ict<NTEMP_REST_REST;++ict) printf("  obsrf - %lf %lf %lf\n",tempfilt_rest[izrest][ict][i_rf_lo],tempfilt_rest[0][ict][nusefilt+i_rf],-2.5*log10(tempfilt_rest[izrest][ict][i_rf_lo]/tempfilt_rest[0][ict][nusefilt+i_rf]));
            
            // printf("obsrf: %lf %lf %lf %lf\n",tempfilt_rest[izrest][icta][i_rf_lo],tempfilt_rest[0][icta][nusefilt+i_rf],tempfilt_rest[izrest][ictb][i_rf_lo],tempfilt_rest[0][ictb][nusefilt+i_rf]);
            // printf("color: %8.3lf %8.3lf %8.3lf -- %lf %lf\n",ctemp[icta],cobs,ctemp[ictb],clza,clzb);
            // printf("fluxes: %8.4lf %8.4lf %8.4lf\n",tempfilt_rest[izrest][icta][i_rf_lo],tempfilt_rest[1][icta][nusefilt+i_rf],tempfilt_rest[izrest][icta][i_rf_hi]);
            
            cfactor=(clzb-clza)/(ctemp[ictb]-ctemp[icta]);
            csum = clza + (cobs-ctemp[icta])*cfactor;
            
            //// I: InterRest/Rudnick interpolated
            mlz = fnu[iobj][i_rf_lo]/pow(10,-.4*csum);
                        
            //// Formal errors (small?)
            ecobs = efnu[iobj][i_rf_lo]*efnu[iobj][i_rf_lo]/fnu[iobj][i_rf_lo]/fnu[iobj][i_rf_lo]+
                    efnu[iobj][i_rf_hi]*efnu[iobj][i_rf_hi]/fnu[iobj][i_rf_hi]/fnu[iobj][i_rf_hi];
                            
            mlze = efnu[iobj][i_rf_lo]*efnu[iobj][i_rf_lo]/fnu[iobj][i_rf_lo]/fnu[iobj][i_rf_lo]+
                            ecobs*cfactor*cfactor;
                    
            mlze = mlz*sqrt(mlze);
            
            *interrest_flux = mlz;
            *flag = rf_flag;
        }
    }
    
}