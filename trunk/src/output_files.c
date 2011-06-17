#include "defs.h"

////////// Dump chisq as a function of redshift to [objid].pz file  
void pz_output(char *filename, long iobj)
{
    
    FILE *pzp;
    long iz;
    
  //  sprintf(pzfile,"%s/%ld.pz",OUTPUT_DIRECTORY,objid[iobj]);
    if (!(pzp = fopen(filename,"w"))) {
        fprintf(stderr,"Oops...can't open the pofz file, %s\n",filename);
        exit(1);
    }
           
    if (APPLY_PRIOR == 1) {
      switch (TEMPLATE_COMBOS) {
        case 1:
            fprintf(pzp,"# z chi2 prior pz t1\n");
            for (iz=0;iz<NZ;++iz)  
                fprintf(pzp," %lf %le %le %le %d\n",ztry[iz],chi2fit[iz],priorkz[Kcolumn][iz],pzout[iz],idtemp1[iz]+1);
            break;
        case 2:
            fprintf(pzp,"# z chi2 prior pz t2a t2b\n");
            for (iz=0;iz<NZ;++iz)  
                fprintf(pzp," %lf %le %le %le %d %d\n",ztry[iz],chi2fit[iz],priorkz[Kcolumn][iz],pzout[iz],idtemp2a[iz]+1,idtemp2b[iz]+1);
            break;
        case 99:
            fprintf(pzp,"# z chi2 prior pz\n");
            for (iz=0;iz<NZ;++iz)  
                fprintf(pzp," %lf %le %le %le\n",ztry[iz],chi2fit[iz],priorkz[Kcolumn][iz],pzout[iz]);
            break;
      }
    } else {
      switch (TEMPLATE_COMBOS) {
        case 1:
            fprintf(pzp,"# z chi2 t1\n");
            for (iz=0;iz<NZ;++iz)  
                fprintf(pzp," %lf %le %d\n",ztry[iz],chi2fit[iz],idtemp1[iz]+1);
            break;
        case 2:
            fprintf(pzp,"# z chi2 t2a t2b\n");
            for (iz=0;iz<NZ;++iz)  
                fprintf(pzp," %lf %le %d %d\n",ztry[iz],chi2fit[iz],idtemp2a[iz]+1,idtemp2b[iz]+1);
            break;
        case 99:
            fprintf(pzp,"# z chi2 prior pz\n");
            for (iz=0;iz<NZ;++iz)  
                fprintf(pzp," %lf %le\n",ztry[iz],chi2fit[iz]);
            break;
      }
    }        
    fclose(pzp);
}
        
////////// Dump full best fit template spectrum to [objid].temp_sed file
void temp_sed_output(char *filename, long iobj, long izbest, long izprior)
{

    FILE *tsp;
    long j,itemp;
    double sedtemp1[NTEMPLMAX], sedtemp2[NTEMPLMAX],flamtofnu,igm_corr;
    
   // sprintf(tsfile,"%s/%ld.temp_sed",OUTPUT_DIRECTORY,objid[iobj]);
    if (!(tsp = fopen(filename,"w"))) {
        fprintf(stderr,"Oops...can't open the .temp_sed file, %s\n",filename);
        exit(1);
    }
            
    fprintf(tsp,"# lambda  tempflux ");
    if (APPLY_PRIOR && FIX_ZSPEC==0) fprintf(tsp,"lambda_zprior tempflux_zprior");
    fprintf(tsp,"\n");
            
    /////// Print template information in header line of temp_sed file
    switch (TEMPLATE_COMBOS) {
        case 1:
            fprintf(tsp,"# z=%6.3lf templ:norm %d: %le\n",ztry[izbest],idtemp1[izbest]+1, atemp1[izbest]/tnorm[idtemp1[izbest]]);
            if (APPLY_PRIOR) fprintf(tsp,"# z_prior=%6.3lf templ:norm %d: %le\n",ztry[izprior],idtemp1[izprior]+1, atemp1[izprior]/tnorm[idtemp1[izprior]]);
            break;
        case 2:
            fprintf(tsp,"# z=%6.3lf templ:norm %d: %le + %d: %le\n",ztry[izbest],idtemp2a[izbest]+1, atemp2a[izbest]/tnorm[idtemp2a[izbest]],
                                                             idtemp2b[izbest]+1, atemp2b[izbest]/tnorm[idtemp2b[izbest]]);
            if (APPLY_PRIOR) fprintf(tsp,"# z=%6.3lf templ:norm %d: %le + %d: %le\n",ztry[izprior],idtemp2a[izprior]+1, atemp2a[izprior]/tnorm[idtemp2a[izprior]],
                                                             idtemp2b[izprior]+1, atemp2b[izprior]/tnorm[idtemp2b[izprior]]);
            break;
        case 99:
            fprintf(tsp,"# z=%6.3lf templ:norm",ztry[izbest]);
            for (itemp=0;itemp<ntemp_all;++itemp)
                    fprintf(tsp," %d: %le",idtempall[izbest][itemp]+1,coeffs[izbest][itemp]/tnorm[itemp]);
            fprintf(tsp,"\n");
            if (APPLY_PRIOR) {
                fprintf(tsp,"# z_prior=%6.3lf templ:norm",ztry[izprior]);
                for (itemp=0;itemp<ntemp_all;++itemp)
                    fprintf(tsp," %d: %le",idtempall[izprior][itemp]+1,coeffs[izprior][itemp]/tnorm[itemp]);
                fprintf(tsp,"\n");
            }
            break;
    }
        
    /////// Change to F_nu and apply IGM if needed, looping over wavelength points in the rebinned templates
        
    for (j=0;j<NTEMPL;++j) {
  
        flamtofnu = templ[j]/5500.;
        flamtofnu = flamtofnu*flamtofnu;
    
        //////// IGM correction following Madau 1995
        igm_corr = 1.0;
        if (APPLY_IGM) {
            if ((templ[j]<912.0)) igm_corr = 0.0 ;
            if ((templ[j]>=912.0)&&(templ[j]<1026.0)) 
                igm_corr = 1.0-dbsum[izbest];
            if ((templ[j]>=1026.0)&&(templ[j]<1216.0))
                igm_corr = 1.0-dasum[izbest];
	    }
     
        igm_corr = igm_corr*flamtofnu;
                    
        switch (TEMPLATE_COMBOS) {
            case 1:
                sedtemp1[j] = atemp1[izbest]*tempf[idtemp1[izbest]][j]*igm_corr;
                break;
            case 2:
                sedtemp1[j] = atemp2a[izbest]*tempf[idtemp2a[izbest]][j]*igm_corr+
                              atemp2b[izbest]*tempf[idtemp2b[izbest]][j]*igm_corr;
                break;
            case 99:
                sedtemp1[j]=0.;
                for (itemp=0;itemp<ntemp_all;++itemp) 
                    sedtemp1[j] += coeffs[izbest][itemp]*tempf[idtempall[izbest][itemp]][j]*igm_corr;
                break;
        }
    
        fprintf(tsp,"%le  %le",templ[j]*(1+ztry[izbest]),sedtemp1[j]); 
    
        if (APPLY_PRIOR && FIX_ZSPEC==0) {
                
            flamtofnu = templ[j]/5500.;
            flamtofnu = flamtofnu*flamtofnu;
        
            igm_corr = 1.0;
            if (APPLY_IGM) {
                if ((templ[j]<912.0)) igm_corr = 0.0 ;
                if ((templ[j]>=912.0)&&(templ[j]<1026.0)) 
                    igm_corr = 1.0-dbsum[izprior];
                if ((templ[j]>=1026.0)&&(templ[j]<1216.0))
                    igm_corr = 1.0-dasum[izprior];
	        }
            igm_corr = igm_corr*flamtofnu;
                
            switch (TEMPLATE_COMBOS) {
                case 1:
                    sedtemp2[j] = atemp1[izprior]*tempf[idtemp1[izprior]][j]*igm_corr;
                    break;
                case 2:
                    sedtemp2[j] = atemp2a[izprior]*tempf[idtemp2a[izprior]][j]*igm_corr+
                                  atemp2b[izprior]*tempf[idtemp2b[izprior]][j]*igm_corr;
                    break;
                case 99:
                    sedtemp2[j]=0.;
                    for (itemp=0;itemp<ntemp_all;++itemp) 
                        sedtemp2[j] += coeffs[izprior][itemp]*tempf[idtempall[izprior][itemp]][j]*igm_corr;
                    break;
            }   
            fprintf(tsp,"  %le  %le",templ[j]*(1+ztry[izprior]),sedtemp2[j]);                    
        }
        fprintf(tsp,"\n");                        
    }
                           
    fclose(tsp);

}  ////// TEMP_SED_FILE

/////// Dump observed SED file to [objid].obs_sed, including columns of 
/////// best-fit templates integrated through the filter curves
void obs_sed_output(char *filename, long iobj, long izbest, long izprior)       
{    
    FILE *obp;
    long j,itemp;
    double sedtemp2[NTEMPLMAX], sigi;
    
 //   sprintf(obfile,"%s/%ld.obs_sed",OUTPUT_DIRECTORY,objid[iobj]);
    if (!(obp = fopen(filename,"w"))) {
        fprintf(stderr,"Oops...can't open the .temp_sed file, %s\n",filename);
        exit(1);
    }
    //////// Header line
    switch (TEMPLATE_COMBOS) {
        case 1:
            fprintf(obp,"# lambda flux_cat err_cat err_full temp1_z ");
            break;
        case 2:
            fprintf(obp,"# lambda flux_cat err_cat err_full temp2_z ");
            break;
        case 99:
            fprintf(obp,"# lambda flux_cat err_cat err_full tempa_z ");
            break;
    }
    if (APPLY_PRIOR && FIX_ZSPEC==0) switch (TEMPLATE_COMBOS) {
        case 1:
            fprintf(obp,"temp1_zprior");
            break;
        case 2:
            fprintf(obp,"temp2_zprior");
            break;
        case 99:
            fprintf(obp,"tempa_zprior");
            break;
    }
    fprintf(obp,"\n");
            
            //////// Add filter/template fluxes IF they satisfy NOT_OBS_THRESHOLD
    for (j=0;j<nusefilt;++j) {
        //////// Add template fluxes to OBS_SED_FILE
        if (fnu[iobj][j] >= NOT_OBS_THRESHOLD && efnu[iobj][j] > 0.) { //// NOT_OBS_THRESHOLD) {
            
            sigi = pow((SYS_ERR*SYS_ERR+temp_errf[izbest][j]*temp_errf[izbest][j]*TEMP_ERR_A2*TEMP_ERR_A2)*fnu[iobj][j]*fnu[iobj][j]+
                 efnu[iobj][j]*efnu[iobj][j],0.5);
            
            switch (TEMPLATE_COMBOS) {
                case 1:
                    fprintf(obp,"%le %e %e %e %e ",lambdac[j],fnu[iobj][j],efnu[iobj][j],sigi,
                                                atemp1[izbest]*tempfilt[izbest][idtemp1[izbest]][j]);
                    break;
                case 2:
                    fprintf(obp,"%le %e %e %e %e ",lambdac[j],fnu[iobj][j],efnu[iobj][j],sigi,
                                    atemp2a[izbest]*tempfilt[izbest][idtemp2a[izbest]][j]+
                                    atemp2b[izbest]*tempfilt[izbest][idtemp2b[izbest]][j]);
                    break;
                case 99:
                    sedtemp2[0]=0.;
                    for (itemp=0;itemp<ntemp_all;++itemp) {
                        sedtemp2[0] += coeffs[izbest][itemp]*tempfilt[izbest][idtempall[izbest][itemp]][j];
                    }
                    fprintf(obp,"%le %e %e %e %e ",lambdac[j],fnu[iobj][j],efnu[iobj][j],sigi,sedtemp2[0]);
                    break;
            }
                    
            if (APPLY_PRIOR && FIX_ZSPEC==0) {

                sigi = pow((SYS_ERR*SYS_ERR+temp_errf[izprior][j]*temp_errf[izprior][j]*TEMP_ERR_A2*TEMP_ERR_A2)*fnu[iobj][j]*fnu[iobj][j]+
                              efnu[iobj][j]*efnu[iobj][j],0.5);
                                
                switch (TEMPLATE_COMBOS) {
                    case 1:
                        fprintf(obp,"%e ",atemp1[izprior]*tempfilt[izprior][idtemp1[izprior]][j]);
                        break;
                    case 2:
                        fprintf(obp,"%e ",atemp2a[izprior]*tempfilt[izprior][idtemp2a[izprior]][j]+
                                          atemp2b[izprior]*tempfilt[izprior][idtemp2b[izprior]][j]);
                        break;
                    case 99:
                        sedtemp2[0]=0.;
                        
                        for (itemp=0;itemp<ntemp_all;++itemp) 
                            sedtemp2[0] += coeffs[izprior][itemp]*tempfilt[izprior][idtempall[izprior][itemp]][j];                
                        
                        fprintf(obp,"%e ",sedtemp2[0]);
                        break;
                }
            }
            fprintf(obp,"\n");
        }
    }
    fclose(obp);
}
