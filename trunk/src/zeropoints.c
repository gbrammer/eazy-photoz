#include "defs.h"

int zeropoints () {

    long i,j,itemp,iobj,izbest,k;

    double *sf,*sfuse, *sfiter, *sweight,*smodel,psiold, psinew, psitol,sigi; 
    long sfcount,MAXITER;
    double save_nmf_tolerance,save_temp_err_a2;
    int save_fix_zspec,nzpfilt;
    
    char zpfile[512];
    FILE *zpf;
    
    izbest = 0;
    
    strcpy(zpfile,OUTPUT_DIRECTORY);
    strcat(zpfile,"/");
    strcat(zpfile,MAIN_OUTPUT_FILE);
    strcat(zpfile,".zeropoint");
    if (!(zpf = fopen(zpfile,"w"))) {
        fprintf(stderr,"\n\nuh,oh. Couldn't open %s.  Does the directory, %s, exist?\n",zpfile,OUTPUT_DIRECTORY);
        exit(1);
    }
    
    sf = malloc(sizeof(double)*nusefilt);
    smodel = malloc(sizeof(double)*nusefilt);
    sfuse = malloc(sizeof(double)*nusefilt);
    sfiter = malloc(sizeof(double)*nusefilt);
    sweight = malloc(sizeof(double)*nusefilt);
    
    for (j=0;j<nusefilt;++j) if (zpfilterid[j] >=0) {
        sfuse[j] = zpfactor[filt_to_col[j]]; 
        sfiter[j] = sfuse[j];
    } else {
        sfuse[j]=1;
        sfiter[j] = sfuse[j];
    }
    
    if (zpcontinue == 0) {
        for (j=0;j<nusefilt;++j) if (zpfilterid[j] >= 0) 
            fprintf(zpf,"F%0d %lf\n",zpfilterid[j],sfuse[j]);
        
        fprintf(stderr,"Applying zeropoint offsets from zphot.zeropoint...");
        
        for (iobj=0; iobj<nobj;++iobj) for (j=0;j<nusefilt;++j) if (zpfilterid[j] >=0 && fnu[iobj][j] >= NOT_OBS_THRESHOLD && efnu[iobj][j] > 0.) {
            fnu[iobj][j]*=sfuse[j];
            efnu[iobj][j]*=sfuse[j];
        }
        fprintf(stderr,"Done\n");
                
        fclose(zpf);
        return(0);
    } else {
        /////// Initial offests
        for (iobj=0; iobj<nobj;++iobj) for (j=0;j<nusefilt;++j) if (zpfilterid[j] >=0 && fnu[iobj][j] >= NOT_OBS_THRESHOLD && efnu[iobj][j] > 0.) {
            fnu[iobj][j]*=sfuse[j];
            efnu[iobj][j]*=sfuse[j];
        }
    }    
    
    fprintf(stderr," Calculating zeropoint offsets...");
            
    nzpfilt = 0;
    for (i=0;i<nusefilt;++i) if (zpfilterid[i] > 0) ++nzpfilt;
        
    if (nzpfilt == 0) {
        fprintf(zpf,"# No filters selected\n");
        fclose(zpf);
        return(0);
    }
    
    psiold = 100; 
    psinew = 10;
    psitol = 1.e-4;
  
  save_fix_zspec = FIX_ZSPEC;
  save_nmf_tolerance = NMF_TOLERANCE;
  save_temp_err_a2 = TEMP_ERR_A2;
  
  TEMP_ERR_A2 = 0.;
  NMF_TOLERANCE /= 100.;
  
  FIX_ZSPEC = 1;

    fprintf(zpf,"# iter ");
    printf("\n# iter ");
    for (j=0;j<nusefilt;++j) if (zpfilterid[j] >= 0) {
        fprintf(zpf," F%0d",zpfilterid[j]);
        printf(" F%0d",zpfilterid[j]);
    }
    fprintf(zpf," tol\n");
    printf(" tol\n");

  sfcount=0;  
  MAXITER=5000;

  while (fabs(psinew) > ZP_OFFSET_TOL && sfcount < MAXITER) {
    ++sfcount;

    fprintf(zpf,"# %8ld ",(sfcount+0)/1);
    printf("# %8ld ",(sfcount+0)/1);
    psiold=psinew;
    psinew=0;

    for (k=0;k<nusefilt;++k) if (zpfilterid[k] > 0) {
      //printf("\nk %ld",k);
      
      sf[k] = 0.;
      sweight[k] = 0.;
      
      for (iobj=0; iobj<nobj; ++iobj) if (zspec[iobj] > 0) { 
                
        //// Downweight filters in turn to calculate residuals, this won't work if
        //// you're trying to fit all filters simultaneously.
        efnu[iobj][k]*=1000.;
            
        ngoodfilters = 0;
      
        for(j=0;j<nusefilt;++j)
            if (fnu[iobj][j] > NOT_OBS_THRESHOLD && efnu[iobj][j] > NOT_OBS_THRESHOLD)
                ++ngoodfilters;
      
        //////// If too few filters, skip  
        //if (ngoodfilters >= N_MIN_COLORS) {
        if (ngoodfilters == nusefilt) {

            /////////////////  Call out to getphotz to get Chisq as a function of z //////////////////
            ////////// Unused arrays are not malloc-ed, is that OK
            getphotz(iobj, pz1, idtemp1, atemp1, 
                       pz2, idtemp2a, idtemp2b, atemp2a, atemp2b,
                       pzall, idtempall, coeffs, &ntemp_all,-1);
            //////////////////////////////////////////////////////////////////////////////////////////

            switch (TEMPLATE_COMBOS) {
                case 1:
                    izbest=0;
                    for (i=0;i<NZ;++i) {
                        chi2fit[i] = pz1[i];
                        pzout[i] = chi2fit[i];
                        if (pz1[i]<chi2fit[izbest]) izbest=i;
                        pz1[i]=1.e30;
                    }
                    break;
                case 2:
                    izbest=0;
                    for (i=0;i<NZ;++i) {
                        chi2fit[i] = pz2[i];
                        pzout[i] = chi2fit[i];
                        if (pz2[i]<chi2fit[izbest]) izbest=i;
                        pz2[i]=1.e30;
                    }
                    break;
                case 99:
                    izbest=0;
                    for (i=0;i<NZ;++i) {
                        chi2fit[i] = pzall[i];
                        pzout[i] = chi2fit[i];
                        if (pzall[i]<chi2fit[izbest]) izbest=i;
                        pzall[i]=1.e30;
                    }
            }
            pzout[izbest] = exp(-0.5*chi2fit[izbest]);
       
            if (fnu[iobj][k] > 0. && efnu[iobj][k] > 0. && fnu[iobj][k]/efnu[iobj][k] > 0.003) { //// 0. rather than flux, err > NOT_OBS_THRESHOLD) 
                    
                    sigi = pow((SYS_ERR*SYS_ERR+temp_errf[izbest][k]*temp_errf[izbest][k]*TEMP_ERR_A2*TEMP_ERR_A2)*fnu[iobj][k]*fnu[iobj][k]+
                                 efnu[iobj][k]*efnu[iobj][k]/1.e6,0.5);
                    
                    //sigi = 1;
                    
                    switch (TEMPLATE_COMBOS) {
                      case 1:
                        smodel[k] = (atemp1[izbest]*tempfilt[izbest][idtemp1[izbest]][k]);
                        break;
                      case 2:
                        smodel[k] = (atemp2a[izbest]*tempfilt[izbest][idtemp2a[izbest]][k]+
                                    atemp2b[izbest]*tempfilt[izbest][idtemp2b[izbest]][k]);
                        break;
                      case 99:
                        smodel[k] = 0.;
                        for (itemp=0;itemp<ntemp_all;++itemp) {
                            smodel[k] += coeffs[izbest][itemp]*tempfilt[izbest][idtempall[izbest][itemp]][k];
                        }
                        break;
                    }
                    /////// Weight by chisq probability // S/N
                    if (smodel[k] > 0) {
                        sf[k] += (fnu[iobj][k])/smodel[k]*pzout[izbest] ; // *fnu[iobj][j]/sigi;
                        sweight[k]+= pzout[izbest]; // fnu[iobj][j]/sigi;
                    }                    
               // }
             
            } 
          
        }
        
        //// Unscale errors
        efnu[iobj][k]/=1000.;
            
      } /// loop over objects in catalog
    
      psinew+=fabs(sfiter[k]-1./(sf[k]/sweight[k]))/nzpfilt;

      sfiter[k] = 1./(sf[k]/sweight[k]);
      sfuse[k] *= sfiter[k];        

      fprintf(zpf,"%lf ",sfuse[k]);
      printf("%lf ",sfuse[k]);
            
    }

      for (iobj=0; iobj<nobj; ++iobj)
        for (k=0;k<nusefilt;++k) if (zpfilterid[k]>=0) 
          if (fnu[iobj][k] >= NOT_OBS_THRESHOLD && efnu[iobj][k] > 0.) {
            fnu[iobj][k]*= sfiter[k];
            efnu[iobj][k]*= sfiter[k];        
      }

    fprintf(zpf," %lf\n",psinew);
    printf(" %lf\n",psinew);
        
        
  } /// WHILE psinew-1 gt tol
  
  for (j=0;j<nusefilt;++j) if (zpfilterid[j] >= 0) 
    fprintf(zpf,"F%0d %lf\n",zpfilterid[j],sfuse[j]);
    
  fclose(zpf);
  
  FIX_ZSPEC = save_fix_zspec;
  NMF_TOLERANCE = save_nmf_tolerance;
  TEMP_ERR_A2 = save_temp_err_a2;
  
  fprintf(stderr,"Done.\n");
  return(0);
      
}
