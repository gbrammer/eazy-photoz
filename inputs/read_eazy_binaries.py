import os
import numpy as np

def readEazyBinary(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same'):
    """
tempfilt, coeffs, temp_sed, pz = readEazyBinary(MAIN_OUTPUT_FILE='photz', \
                                                OUTPUT_DIRECTORY='./OUTPUT', \
                                                CACHE_FILE = 'Same')

    Read Eazy BINARY_OUTPUTS files into structure data.
    
    If the BINARY_OUTPUTS files are not in './OUTPUT', provide either a relative or absolute path
    in the OUTPUT_DIRECTORY keyword.
    
    By default assumes that CACHE_FILE is MAIN_OUTPUT_FILE+'.tempfilt'.
    Specify the full filename if otherwise. 
    """
        
    root = OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE
    
    ###### .tempfilt
    if CACHE_FILE == 'Same':
        CACHE_FILE = root+'.tempfilt'
    
    if os.path.exists(CACHE_FILE) is False:
        print(('File, %s, not found.' %(CACHE_FILE)))
        return -1,-1,-1,-1
    
    f = open(CACHE_FILE,'rb')
    
    s = np.fromfile(file=f,dtype=np.int32, count=4)
    NFILT=s[0]
    NTEMP=s[1]
    NZ=s[2]
    NOBJ=s[3]
    tempfilt = np.fromfile(file=f,dtype=np.double,count=NFILT*NTEMP*NZ).reshape((NZ,NTEMP,NFILT)).transpose()
    lc = np.fromfile(file=f,dtype=np.double,count=NFILT)
    zgrid = np.fromfile(file=f,dtype=np.double,count=NZ)
    fnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
    efnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
    
    f.close()
    
    tempfilt  = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
                 'tempfilt':tempfilt,'lc':lc,'zgrid':zgrid,'fnu':fnu,'efnu':efnu}
    
    ###### .coeff
    f = open(root+'.coeff','rb')
    
    s = np.fromfile(file=f,dtype=np.int32, count=4)
    NFILT=s[0]
    NTEMP=s[1]
    NZ=s[2]
    NOBJ=s[3]
    coeffs = np.fromfile(file=f,dtype=np.double,count=NTEMP*NOBJ).reshape((NOBJ,NTEMP)).transpose()
    izbest = np.fromfile(file=f,dtype=np.int32,count=NOBJ)
    tnorm = np.fromfile(file=f,dtype=np.double,count=NTEMP)
    
    f.close()
    
    coeffs = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
              'coeffs':coeffs,'izbest':izbest,'tnorm':tnorm}
              
    ###### .temp_sed
    f = open(root+'.temp_sed','rb')
    s = np.fromfile(file=f,dtype=np.int32, count=3)
    NTEMP=s[0]
    NTEMPL=s[1]
    NZ=s[2]
    templam = np.fromfile(file=f,dtype=np.double,count=NTEMPL)
    temp_seds = np.fromfile(file=f,dtype=np.double,count=NTEMPL*NTEMP).reshape((NTEMP,NTEMPL)).transpose()
    da = np.fromfile(file=f,dtype=np.double,count=NZ)
    db = np.fromfile(file=f,dtype=np.double,count=NZ)
    
    f.close()
    
    temp_sed = {'NTEMP':NTEMP,'NTEMPL':NTEMPL,'NZ':NZ,\
              'templam':templam,'temp_seds':temp_seds,'da':da,'db':db}
              
    ###### .pz
    if os.path.exists(root+'.pz'):
        f = open(root+'.pz','rb')
        s = np.fromfile(file=f,dtype=np.int32, count=2)
        NZ=s[0]
        NOBJ=s[1]
        chi2fit = np.fromfile(file=f,dtype=np.double,count=NZ*NOBJ).reshape((NOBJ,NZ)).transpose()

        ### This will break if APPLY_PRIOR No
        s = np.fromfile(file=f,dtype=np.int32, count=1)
        
        if len(s) > 0:
            NK = s[0]
            kbins = np.fromfile(file=f,dtype=np.double,count=NK)
            priorzk = np.fromfile(file=f, dtype=np.double, count=NZ*NK).reshape((NK,NZ)).transpose()
            kidx = np.fromfile(file=f,dtype=np.int32,count=NOBJ)
            pz = {'NZ':NZ,'NOBJ':NOBJ,'NK':NK, 'chi2fit':chi2fit, 'kbins':kbins, 'priorzk':priorzk,'kidx':kidx}
        else:
            pz = None
        
        f.close()
        
    else:
        pz = None
    
    if False:
        f = open(root+'.zbin','rb')
        s = np.fromfile(file=f,dtype=np.int32, count=1)
        NOBJ=s[0]
        z_a = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_p = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_m1 = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_m2 = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_peak = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        f.close()
        
    ###### Done.    
    return tempfilt, coeffs, temp_sed, pz


def generate_sed_arrays(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same'):
    """
    Generate full "obs_sed" and "temp_sed" arrays from the stored
    EAZY binary files.
    
    Returns
    -------
    tempfilt : dict
        Dictionary read from the `tempfilt` file.  Keys include the input 
        photometry ('fnu', 'efnu') and the filter central wavelengths
        ('lc').
    
    z_grid : array
        Redshift on the input grid nearest to the best output photo-z.  The
        SEDs are generated here based on the fit coefficients.
    
    obs_sed : (NFILT, NOBJ) array
        Template photometry.  Has fnu units with the same scaling as the 
        input photometry (i.e., AB with zeropoint `PRIOR_ABZP`).
    
    templam : (NTEMPL) array
        Rest-frame wavelengths of the full template spectra
    
    temp_sed : (NTEMPL, NOBJ) array
        Full best-fit template spectra in same units as `obs_sed`.
        
    """
    
    import os
    #import threedhst.eazyPy as eazy
    
    out = readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE=CACHE_FILE)
    tempfilt, coeffs, temp_seds, pz = out
    
    # The redshift gridpoint where the SEDs are evaluated
    z_grid = tempfilt['zgrid'][coeffs['izbest']]
    
    obs_sed = np.zeros_like(tempfilt['fnu'])
    for i in range(tempfilt['NOBJ']):
        obs_sed[:,i] = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][i]],
                         coeffs['coeffs'][:,i])
    
    
    ### Temp_sed
    temp_sed = (np.dot(temp_seds['temp_seds'],coeffs['coeffs']).T*(temp_seds['templam']/5500.)**2).T
    
    ## IGM
    lim1 = np.where(temp_seds['templam'] < 912)
    temp_sed[lim1,:] *= 0
    
    lim2 = np.where((temp_seds['templam'] >= 912) & (temp_seds['templam'] < 1026))
    db = 1.-temp_seds['db'][coeffs['izbest']]
    temp_sed[lim2,:] *= db
    
    lim3 = np.where((temp_seds['templam'] >= 1026) & (temp_seds['templam'] < 1216))
    da = 1.-temp_seds['da'][coeffs['izbest']]
    temp_sed[lim3,:] *= da
    
    templam = temp_seds['templam']
    
    return tempfilt, z_grid, obs_sed, templam, temp_sed
    
def demo():
    """ 
    Demo on the GOODS-N test data
    """
    import matplotlib.pyplot as plt
    
    tempfilt, z_grid, obs_sed, templam, temp_sed = generate_sed_arrays(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same')
    
    idx = 17
    plt.errorbar(tempfilt['lc'], tempfilt['fnu'][:,idx],
                 tempfilt['efnu'][:,idx], marker='.', color='k',
                  linestyle='None', zorder=1000, label='Data')
                  
    plt.plot(tempfilt['lc'], obs_sed[:,idx], marker='.', color='r', alpha=0.8,
             label='Template photometry')
             
    plt.plot(templam*(1+z_grid[idx]), temp_sed[:,idx], color='b', alpha=0.8,
             label='Template spectrum')
    
    plt.loglog(); plt.xlim(2000,3.e4)
    plt.legend(loc='lower right')
    plt.xlabel(r'Observed wavelength, $\mathrm{\AA}$')
    plt.ylabel(r'Observed flux, $f_\nu$ @ PRIOR_ABZP')
    plt.tight_layout()
    