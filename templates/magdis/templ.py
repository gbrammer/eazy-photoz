"""
Fit Georgios' templates with the models from 
https://www.ias.u-psud.fr/DUSTEM/grids.html


http://georgiosmagdis.pbworks.com/w/page/59019974/SED%20Templates

U1 :0-0.025
U2: 0.05-0.275
U3: 0.3-0.6250
U4:0.65 -0.975
U5: 1.0-1.30
U6: 1.325-1.725
U7: 1.75-2.25
U8: 2.27-3.0
(9 z > 3)
(10 GN20 starburst)

"""

def fit_gm():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import nnls
    
    lnorm = 20
    
    nwave = 800
    
    # files = glob.glob('SED_J13/*RES')     
    # files = glob.glob('SED_C11/*RES')     
    # #files = glob.glob('SED_DL07/*RES')     
    # 
    # du = np.zeros((nwave, len(files)))
    # 
    # for i, file in enumerate(files): 
    #     dust = np.loadtxt(file, skiprows=8)
    #     w, f = dust[:,0], dust[:,1]
    #     print(file, len(w))
    #     plt.plot(w, f/np.interp(lnorm, w, f), label=file, alpha=0.5)
    #     du[:,i] = f/np.interp(lnorm, w, f)
    #     du_wave = w
    
    # Magdis
    mag = np.loadtxt('ms.txt')
    mag_wave = mag[:,0]
    mag_flux = mag[:,1:]
    for i in range(10):
        mag_flux[:,i] /= np.interp(lnorm, mag_wave, mag_flux[:,i])
        #plt.plot(mag_wave, mag_flux[:,i], color='k', alpha=0.8)
                
    ### By hand, blackbodies following da Cunha + Magphys
    from astropy.modeling.physical_models import BlackBody
    import astropy.units as u
    from astropy.constants import c
    
    # Eazy for wavelength grid
    ez = np.loadtxt('templates/fsps_full/fsps_QSF_12_v3_001.dat')
    wave_grid = ez[:,0]/1.e4
    nwave = len(wave_grid)
    
    # PAH template.  C11 dies down quickly
    file = 'SED_C11/SED_C11_100.RES'
    #file = 'SED_DL07/SED_DL07_100.RES'
    dust = np.loadtxt(file, skiprows=8)
    du_wave, du_pah = dust[:,0], dust[:,1]
        
    comps = [np.interp(wave_grid, du_wave, du_pah)]
    
    nu = (c/(wave_grid*u.um)).to(u.Hz)
    
    # equilibrium components, da Cunha
    # modified black-bodies, extra factor of 1 turns from Fnu to nu Fnu
    # cold
    for t in np.arange(20,40):
        comps.append(BlackBody(temperature=t*u.K)(wave_grid*u.um)*nu**(1+2.0))
    
    #warm 
    for t in np.arange(30,80):
        comps.append(BlackBody(temperature=t*u.K)(wave_grid*u.um)*nu**(1+1.5))

    # Hot
    for t in [130, 250] :
        comps.append(BlackBody(temperature=t*u.K)(wave_grid*u.um)*nu**(1+1))
    
    _A = np.array(comps).T
    
    nc = _A.shape[1]
    for i in range(nc):
        _A[:,i] /= np.interp(lnorm, wave_grid, _A[:,i])
    
    #######
    mag_int = np.zeros((nwave, 10))
    clip = (wave_grid > 4.) & (wave_grid < 5000)
    
    models_flam = mag_int*0.
    
    for i in range(10):
        mag_int[:,i] = np.interp(wave_grid, mag_wave, mag_flux[:,i])
    
        _a = nnls(_A[clip,:], mag_int[clip,i])
        model = _A.dot(_a[0])

        norm = np.trapz(model/wave_grid, wave_grid)
        
        pl = plt.plot(wave_grid, mag_int[:,i]/norm, linewidth=4, alpha=0.2)
        plt.plot(wave_grid[clip], mag_int[clip,i]/norm, linewidth=4, alpha=0.8, color=pl[0].get_color())
        
        plt.plot(wave_grid, model/norm, linewidth=3, color='w', alpha=0.8)
        plt.plot(wave_grid, model/norm, linewidth=1, color=pl[0].get_color(), alpha=0.8)
        
        mflam = model/wave_grid
        mflam /= np.trapz(mflam, wave_grid)
        
        models_flam[:,i] = mflam
        
        fp = open(f'templates/magdis/magdis_{i+1:02d}.txt','w')
        np.savetxt(fp, np.array([wave_grid*1.e4, mflam]).T, header='wave flam\n wave: A, flux: flam \nnormalized to unit energy', fmt='%.5e')     
        fp.close()
        
        # components
        # plt.plot(wave_grid, (_A*_a[0]), linewidth=1, color='k', alpha=0.2)
        
        
    plt.loglog()
    plt.xlim(0.2, 5000)
    plt.ylim(1.e-6, 10)
    plt.grid()
    plt.savefig('templates/magdis/magdis_fit.png')
        
        
