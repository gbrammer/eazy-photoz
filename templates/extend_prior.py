import numpy as np
import scipy
import matplotlib.pyplot as plt

"""
Use the prior coeff file to extrapoate the eazy prior to arbitrarily high redshift and faint magnitudes.

The coefficients are for Eq. 3 of Brammer et al. (2008).

p(z, mag) ~ z**gamma  exp((-z/z0)**gamma)
"""

def do_both():
    """
    Extend the R and K priors
    """
    go(mags=(15,30.1,0.5), redshifts=(0.01, 12, 0.01), coeff_file='prior_R_zmax7_coeff.dat', outfile='prior_R_extend.dat')
    go(mags=(15,30.1,0.5), redshifts=(0.01, 12, 0.01), coeff_file='prior_K_zmax7_coeff.dat', outfile='prior_K_extend.dat')
    
def go(mags=(15,30.1,0.5), redshifts=(0.01, 12, 0.01), coeff_file='prior_K_zmax7_coeff.dat', outfile='prior_K_extend.dat'):
    
    fp = open(coeff_file)
    lines = fp.readlines()
    fp.close()
    
    mag_list = np.cast[float](lines[0].split()[2:])
    z0 = np.cast[float](lines[1].split()[1:])
    gamma = np.cast[float](lines[2].split()[1:])
    
    z_grid = np.arange(redshifts[0], redshifts[1], redshifts[2])
    NZ = z_grid.shape[0]
    
    mag_grid = np.arange(mags[0], mags[1], mags[2])
    NM = mag_grid.shape[0]
    
    #### Polynomial extrapolation not reliable
    # p_z0 = scipy.polyfit(mag_list, z0, order) 
    # z0_grid = np.maximum(scipy.polyval(p_z0, mag_grid), 0.05)
    # p_gamma = scipy.polyfit(mag_list, gamma, order) 
    # gamma_grid = np.maximum(scipy.polyval(p_gamma, mag_grid), 0.05)
    
    #### Interpolations on the defined grid
    z0_grid = np.interp(mag_grid, mag_list, z0)
    gamma_grid = np.interp(mag_grid, mag_list, gamma)
    
    #### Linear extrapolations of fit coefficients
    p_z0 = scipy.polyfit(mag_list[:3], z0[:3], 1)
    z0_grid[mag_grid < mag_list[0]] = np.maximum(scipy.polyval(p_z0, mag_grid[mag_grid < mag_list[0]]), 0.05)
    p_z0 = scipy.polyfit(mag_list[-3:], z0[-3:], 1)
    z0_grid[mag_grid > mag_list[-1]] = np.maximum(scipy.polyval(p_z0, mag_grid[mag_grid > mag_list[-1]]), 0.05)

    p_gamma = scipy.polyfit(mag_list[:3], gamma[:3], 1)
    gamma_grid[mag_grid < mag_list[0]] = np.maximum(scipy.polyval(p_gamma, mag_grid[mag_grid < mag_list[0]]), 0.05)
    p_gamma = scipy.polyfit(mag_list[-3:], gamma[-3:], 1)
    gamma_grid[mag_grid > mag_list[-1]] = np.maximum(scipy.polyval(p_gamma, mag_grid[mag_grid > mag_list[-1]]), 0.05)
    
    out_matrix = np.zeros((NZ, NM+1))
    out_matrix[:,0] = z_grid
    
    for i in range(NM):
        pz = z_grid**gamma_grid[i] * np.exp(-(z_grid/z0_grid[i])**gamma_grid[i])
        pz /= np.trapz(pz, z_grid)
        plt.plot(z_grid, pz, label=mag_grid[i])
        out_matrix[:,i+1] = pz
    
    plt.legend(ncol=3, loc='upper right')
    
    header = 'z '
    for m in mag_grid: 
        header += '%6.1f' %(m)
    
    # fp = open(outfile,'w')
    # fp.write(header+'\n')
    # np.savetxt(fp, out_matrix, fmt='%6.3e')#['%6.3e']*out_matrix.shape[1])
    # fp.close()
    
    np.savetxt(outfile, out_matrix, fmt='%6.3e', header=header)
    