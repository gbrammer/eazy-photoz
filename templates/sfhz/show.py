"""
Show the redshift dependence of the corr_sfhz_13 template set
"""

import eazy
import glob

import numpy as np
import matplotlib.pyplot as plt

files = glob.glob('*[0-9].fits')
files.sort()

res = eazy.filters.FilterFile(path=None)
templ = [eazy.templates.Template(file=file) for file in files]


vband = res[161]

fig, axes = plt.subplots(3,1, sharex=True, sharey=True, figsize=(8, 8))

for i, z in enumerate([0.1, 3., 8.]):
    
    vz = eazy.filters.FilterDefinition(wave=vband.wave*(1+z), throughput=vband.throughput)
    
    for t in templ:
        tnorm = t.integrate_filter(vz, z=z, flam=True)
        igm = t.igm_absorption(z=z)
        axes[i].plot(t.wave/1.e4, t.flux_flam(z=z)/tnorm*igm, alpha=0.5)
        
    axes[i].set_xlim(0.1, 2)
    axes[i].set_ylim(0.008, 60)
    axes[i].grid()
    axes[i].loglog()
    axes[i].set_ylabel(f'flam, z = {z:.1f}')
    
fig.tight_layout(pad=1)

fig.savefig('corr_sfhz_zdependence.png')

# Parameter variation with redshift
from grizli import utils
from astropy.cosmology import WMAP9

templ = utils.read_catalog('corr_sfhz_13_bin0_av0.01.fits')
zstep = []
for k in templ.meta:
    if k.startswith('Z') & (k not in ['ZMET','ZRED']):
        print(k)
        zstep.append(templ.meta[k])

par = utils.read_catalog('corr_sfhz_13.param.fits')
par[r'$\log_{10}$ sSFR'] = np.log10(par['sfr']/par['mass'])
par[r'$M/L_{V}$'] = par['mass']/par['Lv']
par[r'$L_\mathrm{OIII} / L_{V}$'] = par['LOIII']/par['Lv']
par[r'$A_V$'] = par['Av']

keys = ['lwAgeV',r'$A_V$',r'$\log_{10}$ sSFR',r'$M/L_{V}$',
        r'$L_\mathrm{OIII} / L_{V}$']

fig, axes = plt.subplots(len(keys),1,figsize=(8,8), sharex=True)
for i, k in enumerate(keys):
    axes[i].set_ylabel(k)
    axes[i].grid()
    
    for j in range(len(par)-1):
        axes[i].plot(np.log(1+np.array(zstep)), par[k][j], marker='.', alpha=0.5)

for i in [0, 3, 4]:
    axes[i].semilogy()
        
axes[i].set_xlabel('z')
xt = np.arange(0,12.1,2)

axes[i].set_xticklabels(xt)
axes[i].set_xticks(np.log(1+xt))

axes[0].plot(np.log(1+np.array(zstep)), WMAP9.age(zstep), color='0.8', marker='.', label='WMAP9')
axes[0].legend(loc='upper right')

axes[0].set_ylim(0.008, 14)

axes[0].set_yticks([0.01, 0.1, 1])
axes[0].set_yticklabels(['10 Myr', '100', '1 Gyr'])

fig.tight_layout(pad=1)

fig.savefig('corr_sfhz_parameters.png')

