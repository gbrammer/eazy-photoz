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