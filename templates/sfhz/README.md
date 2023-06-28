Templates with redshift-dependent SFHs that disfavor SFHs that start earlier 
than the age of the universe at a given epoch.  

<img src=https://github.com/gbrammer/eazy-photoz/blob/master/templates/sfhz/corr_sfhz_zdependence.png width=60%/>

The maximum attenuation of the reddened templates also evolves with redshift.

<img src=https://github.com/gbrammer/eazy-photoz/blob/master/templates/sfhz/corr_sfhz_parameters.png width=60%/>

`4590.fits` is the best-fit template to the extreme emission line galaxy
at z=8.5 with a NIRSpec spectrum, provided by A. Carnall, rescaled to roughly
match the normalization of the FSPS templates.

This template is interpreted as a stellar mass / SFR of zero when computing 
those parameters with eazy, so be wary of those numbers for sources where it
dominates the best-fit.

`j0647agn+torus.fits` is a template generated to replicate the remarkable NIRSpec prism spectrum
of the z=4.50 source shown below (M. Killi, in prep), which is fairly flat in the rest UV,
has strong Ly&alpha;, [OIII] and H&alpha; emission lines, and a red continuum, perhaps consistent
with an obscured AGN torus.  The spectrum is from GO-1433 (D. Coe) in the MACS0647 lensing cluster field.

<img src=https://github.com/gbrammer/eazy-photoz/blob/master/templates/sfhz/macsj0647-v1_prism-clear_1433_1045.fnu.png width=60%/>
