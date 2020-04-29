"""
Generate JWST filter throughput curves with Pandeia
"""
def go():
    import os
    import pysynphot as S
    import numpy as np
    import matplotlib.pyplot as plt
    
    import pandeia.engine as etc
    import pandeia.engine.jwst
    import pandeia.engine.calc_utils
    import pandeia.engine.etc3D
    import pandeia.engine.source
    import pandeia.engine.observation
    
    nis = etc.jwst.NIRISS(mode='imaging', config={})
    
    instruments = [etc.jwst.NIRCam, etc.jwst.NIRISS, etc.jwst.MIRI]
        
    modes = [[etc.jwst.NIRISS, 'imaging'], [etc.jwst.NIRCam, 'sw_imaging'], [etc.jwst.NIRCam, 'lw_imaging'], [etc.jwst.MIRI, 'imaging']]

    #modes = [[etc.jwst.NIRCam, 'sw_imaging'], [etc.jwst.NIRCam, 'lw_imaging'], [etc.jwst.MIRI, 'imaging']]

    src = etc.calc_utils.build_default_source()
    #src = etc.source.Source(src_dict)
    src['spectrum']['normalization']['norm_flux'] = 25.0
    src['spectrum']['normalization']['norm_fluxunit'] = 'abmag'
    
    jw = etc.jwst.JWST()
    
    #ote_int = jw.get_ote_eff(wave)
    coll_area = jw.coll_area
    
    ote_file = os.path.join(jw.ref_dir, jw.paths['otefile'])
    ote_bp = S.FileBandpass(ote_file)
    
    filters = {(pandeia.engine.jwst.NIRCam, 'sw_imaging'): ['f070w', 'f090w', 'f115w', 'f150w', 'f200w', 'f150w2', 'f140m', 'f162m', 'f182m', 'f210m', 'f164n', 'f187n', 'f212n'], 
               (pandeia.engine.jwst.NIRCam, 'lw_imaging'): ['f277w', 'f356w', 'f444w', 'f322w2', 'f250m', 'f300m', 'f335m', 'f360m', 'f410m', 'f430m', 'f460m', 'f480m', 'f323n', 'f405n', 'f466n', 'f470n'],
               (pandeia.engine.jwst.NIRISS, 'imaging'): ['f090w', 'f115w', 'f150w', 'f200w', 'f140m', 'f158m', 'f277w', 'f356w', 'f444w', 'f380m', 'f430m', 'f480m'],
               (pandeia.engine.jwst.MIRI, 'imaging'): ['f1065c', 'f1140c', 'f1550c', 'f2300c', 'f560w', 'f770w', 'f1000w', 'f1130w', 'f1280w', 'f1500w', 'f1800w', 'f2100w', 'f2550w']}
    
    fp_info = open('jwst.filters.info','w')
    fp_res = open('jwst.filters.res','w')
    
    # Next number in master FILTER.RES file
    count = 330
    count = 350
    
    today = '030818'
    today = '290419'
    
    for mode in modes:
        ins = mode[0](mode=mode[1])
        ins_str = ins.as_dict()['instrument']['instrument']
        #filters = ins.filters
        for filter in filters[tuple(mode)]:
            label='{0} {1} {2}'.format(mode[0], mode[1], filter)
            
            ins = mode[0](mode=mode[1], config={'instrument':{'filter':filter}})
            
            filter_file = os.path.join(ins.ref_dir, ins.paths[filter])
            filter_bp = S.FileBandpass(filter_file)
            wave = filter_bp.wave
            
            # obs = etc.observation.Observation(instrument=ins)
            # obs.background = 'low'
            # obs.scene.sources[0].spectrum = src['spectrum']
            # det = etc.etc3D.DetectorSignal(observation=obs)
            
            # Manual combine components
            if False:
                ote_int = jw.get_ote_eff(wave)#*coll_area
                filter_eff = ins.get_filter_eff(wave)
                disperser_eff = ins.get_disperser_eff(wave)
                internal_eff = ins.get_internal_eff(wave)
                qe = ins.get_detector_qe(wave)
                #q_yield, fano_factor = ins.get_quantum_yield(wave)
            
                #total_filter = ote_int*filter_eff*internal_eff*qe
            
            total_filter = ins.get_total_eff(wave)
            total_filter *= ins.get_quantum_yield(wave)[0]
                        
            #print(label)
            
            # Resample to R=300
            if filter.endswith('n'):
                R = 1000
            else:
                R = 300
                
            w = np.exp(np.arange(np.log(wave[0]), np.log(wave[-1]), 1./R))
            bp = S.ArrayBandpass(wave*1.e4, total_filter).resample(w*1.e4)

            if total_filter.max() < 0.02:
                continue

            plt.plot(bp.wave, bp.throughput, label=label)
            
            N = len(bp.wave)
            filter_tag = os.path.basename(filter_file).split('.fits')[0]
            filter_label = '{0} lambda_c= {1:.4e}'.format(filter_tag, bp.pivot())
            
            #label_res = '{0:>5d} {1}/{2}/{3}/{4} lambda_c= {5:.4e}'.format(N, ins_str, mode[1], filter, today, bp.pivot())
            #label_info = '{0:>5d} {1}/{2}/{3}/{4} lambda_c= {5:.4e}'.format(count, ins_str, mode[1], filter, today, bp.pivot())
            label_res = '{0:>5d} {1}'.format(N, filter_label)
            label_info = '{0:>5d} {1:>5} {2}'.format(count, N, filter_label)
            
            print(label_info)
            
            fp_res.write(label_res+'\n')
            for i in range(N):
                fp_res.write('{0:<5d} {1:.5e} {2:.5e}\n'.format(i+1, bp.wave[i], bp.throughput[i]))


            fp_info.write(label_info+'\n')
            
            count += 1
            
    fp_info.close()
    fp_res.close()
    
            
            