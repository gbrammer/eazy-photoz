pro graz

    x,2000,8.e4
    y,1.e28,1.e31

    age_use = [10,50,100,200,400,700,1000,1600,2000,3000,5000,11000]
    nuse = n_elements(age_use)

    for i = 1,8 do begin
;        pp = readpeg('graz'+strint(i)+'.dat',/save)
        restore,'graz'+strint(i)+'.dat.SAV'
        pp = pars

        plot,findgen(10),/nodata,/xlog,/ylog,$
            xtit=tid('\lambda'),ytit=tid('f_\lambda'),tit='graz'+strint(i,zero=10)
        
        pp = add_peg_lines(pp,vdisp=300)

        id_use = ainb(age_use,pp.time)

        for j = 0,nuse-1 do $
            forprint,pp.wave, pp.spectra[id_use[j],*],/nocomm,$
             textout='graz'+strint(i,zero=10)+'_'+strint(age_use[j],zero=1.e4)+$
               '.dat'
            
        close,/all

    endfor

    age2 = [100,300,600,1000]
    nuse = n_elements(age_use)

    for i = 9,12 do begin
;        pp = readpeg('graz'+strint(i)+'.dat',/save)
        restore,'graz'+strint(i)+'.dat.SAV'
        pp = pars

        id_use = ainb(age_use,pp.time)
        plot,findgen(10),/nodata,/xlog,/ylog,$
            xtit=tid('\lambda'),ytit=tid('f_\lambda'),tit='graz'+strint(i,zero=10)

        pp = add_peg_lines(pp,vdisp=300)

        for j=0,nuse-1 do $
            if age_use[j] ge age2[i-9] then $
                forprint,pp.wave, pp.spectra[id_use[j],*],/nocomm,$
                  textout='graz'+strint(i,zero=10)+'_'+strint(age_use[j],zero=1.e4)+$
                   '.dat'
            
        close,/all
        
    endfor

end    


pro graz2

    x,800,8.e4
    y,1.e26,1.e31

    age_use = [10,50,100,200,400,700,1000,1600,2000,3000,5000,11000]
    nuse = n_elements(age_use)

    for i = 1,8 do begin
        ;pp = readpeg('graz'+strint(i)+'.dat',/save)
        restore,'graz'+strint(i)+'.dat.SAV'
        pp = pars
        plot,findgen(10),/nodata,/xlog,/ylog,$
            xtit=tid('\lambda'),ytit=tid('f_\lambda'),tit='graz'+strint(i,zero=10)
        
        pp = add_peg_lines(pp,vdisp=300)

        id_use = ainb(age_use,pp.time)

        for j = 0,nuse-1 do $
            oplot,pp.wave, pp.spectra[id_use[j],*],col=j*180./nuse+40,thick=2
            
    endfor

    age2 = [100,300,600,1000]
    nuse = n_elements(age_use)

    for i = 9,12 do begin
        ;pp = readpeg('graz'+strint(i)+'.dat',/save)
        restore,'graz'+strint(i)+'.dat.SAV'
        pp = pars
        id_use = ainb(age_use,pp.time)
        plot,findgen(10),/nodata,/xlog,/ylog,$
            xtit=tid('\lambda'),ytit=tid('f_\lambda'),tit='graz'+strint(i,zero=10)

        pp = add_peg_lines(pp,vdisp=300)

        for j=0,nuse-1 do $
            if age_use[j] ge age2[i-9] then $
                oplot,pp.wave, pp.spectra[id_use[j],*],col=j*180./nuse+40,thick=2
                    
    endfor

end    
