;;; Read PEGASE 2.0 spectra file into struct
;;;
;;; GBB 3/5/07

;;;; 'infile' is a spectrum file output from the PEGASE.2 code.
;;;;
;;;; The function returns a structure with all of the model parameters
;;;; and SEDs at every age in the grid.
;;;;
;;;; EXAMPLE:    IDL> spec = readpeg('spectra1.dat')
;;;;             IDL> help,spec,/struct
;;;;
;;;; You can add emission lines as specified in the PEGASE.2 output
;;;; using the add_peg_lines function below.  Continuing the example
;;;; above, to add lines with a velocity dispersion of 300 km/s, 
;;;;
;;;;             IDL> spec_w_lines = add_peg_lines(spec,vdisp=300)
;;;;
;;;; [requires IDL6.X and strint() readlines() split(), provided below]
function readpeg,infile,savefile=savefile

    if keyword_set(infile) eq 0 then infile = 'spectra1.dat'
    
    nper = 5 ;; items per line
    
    lines = readlines(infile)

    nlines = n_elements(lines)
    
    start=0
    i=-1
    while start eq 0 do if strmid(lines[++i],0,4) eq '****' then start=i+1

    pars = split(lines[start])
    nstep = pars[0]
    nwave = pars[1]
    nlines=pars[2]
    
    wave = fltarr(nwave)
    wave_lines = fltarr(nlines)
    
    spectra = fltarr(nstep,nwave)
    spectra_lines = fltarr(nstep,nlines)
    
    wblock = fix(nwave/nper)
    lblock = fix(nlines/nper)
    
    ;;;;;;;;; Wavelength block
    for i=0,wblock-1 do begin
        j = start+1+i
        wave[i*nper:i*nper+nper-1] = split(lines[j])
    endfor    
    nlast = nwave mod nper
    if (nlast gt 0) then begin
        ++j
        last = split(lines[j])
        wave[nwave-nlast:*] = last
    endif
    
    start = j
    ;;;;;;;; Line wavelengths
    for i=0,lblock-1 do begin
        j = start+1+i
        wave_lines[i*nper:i*nper+nper-1] = split(lines[j])
    endfor    
    nlast = nlines mod nper
    if (nlast gt 0) then begin
        ++j
        last = split(lines[j])
        wave_lines[nlines-nlast:*] = last
    endif

    time = intarr(nstep)
    mgal = fltarr(nstep)
    mstar = fltarr(nstep)
    mwd = fltarr(nstep)
    mbhns = fltarr(nstep)
    msub = fltarr(nstep)
    mgas = fltarr(nstep)
    zgas = fltarr(nstep)
    zstarmass = fltarr(nstep)
    zstarlbol = fltarr(nstep)
    
    lbol = fltarr(nstep)
    tauv = fltarr(nstep)
    ldust_lbol = fltarr(nstep)
    sfr = fltarr(nstep)
    nSNII = fltarr(nstep)
    nSNIa = fltarr(nstep)
    tstarmass = fltarr(nstep)
    tstarlbol = fltarr(nstep)
    
    ;;;;;;;; Grid blocks    
    for step=0,nstep-1 do begin
        print,'Timestep '+strint(step+1)+' of '+strint(nstep)
        
        start=++j
        
        ;;; block 1
        pars = split(lines[j])

        time[step] = pars[0]
        mgal[step] = pars[1]
        mstar[step] = pars[2]
        mwd[step] = pars[3]
        mbhns[step] = pars[4]
        msub[step] = pars[5]
        mgas[step] = pars[6]
        zgas[step] = pars[7]
        zstarmass[step] = pars[8]
        zstarlbol[step] = pars[9]
        
        ;;; block 2
        pars = split(lines[++j])

        lbol[step] = pars[0]
        tauv[step] = pars[1]
        ldust_lbol[step] = pars[2]
        sfr[step] = pars[3]
        nSNII[step] = pars[4]
        nSNIa[step] = pars[5]
        tstarmass[step] = pars[6]
        tstarlbol[step] = pars[7]

        ;;;;;;;;; Spectra block
        start=++j
        for i=0,wblock-1 do begin
            j = start+i
            spectra[step,i*nper:i*nper+nper-1] = split(lines[j])
        endfor    
        nlast = nwave mod nper
        if (nlast gt 0) then begin
            ++j
            last = split(lines[j])
            spectra[step,nwave-nlast:*] = last
        endif
    
        start = ++j
        ;;;;;;;; Line fluxes
        for i=0,lblock-1 do begin
            j = start+i
            spectra_lines[step,i*nper:i*nper+nper-1] = split(lines[j])
        endfor    
        nlast = nlines mod nper
        if (nlast gt 0) then begin
            ++j
            last = split(lines[j])
            spectra_lines[step,nlines-nlast:*] = last
        endif
    endfor
    
    pars = { nstep : nstep, nwave:nwave,nlines:nlines, $
        wave: wave, wave_lines: wave_lines, spectra: spectra,$
        spectra_lines: spectra_lines, $
        time : time , $
        mgal : mgal ,$
        mstar : mstar ,$
        mwd : mwd ,$
        mbhns : mbhns ,$
        msub : msub ,$
        mgas : mgas ,$
        zgas : zgas ,$
        zstarmass : zstarmass ,$
        zstarlbol : zstarlbol ,$
        lbol : lbol ,$
        tauv : tauv ,$
        ldust_lbol : ldust_lbol ,$
        sfr : sfr ,$
        nSNII : nSNII ,$
        nSNIa : nSNIa ,$
        tstarmass : tstarmass ,$
        tstarlbol : tstarlbol }

    if keyword_set(savefile) then save,pars,file=infile+'.SAV'

    return,pars

end    

;;; Add emission lines to PEGASE 2.0 spectra 
function add_peg_lines, pstruct, vdisp=vdisp, plotme=plotme

    ;;; pp is struct from readpeg
    if keyword_set(vdisp) eq 0 then vdisp = 300 ;; km/s
    
    pp = create_struct('vdisp',vdisp,pstruct)
    
    xx = findgen(201)/200.*12.d0-6
    gg = 1.d0/sqrt(2*!pi)*exp(-1.*xx^2/2.d0) ;;; normalized gaussian
    gg[where(gg lt 1.e-7)] = 0.0
    
    sig = vdisp/3.e5
    
    for i=0,pp.nstep-1 do begin
        ;print,'Timestep '+strint(i+1)+' of '+strint(pp.nstep)
        
        lint = pp.wave*0.0
        for j=0,pp.nlines-1 do begin
        
            li = pp.wave_lines[j]
            sigi = sig*li
            
            lint += interpol(gg/sigi*pp.spectra_lines[i,j],xx*sigi+li,pp.wave)
            
        endfor
    
        ;oplot,pp.wave,pp.spectra[i,*];,thick=!p.thick*4>4
        if keyword_set(plotme) then oplot,pp.wave,pp.spectra[i,*]+lint,col=80,thick=1
        
        pp.spectra[i,*]+=lint
    endfor
    
    return,pp
    
end        


;;;;;; For use in /data/dept/brammer/research/drg/PHOTZ/PEGASE/PEGASE.2/READPEG
pro make_spec

    for i = 1,8 do begin
        pp = readpeg('spectra'+strint(i)+'.dat')
        pp = add_peg_lines(pp,vdisp=300)
    
        age_use = [10,50,100,200,400,700,1000,1600,2000,3000,5000,11000]
        nuse = n_elements(age_use)
        id_use = ainb(age_use,pp.time)   
        
        for j = 0,nuse-1 do $
            forprint,pp.wave, pp.spectra[id_use[j],*],/nocomm,$
             textout='spectra'+strint(i)+'_'+strint(age_use[j])+$
               '.dat'
               
    endfor
    
    age2 = [100,300,600,1000]
    nuse = n_elements(age_use)
    id_use = ainb(age_use,pp.time)   
    
    for i = 9,12 do begin
        pp = readpeg('spectra'+strint(i)+'.dat')
        pp = add_peg_lines(pp,vdisp=300)
        
        for j=0,nuse-1 do $
            if age_use[j] ge age2[i-9] then $
                forprint,pp.wave, pp.spectra[id_use[j],*],/nocomm,$
                  textout='spectra'+strint(i)+'_'+strint(age_use[j])+$
                   '.dat'
            
    endfor
    
    
end        

;;;;;;;;;;;;;;        
;;;;;;;;;;;;;; REQUIRED UTILITIES
;;;;;;;;;;;;;;

function split,instring,delim,keep_duplicates=keep_duplicates,count=count

    if keyword_set(delim) eq 0 then delim = ' '
    
    i = 0
    len = strlen(instring)
    dlen = strlen(delim)
    
    outarr = ['remove_later']
    count=0
    
    while i lt len do begin
        sub = strmid(instring,i,len-i)
        idx_1 = strpos(sub,delim)
        ; print, 'first', sub,i,idx_1
        
        if idx_1 eq -1 then begin
            outarr = [outarr,sub]
            count=count+1
            i=len+1
        endif else begin
        
            if keyword_set(keep_duplicates) eq 0 then begin
                if idx_1 ne 0 then begin
                    outarr = [outarr,strmid(instring,i,idx_1)]
                    count=count+1
                endif
            endif else begin
                outarr = [outarr,strmid(instring,i,idx_1)]
                count=count+1
            endelse

            ; print, 'second', sub,i,idx_1
                
            i = idx_1+i+dlen
                     
        endelse
    endwhile
    
    ; print, 'OUTPUT:'
    ; for i = 1,count do help,outarr[i]
    
    return,outarr[1:count]
    
end
                
;;;; Convert int to string with no white space
function strint,int,zeros=zeros,zfill=zfill

    
    if keyword_set(zeros) then begin
        
        if keyword_set(zfill) eq 0 then zfill='0'
        
        zmax = fix(alog10(zeros))
        zdat = fix(alog10(int))
        outstr=string(int,format='(i0)')
        if zdat lt zmax then for i = 0,zmax-zdat-1 do outstr=zfill+outstr
        
    endif else outstr = string(int,format='(i0)')
        
    return,outstr
    
end

;; Read lines of a file into string array
;; requires strint split GBB 2/26/07
function readlines,file,nchar=nchar

    if keyword_set(nchar) eq 0 then nchar=512
    fmt = '(a'+strint(nchar)+')'
    spawn,'wc '+file,wc
    nlines = long((split(wc))[0])
    line=''
    lines=strarr(nlines)

    openr,lun,file,/get_lun

    for i = 0L,nlines-1 do begin
        readf,lun,line,format=fmt
        lines[i] = line
    endfor

    close,lun
    
    return,lines
end        
    
