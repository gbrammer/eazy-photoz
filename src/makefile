# -pg -lc_p for profiling
.c.o:
	gcc -msse2 -c -O3 -Wall $<

all: eazy 

eazy: main.o output_files.o zeropoints.o getpar.o makegrid.o interpol.o getphotz.o kmag_prior.o rest_frame_colors.o inoue_igm.o
	gcc -o eazy main.o output_files.o zeropoints.o getpar.o makegrid.o interpol.o getphotz.o kmag_prior.o inoue_igm.o rest_frame_colors.o -lm 

main.o output_files.o zeropoints.o getpar.o makegrid.o interpol.o getphotz.o temperr.o kmag_prior.o rest_frame_colors.o inoue_igm.o: defs.h

temperr: temperr.o  getpar.o
	gcc -o temperr temperr.o getpar.o -lm

clean:
	rm *.o 
        
