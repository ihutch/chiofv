ifeq ("$(G77)","")     
        G77=gfortran
endif

# The fno-globals (g77 only) flag is necessary to prevent errors because
# of different precision of workspace types in different calls.
G77flags=-O3
ifeq  ("$(G77)","g77")     
  G77flags=-O3 -fno-globals -Wno-globals
endif
# But the test program doesn't compile anyway because trim is not an
# intrinsic in g77. Gfortran works.
G77flags := $(G77flags) -w

libfftpack5.1.a : fftpack5.1.f
	@echo Be patient. This will take a minute.
	$(G77) -c $(G77flags) fftpack5.1.f
	ar cr libfftpack5.1.a fftpack5.1.o

ffttest : fftpack5.1_prb.f
	gfortran -o ffttest fftpack5.1_prb.f -L. -lfftpack5.1

clean :
	rm -f *.o
	rm -f libfftpack5.1.a 
	rm -f ffttest
