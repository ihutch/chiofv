#########################################################################
# This makefile should build on linux provided it has access to the X11 library
# which typically requires the development package, a fortran compiler,
# and git. It won't work on MSWindows.
# Decide the FORTRAN compiler and create the accis graphics routines:
include ACCIS.mk
#########################################################################
LIBRARIES := $(LIBRARIES) -L./fftpack/ -lfftpack5.1
LIBDEPS := $(LIBDEPS) fftpack/libfftpack5.1.a
COMPILE-SWITCHES:=$(COMPILE-SWITCHES) -Wno-unused-dummy-argument  -Wl,--allow-multiple-definition
# -fbounds-check
#########################################################################
OBJECTS=plasmasuscept.o
#########################################################################
% : %.f   # Defeat Fortran Suffix rule f77
# Patterns for compilation etc.
% : %.f makefile $(ACCISX) $(OBJECTS) $(LIBDEPS);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f $(OBJECTS) $(LIBPATH) $(LIBRARIES)

%.o : %.f makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f

% : %.f90  makefile $(ACCISX) $(OBJECTS) $(LIBDEPS);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90 $(OBJECTS) $(LIBPATH) $(LIBRARIES)

%.o : %.f90 makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f90

# Defeat the Modula-2 make booby-trap.
% : %.mod

#########################################################################
# F90 modules and make pattern rules booby-trap.
#
# The heart of the problem: there are lots of suffix rules defined in
# make behind the scenes. One of them is for a Modula-2 program and
# has the effect of a rule of the form % : %.mod. Fortran90 modules
# also use extension .mod so if there's a module in your source file
# 'source.f90' named the same, i.e. 'module source', compiling your
# source will produce a file that is source.mod. Then the Modula-2
# suffix rule, that you knew nothing about, will apply to the
# executable file you are trying to make. If your desired source make
# rule is a pattern rule of the form % : %.f90, as might reasonably be
# the case, it will over-ride the suffix rule ONLY if all the
# prerequisites for it already exist. If one of them happens not to
# exist, because the prerequisite is itself made by some other rule,
# then instead of your rule being used, the Modula-2 rule will be used
# by make. Since the %.mod file was created last time you compiled
# your fortran code, it is older than the executable. The Modula-2
# rule decides therefore that nothing needs to be done, but really it
# does!  This booby-trap can be avoided either by never calling a
# module by the same name as the executable program you are trying to
# make by a pattern rule, or by always defeating the built-in Modula-2
# rule by defining it to be null
#% : %.mod
#

all : chiofv
	./chiofv inputfiles/twogaussLowTe.dat

fftpack/libfftpack5.1.a : fftpack/makefile  fftpack/fftpack5.1.f
	make -C fftpack

clean :
	rm -f *.o *.mod plot*.ps COPTICverif.dat chiofv

mproper :
	make clean
	make -C fftpack clean

