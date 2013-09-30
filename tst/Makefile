#this makefile is for trisrf/trigen

FCOMP = gfortran
FFLAGS =  -O2 -ffixed-line-length-132 -static
BIN = ../bin #relative naming

#rules fortran->object
.f.o:
	$(FCOMP) -c $(FFLAGS) $*.f


all: trigenlarge trisrflarge193 meshsrflarge193 khull
	mv trigenlarge trisrflarge193 meshsrflarge193 khull $(BIN)

#delete object files
clean: 
	rm *.o ../bin/trisrflarge193 ../bin/trigenlarge  ../bin/meshsrflarge193 ../bin/khull

khull: khull.o
	$(FCOMP) $(FFLAGS) -o khull khull.o

trisrflarge193: trisrflarge.o meshsub.o
	$(FCOMP) $(FFLAGS) -o trisrflarge193 trisrflarge.o meshsub.o

trigenlarge: trigenlarge.o meshsub.o
	$(FCOMP) $(FFLAGS) -o trigenlarge trigenlarge.o meshsub.o

meshsrflarge193: meshsrfA.o
	$(FCOMP) $(FFLAGS) -o meshsrflarge193 meshsrfA.o

trisrflarge.o : trisrflarge.h
meshsub.o : tricol.h
meshsrfA.o : meshsrfA.h 

