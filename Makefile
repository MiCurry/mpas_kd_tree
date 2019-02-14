ifeq ($(FC),gfortran)
	FFLAGS= -fdefault-real-8 -fdefault-double-8
	#FFLAGS= -g 
endif

ifeq ($(FC),ifort)
	FFLAGS= -g -fcheck=all
endif

OMP = -fopenmp

default: all

all: cspeed getoptf mpas_kd_tree runner

cspeed: cspeed.c
	gcc -c cspeed.c

getoptf: getoptf.f90
	$(FC) $(FFLAGS) -c getoptf.f90

mpas_kd_tree: mpas_kd_tree.f90
	$(FC) $(FFLAGS) -c mpas_kd_tree.f90

runner: runner.f90 mpas_kd_tree.o getoptf.o cspeed.o
	$(FC) $(FFLAGS) -o mpas_kd_tests runner.f90 getoptf.o mpas_kd_tree.o cspeed.o

clean:
	rm -rf *.mod *.o mpas_kd_tests
