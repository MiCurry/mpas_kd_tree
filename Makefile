ifeq ($(FC),gfortran)
	FFLAGS= -g
endif

ifeq ($(FC),ifort)
	FFLAGS= -g
endif

OMP = -fopenmp

default: all

all: getoptf mpas_kd_tree runner

getoptf: getoptf.f90
	$(FC) $(FFLAGS) -c getoptf.f90

mpas_kd_tree: mpas_kd_tree.f90
	$(FC) $(FFLAGS) -c mpas_kd_tree.f90

runner: runner.f90 mpas_kd_tree.o getoptf.o
	$(FC) $(FFLAGS) -o mpas_kd_tests runner.f90 getoptf.o mpas_kd_tree.o

clean:
	rm -rf *.mod *.o mpas_kd_tests
