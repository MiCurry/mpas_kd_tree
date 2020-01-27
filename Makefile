default:
	@echo "Please provide a compiler name (gnu, intel, pgi)"
	exit 1

gnu:
	( $(MAKE) kdtree \
	 "FC = gfortran" \
	 "FFLAGS = -g -Wall -fcheck=all -pedantic -std=f2003 -fbacktrace")

intel:
	( $(MAKE) kdtree \
	 "FC = ifort " \
	 "FFLAGS = -g -warn all -check all -traceback" )

pgi:
	( $(MAKE) kdtree \
	 "FC = pgfortran" \
	 "FFLAGS = -g -Mbounds -Mchkptr -traceback" )


kdtree:
	$(FC) $(FFLAGS) -c mpas_kd_tree.F90
	$(FC) $(FFLAGS) -o test_kd_tree mpas_kd_tree_tests.F90 mpas_kd_tree.o

clean:
	rm -rf *.o *.mod test_kd_tree
