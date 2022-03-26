FORTRAN = subset.f90

LIBRARIES = $(CMDOBJ)/system.a

subset:
	$(FC) -o $(CMDEXE)/$@ $(FORTRAN) $(LIBRARIES)

