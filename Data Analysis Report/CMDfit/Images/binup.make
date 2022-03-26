FORTRAN = binup.f90

LIBRARIES = $(CMDOBJ)/system.a

binup:
	$(FC) -o $(CMDEXE)/$@ $(FORTRAN) $(LIBRARIES)

