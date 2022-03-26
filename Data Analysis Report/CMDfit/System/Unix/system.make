FORTRAN = *.f *.f90

system:
	$(FC) -c $(FORTRAN)
	ar rs $(CMDOBJ)/system.a *.o
	rm *.o
