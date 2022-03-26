FORTRAN =  *.f90

random:
	$(FC) -c $(FORTRAN)
	ar rs $(CMDOBJ)/random.a *.o
	mv *.mod $(CMDMOD)
	rm *.o
