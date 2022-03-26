FORTRAN =  define_star.f90

define_star:
	$(FC) -c $(FORTRAN)
	mv *.o   $(CMDOBJ)
	mv *.mod $(CMDMOD)
