FORTRAN =  get_header.f90 put_header.f90 addon.f ark_file_io.f90 blknam.f find_free.f locase.f ltcheb.f ripoly.f setbug.f setdt.f upcase.f

system:
	$(FC) -c $(FORTRAN)
	ar rs $(CMDOBJ)/system.a *.o
	rm *.o
	mv *.mod $(CMDMOD)
