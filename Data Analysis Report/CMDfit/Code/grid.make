grid:
	$(FC) -o ${CMDEXE}/$@ cmdfit_system.f90 quadint.f90 bolcor_subs.f90 mass_functions.f90 colteff_subs.f90 bell_ext.f90 likelihood.f90 cmdfit_subs.f90 grid.f90 $(CMDOBJ)/system.a $(CMDOBJ)/define_star.o
	rm *.mod
