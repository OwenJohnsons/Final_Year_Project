monte:
	$(FC) -o ${CMDEXE}/$@ cmdfit_system.f90 quadint.f90 mass_functions.f90 likelihood.f90 bolcor_subs.f90 bell_ext.f90 colteff_subs.f90 cmdfit_subs.f90 monte.f90 $(CMDOBJ)/random.a $(CMDOBJ)/system.a $(CMDOBJ)/define_star.o
	rm *.mod
