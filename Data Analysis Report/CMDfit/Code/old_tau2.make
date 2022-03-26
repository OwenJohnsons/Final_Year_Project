old_tau2:
	$(FC) -o ${CMDEXE}/$@ quadint.f90 cmdfit_system.f90 likelihood.f90 mass_functions.f90 old_tau2.f90  $(CMDOBJ)/random.a $(CMDOBJ)/system.a $(CMDOBJ)/define_star.o
	rm *.mod
