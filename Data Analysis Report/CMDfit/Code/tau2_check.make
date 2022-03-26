tau2_check:
	$(FCC) -o ${CMDEXE}/$@ quadint.f90 cmdfit_system.f90 mass_functions.f90  likelihood.f90 tau2.f90 $(CMDOBJ)/random.a $(CMDOBJ)/system.a $(CMDOBJ)/define_star.o
	rm *.mod
