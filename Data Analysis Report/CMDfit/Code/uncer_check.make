uncer_check:
	$(FCC) -o ${CMDEXE}/$@ quadint.f90 cmdfit_system.f90 mass_functions.f90 likelihood.f90 uncer.f90 $(CMDOBJ)/system.a
