tau2:
	$(FC) -o ${CMDEXE}/$@ cmdfit_system.f90 likelihood.f90 tau2.f90  $(CMDOBJ)/random.a $(CMDOBJ)/system.a $(CMDOBJ)/define_star.o
	rm *.mod
