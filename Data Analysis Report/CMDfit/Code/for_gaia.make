for_gaia:
	$(FC) -o ${CMDEXE}/$@ quadint.f90 for_gaia.f90 $(CMDOBJ)/system.a $(CMDOBJ)/define_star.o
