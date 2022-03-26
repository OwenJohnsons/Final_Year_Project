ascii2cluster:
	$(FCC) -o ${CMDEXE}/$@ ascii2cluster.f90 $(CMDOBJ)/system.a $(CMDOBJ)/define_star.o
