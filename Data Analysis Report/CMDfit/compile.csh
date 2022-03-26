#!/bin/csh

rm $CMDOBJ/*.a $CMDOBJ/*.o $CMDMOD/*.mod $CMDEXE/*

cd Random
make -f random.make

cd ../System/Unix
make -f system.make
cd ..
make -f system.make

cd ../Cluster
make -f define_star.make

cd ../Code
make -f monte.make
make -f grid.make
make -f for_gaia.make
make -f tau2.make
make -f uncer.make
make -f iso.make
make -f old_tau2.make
make -f ascii2cluster.make
make -f sim.make

cd ../Images
make -f binup.make
make -f subset.make

cd ..
