#!csh

setenv FC "g95 -fno-second-underscore -fbounds-check -ftrace=full -O0 -I$CMDMOD"

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
make -f sim.make
make -f uncer.make
make -f iso.make

cd ../Images
make -f binup.make
make -f subset.make

cd ..
