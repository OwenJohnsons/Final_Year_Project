#!/bin/csh

source $CMDDIR/startup.csh

setenv interior_number `grep -n "Lmix=1.9" $CMDDATA/setup.int | awk -F ":" '{ print $1 }'`

setenv bc_number `grep -n "BT-Settl+Bell  Bessell" $CMDDATA/setup.bc | grep baraffe | awk -F ":" '{ print $1 }'`

echo "0.04 3" > monte.in
echo $interior_number >> monte.in
echo $bc_number >> monte.in
echo "I" >> monte.in
echo "R-I" >> monte.in
echo "0.2" >> monte.in
echo "0.5" >> monte.in
echo "6.3 7.5" >> monte.in
echo "0.01" >> monte.in

monte < monte.in > monte.log

echo "hit.cat" > grid.in
echo "0.03" >> grid.in
echo "0.042" >> grid.in
echo "6.3 7.5" >> grid.in
echo "0.01" >> grid.in
echo "9 11" >> grid.in
echo "200" >> grid.in
echo "0.0" >> grid.in
echo "1" >> grid.in

grid < grid.in > grid.log

echo "best_model.fit" > tau2.in
echo "unclipped.cat" >> tau2.in
echo "2" >> tau2.in

tau2 < tau2.in > tau2.log

uncer > uncer.log
