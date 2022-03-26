#!/bin/csh

source $CMDDIR/startup.csh

setenv interior_number `grep -n Geneva $CMDDATA/setup.int | awk -F ":" '{ print $1 }'`

echo "1 16" > monte.in
echo $interior_number >> monte.in
echo 0 >> monte.in
echo V >> monte.in
echo B-V >> monte.in
echo 0 >> monte.in
echo 0.5 >> monte.in
echo 7 7 >> monte.in

monte < monte.in > monte.log

echo "ubv.cat" > grid.in
echo "1" >> grid.in
echo "2" >> grid.in
echo "0.0" >> grid.in
echo "0.0" >> grid.in
echo "7 7" >> grid.in
echo "9.7 10.5" >> grid.in
echo "100" >> grid.in
echo "0.20" >> grid.in
echo "1.0" >> grid.in

grid < grid.in > grid.log

echo "best_model.fit" >> tau2.in
echo "fitted.cat" >> tau2.in
echo "1" >> tau2.in
echo "2" >> tau2.in
echo "1" >> tau2.in

tau2 < tau2.in > tau2.log

uncer > uncer.log
