set -e
#for v in 100 300 1000 3000 10000 30000
#for v in 300 1000 3000 10000 30000 100000 300000 1000000
for v in 1000000
do
  ./power_law_gen2 $v 160 10 160 10 0.5  && nice octave program_solve2.m
done
