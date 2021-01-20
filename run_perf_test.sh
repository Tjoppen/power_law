set -e
g++ power_law_gen.cpp
rm -f perfdata2.csv
for n in 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000
do
  ./a.out $n $n 100 > program.lp
  ./a.out $n $n 100 m > program.m
  wc program.lp
  /usr/bin/time -f "%U" -o time.txt lp_solve -R -S1 -v4 -presolve -presolvel -presolver -presolvem -noint program.lp | tee out.txt
  echo "$n,$(grep nnz program.lp |sed -e 's|// nnz=||'),$(cat time.txt)" >> perfdata2.csv
done
