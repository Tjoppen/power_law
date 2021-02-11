all: power_law_gen power_law_gen2

power_law_gen: power_law_gen.cpp
	g++ -O3 $^ -o $@ -lmatio

power_law_gen2: power_law_gen2.cpp
	g++ -O3 $^ -o $@ -lmatio

