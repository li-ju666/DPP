all: main.c prop.c 
	mpicc main.c prop.c -o simulation -std=c99 -lm -O3
clean: 
	rm -f main.o prop.o simulation
