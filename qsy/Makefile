CC=mpicc

mc: main.c function.c prop.c
	$(CC) function.c main.c prop.c -o mc -O3 -c99
clean: 
	rm -f *.o mc
