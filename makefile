fotosciop: bmp.o maintest.o ip_lib.o
	gcc -Wall --ansi --pedantic bmp.o maintest.o ip_lib.o -ofotosciop -lm

maintest.o: maintest.c
	gcc -c --ansi --pedantic -Wall maintest.c

ip_lib.o: ip_lib.c ip_lib.h
	gcc -c --ansi --pedantic -Wall ip_lib.c

bmp.o: bmp.c bmp.h
	gcc -c -Wall bmp.c

clean:
	@rm -f ip_lib.o maintest.o bmp.o
