main_iplib: bmp.o main_iplib.o ip_lib.o
	gcc -Wall --ansi --pedantic bmp.o main_iplib.o ip_lib.o -omain_iplib -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

main_iplib.o: main_iplib.c
	gcc main_iplib.c -c --ansi --pedantic -Wall -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

ip_lib.o: ip_lib.c ip_lib.h
	gcc ip_lib.c   -c --ansi --pedantic -Wall -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

bmp.o: bmp.c bmp.h
	gcc -c -Wall bmp.c

clean:
	@rm -f ip_lib.o main_iplib.o bmp.o
