#CFLAGS = -O0 -g
CFLAGS = -O2 -ffast-math

all : mul_SSE mul_SSE4_1

mul_SSE : main_SSE.o matrix_multiply_SSE.o
	gcc -o $@ $^

mul_SSE4_1 : main_SSE4_1.o matrix_multiply_SSE4_1.o
	gcc -o $@ $^

.PHONY: clean
clean:
	rm -f *.o

main_SSE.o : main.c
	gcc $(CFLAGS) -DSSE -c -o $@ $^

main_SSE4_1.o : main.c
	gcc $(CFLAGS) -DSSE4_1 -c -o $@ $^

%.o : %.c
	gcc $(CFLAGS) -c -o $@ $^

%.o : %.S
	gcc $(CFLAGS) -c -o $@ $^
