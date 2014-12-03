#CFLAGS = -O0 -g
CFLAGS = -O2 -ffast-math

all : mul_SSE mul_SSE4_1_01 mul_SSE4_1_02

mul_SSE : main_SSE.o matrix_multiply_SSE.o
	gcc -o $@ $^

mul_SSE4_1_01 : main_SSE4_1_01.o matrix_multiply_SSE4_1_01.o
	gcc -o $@ $^

mul_SSE4_1_02 : main_SSE4_1_02.o matrix_multiply_SSE4_1_02.o
	gcc -o $@ $^

.PHONY: clean
clean:
	rm -f *.o

main_SSE.o : main.c
	gcc $(CFLAGS) -DSSE -c -o $@ $^

main_SSE4_1_01.o : main.c
	gcc $(CFLAGS) -DSSE4_1_01 -c -o $@ $^

main_SSE4_1_02.o : main.c
	gcc $(CFLAGS) -DSSE4_1_02 -c -o $@ $^

%.o : %.c
	gcc $(CFLAGS) -c -o $@ $^

%.o : %.S
	gcc $(CFLAGS) -c -o $@ $^
