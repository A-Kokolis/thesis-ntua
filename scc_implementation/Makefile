CC=gcc
CFLAGS=-O3 -lm -Wall

all: main

main: main.c Voltage_Check.c
	$(CC) $(CFLAGS) main.c Voltage_Check.c -o main

clean:
	rm main
