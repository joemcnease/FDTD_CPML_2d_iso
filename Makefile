CC := gcc
CFLAGS := -lm -Ofast #-O3

SRCS := main.c

all: main

main: main.c
	$(CC) $(SRCS) -o main $(CFLAGS)
