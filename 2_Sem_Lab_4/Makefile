CC = gcc
CFLAGS = -Wall -O3 -lm
TARGET = bvp_solver

all: $(TARGET)

$(TARGET): shooting.c newton.c
	$(CC) $(CFLAGS) -o $(TARGET) shooting.c newton.c -lm

clean:
	rm -f $(TARGET)

.PHONY: all clean