CC = gcc
CFLAGS = -Wall -O3
LDFLAGS = -lm
TARGET = wave_eq

all: $(TARGET) run

$(TARGET): wave_eq.c
	$(CC) $(CFLAGS) -o $(TARGET) wave_eq.c $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.dat *.png

.PHONY: all run clean