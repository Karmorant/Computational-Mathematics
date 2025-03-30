CC = gcc
CFLAGS = -Wall -O3
LDFLAGS = -lm
TARGET = method_rk4
SOURCE = irk4_method.c

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

plot:
	gnuplot -e "set terminal png; \
	            set output 'rk4_plot.png'; \
	            plot 'higgins_rk4.dat' using 2:3 with lines title 'y₁(t)', \
	                 '' using 2:3 with lines title 'y₂(t)'"

clean:
	rm -f $(TARGET) *.dat *.png

.PHONY: all run plot clean