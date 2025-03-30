CC = gcc
CFLAGS = -Wall -O3
LDFLAGS = -lm  # Важно: библиотека должна быть в конце
TARGET = adams_results
SOURCE = adams_method.c

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)  # -lm в конце команды

run: $(TARGET)
	./$(TARGET)

plot:
	gnuplot -p -e "set terminal pngcairo; \
	               set output 'adams_plot.png'; \
	               set title 'Higgins System (Adams-Bashforth 4th order)'; \
	               plot 'adams_results.dat' using 1:2 with lines title 'y1(t)', \
	                    '' using 1:3 with lines title 'y2(t)'"

clean:
	rm -f $(TARGET) *.dat *.png

.PHONY: all run plot clean