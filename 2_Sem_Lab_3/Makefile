CC = gcc
CFLAGS = -Wall -O3 -lm
TARGET = periodic_solver
PLOT_SCRIPT = plot_periodic.gp
DATA_FILE = periodic_solution.dat
OUTPUT_IMAGE = periodic_solution.png

all: $(TARGET) plot

$(TARGET): periodic_solver.c
	$(CC) $(CFLAGS) -o $@ $< -lm

plot: $(TARGET)
	./$(TARGET)
	gnuplot $(PLOT_SCRIPT)

clean:
	rm -f $(TARGET) $(DATA_FILE) $(PLOT_SCRIPT) $(OUTPUT_IMAGE)

.PHONY: all plot clean