CC = gcc
CFLAGS = -Wall -O3 -lm
TARGET = rosenbrock_solver
PLOT_SCRIPT = plot.gp
DATA_FILE = solution.dat
OUTPUT_IMAGE = solution.png

all: $(TARGET) plot

$(TARGET): rosenbrock_solver.c
	$(CC) $(CFLAGS) -o $@ $< -lm

plot: $(TARGET)
	./$(TARGET)
	gnuplot $(PLOT_SCRIPT)

clean:
	rm -f $(TARGET) $(DATA_FILE) $(PLOT_SCRIPT) $(OUTPUT_IMAGE)

.PHONY: all plot clean