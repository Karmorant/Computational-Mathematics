CC = gcc
CFLAGS = -Wall -O3 -lm
TARGETS = adams2 adams3 adams4
DATAFILES = adams2.dat adams3.dat adams4.dat
PLOTS = adams2_plot.png adams3_plot.png adams4_plot.png

all: $(TARGETS) data plots

data: $(TARGETS)
	./adams2
	./adams3
	./adams4

plots: $(DATAFILES)
	for order in 2 3 4; do \
		gnuplot -e "set terminal pngcairo enhanced; \
		           set output 'adams$${order}_plot.png'; \
		           set title 'Lorenz XY (Adams-$${order})'; \
		           set xlabel 'X'; set ylabel 'Y'; \
		           plot 'adams$${order}.dat' using 2:3 with lines lc rgb 'red' notitle"; \
	done

adams2: adams_methods.c
	$(CC) $(CFLAGS) -DADAMS_ORDER=2 -o $@ $<

adams3: adams_methods.c
	$(CC) $(CFLAGS) -DADAMS_ORDER=3 -o $@ $<

adams4: adams_methods.c
	$(CC) $(CFLAGS) -DADAMS_ORDER=4 -o $@ $<

clean:
	rm -f $(TARGETS) $(DATAFILES) $(PLOTS)

.PHONY: all data plots clean