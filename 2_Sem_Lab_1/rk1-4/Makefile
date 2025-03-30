CC = gcc
CFLAGS = -Wall -O3 -lm
TARGETS = rk1 rk2 rk3 rk4
DATAFILES = rk1.dat rk2.dat rk3.dat rk4.dat
PLOTS = rk1_plot.png rk2_plot.png rk3_plot.png rk4_plot.png

all: $(TARGETS) data plots

data: $(TARGETS)
	./rk1
	./rk2
	./rk3
	./rk4

plots: $(DATAFILES)
	for order in 1 2 3 4; do \
		gnuplot -e "set terminal pngcairo enhanced; \
		           set output 'rk$${order}_plot.png'; \
		           set title 'Lorenz XY (RK$${order})'; \
		           set xlabel 'X'; set ylabel 'Y'; \
		           plot 'rk$${order}.dat' using 2:3 with lines lc rgb 'blue' notitle"; \
	done

rk1: rk_methods.c
	$(CC) $(CFLAGS) -DRK_ORDER=1 -o $@ $<

rk2: rk_methods.c
	$(CC) $(CFLAGS) -DRK_ORDER=2 -o $@ $<

rk3: rk_methods.c
	$(CC) $(CFLAGS) -DRK_ORDER=3 -o $@ $<

rk4: rk_methods.c
	$(CC) $(CFLAGS) -DRK_ORDER=4 -o $@ $<

clean:
	rm -f $(TARGETS) $(DATAFILES) $(PLOTS)

.PHONY: all data plots clean