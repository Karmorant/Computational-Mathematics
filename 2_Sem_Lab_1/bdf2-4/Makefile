CC = gcc
CFLAGS = -Wall -O3 -lm
TARGETS = bdf2 bdf3 bdf4
DATAFILES = bdf2.dat bdf3.dat bdf4.dat
PLOTS = bdf2_plot.png bdf3_plot.png bdf4_plot.png

all: $(TARGETS) data plots

data: $(TARGETS)
	./bdf2
	./bdf3
	./bdf4

plots: $(DATAFILES)
	for order in 2 3 4; do \
		gnuplot -e "set terminal pngcairo enhanced; \
		           set output 'bdf$${order}_plot.png'; \
		           set title 'Lorenz XY (BDF-$${order})'; \
		           set xlabel 'X'; set ylabel 'Y'; \
		           plot 'bdf$${order}.dat' using 2:3 with lines lc rgb 'green' notitle"; \
	done

bdf2: bdf_methods.c
	$(CC) $(CFLAGS) -DBDF_ORDER=2 -o $@ $<

bdf3: bdf_methods.c
	$(CC) $(CFLAGS) -DBDF_ORDER=3 -o $@ $<

bdf4: bdf_methods.c
	$(CC) $(CFLAGS) -DBDF_ORDER=4 -o $@ $<

clean:
	rm -f $(TARGETS) $(DATAFILES) $(PLOTS)

.PHONY: all data plots clean