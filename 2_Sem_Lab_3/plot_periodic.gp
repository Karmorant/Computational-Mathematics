set terminal pngcairo enhanced font 'Arial,12'
set output 'periodic_solution.png'
set title 'Периодическое решение уравнения y'' - (10 + sin(2πx))y = cos(2πx)'
set xlabel 'x'
set ylabel 'y(x)'
set grid
plot 'periodic_solution.dat' with lines title 'Решение'
