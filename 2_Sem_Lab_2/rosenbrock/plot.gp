set terminal pngcairo enhanced font 'Arial,12'
set output 'solution.png'
set title 'Решение системы ОДУ методом Розенброка-Ваннера 4 порядка (α=1000.0, β=10.0)'
set xlabel 'Время t'
set ylabel 'Значения y1, y2'
plot 'solution.dat' using 1:2 with lines title 'y1(t)', \
     'solution.dat' using 1:3 with lines title 'y2(t)'
