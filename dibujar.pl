reset
set key left top
#set logscale y
#set format y "10^{%L}"
set xrange[:5.0]
set yrange[:0.5]
set grid

#set xlabel 'Distance to the closest BS'
set xlabel '{/Symbol r} {/Symbol l}_d / {/Symbol l}_s'
#set ylabel 'Transmission success probability' offset 1.5,0
set ylabel 'Throughput' offset 1.5,0

set term postscript eps enhanced color 25
set output 'NOMA_small_R.eps'


#plot '4_SIC_dst.txt' using 1:2 with p pt 5 ps 3 title 'SIC {/Symbol l}_d = 1.0', '4_ALOHA_dst.txt' using 1:2 with p pt 7 ps 3 title 'ALOHA {/Symbol l}_d = 1.0'

#plot '5_SIC_dst.txt' using 1:2 with p pt 5 ps 3 title 'SIC {/Symbol l}_d = 5.0', '4_SIC_dst.txt' using 1:2 with p pt 7 ps 3 title 'SIC {/Symbol l}_d = 4.0', '3_SIC_dst.txt' using 1:2 with p pt 9 ps 3 title 'SIC {/Symbol l}_d = 3.0', '2_SIC_dst.txt' using 1:2 with p pt 11 ps 3 title 'SIC {/Symbol l}_d = 2.0', '1_SIC_dst.txt' using 1:2 with p pt 13 ps 3 title 'SIC {/Symbol l}_d = 1.0'

#plot 'NOMA_thp_L5.txt' using 1:2 with p pt 5 ps 3 title 'NOMA Level 5'

#plot '1_PA_L1.txt' using 1:2 with p pt 5 ps 3 title 'NOMA Level 1', '1_PA_L3.txt' using 1:2 with p pt 7 ps 3 title 'NOMA Level 3', '1_PA_L5.txt' using 1:2 with p pt 9 ps 3 title 'NOMA Level 5'

#plot 'PA_thp.txt' using 1:2 with p pt 5 ps 3 title 'NOMA Level 1', 'PA_thp.txt' using 1:3 with p pt 7 ps 3 title 'NOMA Level 3', 'PA_thp.txt' using 1:4 with p pt 9 ps 3 title 'NOMA Level 5',

plot 'NOMA_small_R.txt' using 1:2 with p pt 5 ps 3 title 'NOMA Level 1', 'NOMA_small_R.txt' using 1:3 with p pt 7 ps 3 title 'NOMA Level 3', 'NOMA_small_R.txt' using 1:4 with p pt 9 ps 3 title 'NOMA Level 5',

#plot 'PA_thp.txt' using 1:4 with p pt 11 ps 3 title 'NOMA channel 1', 'NOMA_ch.txt' using 1:2 with p pt 5 ps 3 title 'NOMA channel 5', 'NOMA_ch.txt' using 1:3 with p pt 7 ps 3 title 'NOMA channel 10', 'NOMA_ch.txt' using 1:4 with p pt 9 ps 3 title 'NOMA channel 20',

