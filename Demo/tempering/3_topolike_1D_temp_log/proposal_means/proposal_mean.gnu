unset logs
unset grid
set xtics 1
set term pdfcairo dashed enhanced font 'Arial' fontscale 0.5
set xlabel 'iteration'
set ylabel 'alpha_topo'
set output 'proposal_mean_00.pdf'
pl 'proposal_mean_00' u 1:2 w lp lt 0 lw 3 t '', 'proposal_mean_00' u 1:3 w lp lt 1 lw 3 t '', 'proposal_mean_00' u 1:4 w lp lt 2 lw 3 t '', 'proposal_mean_00' u 1:5 w lp lt 3 lw 3 t '', 'proposal_mean_00' u 1:6 w lp lt 4 lw 3 t '', 'proposal_mean_00' u 1:7 w lp lt 5 lw 3 t '', 'proposal_mean_00' u 1:8 w lp lt 6 lw 3 t '', 'proposal_mean_00' u 1:9 w lp lt 7 lw 3 t '', 'proposal_mean_00' u 1:10 w lp lt 8 lw 3 t '', 'proposal_mean_00' u 1:11 w lp lt 9 lw 3 t ''
set output
