docker build -t coupled_lb_p .
docker rm lbp
docker run --name lbp coupled_lb_p

mkdir -p output
docker cp lbp:/usr/src/coupled-LaplaceBeltrami-Poisson/laplacebeltrami_grid-*.gnuplot ./output/
docker cp lbp:/usr/src/coupled-LaplaceBeltrami-Poisson/laplacebeltrami_solution*d-cycle_*.vtk ./output/
docker cp lbp:/usr/src/coupled-LaplaceBeltrami-Poisson/poisson_grid-*.gnuplot ./output/
docker cp lbp:/usr/src/coupled-LaplaceBeltrami-Poisson/poisson_solution*d-cycle_*.vtk ./output/