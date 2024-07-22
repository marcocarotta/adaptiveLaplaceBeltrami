docker build -t coupled_lb_p .
docker rm lbp
docker run --name lbp coupled_lb_p

mkdir -p output
docker cp lbp:/usr/src/coupled-LaplaceBeltrami-Poisson/laplacebeltrami_grid-0.gnuplot ./output/
docker cp lbp:/usr/src/coupled-LaplaceBeltrami-Poisson/laplacebeltrami_solution2d-cycle_0.vtk ./output/
docker cp lbp:/usr/src/coupled-LaplaceBeltrami-Poisson/poisson_grid-0.gnuplot ./output/
docker cp lbp:/usr/src/coupled-LaplaceBeltrami-Poisson/poisson_solution2d-cycle_0.vtk ./output/