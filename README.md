# coupled Laplace-Beltrami Poisson
The purpose of this project is to solve a coupled Laplace-Beltrami and Poisson problem with a FE method.


# To run 

    docker build -t coupled_lb_p .
    docker run --name lbp coupled_lb_p
    chmod +x copy_output.sh 
    ./copy_output.sh 