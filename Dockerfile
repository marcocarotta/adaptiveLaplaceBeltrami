# Use the deal.II image from Docker Hub
FROM dealii/dealii:latest

# Set the working directory in the container
WORKDIR /usr/src/coupled-LaplaceBeltrami-Poisson

# Copy the current directory contents into the container at /usr/src/myapp
COPY . /usr/src/coupled-LaplaceBeltrami-Poisson

# Compile the C++ program
RUN mkdir build && cd build && cmake .. && make

# Run the executable
# CMD ["./build/coupledLBP"]

# # TODO Maybe need copy back to the original folder some output 
# # Copy output files to a specific directory (ensure the directory exists)
# RUN mkdir -p /usr/src/coupled-LaplaceBeltrami-Poisson/output
# CMD cp solution-*.vtk /usr/src/coupled-LaplaceBeltrami-Poisson/output/ && cp poisson_solution-*.vtk /usr/src/coupled-LaplaceBeltrami-Poisson/output/

# Run the executable and copy output files to a specific directory
CMD ./build/coupledLBP && \
	# mkdir -p /usr/src/coupled-LaplaceBeltrami-Poisson/output && \
    pwd && \
    ls 

