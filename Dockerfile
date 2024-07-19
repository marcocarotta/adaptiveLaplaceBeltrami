# Use the deal.II image from Docker Hub
FROM dealii/dealii:latest

# Set the working directory in the container
WORKDIR /usr/src/coupled-LaplaceBeltrami-Poisson

# Copy the current directory contents into the container at /usr/src/myapp
COPY . /usr/src/coupled-LaplaceBeltrami-Poisson

# Compile the C++ program
RUN mkdir build && cd build && cmake .. && make # seems to be not working

# Run the executable
CMD ["./build/coupledLBP"]

# TODO Maybe need copy back to the original folder some output 
