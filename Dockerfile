# Use the deal.II image from Docker Hub
FROM dealii/dealii:latest

# Set the working directory in the container
WORKDIR /usr/src/adaptiveLaplaceBeltrami

# Copy the current directory contents into the container at /usr/src/myapp
COPY . /usr/src/adaptiveLaplaceBeltrami

# Compile the C++ program
RUN mkdir build && cd build && cmake .. && make

# Run the executable
CMD ["./build/adaptiveLB"]

# TODO Maybe need copy back to the original folder some output 
