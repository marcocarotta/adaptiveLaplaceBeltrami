# Use the deal.II image from Docker Hub
FROM dealii/dealii:latest

# Set the working directory in the container
WORKDIR /usr/src/myapp

# Copy the current directory contents into the container at /usr/src/myapp
COPY . /usr/src/myapp

# Compile the C++ program
RUN mkdir build && cd build && cmake .. && make

# Run the executable
CMD ["./build/my_program"]

