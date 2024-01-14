# Compilers
CC=g++
NVCC=nvcc

# Compiler Flags
CFLAGS=-std=c++20
CFLAGS+=-g

NVFLAGS=-arch=sm_86 # Adjust for your GPU
NVFLAGS+=-G
NVFLAGS+=-g
NVCC_FLAGS=-rdc=true # Allows for CUDA dynamic parallelism

# Source files
CPP_SOURCES := $(wildcard *.cpp)
CU_SOURCES := $(wildcard *.cu)

# Object files
CPP_OBJECTS := $(CPP_SOURCES:.cpp=.o)
CU_OBJECTS := $(CU_SOURCES:.cu=.o)

# Executable name
EXECUTABLE=myApp

# Default rule
all: $(EXECUTABLE)

# Rule for making the executable
$(EXECUTABLE): $(CPP_OBJECTS) $(CU_OBJECTS) device_link.o
	$(NVCC) $(NVFLAGS) $^ -o $@

# Rule for making CPP object files
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Rule for making CUDA object files
%.o: %.cu
	$(NVCC) $(NVFLAGS) $(NVCC_FLAGS) -c $< -o $@

# Device linking rule
device_link.o: $(CU_OBJECTS)
	$(NVCC) $(NVFLAGS) $(NVCC_FLAGS) -dlink $^ -o $@

# Clean rule
clean:
	rm -f $(CPP_OBJECTS) $(CU_OBJECTS) device_link.o $(EXECUTABLE)
