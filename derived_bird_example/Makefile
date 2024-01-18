CC=c++
CFLAGS=-std=c++20 -Wall -Wextra -O2 -fPIC
CFLAGSPY=-shared $(shell python -m pybind11 --includes)
OUTPUT=pyvicsekvision$(shell python-config --extension-suffix)

all: pyvicsekvision

vicsek.o: vicsek.hpp vicsek.cpp
	$(CC) $(CFLAGS) -c vicsek.cpp

vicsek_vision.o: vicsek_vision.hpp vicsek_vision.cpp
	$(CC) $(CFLAGS) -c vicsek_vision.cpp

test: test_vicsek.cpp vicsek.o vicsek_vision.o
	$(CC) $(CFLAGS) test_vicsek.cpp vicsek.o vicsek_vision.o -o test_vicsek

pyvicsekvision: pyvicsekvision.cpp vicsek.o vicsek_vision.o
	$(CC) $(CFLAGS) $(CFLAGSPY) pyvicsekvision.cpp vicsek.o vicsek_vision.o -o $(OUTPUT)

clear:
	rm vicsek.o vicsek_vision.o $(OUTPUT) birds_*_test.* va_*_test.txt
