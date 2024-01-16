CC=c++
CFLAGS=-std=c++20 -Wall -Wextra -O2 -fPIC
CFLAGSPY=-shared $(shell python -m pybind11 --includes)
OUTPUT=pyvicsek$(shell python-config --extension-suffix)

all: pyvicsek

vicsek.o: vicsek.hpp vicsek.cpp
	$(CC) $(CFLAGS) -c vicsek.cpp

test: test_vicsek.cpp vicsek.o get_test_gif.py
	$(CC) $(CFLAGS) vicsek.o test_vicsek.cpp -o test_vicsek
	./test_vicsek 300 7 2
	./test_vicsek 300 25 0.1
	./test_vicsek 300 5 0.1
	python get_test_gif.py

pyvicsek: pyvicsek.cpp vicsek.o
	$(CC) $(CFLAGS) $(CFLAGSPY) pyvicsek.cpp vicsek.o -o $(OUTPUT)
	cp $(OUTPUT) examples/$(OUTPUT)

clear:
	rm vicsek.o test_vicsek $(OUTPUT) examples/$(OUTPUT) birds_*_test.* va_*_test.txt
