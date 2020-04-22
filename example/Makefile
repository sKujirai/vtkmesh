readvtksrc := ../detail/read_vtk_impl.cpp
writevtksrc := ../detail/write_vtk_impl.cpp

all: example

example: example.cpp
	g++ -Wall -Wextra -std=c++17 -o $@ $< $(readvtksrc) $(writevtksrc)

.PHONY: clean

clean:
	rm -rf example out.vtu
