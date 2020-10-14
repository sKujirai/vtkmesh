readvtksrc := ../detail/read_vtk_impl.cpp
writevtksrc := ../detail/write_vtk_impl.cpp

all: example_head example_lib

example_head: example.cpp
	g++ -Wall -Wextra -std=c++17 -o $@ $<

example_lib: example.cpp
	g++ -Wall -Wextra -std=c++17 -DVTKMESH_COMPILED_LIB -o $@ $< $(readvtksrc) $(writevtksrc)

.PHONY: clean

clean:
	rm -rf example_head example_lib out.vtu
