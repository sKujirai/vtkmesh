# vtkmesh
C++ library to handle mesh file written in [VTK XML Formats](https://www.vtk.org/Wiki/VTK_XML_Formats) for FE analysis.

## Requirements
- [Eigen](http://eigen.tuxfamily.org): Template library for linear algebra (v3.3.7)
- [Boost](https://www.boost.org/): Semi-standard C++ library (v1.70.0)

## Example
```cpp
#include <vtkmesh/vtkmesh.hpp>

int main(void){

    // Read mesh
    vtkmesh::VTK_XML_Reader vtk_reader("mesh.vtu");

    // Write
    vtkmesh::VTK_XML_Writer vtkw(vtkmesh::MeshType::ELEMENT_BASED);
    vtkw.setVtkFormat("UnstructuredGrid");
    vtkw.setMeshInfo(vtkr.getCoordinates(), vtkr.getConnectivity());
    vtkw.writeVtkFile("out.vtu");

}
```
