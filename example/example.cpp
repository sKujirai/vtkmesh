#include <iostream>
#include <memory>
#include <cstdlib>
#include <iostream>
#include "../read_vtk.hpp"
#include "../write_vtk.hpp"

int main(void){

    const std::string filename = "./data/mesh.vtu";

    const std::string outputfile = "./out.vtu";

    std::cout << "Reading .vtk file: " << filename << std::endl;

    std::unique_ptr<vtkmesh::VTK_XML_Reader> vtkr;
    try {
        vtkr = std::make_unique<vtkmesh::VTK_XML_Reader>(filename);
    } catch (const std::exception& e){
        std::cout << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "Writing .vtk file: " << outputfile << std::endl;

    std::unique_ptr<vtkmesh::VTK_XML_Writer> vtkw;
    try {
        vtkw = std::make_unique<vtkmesh::VTK_XML_Writer>(vtkmesh::VTK_XML_Writer::ELEMENT_BASED, vtkmesh::VTK_XML_Writer::ZERO_BASED);
    } catch (const std::exception& e){
        std::cout << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    // std::cout << vtkw->getOffset() << std::endl;
    vtkw->setVtkFormat("UnstructuredGrid");
    vtkw->setMeshInfo(vtkr->getCoordinates(), vtkr->getConnectivity());
    vtkw->writeVtkFile(outputfile);

    std::cout << "...Program end" << std::endl;
}