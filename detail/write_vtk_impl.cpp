/**
 * @file  write_vtk_impl.cpp
 * @brief Implementation of class: VTK_XML_Writer
 * @date  2018/03/04
 */

#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "../write_vtk.hpp"

namespace VtkMesh {

    //! Check tag
    bool VTK_XML_Writer::checkType(const std::vector<std::string>& strVec, const std::string& inputStr){
        for(const auto& s : strVec){
            if(s == inputStr){
                return true;
            }
        }
        return false;
    }

    //! Write tag [begin]
    void VTK_XML_Writer::writeTagBegin(const std::string& strbgn){
        buffer += R"(<)" + strbgn + R"(>)" + NEWLINE;
        tagEnd.push(strbgn);
        return;
    }

    //! Write tag [end]
    void VTK_XML_Writer::writeTagEnd(const std::string& strend){
        if(strend != tagEnd.top()){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid order of tag: " << strend << " " << tagEnd.top() << std::endl;
            exit(EXIT_FAILURE);
        }
        buffer += R"(</)" + strend + R"(>)" + NEWLINE;
        tagEnd.pop();

        return;
    }

    //! Write tag [begin dataArray]
    void VTK_XML_Writer::dataArrayBegin(const int numCmp,
                                        const std::string& datTyp,
                                        const std::string& dtName,
                                        const std::string& format){

        std::string numstr = "";
        if(numCmp > 0){
            numstr = R"(NumberOfComponents=')"
                + std::to_string(numCmp)
                + R"(' )";
        }

        if(!checkType(vtkDataType, datTyp)){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Data Type: " << datTyp << std::endl;
            exit(EXIT_FAILURE);
        }
        if(!checkType(vtkDataFormat, format)){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Data Format: " << format << std::endl;
            exit(EXIT_FAILURE);
        }

        std::string tmptag;
        tmptag = R"(<DataArray )"
            + numstr
            + R"(type=')" + datTyp +
            + R"(' Name=')" + dtName + R"(' format=')"
            + format + R"('>)" + NEWLINE;

        buffer += tmptag;
        tagEnd.push("DataArray");

        return;
    }

    //! Write tag [end dataArray]
    void VTK_XML_Writer::dataArrayEnd(){ writeTagEnd("DataArray"); }

    //! Write Vector (std::vector) (Val - offset)
    void VTK_XML_Writer::writeVector(const std::vector<int>& vec, const int offset){
        int nsize = (int)(vec.size());
        for(int i=0;i<nsize;i++){
            buffer += std::to_string(vec[i]-offset) + " ";
        }
        return;
    }

    //! Write Vector of Vector (Val - 1)
    void VTK_XML_Writer::writeVec2nd(const std::vector<std::vector<int>>& vec2nd, const int offset){
        int ncol = (int)(vec2nd.size());
        for(int j=0;j<ncol;j++){
            int nrow = (int)(vec2nd[j].size());
            for(int i=0;i<nrow;i++){
                buffer += std::to_string(vec2nd[j][i]-offset) + " ";
            }
            buffer += NEWLINE;
        }
        return;
    }

    //! Set Value
    void VTK_XML_Writer::setValue(const Eigen::MatrixXd& Val,
                                  const std::string& dlabel,
                                  const std::string& datTyp,
                                  const std::string& format){

        if(!checkType(vtkDataType, datTyp)){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Data Type: " << datTyp << std::endl;
            exit(EXIT_FAILURE);
        }
        if(!checkType(vtkDataFormat, format)){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Data Format: " << format << std::endl;
            exit(EXIT_FAILURE);
        }

        dataArrayBegin(Val.cols(), datTyp, dlabel, format);
        writeMatXd(Val);
        dataArrayEnd();

        return;
    }

    //! Set Value (std::vector)
    void VTK_XML_Writer::setValue(const std::vector<int>& Val,
                                  const std::string& dlabel,
                                  const std::string& datTyp,
                                  const std::string& format){

        if(!checkType(vtkDataType, datTyp)){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Data Type: " << datTyp << std::endl;
            exit(EXIT_FAILURE);
        }
        if(!checkType(vtkDataFormat, format)){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Data Format: " << format << std::endl;
            exit(EXIT_FAILURE);
        }

        dataArrayBegin(1, datTyp, dlabel, format);
        writeVector(Val);
        dataArrayEnd();

        return;
    }

    //! Begin Tag
    void VTK_XML_Writer::begin(const std::string& vtkUnit){
        if(!checkType(vtkDataUnit, vtkUnit)){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid VTK Tag: " << vtkUnit << std::endl;
            exit(EXIT_FAILURE);
        }
        writeTagBegin(vtkUnit);
    }

    //! End Tag
    void VTK_XML_Writer::end(const std::string& vtkUnit){ writeTagEnd(vtkUnit); }

    //! Set nodal value
    void VTK_XML_Writer::setNodalValue(const Eigen::MatrixXd& Val,
                                       const std::string& dlabel,
                                       const std::string& datTyp,
                                       const std::string& format){

        if((Val.rows()) != ntpoin){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Data Size: nPoint: " << ntpoin << " data size : " << Val.rows() << std::endl;
            exit(EXIT_FAILURE);
        }
        setValue(Val, dlabel, datTyp, format);

        return;
    }

    //! Set nodal value (std::vector)
    void VTK_XML_Writer::setNodalValue(const std::vector<int>& Val,
                                       const std::string& dlabel,
                                       const std::string& datTyp,
                                       const std::string& format){

        if((int)Val.size() != ntpoin){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Data Size: nPoint: " << ntpoin << " data size : " << Val.size() << std::endl;
            exit(EXIT_FAILURE);
        }
        setValue(Val, dlabel, datTyp, format);

        return;
    }

    //! Set value on each element
    void VTK_XML_Writer::setElmValue(const Eigen::MatrixXd& Val,
                                     const std::string& dlabel,
                                     const std::string& datTyp,
                                     const std::string& format){

        if((Val.rows()) != nelemn){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Data Size: nElemn: " << nelemn << " data size : " << Val.rows() << std::endl;
            exit(EXIT_FAILURE);
        }
        setValue(Val, dlabel, datTyp, format);

        return;
    }

    //! Set value on each element (std::vector)
    void VTK_XML_Writer::setElmValue(const std::vector<int>& Val,
                                     const std::string& dlabel,
                                     const std::string& datTyp,
                                     const std::string& format){

        if((int)Val.size() != nelemn){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Data Size: nElemn: " << nelemn << " data size : " << Val.size() << std::endl;
            exit(EXIT_FAILURE);
        }
        setValue(Val, dlabel, datTyp, format);

        return;
    }

    //! Set format
    void VTK_XML_Writer::setVtkFormat(const std::string& _vtkType, const std::string& _discretization_method, const std::string& _byteOrder){

        if(!checkType(vtkType, _vtkType)){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid VTK Type: " << _vtkType << std::endl;
            exit(EXIT_FAILURE);
        }
        if(!checkType(vtkByteOrder, _byteOrder)){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid Byte Order: " << _byteOrder << std::endl;
            exit(EXIT_FAILURE);
        }

        buffer += vtkHeader + NEWLINE;
        buffer += R"(<VTKFile xmlns='VTK' byte_order=')" +  _byteOrder
            + R"(' version='0.1' type=')" + _vtkType + R"(' name=')"
            + _discretization_method + R"('>)" + NEWLINE;

        tagEnd.push("VTKFile");

        writeTagBegin(_vtkType);

        return;
    }

    //! Set dummy connectivity for Meshfree
    std::vector<std::vector<int>> VTK_XML_Writer::getDummyConnectivity(const Eigen::MatrixXd& Coords){

        std::vector<std::vector<int>> DummyLnodes(Coords.cols());
        for(int i=0;i<(int)DummyLnodes.size();++i){
            DummyLnodes[i].resize(1, i);
        }

        return DummyLnodes;
    }

    //! Get cell type
    int VTK_XML_Writer::getCelltype(const int ndofnd, const int nnodes){

        switch(ndofnd){
        case 2:
            switch(nnodes){
            case 1:
                return 1;
                break;
            case 2:
                return 3;
                break;
            case 3:
            case 4:
                return 7;
                break;
            case 6:
                return 22;
            case 8:
            case 9:
                return 23;
                break;           
            default:
                return 7;  // for voronoi
                break;
            }
            break;
        case 3:
            switch(nnodes){
            case 1:
                return 1;
                break;
            case 2:
                return 3;
                break;
            case 4:
                return 10;
                break;
            case 8:
                return 12;
                break;            
            default:
                std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
                std::cerr << "Invalid Element: number of nodes: " << nnodes << std::endl;
                exit(EXIT_FAILURE);
                break;
            }
            break;       
        default:
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid DoF" << std::endl;
            exit(EXIT_FAILURE);
            break;
        }
    }

    //! Set mesh info for FEM
    void VTK_XML_Writer::setMeshInfoELM(const Eigen::MatrixXd& Coords,
                                        const std::vector<std::vector<int>>& Lnodes){

        // Set Piece
        ntpoin = Coords.cols();
        ndofnd = Coords.rows();
        nelemn = Lnodes.size();

        buffer += R"(<Piece NumberOfCells=')" + std::to_string(nelemn)
            + R"(' NumberOfPoints=')" + std::to_string(ntpoin)
            + R"(' DegreesOfFreedom=')" + std::to_string(ndofnd)
            + R"('>)" + NEWLINE;
        tagEnd.push("Piece");

        // Set Coordinates
        writeTagBegin("Points");
        buffer += R"(<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>)" + NEWLINE;
        tagEnd.push("DataArray");

        writeCoordinates(Coords);
        dataArrayEnd();
        writeTagEnd("Points");

        // Set Connectivity
        writeTagBegin("Cells");
        dataArrayBegin(0, "Int32", "connectivity", "ascii");
        writeVec2nd(Lnodes, offset);
        dataArrayEnd();
        dataArrayBegin(0, "Int32", "offsets", "ascii");

        int icnt = 0;
        for(int ielm=0;ielm<nelemn;ielm++){
            icnt += Lnodes[ielm].size();
            buffer += std::to_string(icnt) + NEWLINE;
        }

        dataArrayEnd();
        dataArrayBegin(0, "UInt8", "types", "ascii");

        for(int ielm=0;ielm<nelemn;ielm++){
            buffer += std::to_string(getCelltype(ndofnd, Lnodes[ielm].size())) + NEWLINE;
        }

        dataArrayEnd();
        writeTagEnd("Cells");

        return;
    }

    //! Set Coordinates
    void VTK_XML_Writer::setMeshInfo(const Eigen::MatrixXd& Coords,
                                     const std::vector<std::vector<int>>& Lnodes){

        switch(MESH_TYPE){
            case NODE_BASED:
                setMeshInfoELM(Coords, getDummyConnectivity(Coords));
                break;
            case ELEMENT_BASED:
                setMeshInfoELM(Coords, Lnodes);
                break;
            default:
                std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
                std::cerr << "Invalid mesh type: " << MESH_TYPE << std::endl;
                exit(EXIT_FAILURE);
                break;
        }

    }

    //! Write VTK
    void VTK_XML_Writer::writeVtkFile(const std::string& filename){

        while(!tagEnd.empty()){
            writeTagEnd(tagEnd.top());
        }

        // Write data to .vtk file
        if((fp = fopen(filename.c_str(), "w")) == nullptr){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Cannot open output .vtk file: " << filename << std::endl;
            exit(EXIT_FAILURE);
        }
        try {
            fprintf(fp, "%s", buffer.c_str());
        } catch(...){
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Error occurred while writing .vtk file: " << filename << std::endl;
            exit(EXIT_FAILURE);
        }
        fclose(fp);
        buffer.clear();
        fp = nullptr;

        return;
    }

}  // namespace VtkMesh
