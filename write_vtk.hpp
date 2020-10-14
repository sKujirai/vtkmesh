/**
 * @file  write_vtk.hpp
 * @brief Output results in .vtk format
 * @date  2018/03/04
 */

#ifndef WRITE_VTK_H
#define WRITE_VTK_H

#include <cstdio>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <mutex>
#include <type_traits>
#include <Eigen/Core>
#include "config.hpp"

namespace vtkmesh {

    //! Numbering style
    enum class NumStyle : int {
        ZERO_BASED = 0,
        ONE_BASED = 1
    };

    //! Mesh type
    enum class MeshType : int {
        NODE_BASED = 0,
        ELEMENT_BASED = 1
    };

    /**
     * @class VTK_XML_Writer
     * @brief Class to write VTK file
     */
    class VTK_XML_Writer {

    private:
        FILE *fp;
        int ntpoin, nelemn, ndofnd;
        std::string buffer = "";
        std::queue<std::string> tagBegin;
        std::stack<std::string> tagEnd;
        const MeshType MESH_TYPE;
        const int offset;
        const std::string NEWLINE = "\n";
        const std::string vtkHeader = R"(<?xml version='1.0' encoding='UTF-8'?>)";
        const std::vector<std::string> vtkDataFormat {"ascii", "binary", "append"};
        const std::vector<std::string> vtkByteOrder {"LittleEndian", "BigEndian"};
        const std::vector<std::string> vtkDataType {"Float32", "Int32", "UInt8"};
        const std::vector<std::string> vtkType {"UnstructuredGrid", "PolyData"};
        const std::vector<std::string> vtkDataUnit {"PointData", "CellData"};

        //! Check tag
        VTKMESH_API bool checkType(const std::vector<std::string>& strVec, const std::string& inputStr);

        //! Write tag [begin]
        VTKMESH_API void writeTagBegin(const std::string& strbgn);

        //! Write tag [end]
        VTKMESH_API void writeTagEnd(const std::string& strend);

        //! Write tag [begin dataArray]
        VTKMESH_API void dataArrayBegin(const int numCmp, const std::string& datTyp, const std::string& dtName, const std::string& format = "ascii" );

        //! Write tag [end dataArray]
        VTKMESH_API void dataArrayEnd();

        //! Clamp value
        template <typename T>
        float clampValueFloat(const T& val);

        //! Write Vector
        template <typename TVector = Eigen::VectorXd>
        void writeVecXd(const TVector& Vec);

        //! Write Matrix
        template <typename TMatrix = Eigen::MatrixXd>
        void writeMatXd(const TMatrix& Mat);

        //! Write Coordinates
        template <typename TMatrix = Eigen::MatrixXd>
        void writeCoordinates(const TMatrix& Mat);

        //! Write Vector (std::vector) (Val - offset)
        VTKMESH_API void writeVector(const std::vector<int>& vec, const int offset=0);

        //! Write Vector of Vector (Val - offset)
        VTKMESH_API void writeVec2nd(const std::vector<std::vector<int>>& vec2nd, const int offset=0);

        //! Set Value
        VTKMESH_API void setValue(const Eigen::MatrixXd& Val, const std::string& dlabel, const std::string& datTyp = "Float32", const std::string& format = "ascii" );

        //! Set Value (std::vector)
        VTKMESH_API void setValue(const std::vector<int>& Val, const std::string& dlabel, const std::string& datTyp = "Int32", const std::string& format = "ascii" );

        //! Set dummy connectivity for Meshfree
        VTKMESH_API std::vector<std::vector<int>> getDummyConnectivity(const Eigen::MatrixXd& Coords);

        //! Get cell type
        VTKMESH_API int getCelltype(const int ndofnd, const int nnodes);

        //! Set mesh info for FEM
        VTKMESH_API void setMeshInfoELM(const Eigen::MatrixXd& Coords, const std::vector<std::vector<int>>& Lnodes);

        // Prohibit copy & move
        VTK_XML_Writer() = delete;
        VTK_XML_Writer(VTK_XML_Writer const&) = delete;
        VTK_XML_Writer(VTK_XML_Writer&&) = delete;
        VTK_XML_Writer& operator =(VTK_XML_Writer const&) = delete;
        VTK_XML_Writer& operator =(VTK_XML_Writer&&) = delete;

        // Mutex
        std::mutex mtx;

    public:
        VTK_XML_Writer(const MeshType mesh_type, const NumStyle offs = NumStyle::ZERO_BASED)
        : MESH_TYPE(mesh_type), offset(static_cast<int>(offs))
        {}
        virtual ~VTK_XML_Writer(){}

        //! Begin Tag
        VTKMESH_API void begin(const std::string& vtkUnit);

        //! End Tag
        VTKMESH_API void end(const std::string& vtkUnit);

        //! Set nodal value
        VTKMESH_API void setNodalValue(const Eigen::MatrixXd& Val, const std::string& dlabel, const std::string& datTyp = "Float32", const std::string& format = "ascii" );

        //! Set nodal value (std::vector)
        VTKMESH_API void setNodalValue(const std::vector<int>& Val, const std::string& dlabel, const std::string& datTyp = "Int32", const std::string& format = "ascii" );

        //! Set value on each element
        VTKMESH_API void setElmValue(const Eigen::MatrixXd& Val, const std::string& dlabel, const std::string& datTyp = "Float32", const std::string& format = "ascii" );

        //! Set value on each element (std::vector)
        VTKMESH_API void setElmValue(const std::vector<int>& Val, const std::string& dlabel, const std::string& datTyp = "Int32", const std::string& format = "ascii" );

        //! Set format
        VTKMESH_API void setVtkFormat(const std::string& _vtkType, const std::string& _discretization_method = "FiniteElement", const std::string& _byteOrder = "LittleEndian");

        //! Set Coordinates
        VTKMESH_API void setMeshInfo(const Eigen::MatrixXd& Coords, const std::vector<std::vector<int>>& Lnodes=std::vector<std::vector<int>>(0, std::vector<int>(0)));
        //! Write VTK
        VTKMESH_API void writeVtkFile(const std::string& filename);

        // Getter
        inline int getOffset() const { return offset; }

    };

}  // namespace vtkmesh

#include "detail/write_vtk_impl.hpp"
#ifdef VTKMESH_HEADER_ONLY
    #include "detail/write_vtk_impl.cpp"
#endif

#endif
