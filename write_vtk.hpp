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
#include <Eigen/Core>

namespace vtkmesh {

    /**
     * @class VTK_XML_Writer
     * @brief Class to write VTK file
     */
    class VTK_XML_Writer {

    public:
        // Numbering style
        static constexpr int ZERO_BASED = 0;
        static constexpr int ONE_BASED = 1;

        // Mesh type
        static constexpr int NODE_BASED = 0;
        static constexpr int ELEMENT_BASED = 1;

    private:
        FILE *fp;
        int ntpoin, nelemn, ndofnd;
        std::string buffer = "";
        std::queue<std::string> tagBegin;
        std::stack<std::string> tagEnd;
        const int MESH_TYPE;
        const int offset;
        const std::string NEWLINE = "\n";
        const std::string vtkHeader = R"(<?xml version='1.0' encoding='UTF-8'?>)";
        const std::vector<std::string> vtkDataFormat {"ascii", "binary", "append"};
        const std::vector<std::string> vtkByteOrder {"LittleEndian", "BigEndian"};
        const std::vector<std::string> vtkDataType {"Float32", "Int32", "UInt8"};
        const std::vector<std::string> vtkType {"UnstructuredGrid", "PolyData"};
        const std::vector<std::string> vtkDataUnit {"PointData", "CellData"};

        //! Check tag
        bool checkType(const std::vector<std::string>& strVec, const std::string& inputStr);

        //! Write tag [begin]
        void writeTagBegin(const std::string& strbgn);

        //! Write tag [end]
        inline void writeTagEnd(const std::string& strend);

        //! Write tag [begin dataArray]
        void dataArrayBegin(const int numCmp,
                            const std::string& datTyp,
                            const std::string& dtName,
                            const std::string& format = "ascii" );

        //! Write tag [end dataArray]
        void dataArrayEnd();

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
        void writeVector(const std::vector<int>& vec, const int offset=0);

        //! Write Vector of Vector (Val - offset)
        void writeVec2nd(const std::vector<std::vector<int>>& vec2nd, const int offset=0);

        //! Set Value
        void setValue(const Eigen::MatrixXd& Val,
                      const std::string& dlabel,
                      const std::string& datTyp = "Float32",
                      const std::string& format = "ascii" );

        //! Set Value (std::vector)
        void setValue(const std::vector<int>& Val,
                      const std::string& dlabel,
                      const std::string& datTyp = "Int32",
                      const std::string& format = "ascii" );

        //! Set dummy connectivity for Meshfree
        std::vector<std::vector<int>> getDummyConnectivity(const Eigen::MatrixXd& Coords);

        //! Get cell type
        int getCelltype(const int ndofnd, const int nnodes);

        //! Set mesh info for FEM
        void setMeshInfoELM(const Eigen::MatrixXd& Coords,
                            const std::vector<std::vector<int>>& Lnodes);

        // Prohibit copy & move
        VTK_XML_Writer() = delete;
        VTK_XML_Writer(VTK_XML_Writer const&) = delete;
        VTK_XML_Writer(VTK_XML_Writer&&) = delete;
        VTK_XML_Writer& operator =(VTK_XML_Writer const&) = delete;
        VTK_XML_Writer& operator =(VTK_XML_Writer&&) = delete;

        // Mutex
        std::mutex mtx;

    public:
        VTK_XML_Writer(const int mesh_type, const int offs = 0)
        : MESH_TYPE(mesh_type), offset(offs)
        {}
        virtual ~VTK_XML_Writer(){}

        //! Begin Tag
        void begin(const std::string& vtkUnit);

        //! End Tag
        void end(const std::string& vtkUnit);

        //! Set nodal value
        void setNodalValue(const Eigen::MatrixXd& Val,
                           const std::string& dlabel,
                           const std::string& datTyp = "Float32",
                           const std::string& format = "ascii" );

        //! Set nodal value (std::vector)
        void setNodalValue(const std::vector<int>& Val,
                           const std::string& dlabel,
                           const std::string& datTyp = "Int32",
                           const std::string& format = "ascii" );

        //! Set value on each element
        void setElmValue(const Eigen::MatrixXd& Val,
                         const std::string& dlabel,
                         const std::string& datTyp = "Float32",
                         const std::string& format = "ascii" );

        //! Set value on each element (std::vector)
        void setElmValue(const std::vector<int>& Val,
                         const std::string& dlabel,
                         const std::string& datTyp = "Int32",
                         const std::string& format = "ascii" );

        //! Set format
        void setVtkFormat(const std::string& _vtkType,
                          const std::string& _discretization_method = "FiniteElement",
                          const std::string& _byteOrder = "LittleEndian");

        //! Set Coordinates
        void setMeshInfo(const Eigen::MatrixXd& Coords,
                         const std::vector<std::vector<int>>& Lnodes=std::vector<std::vector<int>>(0, std::vector<int>(0)));
        //! Write VTK
        void writeVtkFile(const std::string& filename);

        // Getter
        inline int getOffset() const { return offset; }

    };

}  // namespace vtkmesh

#include "detail/write_vtk_impl.hpp"

#endif
