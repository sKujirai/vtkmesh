/**
 * @file  read_vtk.hpp
 * @brief Read mesh file written in .vtk format
 * @date 2019/02/22
 */

#ifndef READ_VTK_H
#define READ_VTK_H

#include <Eigen/Core>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <mutex>
#include <boost/property_tree/xml_parser.hpp>
#include "config.hpp"

namespace vtkmesh {

    /**
     * @class VTK_XML_Reader
     * @brief VTK (XML format) reader
     */
    class VTK_XML_Reader {

    private:
        boost::property_tree::ptree pt_;
        const std::string filename_;
        std::string method_;

        // VTK tags
        std::vector<std::string> tag_list;
        const std::string VTK_Root = R"(VTKFile)";
        const std::string VTK_UnstructuredGrid = R"(VTKFile.UnstructuredGrid)";
        const std::string VTK_Piece = R"(VTKFile.UnstructuredGrid.Piece)";
        const std::vector<std::string> DataTypeList = {R"(Points)", R"(Cells)", R"(PointData)", R"(CellData)"};
        const std::vector<std::string> MethodList = {R"(FiniteElement)", R"(Meshfree)", R"(MeshfreeELM)"};

        const std::map<std::string, std::map<std::string, bool>> config = {
            {R"(FiniteElement)", {{R"(use mesh)", true}}},
            {R"(Meshfree)", {{R"(use mesh)", false}}},
            {R"(MeshfreeELM)", {{R"(use mesh)", true}}}
        };

        //! Constraints
        enum class Constraint : int {
            FIXED = 1,
            DISPLACEMENT = -1,
            LOAD = -2
        };

        // Joint tags
        std::string jointTags(){
            return std::string("");
        }
        template<typename First, typename... Rest>
        std::string jointTags(First&& Tf, Rest&&... Tr){
            return std::string(".") + static_cast<std::string>(Tf) + jointTags(std::forward<Rest>(Tr)...);
        }

        //! Read value from xml file
        template <typename T>
        T readComponent(std::string&& tagRoot, std::string&& tagChild="");

        //! Read array from xml file
        template <typename T>
        T readDataArray(std::string&& dataType, std::string&& name);

        //! Read array from xml file (return vector value)
        template <typename T>
        std::vector<T> readDataArrayValue(std::string&& dataType, std::string&& tag);

        //! Check whether item exists from tag
        VTKMESH_API bool doesExistItems(std::string&& tagRoot, std::string&& tagChild);

        //! Read connectivity and value
        VTKMESH_API std::tuple<std::vector<std::vector<int>>, Eigen::Matrix3Xd> readPointTableValue(std::string&& tagTable, std::string&& tagValue);

        //! Read analysis method
        VTKMESH_API void readMethod();

        //! Read coordinates
        VTKMESH_API void readCoordinates();

        //! Read connectivity
        VTKMESH_API void readConnectivity();

        //! Read boundary condition
        VTKMESH_API void readBoundaryCondition();

        //! Read connectivity
        VTKMESH_API void readMultiPointConstraints();

        //! Set specimen size, cross-sectional area and length of axis loading
        VTKMESH_API void setGeometricalInformation();

        //! Read material, crystal, ...
        VTKMESH_API void readSubstructure(std::string&& dataType, int valDataArraySize);

        //! Read mesh info.
        VTKMESH_API void readMeshInfo();

        //! Check size of vector
        template <typename T>
        void checkVectorSize(const std::vector<T>& vec, const int trueSize, std::string&& valname);

        //! Mesh info
        struct Mesh {
            int nDoFnd, nPoint, nElemn, mNodes;
            std::vector<int> nNodes;
            std::vector<std::vector<int>> Lnodes;
            Eigen::MatrixXd Coords;
            Eigen::VectorXd specimenSize;
            std::vector<int> isBoundary;
            int idxAxisLoad;
            double areaCrossSect;
            double lenAxisLoad;
        } Mesh_;

        //! Boundary condition
        struct BoundaryCondition {
            int nBcgvn;
            int nBcnzr;
            std::vector<int> iBcgvn;
            std::vector<int> iBcnzr;
            std::vector<int> iDirec;
            std::vector<double> vBcnzr;
        } BC_;

        //! Multi-Point constraits
        struct MultiPointConstraints {
            int nSlMpc;
            std::vector<int> slvMpc;
            std::vector<int> mstMpc;
            std::vector<double> valMpc;
        } MPC_;

        //! Config
        struct Material {
            std::vector<int> grainNumbers;
            std::vector<int> materialNumbers;
            Eigen::Matrix3Xd crystalOrientation;
            Eigen::VectorXd phaseField;
            std::vector<int> newPhaseFlags;
            bool isSetGrainNumber = false;
            bool isSetMaterialNumber = false;
            bool isSetCrystalOrientation = false;
            bool isSetPhaseField = false;
        } Material_;

        // Prohibit copy & move
        VTK_XML_Reader() = delete;
        VTK_XML_Reader(VTK_XML_Reader const&) = delete;
        VTK_XML_Reader(VTK_XML_Reader&&) = delete;
        VTK_XML_Reader& operator =(VTK_XML_Reader const&) = delete;
        VTK_XML_Reader& operator =(VTK_XML_Reader&&) = delete;

        // Mutex
        std::mutex mtx;

    public:
        VTK_XML_Reader(const std::string& file)
        : filename_(file)
        {
            boost::property_tree::xml_parser::read_xml(file, pt_);
            readMeshInfo();
        }
        virtual ~VTK_XML_Reader(){}

        // Getter
        std::string getMethod() const { return method_; }
        // Mesh
        inline int getNumDoF() const { return Mesh_.nDoFnd; }
        inline int getNumPoints() const { return Mesh_.nPoint; }
        inline int getNumElements() const { return Mesh_.nElemn; }
        inline int getMaxNumNodes() const { return Mesh_.mNodes; }
        inline std::vector<int> getNumNodes() const { return Mesh_.nNodes; }
        inline Eigen::MatrixXd getCoordinates() const { return Mesh_.Coords; }
        inline std::vector<std::vector<int>> getConnectivity() const { return Mesh_.Lnodes; }
        inline Eigen::VectorXd getSpecimenSize() const { return Mesh_.specimenSize; }
        inline std::vector<int> getIsBoundary() const { return Mesh_.isBoundary; }
        inline int getIdxAxisLoad() const { return Mesh_.idxAxisLoad; }
        inline double getAreaCrossSect() const { return Mesh_.areaCrossSect; }
        inline double getLenAxisLoad() const { return Mesh_.lenAxisLoad; }
        // Material
        inline std::vector<int> getGrainNumbers() const { return Material_.grainNumbers; }
        inline std::vector<int> getMaterialNumbers() const { return Material_.materialNumbers; }
        inline Eigen::Matrix3Xd getCrystalOrientation() const { return Material_.crystalOrientation; }
        inline Eigen::VectorXd getPhaseField() const { return Material_.phaseField; }
        inline std::vector<int> getNewPhaseFlags() const { return Material_.newPhaseFlags; }
        // B.C.
        inline int getNumBcgvn() const { return BC_.nBcgvn; }
        inline int getNumBcnzr() const { return BC_.nBcnzr; }
        inline std::vector<int> getiBcgvn() const { return BC_.iBcgvn; }
        inline std::vector<int> getiBcnzr() const { return BC_.iBcnzr; }
        inline std::vector<int> getDirection() const { return BC_.iDirec; }
        inline std::vector<double> getvBcnzr() const {return BC_.vBcnzr;}
        // MPC
        inline int getNumSlMpc() const { return MPC_.nSlMpc; }
        inline std::vector<int> getSlvMpc() const { return MPC_.slvMpc; }
        inline std::vector<int> getMstMpc() const { return MPC_.mstMpc; }
        inline std::vector<double> getValMpc() const { return MPC_.valMpc; }

    };

}  // namespace vtkmesh

#include "detail/read_vtk_impl.hpp"
#ifdef VTKMESH_HEADER_ONLY
    #include "detail/read_vtk_impl.cpp"
#endif

#endif
