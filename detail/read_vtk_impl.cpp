/**
 * @file  read_vtk_impl.cpp
 * @brief Implementation of class: VTK_XML_Reader
 * @date 2019/02/22
 */

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <limits>
#include "../read_vtk.hpp"

namespace vtkmesh {

    //! Read analysis method
    void VTK_XML_Reader::readMethod(){

        const std::string tagMethod = VTK_Root + jointTags(R"(<xmlattr>)", R"(Name)");
        method_ = pt_.get<std::string>(tagMethod);

        auto itr = config.find(method_);
        if(itr == config.end()){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Invalid analysis method: " << method_ << " in " << filename_ << std::endl;
            exit(EXIT_FAILURE);
        }

    }

    //! Read coordinates
    void VTK_XML_Reader::readCoordinates(){

        // Read Coordinates
        Eigen::Matrix3Xd Coords3d = Eigen::Matrix3Xd::Zero(3, Mesh_.nPoint);
        std::string tagCoords = readDataArray<std::string>(R"(Points)", R"(Position)");
        std::stringstream ss1(tagCoords);
        std::string item1;
        int ipnt = 0;
        try {
            while(std::getline(ss1, item1)){
                std::stringstream ss2(item1);
                std::string item2;
                int idof = 0;
                while(std::getline(ss2, item2, ' ') && !item1.empty()){
                    if(item2.empty()){ continue; }
                    Coords3d(idof,ipnt) = std::stod(item2);
                    idof++;
                }
                if(idof == 3){
                    ipnt++;
                }
            }
        } catch(std::exception& e){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Error occurred while reading coordinates in " << filename_ << std::endl;
            std::cerr << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }

        if(ipnt != Mesh_.nPoint){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Cannot read coordinates properly in " << filename_ << std::endl;
            exit(EXIT_FAILURE);
        }

        Mesh_.nDoFnd = (Coords3d.row(2).array() < std::numeric_limits<float>::min()).all() ? 2 : 3;
        Mesh_.Coords = Eigen::MatrixXd::Zero(Mesh_.nDoFnd, Mesh_.nPoint);
        Mesh_.Coords = Coords3d.topRows(Mesh_.nDoFnd);

    }

    //! Read connectivity
    void VTK_XML_Reader::readConnectivity(){

        // Read connectivity
        Mesh_.Lnodes.clear();
        Mesh_.Lnodes.resize(Mesh_.nElemn);
        Mesh_.nNodes.clear();
        Mesh_.nNodes.resize(Mesh_.nElemn);
        std::string LndStr = readDataArray<std::string>(R"(Cells)", R"(connectivity)");
        std::stringstream sl1(LndStr);
        std::string cell1;
        int ielm = 0;
        try {
            while(std::getline(sl1, cell1)){
                std::stringstream sl2(cell1);
                std::string icell2;
                int inod = 0;
                bool flag = false;
                while(std::getline(sl2, icell2, ' ') && !cell1.empty()){
                    if(icell2.empty()){ continue; }
                    Mesh_.Lnodes[ielm].push_back(std::stoi(icell2));
                    inod++;
                    flag = true;
                }
                Mesh_.nNodes[ielm] = inod;
                if(flag){ielm++;}
            }
            Mesh_.mNodes = *std::max_element(Mesh_.nNodes.begin(), Mesh_.nNodes.end());
        } catch (std::exception& e){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Error occurred while reading connectivity in " << filename_ << std::endl;
            std::cerr << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }

        if(ielm != Mesh_.nElemn){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Cannot read connectivity properly in " << filename_ << std::endl;
            exit(EXIT_FAILURE);
        }

    }

    //! Read table and corresponding value (nodal value)
    std::tuple<std::vector<std::vector<int>>, Eigen::Matrix3Xd> VTK_XML_Reader::readPointTableValue(std::string&& tagTable, std::string&& tagValue){

        std::vector<std::vector<int>> Table(Mesh_.nPoint, std::vector<int>(3, 0));
        Eigen::Matrix3Xd Value = Eigen::Matrix3Xd::Zero(3, Mesh_.nPoint);

        try {
            std::string TableStr = readDataArray<std::string>(R"(PointData)", std::forward<std::string>(tagTable));
            std::stringstream sl1(TableStr);
            std::string cell1;
            int ipnt = 0;
            while(std::getline(sl1, cell1)){
                std::stringstream sl2(cell1);
                std::string icell2;
                int idof = 0;
                bool flag = false;
                while(std::getline(sl2, icell2, ' ') && !cell1.empty()){
                    if(icell2.empty()){ continue; }
                    Table[ipnt][idof] = std::stoi(icell2);
                    idof++;
                    flag = true;
                }
                if(flag){ipnt++;}
            }
        } catch (std::exception& e){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Error occurred while reading " + tagTable << " " << filename_ << std::endl;
            std::cerr << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }

        try {
            std::string prsdStr = readDataArray<std::string>(R"(PointData)", std::forward<std::string>(tagValue));
            std::stringstream sl1(prsdStr);
            std::string cell1;
            int ipnt = 0;
            while(std::getline(sl1, cell1)){
                std::stringstream sl2(cell1);
                std::string icell2;
                int idof = 0;
                bool flag = false;
                while(std::getline(sl2, icell2, ' ') && !cell1.empty()){
                    if(icell2.empty()){ continue; }
                    Value(idof, ipnt) = std::stod(icell2);
                    idof++;
                    flag = true;
                }
                if(flag){ipnt++;}
            }
        } catch (std::exception& e){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Error occurred while reading " << tagValue << " " << filename_ << std::endl;
            std::cerr << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }

        return std::make_tuple(Table, Value);

    }

    //! Read boundary condition
    void VTK_XML_Reader::readBoundaryCondition(){

        BC_.nBcgvn = 0;
        BC_.nBcnzr = 0;
        BC_.iBcgvn.clear();
        BC_.iBcnzr.clear();
        BC_.iDirec.clear();
        BC_.vBcnzr.clear();

        auto [BC_Table, PrescribedDisplacement] = readPointTableValue(R"(Boundary Condition)", R"(Prescribed Displacement)");

        for(auto ipnt=0;ipnt<Mesh_.nPoint;++ipnt){
            for(auto idof=0;idof<Mesh_.nDoFnd;++idof){
                int itdf = Mesh_.nDoFnd*ipnt + idof;
                switch(BC_Table[ipnt][idof]){
                    case static_cast<int>(Constraint::FIXED):
                        BC_.iBcgvn.push_back(itdf);
                        BC_.nBcgvn++;
                        break;
                    case static_cast<int>(Constraint::DISPLACEMENT):
                        BC_.iBcgvn.push_back(itdf);
                        BC_.nBcgvn++;
                        BC_.iBcnzr.push_back(itdf);
                        BC_.iDirec.push_back(idof);
                        BC_.vBcnzr.push_back(PrescribedDisplacement(idof, ipnt));
                        BC_.nBcnzr++;
                        break;
                    /*case static_cast<int>(Constraint::LOAD):
                        BC_.iBcgvn.push_back(itdf);
                        BC_.nBcgvn++;

                        break;*/
                    default:
                        break;
                }
            }
        }

        if(doesExistItems(R"(PointData)", R"(isBoundary)")){
            Mesh_.isBoundary = readDataArrayValue<int>(R"(PointData)", R"(isBoundary)");
            checkVectorSize(Mesh_.isBoundary, Mesh_.nPoint, R"(isBoundary)");
        }

    }

    //! Read multi-point constraints
    void VTK_XML_Reader::readMultiPointConstraints(){

        MPC_.nSlMpc = 0;
        MPC_.slvMpc.clear();
        MPC_.mstMpc.clear();
        MPC_.valMpc.clear();

        if(doesExistItems(R"(PointData)", R"(Multi-Point Constraints)") && doesExistItems(R"(PointData)", R"(MPC Ratio)")){
            auto [MPC_Table, MPC_Ratio] = readPointTableValue(R"(Multi-Point Constraints)", R"(MPC Ratio)");

            std::vector<std::vector<int>> mpcList(Mesh_.nDoFnd*Mesh_.nPoint);
            std::vector<double> mpcRatioList(Mesh_.nDoFnd*Mesh_.nPoint, 0.);
            for(auto ipnt=0;ipnt<Mesh_.nPoint;++ipnt){
                for(auto idof=0;idof<Mesh_.nDoFnd;++idof){
                    if(MPC_Table[ipnt][idof] > 0){
                        mpcList[MPC_Table[ipnt][idof]].push_back(Mesh_.nDoFnd*ipnt + idof);
                    }
                    mpcRatioList[Mesh_.nDoFnd*ipnt + idof] = MPC_Ratio(idof, ipnt);
                }
            }

            for(const auto& lst : mpcList){
                if(lst.size() > 1){
                    int idxmst = -1;
                    std::vector<int> idxslv;
                    for(const auto& c : lst){
                        if(abs(mpcRatioList[c]) > std::numeric_limits<float>::min()){
                            idxslv.push_back(c);
                        } else {
                            idxmst = c;
                        }
                    }
                    if(idxmst < 0 || idxslv.size() != lst.size() - 1){
                        std::lock_guard<std::mutex> lock(mtx);
                        std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
                        std::cerr << "Error occurred while reading MPC: Invalid MPC Format: " << filename_ << std::endl;
                    }
                    for(int i=0;i<(int)idxslv.size();++i){
                        MPC_.slvMpc.push_back(idxslv[i]);
                        MPC_.mstMpc.push_back(idxmst);
                        MPC_.valMpc.push_back(mpcRatioList[idxslv[i]]);
                        MPC_.nSlMpc++;
                    }
                }
            }
        }

        if(MPC_.nSlMpc == 0){
            MPC_.slvMpc.resize(1);
            MPC_.mstMpc.resize(1);
            MPC_.valMpc.resize(1);
        }

    }

    // Set specimen size, cross-sectional area and length of axis loading
    void VTK_XML_Reader::setGeometricalInformation(){

        // Specimen size
        Mesh_.specimenSize = Mesh_.Coords.rowwise().maxCoeff() - Mesh_.Coords.rowwise().minCoeff();

        // No prescribed displace
        if(BC_.nBcnzr == 0){
            Mesh_.idxAxisLoad = -1;
            Mesh_.lenAxisLoad = 0.;
            Mesh_.areaCrossSect = 0.;
            return;
        }

        Mesh_.idxAxisLoad = BC_.iBcnzr[0];
        int idofAxisLoad = Mesh_.idxAxisLoad%Mesh_.nDoFnd;
        Mesh_.lenAxisLoad = Mesh_.specimenSize(idofAxisLoad);
        Mesh_.areaCrossSect = 1.;
        for(int idof=0;idof<Mesh_.nDoFnd;++idof){
            if(idof != idofAxisLoad){
                Mesh_.areaCrossSect *= Mesh_.specimenSize(idof);
            }
        }

        return;
    }

    //! Check whether items exist or not
    bool VTK_XML_Reader::doesExistItems(std::string&& tagRoot, std::string&& tagChild){

        std::string strTag = VTK_Piece;
        if(tagRoot.length() > 0){
            strTag += jointTags(tagRoot);
        }

        if(tagRoot.length() > 0 && tagChild.length() > 0){
            for(const auto& child : pt_.get_child(strTag.c_str())){
                const auto& childContents = child.second;
                std::string tagName = childContents.get<std::string>(R"(<xmlattr>.Name)");
                if(tagName == tagChild){
                    return true;
                }
            }
        } else {
            if(tagChild.length() > 0){
                strTag += jointTags("<xmlattr>", tagChild);
            }
            auto child = pt_.get_child_optional(strTag);
            if(child){
                return true;
            }
        }
        return false;
    }

    //! Read material, crystal, ...
    void VTK_XML_Reader::readSubstructure(std::string&& dataType, int valDataArraySize){

        // Grain number
        if(doesExistItems(std::forward<std::string>(dataType), R"(Cryst)")){
            Material_.grainNumbers = readDataArrayValue<int>(std::forward<std::string>(dataType), R"(Cryst)");
            checkVectorSize(Material_.grainNumbers, valDataArraySize, R"(Crystal)");
            Material_.isSetGrainNumber = true;
        }

        // Material number
        if(doesExistItems(std::forward<std::string>(dataType), R"(Mater)")){
            Material_.materialNumbers = readDataArrayValue<int>(std::forward<std::string>(dataType), R"(Mater)");
            checkVectorSize(Material_.materialNumbers, valDataArraySize, R"(Material)");
            Material_.isSetMaterialNumber = true;
        }

        if(doesExistItems(std::forward<std::string>(dataType), R"(Crystal Orientation)") &&
        doesExistItems(std::forward<std::string>(dataType), R"(Phase-field)") &&
        doesExistItems(std::forward<std::string>(dataType), R"(Phase Flag)")){

            // Crystal orientation
            Material_.crystalOrientation = Eigen::Matrix3Xd::Zero(3, valDataArraySize);
            std::vector<std::string> orientationVec = readDataArrayValue<std::string>(std::forward<std::string>(dataType), R"(Crystal Orientation)");
            const double V_PI = 3.14159265358979323846;
            try {
                int ipnt = 0;
                for(const auto& c : orientationVec){
                    std::stringstream sl(c);
                    std::string cell;
                    int idof = 0;
                    if(!c.empty()){
                        while(std::getline(sl, cell, ' ')){
                            if(cell.empty()){ continue; }
                            Material_.crystalOrientation(idof, ipnt) = std::stod(cell)*V_PI/180.;
                            idof++;
                        }
                        ipnt++;
                    }
                }
            } catch(std::exception& e){
                std::lock_guard<std::mutex> lock(mtx);
                std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
                std::cerr << "Error occurred while reading crystal orientation: " << filename_ << std::endl;
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }

            // Phase field
            std::vector<double> phaseVec = readDataArrayValue<double>(std::forward<std::string>(dataType), R"(Phase-field)");
            checkVectorSize(phaseVec, valDataArraySize, R"(Phase-field)");
            Material_.phaseField = Eigen::Map<Eigen::VectorXd>(&phaseVec[0], phaseVec.size());

            // Phase flag
            Material_.newPhaseFlags = readDataArrayValue<int>(std::forward<std::string>(dataType), R"(Phase Flag)");
            checkVectorSize(Material_.newPhaseFlags, valDataArraySize, R"(New Phase Flags)");

            Material_.isSetPhaseField = true;
        }

    }


    //! Read mesh info
    void VTK_XML_Reader::readMeshInfo(){

        // Analysis method
        readMethod();

        // Number of nodes, elements, degrees of freedom
        Mesh_.nElemn = readComponent<int>("", "NumberOfCells");
        Mesh_.nPoint = readComponent<int>("", "NumberOfPoints");

        // Coordinates
        readCoordinates();

        // Connectivity
        readConnectivity();

        // Boundary condition
        readBoundaryCondition();

        // Mult-Point constraints
        readMultiPointConstraints();

        // Set specimen size, cross-sectional area and length of axis loading
        setGeometricalInformation();

        // Select data type
        std::string dataType = config.at(method_).at("use mesh") ? R"(CellData)" : R"(PointData)";
        int valDataArraySize = config.at(method_).at("use mesh") ? Mesh_.nElemn : Mesh_.nPoint;
        readSubstructure(std::forward<std::string>(dataType), valDataArraySize);

    }

}  // namespace vtkmesh
