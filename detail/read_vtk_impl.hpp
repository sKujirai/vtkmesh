/**
 * @file  read_vtk_impl.hpp
 * @brief Implementation of class: VTK_XML_Reader
 * @date 2019/02/22
 */

#ifndef READ_VTK_IMPL_H
#define READ_VTK_IMPL_H

#include <Eigen/Core>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <mutex>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/lexical_cast.hpp>

namespace vtkmesh {

    //! Read components
    template <typename T>
    T VTK_XML_Reader::readComponent(std::string&& tagRoot, std::string&& tagChild){

        std::string strTag = VTK_Piece;
        if(tagRoot.length() > 0){
            strTag += jointTags(tagRoot);
        }
        if(tagChild.length() > 0){
            strTag += jointTags("<xmlattr>", tagChild);
        }

        if(boost::optional<T> sVal =
            pt_.get_optional<T>(strTag.c_str())){
            return *sVal;
        } else {
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Cannot read " << strTag << " in " << filename_ << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    //! Read DataArray
    template <typename T>
    T VTK_XML_Reader::readDataArray(std::string&& dataType, std::string&& name){

        // Check whether data type exists in unstructured grid
        auto itr = std::find(DataTypeList.begin(), DataTypeList.end(), dataType);
        if(itr == DataTypeList.end()){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Cannot find data type: " << dataType << " in " << filename_ << std::endl;
            exit(EXIT_FAILURE);
        }

        // Get value
        const std::string tagPath = VTK_Piece + jointTags(dataType);
        try{
            for(const auto& child : pt_.get_child(tagPath.c_str())){
                const auto& childContents = child.second;
                std::string tagName = childContents.get<std::string>(R"(<xmlattr>.Name)");
                if(tagName == name){
                    return childContents.get<T>("");
                }
            }
        } catch(std::exception& e){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Error occurred while reading DataArray: " << tagPath << " in " << filename_ << std::endl;
            exit(EXIT_FAILURE);
        }

        std::lock_guard<std::mutex> lock(mtx);
        std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
        std::cerr << "Cannot find tag: " << name << " in " << filename_ << std::endl;
        exit(EXIT_FAILURE);

    }

    //! Read DataArray as vector
    template <typename T>
    std::vector<T> VTK_XML_Reader::readDataArrayValue(std::string&& dataType, std::string&& tag){

        std::vector<T> valArray;
        try {
            std::string valStr = readDataArray<std::string>(std::forward<std::string>(dataType), std::forward<std::string>(tag));
            std::stringstream sl(valStr);
            std::string cell;
            while(std::getline(sl, cell)){
                if(!cell.empty()){
                    unsigned int idx;
                    while((idx = cell.find_first_of(" \t")) == 0){
                        cell.erase(cell.begin());
                        if(cell.empty()){ break; }
                    }
                    valArray.push_back(boost::lexical_cast<T>(cell));
                }
            }
        } catch(boost::bad_lexical_cast& e){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << "Error occurred while reading " << tag << " " << filename_ << std::endl;
            std::cerr << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }

        return valArray;
    }

    //! Check size of vector
    template <typename T>
    void VTK_XML_Reader::checkVectorSize(const std::vector<T>& vec, const int trueSize, std::string&& valname){

        if((int)vec.size() != trueSize){
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ": ";
            std::cerr << valname << "'s size: " << vec.size() << " must be equal to " << trueSize << " : " << filename_ << std::endl;
            exit(EXIT_FAILURE);
        }

    }

}  // namespace vtkmesh

#endif
