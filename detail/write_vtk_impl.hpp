/**
 * @file  write_vtk_impl.cpp
 * @brief Implementation of class: VTK_XML_Writer
 * @date  2018/03/04
 */

#ifndef WRITE_VTK_IMPL_H
#define WRITE_VTK_IMPL_H

#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <queue>
#include <stack>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <limits>
#include <Eigen/Core>

namespace vtkmesh {

    //! Clamp value
    template <typename T>
    float VTK_XML_Writer::clampValueFloat(const T& val){
        float vclamp = std::abs(val) < std::numeric_limits<float>::min() ? static_cast<float>(0) : static_cast<float>(val);
        return std::clamp(vclamp, -std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    }

    //! Write Vector
    template <typename TVector>
    void VTK_XML_Writer::writeVecXd(const TVector& Vec){
        std::stringstream ss;
        for(int i=0;i<Vec.size();i++){
            ss << std::scientific << clampValueFloat(Vec(i));
            buffer += ss.str() + NEWLINE;
            ss.clear(); ss.str("");
        }
        return;
    }

    //! Write Matrix
    template <typename TMatrix>
    void VTK_XML_Writer::writeMatXd(const TMatrix& Mat){
        std::stringstream ss;
        for(int i=0;i<Mat.rows();i++){
            for(int j=0;j<Mat.cols();j++){
                ss << std::scientific << clampValueFloat(Mat(i,j));
                buffer += ss.str() + " ";
                ss.clear(); ss.str("");
            }
            buffer += NEWLINE;
        }
        return;
    }

    //! Write Coordinates
    template <typename TMatrix>
    void VTK_XML_Writer::writeCoordinates(const TMatrix& Mat){
        std::lock_guard<std::mutex> lock(mtx);
        std::stringstream ss;
        for(int j=0;j<Mat.cols();j++){
            int idof = 0;
            for(int i=0;i<Mat.rows();i++){
                ss << std::scientific << clampValueFloat(Mat(i,j));
                buffer += ss.str() + " ";
                ss.clear(); ss.str("");
                idof++;
            }
            for(int k=idof;k<3;k++){
                buffer += std::to_string(0.) + " ";
            }
            buffer += NEWLINE;
        }
        return;
    }

}  // namespace vtkmesh

#endif
