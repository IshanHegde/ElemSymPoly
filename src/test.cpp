#include "Eigen/Dense"
#include <iostream>


int main()
    {


        Eigen::MatrixXf mat(4,4);
        mat << 1, 2, 6, 9,
         3, 1, 7, 2,
         5, 1, 9.02, 2,
         3, 1, 75, 17;

        mat=mat.array()-1.0;

        std::cout<<mat<<std::endl;

    }