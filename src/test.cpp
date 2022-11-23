#include "Newton.h"
#include <iostream>
#include <vector>
#include "Eigen/Dense"

struct function
    {
        Eigen::VectorXd operator()(Eigen::VectorXd in )
            {
                Eigen::VectorXd out(in.rows());

                for (int i =0;i<in.rows();i++)
                    {
                        out(i)=in(i)*in(i)+87*in(i)+0.188;
                    }
                return out;
            }
    };

int main()
    {

        Eigen::VectorXd x(10);
        Eigen::VectorXd y(10);

        for (int i=0;i<10;i++)
            {
                x(i)=i;
                y(i)=2;
            }
        
        std::cout<<y.colwise().sum()<<std::endl;

        

    }