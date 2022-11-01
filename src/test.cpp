#include "Jmle.h"
#include <iostream>


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

        Eigen::MatrixXd y = Eigen::MatrixXd::Identity(10,10);

        for (int i=0;i<10;i++)
            {
                x(i)=i/18+i*i*0.189+89;
            }

        function f;
        bool check;
        check=true;
        newt(x,check,f);

        //LUdcmp obj = LUdcmp(y);



        std::cout<<x<<std::endl;

    }