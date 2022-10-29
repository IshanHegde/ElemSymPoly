#include "RunningStats.h"
#include <iostream>

int main()
    {
        RunningStats rs;

        rs.Push(1000.0);
        rs.Push(-0.0000983);


        std::cout<<rs.Mean()<<std::endl;

        std::cout<<rs.StandardDeviation()<<std::endl;



    }