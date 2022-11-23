#ifndef RUNNINGSTATS_H
#define RUNNINGSTATS_H
#include <cmath>
#include <iostream>
#include <pybind11/eigen.h>
//#include "Eigen/Dense"
// creddit to Dr. John D. Cook 
// link: https://www.johndcook.com/blog/skewness_kurtosis/

class RunningStats
{
    public:
        RunningStats();
        void Clear();
        void Push(double x,double w);
        long long NumDataValues() const;
        double Mean() const;
        double Variance() const;
        double StandardDeviation() const;

        friend RunningStats operator+(const RunningStats a, const RunningStats b);
        RunningStats& operator+=(const RunningStats &rhs);

    private:
        long long n;
        double w_sum, w_sum2, mean, S;
};

RunningStats::RunningStats() 
{
    Clear();
}

void RunningStats::Clear()
{
    n = 0;
    w_sum = w_sum2 = mean = S = 0.0;
}

void RunningStats::Push(double x,double w)
{
    double mean_old;

    w_sum = w_sum + w;
    w_sum2 = w_sum2 + w*w;
    mean_old = mean;
    mean = mean_old + (w / w_sum) * (x - mean_old);
    n++;
    
    S = S + w * (x - mean_old) * (x - mean);
    
}

long long RunningStats::NumDataValues() const
{
    return n;
}

double RunningStats::Mean() const
{
    return mean;
}

double RunningStats::Variance() const
{
    return S / (w_sum - 1);
}

double RunningStats::StandardDeviation() const
{
    return sqrt( Variance() );
}

template <typename T>
RunningStats get_stats_obj(const T & array, const Eigen::VectorXd & weights) 
    {
        RunningStats obj;
        u_int64_t size = array.size();
        u_int64_t i;

        for (i=0;i<size;i++)
            {
                obj.Push(array[i],weights.coeff(i));
            }
                
        return obj;
    }

#endif