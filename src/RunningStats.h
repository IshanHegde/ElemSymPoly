#ifndef RUNNINGSTATS_H
#define RUNNINGSTATS_H
#include <cmath>
#include <pybind11/eigen.h>
// creddit to Dr. John D. Cook 
// link: https://www.johndcook.com/blog/skewness_kurtosis/

class RunningStats
{
    public:
        RunningStats();
        void Clear();
        void Push(double x);
        long long NumDataValues() const;
        double Mean() const;
        double Variance() const;
        double StandardDeviation() const;

        friend RunningStats operator+(const RunningStats a, const RunningStats b);
        RunningStats& operator+=(const RunningStats &rhs);

    private:
        long long n;
        double M1, M2;
};

RunningStats::RunningStats() 
{
    Clear();
}

void RunningStats::Clear()
{
    n = 0;
    M1 = M2 = 0.0;
}

void RunningStats::Push(double x)
{
    double delta, delta_n, term1;

    long long n1 = n;
    n++;
    delta = x - M1;
    delta_n = delta / n;
    term1 = delta * delta_n * n1;
    M1 += delta_n;
    M2 += term1;
}

long long RunningStats::NumDataValues() const
{
    return n;
}

double RunningStats::Mean() const
{
    return M1;
}

double RunningStats::Variance() const
{
    return M2/(n-1.0);
}

double RunningStats::StandardDeviation() const
{
    return sqrt( Variance() );
}

template <typename T>
RunningStats get_stats_obj(const T & array) 
    {
        RunningStats obj;

        for (auto x: array)
            {
                obj.Push(x);
            }
                
        return obj;
    }

#endif