#ifndef RUNNINGSTATS_H
#define RUNNINGSTATS_H
#include <cmath>
#include <vector>

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
        double M1, M2, M3, M4;
};

RunningStats::RunningStats() 
{
    Clear();
}

void RunningStats::Clear()
{
    n = 0;
    M1 = M2 = M3 = M4 = 0.0;
}

void RunningStats::Push(double x)
{
    double delta, delta_n, delta_n2, term1;

    long long n1 = n;
    n++;
    delta = x - M1;
    delta_n = delta / n;
    delta_n2 = delta_n * delta_n;
    term1 = delta * delta_n * n1;
    M1 += delta_n;
    M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
    M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
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

RunningStats get_stats_obj(const std::vector<double> & array) 
    {
        RunningStats obj;

        const uint64_t size = array.size();

        for(int i=0;i<size;i++)
            {
                obj.Push(array.at(i));
            }
        
        return obj;
    }

#endif