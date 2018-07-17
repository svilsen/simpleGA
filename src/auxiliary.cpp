#include <RcppEigen.h>
#include <algorithm>

Eigen::VectorXd proportionalSum(const Eigen::VectorXd & x)
{
    const std::size_t N = x.size();
    Eigen::VectorXd proportionalSum(N + 1);
    
    proportionalSum[0] = 0;
    for (std::size_t i = 1; i < N + 1; i++)
    {
        proportionalSum[i] = std::exp(x[i - 1]) + proportionalSum[i - 1];
    }
    
    proportionalSum = proportionalSum / proportionalSum[N];
    return proportionalSum;
}

std::vector<int> sorted_index(const Eigen::VectorXd & x)
{
    std::vector<int> x_sorted(x.size());
    std::iota(x_sorted.begin(), x_sorted.end(), 0);
    auto comparator = [&x](int i, int j){ return x[i] > x[j]; };
    
    std::sort(x_sorted.begin(), x_sorted.end(), comparator);
    
    return x_sorted;
}