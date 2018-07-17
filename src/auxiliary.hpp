#ifndef auxiliary
#define auxiliary

#include <RcppEigen.h>

Eigen::VectorXd proportionalSum(const Eigen::VectorXd & x);

std::vector<int> sorted_index(const Eigen::VectorXd & x);

#endif