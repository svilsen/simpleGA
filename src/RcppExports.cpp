// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// GA_Cpp
Eigen::Vector2d GA_Cpp(const int& population_size, const int& bitstring_size, const double& lower, const double& upper, const int& maximum_number_of_iterations, const double& tolerance, const int& maximum_number_of_iterations_equal, const bool& trace);
RcppExport SEXP _simpleGA_GA_Cpp(SEXP population_sizeSEXP, SEXP bitstring_sizeSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP maximum_number_of_iterationsSEXP, SEXP toleranceSEXP, SEXP maximum_number_of_iterations_equalSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type population_size(population_sizeSEXP);
    Rcpp::traits::input_parameter< const int& >::type bitstring_size(bitstring_sizeSEXP);
    Rcpp::traits::input_parameter< const double& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const double& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const int& >::type maximum_number_of_iterations(maximum_number_of_iterationsSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const int& >::type maximum_number_of_iterations_equal(maximum_number_of_iterations_equalSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(GA_Cpp(population_size, bitstring_size, lower, upper, maximum_number_of_iterations, tolerance, maximum_number_of_iterations_equal, trace));
    return rcpp_result_gen;
END_RCPP
}
// GA_Cpp_with_R
Eigen::Vector2d GA_Cpp_with_R(const int& population_size, const int& bitstring_size, Rcpp::Function fitness, const double& lower, const double& upper, const int& maximum_number_of_iterations, const double& tolerance, const int& maximum_number_of_iterations_equal, const bool& trace);
RcppExport SEXP _simpleGA_GA_Cpp_with_R(SEXP population_sizeSEXP, SEXP bitstring_sizeSEXP, SEXP fitnessSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP maximum_number_of_iterationsSEXP, SEXP toleranceSEXP, SEXP maximum_number_of_iterations_equalSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type population_size(population_sizeSEXP);
    Rcpp::traits::input_parameter< const int& >::type bitstring_size(bitstring_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type fitness(fitnessSEXP);
    Rcpp::traits::input_parameter< const double& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const double& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const int& >::type maximum_number_of_iterations(maximum_number_of_iterationsSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const int& >::type maximum_number_of_iterations_equal(maximum_number_of_iterations_equalSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(GA_Cpp_with_R(population_size, bitstring_size, fitness, lower, upper, maximum_number_of_iterations, tolerance, maximum_number_of_iterations_equal, trace));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_simpleGA_GA_Cpp", (DL_FUNC) &_simpleGA_GA_Cpp, 8},
    {"_simpleGA_GA_Cpp_with_R", (DL_FUNC) &_simpleGA_GA_Cpp_with_R, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_simpleGA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
