#include <Rcpp.h>
#include <RcppEigen.h>

#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

struct RandomVariates 
{
    boost::mt19937 rng;
    
    boost::random::uniform_01<> uniform_real;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_01<> > call_uniform_real;
    
    boost::random::uniform_int_distribution<> uniform_int;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > call_uniform_int;
    
    boost::random::uniform_int_distribution<> uniform_binary;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > call_uniform_binary;
    
    RandomVariates();
    RandomVariates(int seed, int N);
};

RandomVariates::RandomVariates(int seed, int N) : 
    rng(seed), 
    uniform_real(), call_uniform_real(rng, uniform_real),
    uniform_int(0, N - 1), call_uniform_int(rng, uniform_int),
    uniform_binary(0, 1), call_uniform_binary(rng, uniform_binary) {};

struct Interval 
{
    double lower;
    double upper;
    
    Interval();
    Interval(const double & _lower, const double & _upper);
};

Interval::Interval(const double & _lower, const double & _upper) : lower(_lower), upper(_upper) {}; 

struct Population 
{
    int population_size;
    int bitstring_size;
    RandomVariates random_variate;
    Interval interval;
    
    Eigen::MatrixXd individuals;
    Eigen::Vector2d fittest_individual;
    
    double fitness(const double & x) 
    {
        return x + std::abs(std::sin(4 * x));
    }
    
    Population(const int & _population_size, const int & _bitstring_size, 
               const double & _lower, const double & _upper, const int & _seed) 
    {
        population_size = _population_size;
        bitstring_size = _bitstring_size;
        
        random_variate = RandomVariates(_bitstring_size, _seed);
        interval = Interval(_lower, _upper);
        
        fittest_individual[0] = 0.0;
        fittest_individual[1] = fitness(0.0);
        
        individuals = Eigen::MatrixXd(population_size, bitstring_size);
        for (int i = 0; i < population_size; i++) 
        {
            for (int j = 0; j < bitstring_size; j++) 
            {
                individuals(i, j) = random_variate.call_uniform_binary();
            }
        }
    }
    
    
};

//' @title GA (Eigen)
//' 
//' @description Simple genetic algorithm for maximising function. Implemented using the boost and Eigen libraries.
//' 
//' @param population_size The size of the population.
//' @param bitstring_size The size of the individual.
//' @param lower The lower bound of the interval.
//' @param upper The upper bound of the interval.
//' @param maximum_number_of_iterations The maximum allowed number of interations.
//' @param tolerance The convergence tolerance.
//' @param seed A random initial seed.
//' 
//' @return A list with two elements: the fittest decoded individual found in the entire run and its fitness.
//' @export
//[[Rcpp::export()]]
Eigen::Vector2d GACppEigen(const int & population_size, const int & bitstring_size,
                           const double & lower, const double & upper,
                           const int & maximum_number_of_iterations, const double & tolerance,
                           const int & seed) 
{
    Population population(population_size, bitstring_size, lower, upper, seed);
    
    Rcpp::Rcout << "Random binary matrix:" << population.individuals;
    
    int i = 0;
    bool converged = true; // false;
    while (!converged) 
    {
        converged = (i < maximum_number_of_iterations);
    }
    
    return population.fittest_individual;
}

