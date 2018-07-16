#include <Rcpp.h>
#include <RcppEigen.h>

#include <ctime>
#include <algorithm>

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
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_01<> > generate_uniform_real;
    
    boost::random::uniform_int_distribution<> uniform_int;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > generate_uniform_int;
    
    boost::random::uniform_int_distribution<> uniform_binary;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > generate_uniform_binary;
    
    RandomVariates(int N) : 
        rng(std::time(0)),
        uniform_real(), generate_uniform_real(rng, uniform_real),
        uniform_int(0, N - 1), generate_uniform_int(rng, uniform_int),
        uniform_binary(0, 1), generate_uniform_binary(rng, uniform_binary) { };
};



struct Interval 
{
    double lower;
    double upper;
    
    Interval(const double & _lower, const double & _upper) : lower(_lower), upper(_upper) { }; 
};

double fitness(const double & x) 
{
    return x + std::abs(std::sin(4.0 * x));
}

struct Population 
{
    // Objects
    int population_size;
    int bitstring_size;
    RandomVariates random_variate;
    Interval interval;
    
    Eigen::MatrixXd individuals;
    Eigen::VectorXd individuals_decoded;
    Eigen::VectorXd individuals_fitness;
    
    Eigen::Vector2d fittest_individual;
    
    // Constructors
    Population(const int & _population_size, const int & _bitstring_size, 
               RandomVariates & _random_variate, Interval & _interval) :
        population_size(_population_size), bitstring_size(_bitstring_size), random_variate(_random_variate), interval(_interval)
    {
        individuals = Eigen::MatrixXd(population_size, bitstring_size);
        individuals_decoded = Eigen::VectorXd(population_size);
        individuals_fitness = Eigen::VectorXd(population_size);
        fittest_individual = Eigen::Vector2d();
        for (int i = 0; i < population_size; i++) 
        {
            for (int j = 0; j < bitstring_size; j++) 
            {
                individuals(i, j) = random_variate.generate_uniform_binary();
            }
            
            const double & individuals_decoded_i = decode(individuals.row(i));
            individuals_decoded[i] = individuals_decoded_i;
            
            const double & individuals_fitness_i = fitness(individuals_decoded[i]);
            individuals_fitness[i] = individuals_fitness_i;
            if (i == 0) 
            {
                fittest_individual[0] = individuals_decoded[0];
                fittest_individual[1] = individuals_fitness[0];
            }
            else if (individuals_fitness[i] > fittest_individual[1]) 
            {
                fittest_individual[0] = individuals_decoded[i];
                fittest_individual[1] = individuals_fitness[i];
            }
        }
    }
    
    Population(const Eigen::MatrixXd & _individuals, 
               RandomVariates & _random_variate, Interval & _interval) :
        population_size(_individuals.rows()), bitstring_size(_individuals.cols()), 
        random_variate(_random_variate), interval(_interval)
    {
        crossover_mutation(_individuals);
        decode_fitness_population();
    }
    
    Population(const Eigen::MatrixXd & _individuals, const Eigen::VectorXd & _individuals_decoded,
               const Eigen::VectorXd & _individuals_fitness, const Eigen::Vector2d _fittest_individual,
               RandomVariates & _random_variate, const Interval & _interval) :
        population_size(_individuals.rows()), bitstring_size(_individuals.cols()),
        random_variate(_random_variate), interval(_interval), individuals(_individuals),
        individuals_decoded(_individuals_decoded), individuals_fitness(_individuals_fitness),
        fittest_individual(_fittest_individual) { }
    
    // Fucntions
    double decode(const Eigen::VectorXd & x) 
    {
        double bit_sum = 0.0;
        for (int i = 0; i < bitstring_size; i++) 
        {
            if (x[i] > 0.0) 
            {
                double bit = std::pow(2.0, i);
                bit_sum += bit;
            }
        }
        
        const double & scaling = std::pow(2.0, bitstring_size) - 1.0;
        const double & decoded = interval.lower + (interval.upper - interval.lower) / scaling * bit_sum;
        return decoded;
    }
    
    void decode_fitness_population() 
    {
        fittest_individual = Eigen::Vector2d();
        fittest_individual << -HUGE_VAL, -HUGE_VAL;
        
        individuals_decoded = Eigen::VectorXd(population_size);
        individuals_fitness = Eigen::VectorXd(population_size);
        for (int i = 0; i < population_size; i++) 
        {
            const double & individuals_decoded_i = decode(individuals.row(i));
            individuals_decoded[i] = individuals_decoded_i;
            
            const double & individuals_fitness_i = fitness(individuals_decoded_i);
            individuals_fitness[i] = individuals_fitness_i;
            
            if (fittest_individual[1] < individuals_fitness_i) 
            {
                fittest_individual[0] = individuals_decoded_i;
                fittest_individual[1] = individuals_fitness_i;
            }
            
        }
    }
    
    void crossover_mutation(const Eigen::MatrixXd & _individuals) 
    {
        const double probability_of_mutation = 1.0 / bitstring_size;
        individuals = Eigen::MatrixXd::Zero(_individuals.rows(), _individuals.cols());

        const int & half_population_size = population_size / 2.0;
        for (int i = 0; i < half_population_size; i++) 
        {
            int c = random_variate.generate_uniform_int();
            for (int j = 0; j < bitstring_size; j++) 
            {
                // Crossover
                if (j <= c)
                {
                    individuals(2 * i, j) = _individuals(2 * i, j);
                    individuals(2 * i + 1, j) = _individuals(2 * i + 1, j);
                }
                else
                {
                    individuals(2 * i, j) = _individuals(2 * i + 1, j);
                    individuals(2 * i + 1, j) = _individuals(2 * i, j);
                }
                
                // Mutation
                for (int k = 0; k < 2; k++) 
                {
                    double u = random_variate.generate_uniform_real();
                    if (u < probability_of_mutation)
                    {
                        const int & individuals_kj = individuals(2 * i + k, j);
                        individuals(2 * i + k, j) = (individuals_kj + 1) % 2;
                    }
                }
            }
        }
    }
};


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


void select_parents(const Population & current_population, Eigen::MatrixXd & parents, RandomVariates & random_variate) 
{
    const int & N = current_population.population_size;
    const Eigen::VectorXd & proportionalSumFitness = proportionalSum(current_population.individuals_fitness);
    for (int i = 0; i < 2 * N; i++) 
    {
        double u = random_variate.generate_uniform_real();
        int index = 0;
        while ((proportionalSumFitness[index + 1] < u) & ((index + 1) < N)) 
            index++;
        
        parents.row(i) = current_population.individuals.row(index);
    }
}

std::vector<int> sorted_index(const Eigen::VectorXd & x)
{
    std::vector<int> x_sorted(x.size());
    std::iota(x_sorted.begin(), x_sorted.end(), 0);
    auto comparator = [&x](int i, int j){ return x[i] > x[j]; };
    
    std::sort(x_sorted.begin(), x_sorted.end(), comparator);
    
    return x_sorted;
}

void survivor_selection(Population & current_population, const Population & children, RandomVariates & random_variate) 
{
    Eigen::VectorXd total_fitness(current_population.population_size + children.population_size);
    total_fitness << current_population.individuals_fitness, 
                     children.individuals_fitness;
    
    std::vector<int> sorted_fitness = sorted_index(total_fitness);
    
    Eigen::MatrixXd new_individuals(current_population.population_size, current_population.bitstring_size);
    Eigen::VectorXd new_individuals_decoded(current_population.population_size);
    Eigen::VectorXd new_individuals_fitness(current_population.population_size);
    for (int i = 0; i < current_population.population_size; i++) 
    {
        int index = sorted_fitness[i];
        if (index < current_population.population_size) 
        {
            new_individuals.row(i) = current_population.individuals.row(index);
            new_individuals_decoded[i] = current_population.individuals_decoded[index];
            new_individuals_fitness[i] = current_population.individuals_fitness[index];
        }
        else 
        {
            index -= current_population.population_size;
            new_individuals.row(i) = children.individuals.row(index);
            new_individuals_decoded[i] = children.individuals_decoded[index];
            new_individuals_fitness[i] = children.individuals_fitness[index];
        }
    }
    
    if (current_population.fittest_individual[1] < new_individuals_fitness[0]) 
    {
        current_population.fittest_individual[0] = new_individuals_decoded[0];
        current_population.fittest_individual[1] = new_individuals_fitness[0];
    }
    
    current_population.individuals = new_individuals;
    current_population.individuals_decoded = new_individuals_decoded;
    current_population.individuals_fitness = new_individuals_fitness;
}

//' @title GA (Cpp)
//' 
//' @description Simple genetic algorithm for maximising function. Implemented using the boost and Eigen libraries.
//' 
//' @param population_size The size of the population.
//' @param bitstring_size The size of the individual.
//' @param lower The lower bound of the interval.
//' @param upper The upper bound of the interval.
//' @param maximum_number_of_iterations The maximum allowed number of interations.
//' @param tolerance The convergence tolerance.
//' @param maximum_number_of_iterations_equal The number of iterations without changing the fittest individual before forced convergence.
//' @param trace TRUE/FALSE: Show trace?
//' 
//' @return A list with two elements: the fittest decoded individual found in the entire run and its fitness.
//' @export
//' @example inst/examples/GACpp.R
//[[Rcpp::export()]]
Eigen::Vector2d GACpp(const int & population_size, const int & bitstring_size,
                      const double & lower, const double & upper,
                      const int & maximum_number_of_iterations, const double & tolerance,
                      const int & maximum_number_of_iterations_equal, const bool & trace = false) 
{
    Interval interval(lower, upper);
    RandomVariates random_variates(bitstring_size);
    Population population(population_size, bitstring_size, random_variates, interval);
    
    int i = 0;
    int j = 0;
    bool converged = false;
    while (!converged) 
    {
        double fitness_old = population.fittest_individual[1];
        
        // Parent selection
        Eigen::MatrixXd parents(2 * population_size, bitstring_size);
        select_parents(population, parents, random_variates);
        
        // Crossover and mutation
        Population children(parents, random_variates, interval);
        
        // Survivor selection
        survivor_selection(population, children, random_variates);
        
        double fitness_change = std::abs(population.fittest_individual[1] - fitness_old);
        if (fitness_change < tolerance) 
            j++;
        else 
            j = 0;
        
        i++;
        converged = (i > maximum_number_of_iterations) | (j > maximum_number_of_iterations_equal);
        
        if (trace)
            Rcpp::Rcout << "Iterations (total): " << i << "\n"
                        << "Iterations (equal): " << j << "\n" 
                        << "Fitness difference: " << fitness_change << "\n"
                        << "-----------------------------\n";
    }
    
    return population.fittest_individual;
}

