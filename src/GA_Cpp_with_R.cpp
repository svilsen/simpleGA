#include <Rcpp.h>
#include <RcppEigen.h>

#include "auxiliary.hpp"
#include "random_variate_generators.hpp"
#include "bounded_interval.hpp"

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
    Rcpp::Function fitness;
    
    // Constructors
    Population(const int & _population_size, const int & _bitstring_size, 
               RandomVariates & _random_variate, Interval & _interval, Rcpp::Function f) :
        population_size(_population_size), bitstring_size(_bitstring_size), random_variate(_random_variate), interval(_interval), 
        fitness(f)
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
            
            Rcpp::NumericVector individuals_fitness_i = fitness(individuals_decoded[i]);
            individuals_fitness[i] = individuals_fitness_i[0];
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
               RandomVariates & _random_variate, Interval & _interval, Rcpp::Function f) :
        population_size(_individuals.rows()), bitstring_size(_individuals.cols()), 
        random_variate(_random_variate), interval(_interval), fitness(f)
    {
        crossover_mutation(_individuals);
        decode_fitness_population();
    }
    
    Population(const Eigen::MatrixXd & _individuals, const Eigen::VectorXd & _individuals_decoded,
               const Eigen::VectorXd & _individuals_fitness, const Eigen::Vector2d _fittest_individual,
               RandomVariates & _random_variate, const Interval & _interval, Rcpp::Function f) :
        population_size(_individuals.rows()), bitstring_size(_individuals.cols()),
        random_variate(_random_variate), interval(_interval), individuals(_individuals),
        individuals_decoded(_individuals_decoded), individuals_fitness(_individuals_fitness),
        fittest_individual(_fittest_individual), fitness(f) { }
    
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
            
            Rcpp::NumericVector individuals_fitness_i = fitness(individuals_decoded_i);
            individuals_fitness[i] = individuals_fitness_i[0];
            
            if (fittest_individual[1] < individuals_fitness_i[0]) 
            {
                fittest_individual[0] = individuals_decoded_i;
                fittest_individual[1] = individuals_fitness_i[0];
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
    
    void select_parents(Eigen::MatrixXd & parents) 
    {
        const int & N = population_size;
        const Eigen::VectorXd & proportionalSumFitness = proportionalSum(individuals_fitness);
        for (int i = 0; i < 2 * N; i++) 
        {
            double u = random_variate.generate_uniform_real();
            int index = 0;
            while ((proportionalSumFitness[index + 1] < u) & ((index + 1) < N)) 
                index++;
            
            parents.row(i) = individuals.row(index);
        }
    }
    
    void survivor_selection(const Population & children) 
    {
        Eigen::VectorXd total_fitness(population_size + children.population_size);
        total_fitness << individuals_fitness, 
                         children.individuals_fitness;
        
        std::vector<int> sorted_fitness = sorted_index(total_fitness);
        
        Eigen::MatrixXd new_individuals(population_size, bitstring_size);
        Eigen::VectorXd new_individuals_decoded(population_size);
        Eigen::VectorXd new_individuals_fitness(population_size);
        for (int i = 0; i < population_size; i++) 
        {
            int index = sorted_fitness[i];
            if (index < population_size) 
            {
                new_individuals.row(i) = individuals.row(index);
                new_individuals_decoded[i] = individuals_decoded[index];
                new_individuals_fitness[i] = individuals_fitness[index];
            }
            else 
            {
                index -= population_size;
                new_individuals.row(i) = children.individuals.row(index);
                new_individuals_decoded[i] = children.individuals_decoded[index];
                new_individuals_fitness[i] = children.individuals_fitness[index];
            }
        }
        
        if (fittest_individual[1] < new_individuals_fitness[0]) 
        {
            fittest_individual[0] = new_individuals_decoded[0];
            fittest_individual[1] = new_individuals_fitness[0];
        }
        
        individuals = new_individuals;
        individuals_decoded = new_individuals_decoded;
        individuals_fitness = new_individuals_fitness;
    }
    
};

//' @title GA (Cpp -- extended)
//' 
//' @description Simple genetic algorithm for maximising function. Implemented using the boost and Eigen libraries.
//' 
//' @param population_size The size of the population.
//' @param bitstring_size The size of the individual.
//' @param fitness A univariate function to be maximised.
//' @param lower The lower bound of the interval.
//' @param upper The upper bound of the interval.
//' @param maximum_number_of_iterations The maximum allowed number of interations.
//' @param tolerance The convergence tolerance.
//' @param maximum_number_of_iterations_equal The number of iterations without changing the fittest individual before forced convergence.
//' @param trace TRUE/FALSE: Show trace?
//' 
//' @return A vector with two elements: the fittest decoded individual found in the entire run and its fitness.
//' @export
//' @example inst/examples/GACppwithR.R
//[[Rcpp::export()]]
Eigen::Vector2d GA_Cpp_with_R(const int & population_size, const int & bitstring_size,
                              Rcpp::Function fitness, const double & lower, const double & upper,
                              const int & maximum_number_of_iterations, const double & tolerance,
                              const int & maximum_number_of_iterations_equal, const bool & trace = false) 
{
    Interval interval(lower, upper);
    RandomVariates random_variates(bitstring_size);
    Population population(population_size, bitstring_size, random_variates, interval, fitness);
    
    int i = 0;
    int j = 0;
    bool converged = false;
    while (!converged) 
    {
        double fitness_old = population.fittest_individual[1];
        
        // Parent selection
        Eigen::MatrixXd parents(2 * population_size, bitstring_size);
        population.select_parents(parents);
        
        // Crossover and mutation
        Population children(parents, random_variates, interval, fitness);
        
        // Survivor selection
        population.survivor_selection(children);
        
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

