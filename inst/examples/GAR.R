###############################################################
############ Simple example using R implementation ############
###############################################################
population_size = 10
bitstring_size = 10

lower = -pi
upper = pi

maximum_number_of_iterations = 10
maximum_number_of_iterations_equal = 3

tolerance = (upper - lower) / (2^bitstring_size - 2)
trace = FALSE

GAR(population_size, bitstring_size, lower, upper,
    maximum_number_of_iterations, tolerance,
    maximum_number_of_iterations_equal, trace)

###############################################################
############# Comparing R and Rcpp implementation #############
###############################################################
library("microbenchmark")

microbenchmark(
    GA_R <- GAR(population_size, bitstring_size, lower, upper,
                maximum_number_of_iterations, tolerance,
                maximum_number_of_iterations_equal, trace),
    GA_Cpp <- GACpp(population_size, bitstring_size, lower, upper,
                    maximum_number_of_iterations, tolerance,
                    maximum_number_of_iterations_equal, trace),
    Optim_R <- optim(runif(1, -pi, pi), simpleGA:::negative_fitness, 
                     method = "L-BFGS-B", lower = lower, upper = upper),
    times = 100
)
