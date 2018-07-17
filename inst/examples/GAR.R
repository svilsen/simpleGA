###############################################################
############ Simple example using R implementation ############
###############################################################
population_size = 3
bitstring_size = 10

lower = -pi
upper = pi

maximum_number_of_iterations = 10
maximum_number_of_iterations_equal = 3

tolerance = 1e-3
trace = FALSE

f = function(x) x + abs(sin(4L * x))

GA_R(population_size, bitstring_size, 
     f, lower, upper,
     maximum_number_of_iterations, tolerance,
     maximum_number_of_iterations_equal, trace)

###############################################################
############# Comparing R and Rcpp implementation #############
###############################################################
library("microbenchmark")

population_size = 3
bitstring_size = 10

lower = -pi
upper = pi

maximum_number_of_iterations = 10
maximum_number_of_iterations_equal = 3

tolerance = 1e-3
trace = FALSE

f = function(x) x + abs(sin(4L * x))

negative_f = function(x) -x - abs(sin(4L * x))

microbenchmark(
    GAR <- GA_R(population_size, bitstring_size, 
                f, lower, upper,
                maximum_number_of_iterations, tolerance,
                maximum_number_of_iterations_equal, trace),
    GACpp <- GA_Cpp(population_size, bitstring_size, 
                    lower, upper,
                    maximum_number_of_iterations, tolerance,
                    maximum_number_of_iterations_equal, trace),
    GACppwithR <- GA_Cpp_with_R(population_size, bitstring_size, 
                                f, lower, upper,
                                maximum_number_of_iterations, tolerance,
                                maximum_number_of_iterations_equal, trace),
    OptimR <- optim(runif(1, -pi, pi), negative_f, 
                    method = "L-BFGS-B", lower = lower, upper = upper),
    times = 500
)
