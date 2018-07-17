##############################################################
########## Simple example using RCpp implementation ##########
##############################################################
population_size = 3
bitstring_size = 10

lower = -pi
upper = pi

maximum_number_of_iterations = 10
maximum_number_of_iterations_equal = 3

tolerance = 1e-3
trace = FALSE

f = function(x) x + abs(sin(4L * x))

GA_Cpp_with_R(population_size, bitstring_size, 
              f, lower, upper,
              maximum_number_of_iterations, tolerance,
              maximum_number_of_iterations_equal, trace)

#################################################################
########## Secondary example using RCpp implementation ##########
#################################################################
population_size = 5
bitstring_size = 12

lower = -2.0
upper = 2.0

maximum_number_of_iterations = 20
maximum_number_of_iterations_equal = 5

tolerance = 1e-3
trace = FALSE

g = function(x) -x * x

parabola_maximum <- GA_Cpp_with_R(population_size, bitstring_size, 
                                  g, lower, upper,
                                  maximum_number_of_iterations, tolerance,
                                  maximum_number_of_iterations_equal, trace)

cat("argmax:", round(parabola_maximum[1], 4), ":: max:", g(round(parabola_maximum[1], 4)))

#####################################################################
########## Comparing the speed of the Rcpp implementations ##########
#####################################################################
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

microbenchmark(PureCpp <- GA_Cpp(population_size, bitstring_size,
                                 lower, upper,
                                 maximum_number_of_iterations, tolerance,
                                 maximum_number_of_iterations_equal, trace), 
               RtoCpp <- GA_Cpp_with_R(population_size, bitstring_size, 
                                       f, lower, upper,
                                       maximum_number_of_iterations, tolerance,
                                       maximum_number_of_iterations_equal, trace), 
               PureR <- GA_R(population_size, bitstring_size, 
                             f, lower, upper,
                             maximum_number_of_iterations, tolerance,
                             maximum_number_of_iterations_equal, trace),
               times = 400)

