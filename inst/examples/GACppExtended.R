##############################################################
########## Simple example using RCpp implementation ##########
##############################################################
population_size = 2
bitstring_size = 10

lower = -pi
upper = pi

maximum_number_of_iterations = 10
maximum_number_of_iterations_equal = 3

tolerance = 1e-3
trace = FALSE

f = function(x) x + abs(sin(4L * x))

GACppExtended(population_size, bitstring_size, 
              f, lower, upper,
              maximum_number_of_iterations, tolerance,
              maximum_number_of_iterations_equal, trace)
