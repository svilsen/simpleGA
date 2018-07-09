##############################################################
########## Simple example using RCpp implementation ##########
##############################################################
population_size = 10
bitstring_size = 10

lower = -pi
upper = pi

maximum_number_of_iterations = 10
maximum_number_of_iterations_equal = 3

tolerance = 1e-3
trace = FALSE

GACpp(population_size, bitstring_size, lower, upper,
      maximum_number_of_iterations, tolerance, 
      maximum_number_of_iterations_equal, trace)
