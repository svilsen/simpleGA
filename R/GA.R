fitness <- function(x) 
    x + abs(sin(4L * x))

negative_fitness <- function(x) 
    -x - abs(sin(4L * x))

decoder <- function(x, lower, upper) {
    N = length(x)
    res = lower + (upper - lower) / (2^N - 1) * sum(2^seq(0, N - 1) * x)
    return(res)
}

initialise_population <- function(population_size, bitstring_size, lower, upper) {
    individuals <-  matrix(sample(c(0, 1), population_size * bitstring_size, replace = TRUE), 
                           nrow = population_size, ncol = bitstring_size)
    decoded_individuals <- apply(individuals, 1, decoder, lower = lower, upper = upper)
    fitness_population <- fitness(decoded_individuals)
    
    res <- list(Individuals = individuals, Fitness = fitness_population,
                Fittest = list(CandidateSolution = decoded_individuals[which.max(fitness_population)],
                               Value = fitness_population[which.max(fitness_population)]))
    return(res)
}

select_parents <- function(current_population) {
    N = length(current_population$Fitness)
    non_zero_fitness = exp(current_population$Fitness)
    sample_probabilities = rmultinom(2 * N, 1, prob = non_zero_fitness / sum(non_zero_fitness))
    
    return(current_population$Individuals[apply(sample_probabilities, 2, which.max), ])
}

crossover_mutation <- function(parents, lower, upper) {
    M = dim(parents)[1] / 2
    N = dim(parents)[2]
    
    res = matrix(0, nrow = 2 * M, ncol = N)
    for (i in 1:M) {
        cp = sample(N, 1)
        c_1 = parents[2 * i - 1, ] * c(rep(1, cp), rep(0, N - cp)) + 
            parents[2 * i, ] * c(rep(0, cp), rep(1, N - cp))
        c_2 = parents[2 * i - 1, ] * c(rep(0, cp), rep(1, N - cp)) + 
            parents[2 * i, ] * c(rep(1, cp), rep(0, N - cp))
        
        pi_m = rbinom(n = N, size = 1, prob = 1 / N)
        m_1 = (c_1 + pi_m) %% 2
        
        pi_m = rbinom(n = N, size = 1, prob = 1 / N)
        m_2 = (c_2 + pi_m) %% 2
        
        res[2 * i - 1, ] = m_1
        res[2 * i, ] = m_2
    }
    
    return(res)
}

survivor_selection <- function(population, children, lower, upper) {
    children_decoded = apply(children, 1, decoder, lower = lower, upper = upper)
    children_fitness = fitness(children_decoded)
    
    total_population = c(population$Fitness, children_fitness)
    
    N = length(population$Fitness)
    sorted_index = order(total_population, decreasing = TRUE)[1:N]
    
    individuals = rbind(population$Individuals, children)[sorted_index, ]
    
    fitness_population = total_population[sorted_index]
    
    fittest_individual = list(CandidateSolution = decoder(individuals[which.max(fitness_population), ], lower, upper),
                              Value = fitness_population[which.max(fitness_population)])
    
    res = list(Individuals = individuals, Fitness = fitness_population,
               Fittest = fittest_individual)
    return(res)
}


#' @title GA (R)
#' 
#' @description Simple genetic algorithm for maximising function. Implemented using the boost and Eigen libraries.
#' 
#' @param population_size The size of the population.
#' @param bitstring_size The size of the individual.
#' @param lower The lower bound of the interval.
#' @param upper The upper bound of the interval.
#' @param maximum_number_of_iterations The maximum allowed number of interations.
#' @param tolerance The convergence tolerance.
#' @param maximum_number_of_iterations_equal The number of iterations without changing the fittest individual before forced convergence.
#' @param trace TRUE/FALSE: Show trace?
#' 
#' @return A list with two elements: the fittest decoded individual found in the entire run and its fitness.
#' @export
#' @example inst/examples/GAR.R
GAR <- function(population_size, bitstring_size,
                lower, upper, maximum_number_of_iterations, tolerance,
                maximum_number_of_iterations_equal, trace = FALSE) {
    population = initialise_population(population_size, bitstring_size, lower, upper)
    
    i = 1
    j = 1
    converged = F
    while (!converged) {
        fitness_old = population$Fittest$Value
        
        ## Parent selection
        parents = select_parents(population)
        
        ## Crossover and mutation
        children = crossover_mutation(parents, lower, upper)
        
        ## Survivor selection
        population = survivor_selection(population, children, lower, upper)
        
        fitness_change = abs(population$Fittest$Value - fitness_old)
        if (fitness_change < tolerance) 
            j = j + 1
        else 
            j = 0
        
        i = i + 1
        converged = (i > maximum_number_of_iterations) | (j > maximum_number_of_iterations_equal)
        
        if (trace)
            cat("Iterations (total):", i, "\n", 
                "Iterations (equal):", j, "\n",
                "Fitness difference:", fitness_change, "\n",
                "-----------------------------\n")
    }
    
    res = unlist(population$Fittest)
    return(res)
}


