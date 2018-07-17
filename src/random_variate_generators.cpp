#include <ctime>

#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "random_variate_generators.hpp"

RandomVariates::RandomVariates(int N) : 
    rng(std::time(0)),
    uniform_real(), generate_uniform_real(rng, uniform_real),
    uniform_int(0, N - 1), generate_uniform_int(rng, uniform_int),
    uniform_binary(0, 1), generate_uniform_binary(rng, uniform_binary) { };