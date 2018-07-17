#ifndef random_variate_generators
#define random_variate_generators

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
    
    RandomVariates(int N);
};

#endif
