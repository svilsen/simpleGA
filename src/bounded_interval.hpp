#ifndef boundedInterval
#define boundedInterval

struct Interval 
{
    double lower;
    double upper;
    
    Interval(const double & _lower, const double & _upper);
};

#endif