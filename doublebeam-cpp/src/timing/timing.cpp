#include <iostream>
#include "timing.h"


std::ostream& operator<<(std::ostream& os, const TimingResults& results){
    return os << results.mean << " ns" << " +/- " << results.standard_deviation << "ns ("
    << results.number_of_runs << " runs)";
}
