#ifndef UTILS_HIST_H_
#define UTILS_HIST_H_


#include <vector>
#include <numeric>
#include <fstream>
#include <algorithm>

#include "histogram.hpp"


void inline histogram(std::vector<double> data, double vmin, double vmax, int num_bins, const char *fname) {
    Histogram hist = Histogram(data, vmin, vmax, num_bins, false);
    hist.WriteToFile(fname);
}


#endif  // UTILS_HIST_H_
