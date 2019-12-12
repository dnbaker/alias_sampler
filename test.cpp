#include "alias_sampler.h"
#include <iostream>

int main() {
    std::vector<double> v{0.1, 0.4, 0.5, 0.234};
    alias::AliasSampler<float, int> as(v.begin(), v.end());
    std::cout << as.sample() << '\n';
}
