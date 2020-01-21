#include "alias_sampler.h"
#include <iostream>

int main() {
    std::vector<double> v{0.10, 0.225, 0.45, 0.10, .125};
    for(auto &i: v)
        i *= 50;
    alias::AliasSampler<float, std::mt19937> as(v.begin(), v.end());
    alias::MaskedAliasSampler<float, std::mt19937> mas(v.begin(), v.end());
    std::vector<int> counts(v.size()), mcounts(v.size());
    unsigned nrounds = 1500;
    for(size_t i = 0; i < nrounds; ++i) {
        ++counts.at(as());
        ++mcounts.at(mas());
    }
    for(auto &i: v)
        i /= 50;
    for(size_t i = 0; i < v.size(); ++i) {
#if VERBOSE_AF
        std::fprintf(stderr, "Expected [%f/%f/%f/%f]. Got [%f/%f/%f/%f].\n",
                     v[0], v[1], v[2], v[3], counts[0] / double(nrounds), counts[1] / double(nrounds), counts[2] / double(nrounds), counts[3] / double(nrounds));
#endif
        assert(std::abs(v[i] - counts[i] / double(nrounds)) / v[i] < 0.15);
        assert(std::abs(v[i] - mcounts[i] / double(nrounds)) / v[i] < 0.15);
    }
    std::vector<int> zomg(1000);
    as(zomg.begin(),zomg.end());
    mas(zomg.begin(),zomg.end());
}
