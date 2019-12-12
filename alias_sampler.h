#ifndef ALIAS_SAMPLER_H__
#define ALIAS_SAMPLER_H__
#include <type_traits>
#include <algorithm>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <memory>
#include "./div.h"

namespace alias {
template<typename FT, typename IT=std::uint32_t, typename RNG=std::mt19937_64, bool mutable_rng=true>
struct AliasSampler {
    size_t n_;
    mutable RNG rng_;
    mutable std::uniform_real_distribution<FT> urd_;
    std::unique_ptr<FT []> prob_;
    std::unique_ptr<IT []> alias_;
    schism::Schismatic<IT> div_;
    template<typename It>
    AliasSampler(It start, It end, uint64_t seed=13): n_(std::distance(start, end)), div_(n_)
    {
        if(n_ > std::numeric_limits<IT>::max()) throw std::out_of_range("n_ > max of int");
        alias_ = std::make_unique<IT[]>(n_);
        auto cp(std::make_unique<FT[]>(n_));
        prob_ = std::make_unique<FT[]>(n_);
        FT sum = 0.;
        for(IT i = 0; i < n_; ++i) {
            auto v = *start++;
            sum += v;
            cp[i] = v;
        }
        sum = double(n_) / sum;
        
        for(IT i = 0; i < n_; ++i)
            cp[i] *= sum;
        auto lb = std::make_unique<IT []>(n_), sb = std::make_unique<IT []>(n_);
        IT nsb = 0, nlb = 0;
        for(IT k = n_; k--;) {
            if(cp[k] < 1)
                sb[nsb++] = k;
            else
                lb[nlb++] = k;
        }
        while(nsb > 0&& nlb > 0) {
            auto csb = sb[--nsb];
            auto clb = lb[--nlb];
            assert(prob_.get());
            prob_[csb] = cp[csb];
            assert(alias_.get());
            alias_[csb] = clb;
            cp[clb] += cp[csb] - 1.;
            if(cp[clb] < 1) {
                sb[nsb++] = clb;
            } else {
                lb[nlb++] = clb;
            }
        }
        while(nlb) {
            prob_[lb[--nlb]] = 1.;
        }
        while(nsb) {
            prob_[sb[--nsb]] = 1.;
        }
    }
    IT operator()() const {return sample();}
    IT operator()()       {return sample();}
    IT sample() {
        const auto ind = div_.mod(rng_());
        return urd_(rng_) < prob_[ind] ? ind : alias_[ind];
    }
    IT sample() const {
        if constexpr(!mutable_rng) throw std::runtime_error("Not permitted.");
        return const_cast<AliasSampler *>(this)->sample();
    }
};
} // alias
#endif
