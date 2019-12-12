#ifndef ALIAS_SAMPLER_H__
#define ALIAS_SAMPLER_H__
#include <type_traits>
#include <cstdint>
#include <random>
#include <cstdexcept>
#include <memory>

namespace alias {
template<typename FT, typename IT=std::uint32_t, typename RNG=std::mt19937_64, bool mutable_rng=true>
struct alias_sampler {
    size_t n_;
    std::conditional<mutable_rng, RNG, mutable RNG> rng_;
    mutable std::uniform_real_distribution<FT> urd_;
    std::unique<FT []> prob_;
    std::unique<IT []> alias_;
    template<typename It>
    alias_sampler(It start, It end, uint64_t seed=13): n_(std::distance(start, end))
    {
        alias_ = std::make_unique<IT[]>(n);
        if(n_ > std::numeric_limits<IT>::max()) throw std::out_of_range();
        auto cp(std::make_unique<FT[]>(n));
        FT sum = 0.;
        for(auto it = cp.get(); start != end; ++it, ++start) {
            auto nv = *start++;
            sum += nv;
            *it = nv;
        }
        sum = double(n_) / sum;
        std::for_each(cp.get(), cp.get() + n_, [](auto &x) {x *= sum;});
        auto lb = std::make_unique<IT []>(n_), sb = std::make_unique<IT []>(n_);
        IT nsb = 0, nlb = 0;
        for(IT k = n_; k--;)
            if(cp[k] < 1)
                sb[nsb++] = k;
            else
                lb[nlb++] = k;

        while(nsb && nlb) {
            auto csb = sb[--nsb];
            auto clb = lb[--nlb];
            prob_[csb] = normp(csb);
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
    IT operator()() const {return sample();}
    IT sample() {
        const auto ind = rng_() % n_; // Accelerate with fastrange later
        return urd_(rng_) < prob_[ind] ? ind : alias_[ind];
    }
    IT sample() const {
        if constexpr(mutable_rng) {
            const auto ind = rng_() % n_; // Accelerate with fastrange later
            return urd_(rng_) < prob_[ind] ? ind : alias_[ind];
        } else {
            throw std::runtime_error("Not permitted.");
        }
    }
};
} // alias
#endif
