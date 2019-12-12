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
template<typename FT=float, typename IT=std::uint32_t, typename RNG=std::mt19937_64, bool mutable_rng=true>
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
        for(auto it = cp.get(), e = it + n_; it < e; ++it) {
            auto v = *start++;
            sum += v;
            *it = v;
        }
        sum = double(n_) / sum;
        
        std::for_each(cp.get(), cp.get() + n_, [sum](auto &x) {x *= sum;});
        std::vector<IT> sb, lb;
        for(IT k = n_; k--;cp[k] < 1 ? sb.push_back(k): lb.push_back(k));

        while(sb.size() && lb.size()) {
            auto csb = sb.back(), clb = lb.back();
            sb.pop_back();
            lb.pop_back();
            prob_[csb] = cp[csb];
            alias_[csb] = clb;
            cp[clb] += cp[csb] - 1.;
            if(cp[clb] < 1)
                sb.push_back(clb);
            else
                lb.push_back(clb);
        }
        for(const auto v: lb) prob_[v] = 1.;
        for(const auto v: sb) prob_[v] = 1.;
    }
    IT operator()() const {return sample();}
    IT operator()()       {return sample();}
    IT sample() {
        const auto ind = div_.mod(rng_());
        return urd_(rng_) < prob_[ind] ? ind : alias_[ind];
    }
    IT sample() const {
        CONST_IF(!mutable_rng) throw std::runtime_error("Not permitted.");
        return const_cast<AliasSampler *>(this)->sample();
    }
};
} // alias
#endif
