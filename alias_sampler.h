#ifndef ALIAS_SAMPLER_H__
#define ALIAS_SAMPLER_H__
#include <type_traits>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <memory>
#include "./div.h"

namespace alias {
using std::size_t;
inline constexpr uint64_t roundup(uint64_t i) {
    --i;
    i |= i >> 1;
    i |= i >> 2;
    i |= i >> 4;
    i |= i >> 8;
    i |= i >> 16;
    i |= i >> 32;
    ++i;
    return i;
}
inline constexpr uint32_t roundup(uint32_t i) {
    --i;
    i |= i >> 1;
    i |= i >> 2;
    i |= i >> 4;
    i |= i >> 8;
    i |= i >> 16;
    ++i;
    return i;
}
inline constexpr int32_t roundup(int32_t i) {return roundup(uint32_t(i));}
inline constexpr int64_t roundup(int64_t i) {return roundup(uint64_t(i));}
template<typename T, typename=std::enable_if_t<!std::is_same<T, uint32_t>::value &&
                                               !std::is_same<T, int32_t>::value &&
                                               !std::is_same<T, uint64_t>::value &&
                                               !std::is_same<T, int64_t>::value &&
                                               std::is_integral<T>::value>
        >
inline constexpr T roundup(T i) {
    return roundup(uint64_t(i));
}

template<typename FT=float,
         typename RNG=std::mt19937_64,
         typename IT=std::uint32_t,
         bool mutable_rng=true>
class AliasSampler {
protected:
    size_t n_;
    mutable RNG rng_;
    mutable std::uniform_real_distribution<FT> urd_;
    std::unique_ptr<FT []> prob_;
    std::unique_ptr<IT []> alias_;
    schism::Schismatic<IT> div_;
    AliasSampler(): div_(1) {}
public:
    template<typename It>
    AliasSampler(It start, It end, uint64_t seed=13): n_(std::distance(start, end)), rng_(seed), div_(n_)
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
        for(IT k = n_; k--; (cp[k] < 1 ? &sb: &lb)->push_back(k));

        while(sb.size() && lb.size()) {
            auto csb = sb.back(), clb = lb.back();
            sb.pop_back();
            lb.pop_back();
            prob_[csb] = cp[csb];
            alias_[csb] = clb;
            cp[clb] += cp[csb] - 1.;
            (cp[clb] < 1. ? &sb: &lb)->push_back(clb);
        }
        for(const auto v: lb) prob_[v] = 1.;
        for(const auto v: sb) prob_[v] = 1.;
    }
    IT operator()() const noexcept {return sample();}
    IT operator()()       noexcept {return sample();}
    IT sample() noexcept {
        const auto ind = div_.mod(rng_());
        return urd_(rng_) < prob_[ind] ? ind : alias_[ind];
    }
    IT sample() const {
        CONST_IF(!mutable_rng) throw std::runtime_error("Not permitted.");
        return const_cast<AliasSampler *>(this)->sample();
    }
};

template<typename FT=float,
         typename RNG=std::mt19937_64,
         typename IT=std::uint32_t,
         bool mutable_rng=true>
class MaskedAliasSampler: public AliasSampler<FT, RNG, IT, mutable_rng> {
    /*
     * Same as AliasSampler, but pad with events of 0 probability
     * to round up a power of two so you can bitmask.
     */
protected:
    IT bitmask_;
public:
    template<typename It>
    MaskedAliasSampler(It start, It end, uint64_t seed=13)
    {
        this->rng_.seed(seed);
        if((this->n_ = std::distance(start, end)) > std::numeric_limits<IT>::max())
            throw std::out_of_range("n_ > max of int");
        bitmask_ = roundup(this->n_) - 1;
        auto ru = size_t(bitmask_) + 1;
        this->alias_  = std::make_unique<IT[]>(ru);
        auto cp = std::make_unique<FT[]>(ru);
        this->prob_   = std::make_unique<FT[]>(ru);
        FT sum = 0.;
        IT i;
        for(i = 0; i < this->n_; ++i) {
            auto v = *start++;
            sum += v;
            cp[i] = v;
        }
        std::memset(&cp[this->n_], 0, ru - this->n_);

        for(i = 0, sum = double(ru) / sum; i < this->n_; cp[i++] *= sum);

        std::vector<IT> sb, lb;
        for(i = ru; i--; (cp[i] < 1 ? &sb: &lb)->push_back(i));

        while(sb.size() && lb.size()) {
            auto csb = sb.back(), clb = lb.back();
            sb.pop_back();
            lb.pop_back();
            this->prob_[csb] = cp[csb];
            this->alias_[csb] = clb;
            ((cp[clb] += cp[csb] - 1.) < 1. ? &sb: &lb)->push_back(clb);
        }
        for(const auto v: lb) this->prob_[v] = 1.;
        for(const auto v: sb) this->prob_[v] = 1.;
    }
    IT operator()() const noexcept {return sample();}
    IT operator()()       noexcept {return sample();}
    IT sample() noexcept {
        const auto ind = this->rng_() & bitmask_;
        return this->urd_(this->rng_) < this->prob_[ind] ? ind : this->alias_[ind];
    }
    IT sample() const {
        CONST_IF(!mutable_rng) throw std::runtime_error("Not permitted.");
        return const_cast<MaskedAliasSampler *>(this)->sample();
    }
};
} // alias
#endif
