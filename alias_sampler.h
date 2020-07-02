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
#include "./tsg.h"

namespace alias {
using std::size_t;
#ifdef roundup
#undef roundup
#define REDEF_roundup
#endif

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
    bool sort_while_sampling = false;
    template<typename It>
    AliasSampler(It start, It end, uint64_t seed=13): n_(std::distance(start, end)), rng_(seed), div_(n_)
    {
        if(n_ > std::numeric_limits<IT>::max()) throw std::out_of_range("n_ > max of int");
        alias_ = std::make_unique<IT[]>(n_);
        auto cp(std::make_unique<FT[]>(n_));
        prob_ = std::make_unique<FT[]>(n_);
        double sum = 0.;
        for(auto it = cp.get(), e = it + n_; it < e; ++it) {
            auto v = *start++;
            sum += v;
            *it = v;
        }
        sum = double(n_) / sum;

#ifdef _BLAZE_MATH_VECTOR_H_
        // Use blaze if available
        blaze::CustomVector<FT, blaze::unaligned, blaze::unpadded> cv(cp.get(), n_);
        cv *= sum;
#elif defined(_OPENMP)
        // Otherwise openmp
        _Pragma("omp parallel for")
        for(size_t i = 0; i < n_; ++i)
            cp[i] *= sum;
#else
        // Otherwise, nothing
        std::for_each(cp.get(), cp.get() + n_, [sum](auto &x) {x *= sum;});
#endif
        std::vector<IT> sb, lb;
        for(IT k = n_; k--; (cp[k] < 1 ? &sb: &lb)->push_back(k));

        while(sb.size() && lb.size()) {
            auto csb = sb.back(), clb = lb.back();
            sb.pop_back();
            lb.pop_back();
            prob_[csb] = cp[csb];
            alias_[csb] = clb;
            ((cp[clb] += cp[csb] - 1.) < 1. ? &sb: &lb)->push_back(clb);
        }
        for(const auto v: lb) prob_[v] = 1.;
        for(const auto v: sb) prob_[v] = 1.;
    }
    IT operator()() const noexcept {return sample();}
    IT operator()()       noexcept {return sample();}
    template<typename It1, typename It2>
    void conditional_sort(It1 it1, It2 it2) const {
#ifndef PDQSORT_H
            std::sort(it1, it2);
#  else
            pdqsort(it1, it2);
#endif
    }
    std::vector<IT> operator()(size_t n) {
        std::vector<IT> ret(n);
        this->operator()(ret.data(), ret.data() + n);
        return ret;
    }
    template<typename P, typename=std::enable_if_t<std::is_integral<P>::value>>
    void operator()(P *beg, P *end, uint64_t seed) {
        this->seed(seed);
        this->operator()(beg, end);
    }
    template<typename P, typename=std::enable_if_t<std::is_integral<P>::value>>
    void sample(P *beg, P *end, uint64_t seed) {
        this->operator()(beg, end, seed);
    }
    template<typename P, typename=std::enable_if_t<std::is_integral<P>::value>>
    void sample(P *beg, P *end) {
        this->operator()(beg, end);
    }
    template<typename P, typename=std::enable_if_t<std::is_integral<P>::value>>
    void operator()(P *beg, P *end) {
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
        for(size_t i = 0, d = end - beg; i < d; ++i) {
            beg[i] = sample();
        }
    }
    template<typename Iterator>
    void operator()(Iterator beg, Iterator end, uint64_t seed=0) {
        if(seed) this->seed(seed);
        while(beg != end) *beg++ = this->sample(); 
    }
    std::vector<IT> operator()(size_t n) const {
        CONST_IF(!mutable_rng) throw std::runtime_error("Not permitted.");
        return const_cast<AliasSampler *>(this)->operator()(n);
    }
    IT sample() noexcept {
#ifdef ALIAS_THREADSAFE
        static thread_local tsg::ThreadSeededGen<RNG> rng;
#else
        auto &rng(rng_);
#endif
        auto ind = div_.mod(rng());
        if(urd_(rng) >= prob_[ind]) ind = alias_[ind];
        assert(ind < n_);
        return ind;
    }
    IT sample() const {
        CONST_IF(!mutable_rng) throw std::runtime_error("Not permitted.");
        return const_cast<AliasSampler *>(this)->sample();
    }
    auto sample(size_t n) const {return this->operator()(n);}
    auto sample(size_t n)       {return this->operator()(n);}
    void seed(uint64_t seed) {this->rng_.seed(seed);}
    size_t size() const {return n_;}
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
    using super = AliasSampler<FT, RNG, IT, mutable_rng>;
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
    std::vector<IT> operator()(size_t n) {
        std::vector<IT> ret(n);
        this->operator()(ret.begin(), ret.end(), n);
        return ret;
    }
    std::vector<IT> operator()(size_t n) const {
        CONST_IF(!mutable_rng) throw std::runtime_error("Not permitted.");
        return const_cast<MaskedAliasSampler *>(this)->operator()(n);
    }
    template<typename Iterator>
    void operator()(Iterator beg, Iterator end, uint64_t seed=0) {
        if(seed) this->seed(seed);
        do *beg++ = this->sample(); while(beg != end);
    }
    template<typename P, typename=std::enable_if_t<std::is_integral<P>::value>>
    void operator()(P *beg, P *end, uint64_t seed) {
        this->seed(seed);
        this->operator()(beg, end);
    }
    template<typename P, typename=std::enable_if_t<std::is_integral<P>::value>>
    void sample(P *beg, P *end, uint64_t seed) {
        this->operator()(beg, end, seed);
    }
    template<typename P, typename=std::enable_if_t<std::is_integral<P>::value>>
    void sample(P *beg, P *end) {
        this->operator()(beg, end);
    }
    template<typename P, typename=std::enable_if_t<std::is_integral<P>::value>>
    void operator()(P *beg, P *end) {
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
        for(size_t i = 0, d = end - beg; i < d; ++i) {
            beg[i] = sample();
        }
    }
    auto sample(size_t n) const {return this->operator()(n);}
    auto sample(size_t n)       {return this->operator()(n);}
    IT sample() noexcept {
#ifdef ALIAS_THREADSAFE
        static thread_local tsg::ThreadSeededGen<RNG> rng;
#else
        auto &rng(this->rng_);
#endif
        auto ind = rng() & bitmask_;
        if(this->urd_(rng) >= this->prob_[ind]) ind = this->alias_[ind];
        return ind;
    }
    IT sample() const {
        CONST_IF(!mutable_rng) throw std::runtime_error("Not permitted.");
        return const_cast<MaskedAliasSampler *>(this)->sample();
    }
    void seed(uint64_t seed) {
        this->rng_.seed(seed);
    }
    size_t size() const {return super::size();}
};
} // alias
#ifdef REDEF_roundup
#  define roundup(x, y) ((((x) % (y)) == 0) ? (x) : ((x) + ((y) - ((x) % (y)))))
#  undef REDEF_roundup
#endif

#endif /* ALIAS_SAMPLER_H__ */
