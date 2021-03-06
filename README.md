# alias_sampler 🎭 [![Build Status](https://travis-ci.com/dnbaker/alias_sampler.svg?branch=master)](https://travis-ci.com/dnbaker/alias_sampler)

Alias sampling, for sampling from discrete probability distributions in constant time.

libstdc++ uses binary search in a cumulative probability array, which both takes logarithmic time and pollutes the cache. This doesn't.

AliasSampler uses [fastmod](https://lemire.me/blog/2019/02/08/faster-remainders-when-the-divisor-is-a-constant-beating-compilers-and-libdivide/), while MaskedAliasSampler pads to a power of two with items at probability 0 and uses a bitmask.

These samplers provide single samples (`sampler.sample()` or `sampler()`), a batch of samples (`sampler.sample(1000)`/`sampler(1000)`),
or fills a range spanned by iterators with samples (`sampler.sample(vec.begin(), vec.end(), sampling_seed)`).

Similarly, uses blaze for vector multiplications if available, but not otherwise.

In order to make sampling threadsafe, pass `-DALIAS_THREADSAFE`, which will provide a separate RNG for each thread.
