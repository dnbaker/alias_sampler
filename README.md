# alias_sampler 🎭 [![Build Status](https://travis-ci.com/dnbaker/alias_sampler.svg?branch=master)](https://travis-ci.com/dnbaker/alias_sampler)

Alias sampling, for sampling from discrete probability distributions in constant time.

libstdc++ uses binary search in a cumulative probability array, which both takes logarithmic time and pollutes the cache. This doesn't.

Travis CI included.

These samplers provide both single samples (`sampler.sample()` or `sampler()`) and `sampler.sample(1000)`/`sampler(1000)`).
Beyond a certain each successive sample of one item requires a cache miss; a bulk sample samples a set of integers
in sequence and sorts them in order to do one pass through the alias table in order, which allows for both prefetching and minimizing cache misses.
