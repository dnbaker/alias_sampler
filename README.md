# alias_sampler  [![Build Status](https://travis-ci.com/dnbaker/alias_sampler.svg?branch=master)](https://travis-ci.com/dnbaker/alias_sampler)

Alias sampling, for sampling from discrete probability distributions in constant time.

libstdc++ uses binary search in a cumulative probability array, which both takes logarithmic time and pollutes the cache. This doesn't.

Travis CI included.

The only disadvantage here is the issue of cache. Beyond a certain size, each sample requires a cache miss.
For very large arrays and large numbers of samples, it may be suitable
to create several sub-indexes to sample from in sequence to minimize this.
