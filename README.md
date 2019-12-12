# alias_sampler

Alias sampling, for sampling from discrete probability distributions in constant time.

libstdc++ uses binary search in a cumulative probability array, which both takes logarithmic time and pollutes the cache. This doesn't.
