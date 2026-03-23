# Expert Knowledge

## Schnorr lattice factoring (tested, does not scale)

Schnorr (2021) proposed reducing factoring to finding short vectors in a lattice constructed from the prime factorizations of smooth numbers near sqrt(N). The idea: if you can find a short enough vector, it encodes a multiplicative relation that yields a factor.

Tested on small numbers — it works for 30-40 digits. But the lattice dimension grows roughly as the number of primes in the factor base, which for 90-digit numbers means ~50K dimensions. LLL/BKZ reduction in 50K dimensions is completely infeasible. The approach hits a wall around 50 digits.

**Conclusion**: The lattice dimension is tied to the smoothness bound, which grows sub-exponentially with N. So this doesn't escape the smoothness bottleneck — it just reformulates it as a lattice problem that's equally hard. Any lattice-based approach that depends on a large smooth factor base will have the same problem.
