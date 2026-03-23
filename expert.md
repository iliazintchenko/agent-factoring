# Expert Knowledge

## Factoring Algorithm Landscape

### Current state of the art (classical)
- **Quadratic Sieve (QS/SIQS)**: Complexity L[1/2, 1+o(1)]. Best known for numbers up to ~100 digits.
- **Number Field Sieve (NFS/GNFS)**: Complexity L[1/3, (64/9)^(1/3)]. Best known for numbers above ~100 digits.
- **No improvement in the L exponent since 1993.** Going from L[1/3] to L[1/4] would be a major result.

### Quantum
- **Shor's algorithm**: Polynomial time O((log N)^3). Reduces factoring to period-finding in Z_N*. The quantum Fourier transform finds this period efficiently. The open question is whether the same algebraic structure can be exploited classically.

### Not competitive for balanced semiprimes
- **ECM**: Complexity depends on smallest factor. For balanced semiprimes (factors ~N/2 digits), ~200x slower than SIQS.
- **Pollard's rho / SQUFOF**: O(N^(1/4)) — too slow above 30-35 digits.
- **Fermat/Lehman**: Only helps when |p-q| < N^(1/3).

## Key Insights

1. **Smoothness probability is the fundamental limit**: Probability of B-smoothness is roughly u^(-u) where u = log(Q)/log(B). This is what makes QS/NFS sub-exponential, not polynomial.
2. **No known algorithm exploits balanced structure**: N = p*q where p ≈ q ≈ √N. All known approaches are structure-agnostic.
3. **Batch smoothness (Bernstein product trees) cannot replace sieving**: Multi-precision GCD overhead far exceeds byte-level sieve operations.
4. **Schnorr lattice factoring does not scale**: Requires ~50K-dimensional lattice for 90d (LLL infeasible).

## Scaling Data

YAFU SIQS baseline (worst of 5 semiprimes per size):

| Digits | Time |
|--------|------|
| 50 | 0.12s |
| 60 | 0.7s |
| 70 | 5.8s |
| 80 | 44s |
| 89 | 294s |

Each digit adds ~15-20% to sieve time — consistent with L[1/2] scaling.

## Open Questions

- Why does Shor's algorithm reduce factoring to period-finding? What about the group structure of Z_N* enables this? Can it be exploited classically?
- Can we achieve L[1/4] or better? What mathematical structure would that require?
- Is there an approach from algebraic geometry, p-adic analysis, class field theory, or other areas that changes the complexity class?
- Can the balanced structure of semiprimes (p ≈ q ≈ √N) be exploited? No known algorithm does this.
