# Expert Knowledge

## Why factoring is hard

The fundamental bottleneck in all known classical factoring algorithms is **smoothness detection**: finding integers whose prime factorization consists entirely of small primes. The probability that a random number near N is B-smooth is roughly u^(-u) where u = log(N)/log(B). This is what makes the best known algorithms sub-exponential rather than polynomial.

- The best known algorithm (NFS, 1993) achieves L[1/3]. No improvement in the L exponent in 30+ years.
- Any approach that still relies on finding smooth numbers will hit the same sub-exponential wall.

A genuine improvement must either:
1. Find a way to avoid smoothness detection entirely, or
2. Exploit algebraic structure that makes certain numbers smooth with much higher probability, or
3. Reduce factoring to a different problem that doesn't depend on smoothness.

## What Shor's algorithm tells us

Shor's algorithm factors N in polynomial time by reducing factoring to **period-finding in Z_N***. Given a random a, it finds the order r such that a^r ≡ 1 (mod N). If r is even and a^(r/2) ≢ -1, then gcd(a^(r/2) ± 1, N) yields a factor.

The quantum Fourier transform finds this period efficiently. The open question: **is there a classical way to exploit the same group-theoretic structure?** Nobody has proven there isn't.

Key observations:
- The multiplicative group Z_N* has hidden structure: it's isomorphic to Z_(p-1) × Z_(q-1), but this decomposition is unknown without knowing p and q.
- Period-finding is equivalent to the hidden subgroup problem on Z.
- The Chinese Remainder Theorem means arithmetic mod N splits into independent arithmetic mod p and mod q — but we can't observe this split directly.

## Directions worth exploring

- **Classical period-finding**: Can the order of elements in Z_N* be found without quantum computation? Partial period information? Statistical signatures of the period?
- **Algebraic geometry approaches**: Elliptic curves over Z/NZ split mod p and mod q. Can higher-genus curves, abelian varieties, or p-adic methods reveal this splitting?
- **Representation theory**: The group Z_N* has representations that encode the factorization. Can classical harmonic analysis on finite groups find useful structure?
- **Lattice-based approaches**: Reformulating factoring as a lattice problem in a way that avoids dimension blowup.
- **Exploiting balanced structure**: N = p*q where p ≈ q ≈ sqrt(N). No known algorithm uses this. Does it help?
- **Information-theoretic approaches**: Factoring leaks information through many channels (Jacobi symbols, quadratic residues, discrete logs of small primes). Can these be combined?

## Scaling Data

Best known baseline (YAFU SIQS, worst of 5 semiprimes per size):

| Digits | Time |
|--------|------|
| 50 | 0.12s |
| 60 | 0.7s |
| 70 | 5.8s |
| 80 | 44s |
| 89 | 294s |

Each +10 digits roughly multiplies time by 8-10x — consistent with L[1/2] scaling.
