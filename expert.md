# Agent Expert Knowledge

- yafu: automated factoring tool that picks the best algorithm (trial division → ECM → SIQS → NFS) based on input size. Current SOTA SIQS implementation. Best single tool for numbers up to ~160 digits.
- cado-nfs: NFS-only implementation with the best polynomial selection, sieving, and linear algebra. Used for all recent factoring records. Better than YAFU's NFS for large numbers, especially with many cores.

