import random
from sympy import randprime, isprime
import math
import json


def generate_balanced_semiprime(total_digits):
    """
    Generate a balanced semiprime with exactly `total_digits` decimal digits.
    Both prime factors have roughly total_digits/2 digits.
    Factors are guaranteed not to be too close together (defeats Fermat's method).
    """
    # each factor gets half the digits; if odd, one factor gets the extra digit
    d_small = total_digits // 2
    d_large = total_digits - d_small

    lo_small = 10 ** (d_small - 1)
    hi_small = 10 ** d_small
    lo_large = 10 ** (d_large - 1)
    hi_large = 10 ** d_large

    while True:
        p = randprime(lo_small, hi_small)
        q = randprime(lo_large, hi_large)

        n = p * q

        # check it has the right number of digits
        if len(str(n)) != total_digits:
            continue

        # ensure factors aren't too close (anti-Fermat)
        if abs(p - q) < math.isqrt(n):
            continue

        # return smaller factor first
        if p > q:
            p, q = q, p

        return n, p, q


def generate_test_suite(digit_counts, count_per_size=5, seed=42):
    """
    Generate a frozen test suite of balanced semiprimes.

    Args:
        digit_counts: list of digit lengths to generate, e.g. [30, 40, 50, ...]
        count_per_size: how many semiprimes per digit length
        seed: random seed for reproducibility

    Returns:
        dict mapping digit count -> list of (n, p, q) tuples
    """
    random.seed(seed)
    suite = {}

    for d in sorted(digit_counts):
        print(f"Generating {count_per_size} balanced semiprimes with {d} digits...")
        suite[d] = []
        for i in range(count_per_size):
            n, p, q = generate_balanced_semiprime(d)
            suite[d].append(str(n))
            print(f"  [{i+1}/{count_per_size}] done")

    return suite


def save_suite(suite, filename="test_suite.json"):
    with open(filename, "w") as f:
        json.dump(suite, f, indent=2)
    print(f"\nSaved to {filename}")


def load_suite(filename="test_suite.json"):
    with open(filename) as f:
        return json.load(f)


if __name__ == "__main__":
    digit_counts = list(range(30, 101))
    suite = generate_test_suite(digit_counts, count_per_size=5)
    save_suite(suite)

    # print summary
    print("\n=== Test Suite Summary ===")
    for d in sorted(suite.keys(), key=int):
        entries = suite[d]
        print(f"\n{d} digits:")
        for e in entries:
            print(f"  {e}")
