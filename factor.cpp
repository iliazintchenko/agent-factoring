#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

// Simple trial-division factoring for numbers that fit in unsigned __int128
typedef unsigned __int128 u128;

static u128 parse(const char *s) {
    u128 r = 0;
    for (; *s >= '0' && *s <= '9'; s++)
        r = r * 10 + (*s - '0');
    return r;
}

static void print_u128(u128 v) {
    if (v == 0) { putchar('0'); return; }
    char buf[40];
    int i = 0;
    while (v > 0) { buf[i++] = '0' + (int)(v % 10); v /= 10; }
    while (i--) putchar(buf[i]);
}

static u128 isqrt128(u128 n) {
    if (n == 0) return 0;
    u128 x = (u128)sqrtl((long double)n);
    while (x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

static u128 trial_divide(u128 n) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    u128 lim = isqrt128(n);
    for (u128 i = 5; i <= lim; i += 6) {
        if (n % i == 0) return i;
        if (n % (i + 2) == 0) return i + 2;
    }
    return n;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s <number>\n", argv[0]); return 1; }
    u128 n = parse(argv[1]);
    printf("Factoring: %s\n", argv[1]);
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    u128 p = trial_divide(n);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    if (p == n) {
        printf("  prime (or 1)\n");
    } else {
        printf("  ");  print_u128(p);
        printf(" * ");  print_u128(n / p);
        printf("\n");
    }
    printf("  time: %.6f s\n", elapsed);
    return 0;
}
