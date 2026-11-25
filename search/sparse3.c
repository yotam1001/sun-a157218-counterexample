/*
 * sparse3.c - find all N in [lo, hi) that are NOT of the form  N = p + s,
 * where p is an odd prime (optionally p == r mod 6) and s runs over a fixed
 * sparse list of "allowed sums" (e.g. s = 2^x + 7*2^y, or s = F_i + F_j^2).
 *
 * Usage: sparse3 SPEC LO HI [HLOG=26] [K=256]
 *
 * SPEC file:
 *   line 1: r        0 = any odd prime, 1 = primes == 1 mod 6, 5 = primes == 5 mod 6
 *   line 2: parity   0 = test all N, 1 = odd N only, 2 = even N only
 *   rest  : the allowed sums s, one per line, strictly increasing, covering all s < HI
 *
 * Output (stdout): one line "X <N>" per N in [lo, hi) with no representation.
 * Stderr: progress and a "hardest cases" profile.
 *
 * Method. Odd numbers 2g+1 are stored one bit each. For odd N = 2j+1 and even
 * s, p = N - s has index g = j - s/2; for even N = 2j and odd s, g = j - (s+1)/2.
 * So "N is represented" is an OR of shifted copies of the prime bitmap. The
 * first K shifts are done 64 numbers at a time, in 256 KB segments that stay in
 * the L2 cache; the few N still uncovered are finished one by one: bitmap
 * lookups while the shift fits in the history window (2^HLOG bits),
 * deterministic Miller-Rabin beyond it.
 *
 * Build: gcc -O3 -march=native -o sparse3 sparse3.c -lm
 */
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint64_t u64;
typedef int64_t i64;
typedef unsigned __int128 u128;

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}
/* deterministic for all 64-bit n (first 12 primes as bases) */
static int is_prime64(u64 n) {
    static const u64 bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return 0;
    for (int i = 0; i < 12; i++) {
        if (n == bases[i]) return 1;
        if (n % bases[i] == 0) return 0;
    }
    u64 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        u64 x = powmod(bases[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int comp = 1;
        for (int k = 1; k < s; k++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { comp = 0; break; }
        }
        if (comp) return 0;
    }
    return 1;
}

static uint32_t *base_primes;
static size_t n_base;

static void make_base_primes(u64 limit) {
    size_t n = (size_t)limit + 2;
    unsigned char *c = calloc(n, 1);
    size_t cnt = 0;
    for (size_t i = 2; i < n; i++) {
        if (!c[i]) {
            cnt++;
            for (size_t j = i * i; j < n; j += i) c[j] = 1;
        }
    }
    base_primes = malloc(cnt * sizeof(uint32_t));
    n_base = 0;
    for (size_t i = 3; i < n; i++)
        if (!c[i]) base_primes[n_base++] = (uint32_t)i;
    free(c);
}

/* Fill bits [0, nbits) of dst: bit b <-> odd number 2*(G+b)+1, set iff that
 * number is an admissible prime (odd prime, and == r mod 6 if r != 0).
 * G may be negative; bits with G+b < 0 are 0. nbits must be a multiple of 64. */
static void sieve_range(i64 G, size_t nbits, u64 *dst, int r) {
    size_t nw = nbits / 64;
    memset(dst, 0xff, nw * 8);
    i64 gend = G + (i64)nbits;
    if (gend <= 0) { memset(dst, 0, nw * 8); return; }
    u64 nmax = 2 * (u64)(gend - 1) + 1;
    i64 gstart = G < 0 ? 0 : G;
    for (size_t i = 0; i < n_base; i++) {
        u64 p = base_primes[i];
        if (p * p > nmax) break;
        u64 first = p * p; /* odd */
        u64 lo_n = 2 * (u64)gstart + 1;
        if (first < lo_n) {
            first = ((lo_n + p - 1) / p) * p;
            if ((first & 1) == 0) first += p;
        }
        i64 g = (i64)((first - 1) / 2);
        for (i64 b = g - G; b < (i64)nbits; b += (i64)p) dst[b >> 6] &= ~(1ULL << (b & 63));
    }
    /* clear negative g and g == 0 (the number 1) */
    for (i64 b = 0; b < (i64)nbits && G + b <= 0; b++) dst[b >> 6] &= ~(1ULL << (b & 63));
    if (r == 1 || r == 5) {
        /* keep g == 0 mod 3 (n == 1 mod 6) or g == 2 mod 3 (n == 5 mod 6) */
        int want = (r == 1) ? 0 : 2;
        u64 pat[3];
        for (int ph = 0; ph < 3; ph++) { /* ph = g mod 3 of bit 0 */
            u64 m = 0;
            for (int b = 0; b < 64; b++)
                if ((ph + b) % 3 == want) m |= 1ULL << b;
            pat[ph] = m;
        }
        i64 g0m = ((G % 3) + 3) % 3;
        for (size_t w = 0; w < nw; w++) dst[w] &= pat[(g0m + (i64)(w % 3) * 64) % 3];
    }
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s SPEC LO HI [HLOG=26] [K=256]\n", argv[0]);
        return 2;
    }
    u64 lo = strtoull(argv[2], 0, 10), hi = strtoull(argv[3], 0, 10);
    int hlog = argc > 4 ? atoi(argv[4]) : 26;
    int K = argc > 5 ? atoi(argv[5]) : 256;
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("spec"); return 2; }
    int r, parity;
    if (fscanf(f, "%d %d", &r, &parity) != 2) { fprintf(stderr, "bad spec header\n"); return 2; }
    size_t cap = 1 << 16, ns = 0;
    u64 *sums = malloc(cap * 8), v, prev = 0;
    while (fscanf(f, "%" SCNu64, &v) == 1) {
        if (ns && v <= prev) { fprintf(stderr, "sums not strictly increasing\n"); return 2; }
        prev = v;
        if (v >= hi) continue;
        if (ns == cap) { cap *= 2; sums = realloc(sums, cap * 8); }
        sums[ns++] = v;
    }
    fclose(f);

    /* shifts per parity class: cls 0 = odd N (even s, d = s/2), cls 1 = even N (odd s, d = (s+1)/2) */
    u64 *D[2], *SV[2];
    size_t nd[2] = {0, 0};
    for (int c = 0; c < 2; c++) { D[c] = malloc((ns + 1) * 8); SV[c] = malloc((ns + 1) * 8); }
    for (size_t i = 0; i < ns; i++) {
        int c = (sums[i] & 1) ? 1 : 0;
        if (sums[i] == 0) continue;
        SV[c][nd[c]] = sums[i];
        D[c][nd[c]++] = (sums[i] + 1) / 2;
    }

    const size_t H = (size_t)1 << hlog; /* history bits */
    const size_t S = (size_t)1 << 21;   /* segment bits: 256 KB, stays in L2 */
    const size_t B = H;                 /* block bits processed before the history is slid */
    u64 *buf = calloc((H + B) / 64 + 2, 8);
    u64 *cov = malloc(S / 64 * 8);
    if (hlog < 21) { fprintf(stderr, "HLOG must be >= 21\n"); return 2; }
    make_base_primes((u64)sqrt((double)hi + 4.0 * (double)(H + B)) + 16);

    i64 J0 = (i64)((lo / 2) & ~63ULL);
    i64 Jend = (i64)(hi / 2 + 1);
    sieve_range(J0 - (i64)H, H, buf, r);

    u64 n_exc = 0, n_scalar = 0, n_mr = 0;
    /* hardest cases: N whose smallest working sum has the largest rank */
    u64 hard_rank = 0, hard_N = 0, hard_s = 0;
    size_t off = 0; /* bits of the current block already processed */
    for (; J0 < Jend; J0 += (i64)S) {
        /* buffer bit b <-> prime index g = (J0 - off - H) + b; this segment is bits [H+off, H+off+S) */
        u64 *seg = buf + (H + off) / 64;
        sieve_range(J0, S, seg, r);
        seg[S / 64] = 0;
        for (int c = 0; c < 2; c++) {
            if (parity == 1 && c == 1) continue;
            if (parity == 2 && c == 0) continue;
            memset(cov, 0, S / 64 * 8);
            size_t k = 0;
            for (; k < nd[c] && k < (size_t)K && D[c][k] <= H; k++) {
                size_t q0 = H + off - (size_t)D[c][k];
                size_t wi = q0 >> 6;
                int sh = (int)(q0 & 63);
                const u64 *src = buf + wi;
                if (sh == 0) {
                    for (size_t w = 0; w < S / 64; w++) cov[w] |= src[w];
                } else {
                    for (size_t w = 0; w < S / 64; w++)
                        cov[w] |= (src[w] >> sh) | (src[w + 1] << (64 - sh));
                }
            }
            size_t kbit = k;
            for (size_t w = 0; w < S / 64; w++) {
                u64 unc = ~cov[w];
                while (unc) {
                    int t = __builtin_ctzll(unc);
                    unc &= unc - 1;
                    i64 j = J0 + (i64)(w * 64 + t);
                    u64 N = c == 0 ? 2 * (u64)j + 1 : 2 * (u64)j;
                    if (N < lo || N >= hi) continue;
                    n_scalar++;
                    int found = 0;
                    size_t kk = kbit;
                    for (; kk < nd[c]; kk++) {
                        u64 d = D[c][kk];
                        if (SV[c][kk] >= N) break;
                        if (d <= H) {
                            size_t b = H + off + (size_t)(j - J0) - (size_t)d;
                            if (buf[b >> 6] >> (b & 63) & 1) { found = 1; break; }
                        } else {
                            u64 p = N - SV[c][kk];
                            if (p < 3 || (p & 1) == 0) continue;
                            if (r && p % 6 != (u64)r) continue;
                            n_mr++;
                            if (is_prime64(p)) { found = 1; break; }
                        }
                    }
                    if (!found) {
                        printf("X %" PRIu64 "\n", N);
                        fflush(stdout);
                        n_exc++;
                    } else if (kk > hard_rank) {
                        hard_rank = kk; hard_N = N; hard_s = SV[c][kk];
                    }
                }
            }
        }
        off += S;
        if (off == B) { /* slide: keep the last H bits as history */
            memmove(buf, buf + B / 64, H / 64 * 8);
            off = 0;
        }
    }
    fprintf(stderr,
            "done [%" PRIu64 ",%" PRIu64 ") exceptions=%" PRIu64 " scalar=%" PRIu64 " mr=%" PRIu64
            " hardest: N=%" PRIu64 " first_sum=%" PRIu64 " rank=%" PRIu64 "\n",
            lo, hi, n_exc, n_scalar, n_mr, hard_N, hard_s, hard_rank);
    return 0;
}
