"""Why is 445466268731 a counterexample to A157218, and how many more should there be?

    python structure.py table        local densities of the counterexample, prime by prime
    python structure.py rank         how extreme its local factor is among random m
    python structure.py predict      heuristic expected number of counterexamples below X

Setting. For odd m we need a prime p == 1 (mod 6) with p = m - 2^x - 7*2^y, x, y >= 1.
Reducing mod 3 fixes the parities of (x, y):

    m == 5 (mod 6)   x, y both odd     candidates  m - 2*4^a - 14*4^b   (a, b >= 0)
    m == 3 (mod 6)   x, y both even    candidates  m -   4^a -  7*4^b   (a, b >= 1)
    m == 1 (mod 6)   x, y of opposite parity (twice as many candidates)

Local model. For a prime q >= 5 let d = ord_q(4). The residue of a candidate mod q depends
only on (a mod d, b mod d), so the proportion of candidates divisible by q is

    rho_q(r) = #{(a, b) mod d : A*4^a + B*4^b == r (mod q)} / d^2,     r = m mod q,

with (A, B) = (2, 14) or (1, 7). A random integer is divisible by q with probability 1/q, so q
changes the expected number of prime candidates by the factor f_q(r) = (1 - rho_q(r)) / (1 - 1/q).
Treating primes as independent gives the local factor S(m) = prod_q f_q(m mod q) and the
heuristic expected number of representations

    E(m) = S(m) * T(m),     T(m) = sum over admissible (x, y) of 3 / ln(m - 2^x - 7*2^y),

(3/ln p is the density of primes among integers == 1 mod 6). The chance that m has no
representation is then about exp(-E(m)).
"""
import math
import random
import sys

M0 = 445466268731


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i * i :: i] = bytearray(len(s[i * i :: i]))
    return [i for i in range(n + 1) if s[i]]


def order4(q):
    d, v = 1, 4 % q
    while v != 1:
        v = v * 4 % q
        d += 1
    return d


def rho_table(q, A, B):
    """rho_q(r) for r = 0..q-1 as a list of fractions (floats), for candidates m - A*4^a - B*4^b."""
    d = order4(q)
    pw = [pow(4, i, q) for i in range(d)]
    cnt = [0] * q
    for u in pw:
        for v in pw:
            cnt[(A * u + B * v) % q] += 1
    return d, [c / (d * d) for c in cnt]


_MODEL_CACHE = {}


def model_primes(Q=500, max_order=40):
    """Primes used in the model: all q <= Q, plus primes below 2e5 whose order of 4 is <= max_order."""
    key = (Q, max_order)
    if key not in _MODEL_CACHE:
        ps = set(p for p in primes_upto(Q) if p >= 5)
        big = [p for p in primes_upto(200000) if p >= 5]
        for d in range(1, max_order + 1):
            n = 4**d - 1
            ps.update(p for p in big if n % p == 0)
        _MODEL_CACHE[key] = sorted(ps)
    return _MODEL_CACHE[key]


CASES = {5: (2, 14), 3: (1, 7)}  # m mod 6 -> (A, B)


def tables(mod6, Q=500):
    A, B = CASES[mod6]
    out = {}
    for q in model_primes(Q):
        d, rho = rho_table(q, A, B)
        out[q] = (d, rho, [(1 - x) / (1 - 1 / q) for x in rho])
    return out


def S_of(m, tabs):
    s = 1.0
    for q, (_, _, f) in tabs.items():
        s *= f[m % q]
    return s


def T_of(m, mod6=5):
    """Sum of 3/ln(p) over admissible (x, y) (parities fixed by m mod 6), p = m - 2^x - 7*2^y > 6."""
    par = {5: (1, 1), 3: (0, 0)}[mod6]
    t, n = 0.0, 0
    x = 1
    while 2**x < m:
        y = 1
        while 2**x + 7 * 2**y < m:
            if x % 2 == par[0] and y % 2 == par[1]:
                p = m - 2**x - 7 * 2**y
                if p > 6:
                    t += 3 / math.log(p)
                    n += 1
            y += 1
        x += 1
    return t, n


def cmd_table():
    m = M0
    tabs = {q: v for q, v in tables(5, Q=500).items() if q <= 100 or v[0] <= 12}
    print(f"m = {m} = 11^3 * 2999 * 111599,  m mod 6 = {m % 6}  ->  x, y both odd")
    print("primes q <= 100, and larger primes for which 4 has order <= 12:")
    print(f"{'q':>6} {'ord(4)':>6} {'m mod q':>8} {'rho_q(m)':>9} {'max rho':>8} {'1/q':>7} {'factor f_q':>10} {'best f_q':>9}")
    S = 1.0
    for q, (d, rho, f) in tabs.items():
        r = m % q
        S *= f[r]
        star = " *" if abs(rho[r] - max(rho)) < 1e-12 else ""
        print(f"{q:>6} {d:>6} {r:>8} {rho[r]:>9.4f} {max(rho):>8.4f} {1 / q:>7.4f} {f[r]:>10.4f} {min(f):>9.4f}{star}")
    print("(* = m sits in the worst residue class for this prime)")
    tabs_full = tables(5, Q=500)
    S_full = S_of(m, tabs_full)
    t, n = T_of(m)
    print(f"\nadmissible pairs (x, y both odd): {n};  T(m) = sum 3/ln p = {t:.2f}  (expected count for a 'typical' m)")
    print(f"local factor S(m): primes in the table: {S:.4f};  full model: {S_full:.4f}")
    print(f"heuristic expected number of representations E(m) = S*T = {S_full * t:.2f};  exp(-E) = {math.exp(-S_full * t):.2e}")
    # exact version: sieve the actual candidates by the model primes
    surv = 0.0
    ns = 0
    corr = 1.0
    for q in tabs_full:
        corr *= 1 / (1 - 1 / q)
    x = 1
    while 2**x < m:
        y = 1
        while 2**x + 7 * 2**y < m:
            if x % 2 and y % 2:
                p = m - 2**x - 7 * 2**y
                if p > 6 and all(p % q for q in tabs_full):
                    surv += 3 / math.log(p) * corr
                    ns += 1
            y += 1
        x += 1
    print(f"exact sieve of the {n} candidates by the model primes: {ns} survive; expected primes among them = {surv:.2f}; exp(-E) = {math.exp(-surv):.2e}")


def cmd_rank(samples=200000, seed=1):
    random.seed(seed)
    tabs = tables(5, Q=500)
    s0 = S_of(M0, tabs)
    vals = []
    for _ in range(samples):
        m = random.randrange(10**11, 10**12)
        m += (5 - m) % 6
        vals.append(S_of(m, tabs))
    vals.sort()
    below = sum(1 for v in vals if v <= s0)
    q = lambda p: vals[int(p * (len(vals) - 1))]  # noqa: E731
    print(f"S(m) over {samples} random m == 5 (mod 6):  min {vals[0]:.3f}  1% {q(0.01):.3f}  10% {q(0.1):.3f}  "
          f"median {q(0.5):.3f}  90% {q(0.9):.3f}  max {vals[-1]:.3f}  mean {sum(vals) / len(vals):.3f}")
    print(f"S({M0}) = {s0:.4f}: {below} of {samples} samples are this small or smaller ({100 * below / samples:.4f}%)")
    inf = 1.0
    for q_, (_, _, f) in tabs.items():
        inf *= min(f)
    print(f"smallest possible S over all residue classes (model primes): {inf:.4f}")


def log_dist(tabs, step=0.002):
    """Distribution of ln S(m) for m uniform over residues (independent primes): dict bin -> prob."""
    dist = {0: 1.0}
    for q, (_, _, f) in tabs.items():
        vals = {}
        for r in range(q):
            k = round(math.log(f[r]) / step) if f[r] > 0 else None
            if k is None:
                continue  # f = 0 cannot happen for q >= 5 here
            vals[k] = vals.get(k, 0) + 1 / q
        new = {}
        for k1, p1 in dist.items():
            for k2, p2 in vals.items():
                new[k1 + k2] = new.get(k1 + k2, 0) + p1 * p2
        # prune negligible mass in the upper tail only (large S never matters for exp(-S*T))
        dist = new
    return dist, step


def kappa(Q, m):
    """Second-order correction. If candidates are first sieved by all primes <= Q, each survivor is
    prime with probability pi_Q = (3/ln m) * prod_{5<=q<=Q} (1-1/q)^-1, and
    P(no prime) = (1-pi_Q)^(#survivors) = exp(-kappa * E) with kappa = -ln(1-pi_Q)/pi_Q >= 1.
    Q = 1 (no sieving) gives the plain Poisson value kappa ~ 1.06."""
    pi = 3 / math.log(m)
    for q in primes_upto(Q):
        if q >= 5:
            pi /= 1 - 1 / q
    return -math.log(1 - pi) / pi


def cmd_predict():
    print("Heuristic expected number of counterexamples (independent-primes model).")
    print("P(m has no representation) = exp(-kappa * S(m) * T(m)); the model has one loose end, the sieving")
    print("level Q behind kappa (see kappa()). Rows: expected number of counterexamples in each range of m.")
    ranges = [(1.5e8, 1e12), (1e12, 1e13), (1e13, 1e14), (1e14, 1e16), (1e16, 1e18)]
    grid = [1.5e8 * 10 ** (i / 40) for i in range(0, 404)]
    grid = [g for g in grid if g <= 1e18 * 1.0001]
    pts = sorted(set(grid + [r for rr in ranges for r in rr]))
    for mod6, label in ((5, "x, y both odd  (m == 5 mod 6 for A157218, m == 3 mod 6 for A157225)"),
                        (3, "x, y both even (m == 3 mod 6 for A157218, m == 1 mod 6 for A157225)")):
        tabs = tables(mod6, Q=500)
        dist, step = log_dist(tabs)
        items = [(math.exp(k * step), p) for k, p in sorted(dist.items())]
        print()
        print(f"class: {label}")
        print(f"  {'model':<22}" + "".join(f"{f'[{lo:.1e},{hi:.0e})':>20}" for lo, hi in ranges))
        for name, Q in (("Poisson, kappa=1", None), ("sieved to Q=1", 1), ("sieved to Q=100", 100), ("sieved to Q=500", 500)):
            tot = {r: 0.0 for r in ranges}
            for lo, hi in zip(pts[:-1], pts[1:]):
                mid = int(math.sqrt(lo * hi))
                mid += (mod6 - mid) % 6
                t, _ = T_of(mid, mod6)
                k = 1.0 if Q is None else kappa(Q, mid)
                pfail = sum(p * math.exp(-k * sv * t) for sv, p in items)
                for r in ranges:
                    if lo >= r[0] * 0.999999 and hi <= r[1] * 1.000001:
                        tot[r] += pfail * (hi - lo) / 6
            print(f"  {name:<22}" + "".join(f"{tot[r]:>20.3g}" for r in ranges))
    print()
    print("Observed: see README (counts from the exhaustive search).")


if __name__ == "__main__":
    {"table": cmd_table, "rank": cmd_rank, "predict": cmd_predict}[sys.argv[1]]()
