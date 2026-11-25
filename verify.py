"""Standalone check of a claimed counterexample. No dependencies beyond Python 3.

    python verify.py M C R                    e.g.  python verify.py 445466268731 7 1
    python verify.py --certificate M C R      prints a factor of every candidate p (a checkable certificate)

Claim being checked: the odd number M has NO representation
    M = p + 2^x + C*2^y,   x, y >= 1,   p prime,   p == R (mod 6)      (R = 0 means: any odd prime)

Everything is written straight from the definition in the OEIS entry: a double
loop over (x, y) and a deterministic Miller-Rabin test (exact for p < 3.3e24).
It deliberately shares no code with the search engine.
"""
import sys


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:  # these 13 bases are a proven deterministic set below 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def representations(m: int, c: int, r: int):
    out = []
    x = 1
    while 2**x < m:
        y = 1
        while 2**x + c * 2**y < m:
            p = m - 2**x - c * 2**y
            if p % 2 == 1 and p > 2 and (r == 0 or p % 6 == r) and is_prime(p):
                out.append((p, x, y))
            y += 1
        x += 1
    return out


def least_factor(n: int) -> int:
    """Smallest prime factor by trial division (rigorous; fine for n < 1e13)."""
    if n % 2 == 0:
        return 2
    f = 3
    while f * f <= n:
        if n % f == 0:
            return f
        f += 2
    return n


def certificate(m: int, c: int, r: int):
    """Print every candidate p = m - 2^x - c*2^y > 1 in the right class with a proper factor."""
    print(f"# certificate that {m} is not p + 2^x + {c}*2^y with p prime"
          + (f", p == {r} (mod 6)" if r else ""))
    print("# every (x, y) with x, y >= 1 and p > 1 falls in exactly one of these cases:")
    print("#   p not in the required residue class mod 6 (not listed), or")
    print("#   p listed below with its least prime factor q < p, so p is composite.")
    print("# x y p least_prime_factor")
    total = skipped = listed = 0
    x = 1
    while 2**x < m:
        y = 1
        while 2**x + c * 2**y < m:
            p = m - 2**x - c * 2**y
            y += 1
            if p <= 1:
                continue
            total += 1
            if (r and p % 6 != r) or (not r and p % 2 == 0):
                skipped += 1
                continue
            q = least_factor(p)
            assert q < p, f"p = {p} is PRIME: not a counterexample"
            listed += 1
            print(x, y - 1, p, q)
        x += 1
    print(f"# total candidates {total}; wrong residue class {skipped}; composite (listed) {listed}; prime 0")


def main():
    if sys.argv[1] == "--certificate":
        certificate(*(int(a) for a in sys.argv[2:5]))
        return 0
    m, c, r = (int(a) for a in sys.argv[1:4])
    assert m % 2 == 1, "M must be odd"
    # self-test against the first terms of A157218 (c=7, r=1), copied from the OEIS data section
    a157218 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 2, 1, 0, 2, 3, 1, 1, 3, 1, 1, 4, 2, 3, 2, 1, 3,
               3, 2, 3, 5, 1, 2, 5, 2, 4, 5, 1, 4, 3, 1, 4, 7, 1, 5, 7, 2]
    assert [len(representations(2 * n - 1, 7, 1)) for n in range(1, 51)] == a157218, "self-test failed"
    reps = representations(m, c, r)
    cls = f"p == {r} (mod 6)" if r else "p an odd prime"
    print(f"M = {m} (n = {(m + 1) // 2}):  number of representations p + 2^x + {c}*2^y with {cls}: {len(reps)}")
    for p, x, y in reps[:5]:
        print(f"   {m} = {p} + 2^{x} + {c}*2^{y}")
    for d in (-2, 2):
        print(f"   neighbour M{d:+d}: {len(representations(m + d, c, r))} representations")
    print("COUNTEREXAMPLE CONFIRMED" if not reps else "not a counterexample")
    return 0 if not reps else 1


if __name__ == "__main__":
    sys.exit(main())
