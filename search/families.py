"""Definitions of four sister conjectures of Zhi-Wei Sun (OEIS, February 2009):

    every odd m beyond a small bound is  p + 2^x + c*2^y,   x, y >= 1,
    with p a prime in a fixed residue class mod 6            (c = 7 or 11).

This file is the single source of truth for the definitions. It
  * writes SPEC files for the C engine (sparse3.c):   python families.py spec NAME LIMIT OUT
  * brute-forces straight from the definition:        python families.py brute NAME LIMIT
                                                      python families.py rep NAME N
"""
import sys
from itertools import count


def pow2_upto(limit, start=1):
    out = []
    for x in count(start):
        if 2**x > limit:
            break
        out.append((x, 2**x))
    return out


def _pairs_pow2(c):
    def pairs(limit):
        for _, u in pow2_upto(limit):
            for _, v in pow2_upto(limit // c):
                if u + c * v <= limit:
                    yield (u, c * v)

    return pairs


FAMILIES = {}


def fam(name, **kw):
    kw.setdefault("r", 0)
    kw.setdefault("parity", 0)
    kw.setdefault("index", lambda n: n)
    kw.setdefault("exceptions", [])
    FAMILIES[name] = kw


odd_index = lambda n: 2 * n - 1  # noqa: E731

# r = required residue of p mod 6; parity = 1: only odd m are tested; minN = first m covered by the claim;
# exceptions = the exceptional m >= minN listed in the conjecture; verified = bound published in OEIS.
fam("A157218", oeis="A157218", pairs=_pairs_pow2(7), r=1, parity=1, index=odd_index, minN=35, verified=15 * 10**7,
    conj="every odd m > 34 is p + 2^x + 7*2^y with p == 1 mod 6  [FALSE: m = 445466268731]")
fam("A157225", oeis="A157225", pairs=_pairs_pow2(7), r=5, parity=1, index=odd_index, minN=21, verified=15 * 10**7,
    exceptions=[25, 31, 1983],
    conj="every odd m > 20 except 25, 31, 1983 is p + 2^x + 7*2^y with p == 5 mod 6  "
         "[FALSE: m = 1433987797, Adamczewski 2026]")
fam("A157237", oeis="A157237", pairs=_pairs_pow2(11), r=1, parity=1, index=odd_index, minN=31, verified=15 * 10**7,
    exceptions=[35, 41, 47, 101, 167, 2021, 119171],
    conj="every odd m > 30 except 35,41,47,101,167,2021,119171 is p + 2^x + 11*2^y with p == 1 mod 6")
fam("A157242", oeis="A157242", pairs=_pairs_pow2(11), r=5, parity=1, index=odd_index, minN=29, verified=10**12,
    exceptions=[33, 39, 45, 171, 247],
    conj="every odd m > 28 except 33,39,45,171,247 is p + 2^x + 11*2^y with p == 5 mod 6")

# ----------------------------------------------------------------------------- tools


def allowed_sums(name, limit):
    return sorted({u + v for u, v in FAMILIES[name]["pairs"](limit)})


def write_spec(name, limit, path):
    f = FAMILIES[name]
    with open(path, "w", newline="\n") as fh:
        fh.write(f"{f['r']} {f['parity']}\n")
        for s in allowed_sums(name, limit):
            fh.write(f"{s}\n")


def brute_exceptions(name, limit):
    """All N <= limit with no representation, straight from the definition."""
    from sympy import isprime

    f = FAMILIES[name]
    sums = allowed_sums(name, limit)
    out = []
    for N in range(1, limit + 1):
        if f["parity"] == 1 and N % 2 == 0:
            continue
        if f["parity"] == 2 and N % 2 == 1:
            continue
        ok = False
        for s in sums:
            p = N - s
            if p < 3:
                break
            if p % 2 == 0 or (f["r"] and p % 6 != f["r"]):
                continue
            if isprime(p):
                ok = True
                break
        if not ok:
            out.append(N)
    return out


def has_rep(name, N):
    """Witness (p, u, v) for a single N, or None. Straight from the definition."""
    from sympy import isprime

    f = FAMILIES[name]
    for u, v in f["pairs"](N):
        p = N - u - v
        if p >= 3 and p % 2 == 1 and (not f["r"] or p % 6 == f["r"]) and isprime(p):
            return (p, u, v)
    return None


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "spec":
        write_spec(sys.argv[2], int(sys.argv[3]), sys.argv[4])
    elif cmd == "list":
        for k, f in FAMILIES.items():
            print(k, f["oeis"], f["verified"], f["conj"], sep="\t")
    elif cmd == "brute":
        print(*brute_exceptions(sys.argv[2], int(sys.argv[3])))
    elif cmd == "rep":
        print(has_rep(sys.argv[2], int(sys.argv[3])))
