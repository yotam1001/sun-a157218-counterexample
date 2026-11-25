"""Validate the C engine against (1) brute force from the definitions,
(2) the OEIS data sections, (3) published counterexamples far out.

    python test_engine.py            # all checks (a few minutes)
"""
import json
import sys
import urllib.request

import families
import run

LIMIT = 200_000


def oeis_data(anum):
    url = f"https://oeis.org/search?q=id:{anum}&fmt=json"
    req = urllib.request.Request(url, headers={"User-Agent": "oeis-research-client/0.1"})
    with urllib.request.urlopen(req, timeout=60) as r:
        d = json.loads(r.read().decode())
    d = d if isinstance(d, list) else d["results"]
    return [int(x) for x in d[0]["data"].split(",")], int(d[0]["offset"].split(",")[0])


def count_reps(name, N):
    """a(n) from the definition: number of (p, u, v) with multiplicity as generated."""
    from sympy import isprime

    f = families.FAMILIES[name]
    c = 0
    for u, v in f["pairs"](N):
        p = N - u - v
        if p >= 3 and p % 2 == 1 and (not f["r"] or p % 6 == f["r"]) and isprime(p):
            c += 1
    return c


def main():
    ok = True
    only = sys.argv[1:]
    for name, f in families.FAMILIES.items():
        if only and name not in only:
            continue
        # (1) engine vs brute force on [1, LIMIT)
        run.ensure_spec(name, 10**15)
        brute = families.brute_exceptions(name, LIMIT - 1)
        eng = sorted(run.run_chunk(name, 1, LIMIT)[0])
        eng24 = sorted(run.run_chunk(name, 1, LIMIT, hlog=21, k=3)[0])  # exercise the scalar and Miller-Rabin paths
        status = "ok" if brute == eng == eng24 else "MISMATCH"
        ok &= status == "ok"
        big = [x for x in brute if x >= f["minN"]]
        claimed = sorted(f["exceptions"])
        claim_status = "matches-claim" if big == [x for x in claimed if x < LIMIT] else f"DIFFERS-FROM-CLAIM {big}"
        print(f"{name:12s} engine-vs-brute: {status:8s} exceptions>=minN below {LIMIT}: {big} {claim_status}")
        # (2) OEIS data section
        if f["oeis"]:
            try:
                data, off = oeis_data(f["oeis"])
                mine = [count_reps(name, f["index"](n)) for n in range(off, off + len(data))]
                s = "ok" if mine == data else "MISMATCH"
                ok &= s == "ok"
                print(f"{'':12s} OEIS data ({len(data)} terms from n={off}): {s}")
                if s != "ok":
                    bad = [(off + i, a, b) for i, (a, b) in enumerate(zip(mine, data)) if a != b][:5]
                    print("             first differences (n, ours, oeis):", bad)
            except Exception as e:  # noqa: BLE001
                print(f"{'':12s} OEIS data: could not fetch ({e})")
        # (3) far counterexamples, through the engine and through brute force
        for N in f["exceptions"]:
            if N < LIMIT:
                continue
            eng, _ = run.run_chunk(name, N - 1000, N + 1000)
            w = families.has_rep(name, N)
            s = "ok" if eng == [N] and w is None else f"MISMATCH engine={eng} witness={w}"
            ok &= s == "ok"
            print(f"{'':12s} published exception {N}: {s}")
    print("ALL OK" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
