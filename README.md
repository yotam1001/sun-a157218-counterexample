# A counterexample to a conjecture of Zhi-Wei Sun (OEIS A157218)

**Conjecture** (Zhi-Wei Sun, 24 Feb 2009, [OEIS A157218](https://oeis.org/A157218)).
Every odd integer greater than 34 can be written as

> p + 2ˣ + 7·2ʸ, with p a prime congruent to 1 modulo 6 and x, y positive integers.

It had been verified for odd integers below 1.5·10⁸ (Sun; Qing-Hu Hou).

**Result.** The conjecture is false. The odd number

> **m = 445466268731 = 11³ · 2999 · 111599**

has no such representation. An exhaustive search shows that it is the **only** counterexample below 10¹².

In OEIS terms: a(222733134366) = 0, since 2·222733134366 − 1 = 445466268731.

## Check it yourself in one minute

Plain Python 3, no packages:

```
python verify.py 445466268731 7 1
```

```
M = 445466268731 (n = 222733134366):  number of representations p + 2^x + 7*2^y with p == 1 (mod 6): 0
   neighbour M-2: 42 representations
   neighbour M+2: 87 representations
COUNTEREXAMPLE CONFIRMED
```

`verify.py` is a double loop over (x, y) written directly from the OEIS definition, with its own
deterministic Miller–Rabin test. Before doing anything it re-derives the first 50 terms of the
OEIS data section, so you know it computes the same sequence. It shares no code with the search.

The same check in PARI/GP:

```
cnt(m)=my(c=0,q);for(x=1,logint(m,2),for(y=1,logint(m,2),q=m-2^x-7*2^y;if(q>1&&q%6==1&&isprime(q),c++)));c;
vector(50,n,cnt(2*n-1))    \\ the first 50 terms of A157218
cnt(445466268731)          \\ 0
```

## The certificate

There are 1329 pairs (x, y) with p = m − 2ˣ − 7·2ʸ > 1.

- For 987 of them p is not ≡ 1 (mod 6). Since m ≡ 5 (mod 6), the residue of 2ˣ + 7·2ʸ must be 4,
  which happens exactly when x and y are both odd.
- For the remaining 342, p is composite. [`certificate/445466268731.txt`](certificate/445466268731.txt)
  lists each of them with its least prime factor (`python verify.py --certificate 445466268731 7 1`
  regenerates the file). The largest least prime factor is 566999.

So the whole disproof is 342 divisions that anyone can redo.

## How it was found

[`search/sparse3.c`](search/sparse3.c) lists every odd m in a range that has no representation.

1. Sieve the primes ≡ 1 (mod 6) in a sliding window, one bit per odd number.
2. "m is represented" is the OR of shifted copies of that bitmap, one shift per value
   s = 2ˣ + 7·2ʸ. The 256 smallest shifts are applied 64 numbers at a time, in 256 KB segments that
   stay in the L2 cache.
3. The few m still uncovered are finished one by one: bitmap lookups while the shift fits in the
   window, a deterministic 64-bit Miller–Rabin test beyond it.

It scans about 4·10⁹ integers per second on a 6-core desktop; the range up to 10¹² takes four minutes.

Reproduce the full search (Linux, macOS, or Windows with WSL; needs gcc and Python 3 with `sympy`
for the tests):

```
cd search
bash build.sh
python test_engine.py             # engine vs brute force, vs OEIS data, vs known exceptions
python run.py A157218 1 1e12      # writes data/A157218.exceptions.txt
```

`test_engine.py` checks that the engine's list of unrepresentable m below 200000 equals a
brute-force computation written from the definition, and that brute-force counts reproduce the OEIS
data section, for A157218 and three sister sequences (A157225, A157237, A157242). The search log of
our run is in [`search/data/`](search/data/): the only unrepresentable odd numbers below 10¹² are
1, 3, …, 21, 27, 33 (all below Sun's bound of 34) and 445466268731.

The same engine, run on the sister conjecture A157225 (p ≡ 5 mod 6), rediscovers the counterexample
1433987797 published by Adamczewski in 2026 — an end-to-end test on a case we did not find ourselves.

## Why this number?

The 342 candidates are p = m − 2·4ᵃ − 14·4ᵇ. For a prime q, whether q divides a candidate depends
only on a and b modulo the order of 4 mod q, so each small prime removes a fixed proportion of the
candidates — and that proportion depends on m mod q. The number 445466268731 lies in the *worst*
residue class for 5, 7, 11, 13, 31 and 41 at the same time:

| q | order of 4 | m mod q | candidates divisible by q | for a random integer |
|---|---|---|---|---|
| 5 | 2 | 1 | 25.0 % | 20.0 % |
| 7 | 3 | 1 | 33.3 % | 14.3 % |
| 11 | 5 | 0 | 20.0 % | 9.1 % |
| 13 | 6 | 3 | 8.3 % | 7.7 % |
| 31 | 5 | 15 | 8.0 % | 3.2 % |
| 41 | 10 | 27 | 4.0 % | 2.4 % |

(11 divides m itself, which is why 11 hits so hard: 11 | p exactly when b − a lies in one class mod 5.)
Together the small primes cut the expected number of prime candidates to about 0.62 of normal —
only 1 % of all m are that unlucky — and the 41 candidates that survive the small primes all
happen to be composite.

[`analysis/structure.py`](analysis/structure.py) makes this quantitative (`python structure.py table`,
`rank`, `predict`). A standard independent-primes heuristic then says:

- counterexamples should be **rare but never stop**: the expected number in [10ᵏ, 10ᵏ⁺¹) is of
  order one and grows slowly with k;
- the model predicts about 1.6 counterexamples per conjecture in [1.5·10⁸, 10¹²). We found 1 for
  A157218, and 4 for its sister A157225: 1433987797 (known since 2026) and three new ones,
  45483661477, 85305551767 and 367084269651.

This is a heuristic, not a theorem. The model is exponentially sensitive to second-order choices
(reasonable variants give between 0.1 and 6), so only its qualitative message should be trusted.
**We do not know how to prove that there are infinitely many counterexamples.** The known proofs
for the classical form p + 2ᵃ + 2ᵇ (Crocker 1971; Pan 2011) rely on the factorisation of 2ᵃ + 2ᵇ,
which has no obvious analogue for 2ˣ + 7·2ʸ, and a covering-congruence proof in the style of Erdős
would have to cover all pairs (a, b).

## Files

| Path | What |
|---|---|
| `verify.py` | Standalone verifier and certificate generator |
| `certificate/445466268731.txt` | A prime factor for each of the 342 candidates |
| `search/` | The search engine, its definitions, tests, and the log of the run to 10¹² |
| `analysis/structure.py` | Local-density analysis and the heuristic count of counterexamples |
| `paper/` | A short write-up: `paper.pdf` and its LaTeX source |

## Credits

Found on 25 November 2026 by Yotam Moshe, Idan Gan Levi, Ohav Betito, Boaz Cohen and Shay Gez Akiva.
The search code, the analysis and this write-up were produced with Claude (Anthropic) as an AI
assistant; every claim above is backed by code in this repository that you can run.

The conjecture is due to Zhi-Wei Sun. Background: Z.-W. Sun, *Mixed sums of primes and other terms*
([arXiv:0901.3075](https://arxiv.org/abs/0901.3075)); R. Crocker, *On a sum of a prime and two powers
of two*, Pacific J. Math. 36 (1971), 103–107.

Code: MIT License. OEIS content is quoted under CC BY-SA 4.0.
