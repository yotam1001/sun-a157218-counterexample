"""Run the sparse3 engine over a range, in parallel, and collect the exceptions.

    python run.py NAME LO HI [--jobs 6] [--chunk 1e10]

Writes data/NAME.exceptions.txt (all N in [LO, HI) with no representation) and
appends one line per finished chunk to data/NAME.chunks.tsv, so a run can be
resumed: finished chunks are skipped.

On Windows the engine runs inside WSL; elsewhere it runs directly.
"""
import argparse
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import families

HERE = Path(__file__).parent
DATA = HERE / "data"
SPECS = HERE / "specs"
WSL = platform.system() == "Windows"


def engine_cmd(spec_rel, lo, hi, hlog=26, k=256):
    cmd = ["./bin/sparse3", spec_rel, str(lo), str(hi), str(hlog), str(k)]
    if WSL:
        cmd = ["wsl", "-d", os.environ.get("WSL_DISTRO", "Ubuntu-24.04"), "--"] + cmd
    return cmd


def run_chunk(name, lo, hi, hlog=26, k=256):
    """Return (exceptions, stderr_summary_line)."""
    spec_rel = f"specs/{name}.spec"
    p = subprocess.run(engine_cmd(spec_rel, lo, hi, hlog, k), cwd=HERE, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"engine failed on {name} [{lo},{hi}): {p.stderr[-500:]}")
    exc = [int(line.split()[1]) for line in p.stdout.splitlines() if line.startswith("X ")]
    done = [line for line in p.stderr.replace("\x00", "").splitlines() if line.startswith("done")]
    if not done:
        raise RuntimeError(f"engine gave no summary on {name} [{lo},{hi}): {p.stderr[-500:]}")
    return exc, done[-1]


def ensure_spec(name, hi):
    SPECS.mkdir(exist_ok=True)
    path = SPECS / f"{name}.spec"
    marker = SPECS / f"{name}.limit"
    if not path.exists() or not marker.exists() or int(marker.read_text()) < hi:
        families.write_spec(name, hi, path)
        marker.write_text(str(hi))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("lo", type=float)
    ap.add_argument("hi", type=float)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--chunk", type=float, default=1e10)
    a = ap.parse_args()
    lo, hi, chunk = int(a.lo), int(a.hi), int(a.chunk)
    DATA.mkdir(exist_ok=True)
    ensure_spec(a.name, hi)
    chunks_file = DATA / f"{a.name}.chunks.tsv"
    done = {}
    if chunks_file.exists():
        for line in chunks_file.read_text().splitlines():
            c_lo, c_hi, exc, summary = line.split("\t")
            done[(int(c_lo), int(c_hi))] = [int(x) for x in exc.split(",") if x]
    todo = [(x, min(x + chunk, hi)) for x in range(lo, hi, chunk)]
    todo = [c for c in todo if c not in done]
    t0 = time.time()

    def work(c):
        exc, summary = run_chunk(a.name, c[0], c[1])
        return c, exc, summary

    with ThreadPoolExecutor(a.jobs) as ex:
        for c, exc, summary in ex.map(work, todo):
            done[c] = exc
            with chunks_file.open("a") as fh:
                fh.write(f"{c[0]}\t{c[1]}\t{','.join(map(str, exc))}\t{summary}\n")
            print(f"[{time.time() - t0:7.0f}s] {a.name} [{c[0]:.3e},{c[1]:.3e}) exceptions={exc[:8]} | {summary}",
                  flush=True)
    allx = sorted({x for c, v in done.items() for x in v})
    (DATA / f"{a.name}.exceptions.txt").write_text("\n".join(map(str, allx)) + ("\n" if allx else ""))
    covered = sorted(done)
    print(f"{a.name}: chunks={len(covered)} span=[{covered[0][0]},{covered[-1][1]}) exceptions={allx}")


if __name__ == "__main__":
    sys.exit(main())
