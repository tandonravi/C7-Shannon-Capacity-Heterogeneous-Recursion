#!/usr/bin/env python3
"""Reproduce every bound in the paper from the certificate ledger.

    python3 verify_arithmetic.py            full report
    python3 verify_arithmetic.py --table    the Theorem 1 codebook table only

Exact integer arithmetic throughout; decimal expansions are produced by integer
bisection, so every digit printed is correct as printed.  Runs in well under a
second.  By default the certificate ledger is assumed (derive it with
verify_certificates.py), and `--recomputed` re-derives the whole certificate
layer live inside this run; see also
verify_certificates.py.
"""

import argparse
import json
import sys
from pathlib import Path

from shannon.construction import (baselines, load_certificates,
                                  nth_root_digits, theorem2, warmup)
from shannon.gadget import format_codebook_table

ROOT = Path(__file__).resolve().parent


def _rule(title):
    print(f"\n{title}\n" + "-" * len(title))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", action="store_true",
                    help="print only the per-application codebook table")
    ap.add_argument("--digits", type=int, default=50)
    ap.add_argument("--recomputed", action="store_true",
                    help="derive the five certificates live from "
                         "data/base_c7.json (about half a minute) and use "
                         "those values; raises if they differ from the ledger")
    args = ap.parse_args(argv)
    if args.digits < 0:
        ap.error("--digits must be a nonnegative integer")

    led = load_certificates(recomputed=args.recomputed)
    t2 = theorem2(led)

    if args.table:
        print(format_codebook_table(t2["trace"]))
        return 0

    _rule("Certificate ledger"
          + (" (derived live from data/base_c7.json)" if args.recomputed else ""))
    for name, entry in led["certificates"].items():
        print(f"  {name:<4} {entry['value']:>18}   {entry['definition']}")
    if not args.recomputed:
        print("\n  (run verify_certificates.py, or pass --recomputed, to derive"
              "\n   all five from data/base_c7.json rather than trust the ledger)")
    print(f"\n  base gadget  pi(G5) = "
          f"{tuple(led['base_gadget']['profile'])}")

    _rule("Published baselines, recomputed from scratch")
    base = baselines(led)
    for label, value, n in [
            ("Gao [8]                 ", base["gao_a200"], 200),
            ("BPZ v1 [9]              ", base["bpz1_a200"], 200),
            ("BPZ v2 [10]             ", base["bpz2_M"], 500)]:
        print(f"  {label} {nth_root_digits(value, n, args.digits)}")

    _rule("This paper")
    wu = warmup(led)
    print(f"  Section 3.1 (d = 200)    {nth_root_digits(wu['a200'], 200, args.digits)}")
    print(f"  Theorem 2   (d = 500)    {nth_root_digits(t2['M'], 500, args.digits)}")
    print(f"\n  |I| in C_7^500 has {len(str(t2['M']))} decimal digits")

    _rule("Theorem 1 applications (codebooks used at each step)")
    print(format_codebook_table(t2["trace"]))

    _rule("Checks against data/paper_claims.json")
    claims = json.load(open(ROOT / "data" / "paper_claims.json"))
    failures = []
    if t2["M"] != claims["code_sizes"]["eq106_Mstar"]:
        failures.append("M* disagrees with the manuscript")
    if wu["a200"] != claims["code_sizes"]["eq72_a200_warmup"]:
        failures.append("warmup a200 disagrees with the manuscript")
    printed = claims["bounds"]["eq165_theorem2"]
    if nth_root_digits(t2["M"], 500, len(printed.split(".")[1])) != printed:
        failures.append("Theorem 2 decimal expansion disagrees")
    # Exhaustive profile-by-profile checks live in tests/test_paper_claims.py.

    if failures:
        for f in failures:
            print(f"  FAIL  {f}")
        return 1
    print("  all headline claims reproduce exactly")
    print("  (run `python3 -m pytest tests` for the full check)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
