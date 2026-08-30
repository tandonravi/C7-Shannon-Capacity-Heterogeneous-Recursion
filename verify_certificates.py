#!/usr/bin/env python3
"""Derive the certificates C1-C4 and q15 from the five-dimensional base data.

    python3 verify_certificates.py

This is the second tier of verification.  verify_arithmetic.py *assumes* the
certificate ledger; this script *derives* it, so that no number in the paper
rests on an unchecked count.  Everything it needs is in data/base_c7.json;
there is nothing to download.  Takes about a minute.

The counting machinery used here is validated end to end against brute-force
enumeration by tests/test_c4_bruteforce.py, which runs this same code at base
dimension 1 where the whole pipeline fits in C_7^{box 6}.

    data/base_c7.json  ->  C1, C2, C3, C4, q15  ->  checked against
    data/certificates.json

Method: see the module docstring of shannon/certificates.py.
"""

import json
import sys
from pathlib import Path

from shannon.certificates import compute_all

ROOT = Path(__file__).resolve().parent


def main():
    got = compute_all()
    led = json.load(open(ROOT / "data" / "certificates.json"))

    print("\nAgainst data/certificates.json")
    print("-" * 30)
    failures = []
    for name in ("C1", "C2", "C3", "C4", "q15"):
        printed = led["certificates"][name]["value"]
        mine = got[name]
        ok = mine == printed
        print(f"  {name:<4} derived {mine:>18}   paper {printed:>18}   "
              f"{'ok' if ok else 'MISMATCH'}")
        if not ok:
            failures.append(name)

    expected = tuple(led["base_gadget"]["profile"])
    print(f"\n  base gadget profile {got['base_profile']}  expected {expected}  "
          f"{'ok' if got['base_profile'] == expected else 'MISMATCH'}")
    if got["base_profile"] != expected:
        failures.append("base profile")
    if not all(got["base_axioms"].values()):
        failures.append("base gadget axioms")

    c1 = led["certificates"]["C1"]
    ok = got["C1_pre_exchange"] == c1["pre_exchange_value"]
    print(f"  pre-exchange |{{(TxT)(I_10)}} \\ N(X_10^0)| derived "
          f"{got['C1_pre_exchange']}   ledger {c1['pre_exchange_value']}   "
          f"{'ok' if ok else 'MISMATCH'}")
    if not ok:
        failures.append("pre-exchange C1")

    c4 = led["certificates"]["C4"]
    recorded = {k: v for k, v in c4["component_decomposition"].items()
                if not k.startswith("_")}
    sizes = json.load(open(ROOT / "data" / "paper_claims.json"))["table5_components"]
    sizes = {k: v for k, v in sizes.items() if not k.startswith("_")}
    print("\n  Table 5, per component (size | contribution)")
    for name in recorded:
        sz_ok = got["component_sizes"].get(name) == sizes.get(name)
        ct_ok = got["C4_components"].get(name) == recorded[name]
        print(f"    {name:<16} size {got['component_sizes'].get(name):>19} "
              f"{'ok' if sz_ok else 'MISMATCH':>8}   "
              f"contribution {got['C4_components'].get(name):>18} "
              f"{'ok' if ct_ok else 'MISMATCH':>8}")
        if not sz_ok:
            failures.append(f"size of {name}")
        if not ct_ok:
            failures.append(f"contribution of {name}")
    if got["I30pp"] != c4["size_I30pp"]:
        failures.append("|I_30^++|")
    if got["X30_0L"] != c4["size_X30_0L"]:
        failures.append("|X_30^{0,L}|")
    if sum(got["C4_components"].values()) != got["C4"]:
        failures.append("C4 component sum")
    print(f"\n  |I_30^++|  derived {got['I30pp']}   ledger {c4['size_I30pp']}   "
          f"{'ok' if got['I30pp'] == c4['size_I30pp'] else 'MISMATCH'}")
    print(f"  |X_30^0L|  derived {got['X30_0L']}   ledger {c4['size_X30_0L']}   "
          f"{'ok' if got['X30_0L'] == c4['size_X30_0L'] else 'MISMATCH'}")

    if failures:
        print(f"\n  UNRESOLVED: {', '.join(failures)}")
        print("  See the corresponding 'status' field in data/certificates.json.")
        return 1
    print("\n  all certificates derived and confirmed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
