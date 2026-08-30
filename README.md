# Strengthening Recursive Constructions for Zero-Error Shannon Capacity

Companion verification code and finite certificates for the paper
**“Strengthening Recursive Constructions for Zero-Error Shannon Capacity.”**

The repository verifies the construction of an independent set in
\(C_7^{\boxtimes 500}\) that yields

```text
Theta(C_7) >= 3.2588326203532663091215390518104754376053875943219...
```

All certificate values and construction arithmetic are computed with exact
integers. Decimal root expansions are obtained by integer bisection and are
truncated, not rounded.

The runtime code uses only the Python standard library. Python 3.9 or newer is
recommended; `pytest` is needed for the test suite.

## Quick start

From the repository root, run:

```bash
python3 verify_arithmetic.py
python3 verify_certificates.py
python3 verify_arithmetic.py --recomputed
python3 -m pytest tests -q
```

The commands serve different purposes:

- `verify_arithmetic.py` checks the recursive arithmetic using the recorded
  certificate ledger.
- `verify_certificates.py` derives all five finite neighborhood counts from the
  five-dimensional base data.
- `verify_arithmetic.py --recomputed` derives the finite counts and then uses
  those live-derived values in the complete arithmetic chain.
- `python3 -m pytest tests -q` runs the full set-level, arithmetic, provenance,
  and regression test suite.

To print the complete table of heterogeneous codebook parameters used in the
applications of Theorem 1, run:

```bash
python3 verify_arithmetic.py --table
```

The full certificate derivation and test suite may use roughly 2 GiB of memory;
runtime depends on the machine.

## What the repository verifies

The verification has two layers.

### Exact arithmetic layer

Starting from the five recorded finite counts, the code reconstructs:

- the Section 3.1 heterogeneous Gao recursion in dimension 200;
- the role-specific gadgets used in the BPZ recursion;
- the ordered BPZ combining-rule applications;
- the final 257-digit independent-set cardinality in dimension 500;
- the displayed lower bounds for the new construction and the Gao/BPZ
  baselines.

The numerical values printed in the paper are transcribed in
`data/paper_claims.json`, and the tests compare the computed profiles, vectors,
table entries, final integers, and decimal prefixes against that ledger.

### Finite-certificate layer

Starting from `data/base_c7.json`, the code verifies the five-dimensional base
gadget and derives all five nontrivial neighborhood counts:

| Certificate | Value |
|---|---:|
| \(C_1=|J^+\setminus N(X_{10}^0)|\) | 27,488 |
| \(q_{15}=|J_{15}\setminus N(X_{15}^0)|\) | 12,872,271 |
| \(C_2=|J_{15}\setminus N(X_{15}^{0,X})|\) | 12,839,823 |
| \(C_3=|X_{15}\setminus N(X_{15}^{0,X})|\) | 14,045,805 |
| \(C_4=|J_{30}^{++}\setminus N(X_{30}^{0,L})|\) | 841,760,069,965,664 |

The seven individual component contributions to \(C_4\) are also recomputed.

## Set-level checks

The arithmetic implementation is supplemented by finite set-level tests:

- genuine Gao gadgets are built in small powers of \(C_7\), and every clause of
  the gadget definition is checked after ordinary and heterogeneous products;
- Theorem 1 is checked against explicit sets, including its reduction to Gao’s
  ordinary product when all three codebooks coincide;
- all BPZ combining rules and terminal codes are checked for their required
  separation properties;
- the binary BPZ rule `S2a` is checked to produce the same sets as Gao’s product;
- the atom-block counting method used for the large neighborhood certificates
  is compared with brute-force enumeration at toy scale;
- the ordered wiring of the Section 4.2 construction is rebuilt at toy scale.

These finite checks validate the implementation and the particular
constructions used in the paper. The general proof of Theorem 1 is mathematical
and appears in the paper.

## Repository layout

```text
.github/workflows/verify.yml   GitHub Actions verification workflow
data/
  base_c7.json                 five-dimensional BPZ base gadget
  certificates.json            finite certificate ledger
  paper_claims.json            numerical claims transcribed from the paper
scripts/
  make_manifest.py             generate file hashes for a release
  verify_upstream_online.py    optional comparison with the pinned BPZ commit
shannon/
  graph.py                     confusability and neighborhoods in Z_7^d
  counting.py                  product-set neighborhood counting
  certificates.py              derivation of C1-C4 and q15
  gadget.py                    Gao and heterogeneous propagation formulas
  bpz.py                       seven-family rules and terminal codes
  construction.py              the paper’s recursions and published baselines
tests/                         set-level, arithmetic, provenance, and regression tests
upstream/                      pinned BPZ Lean source snapshots
verify_arithmetic.py
verify_certificates.py
```

## Provenance of the base data and BPZ rules

The base gadget and BPZ combining rules are derived from the BPZ Lean
repository at commit

```text
aa21eeb12b75b0413d3fa9fb4208b5d0bf2c4d65
```

Verbatim snapshots of the relevant Lean files are included in `upstream/`.
`tests/test_upstream_pinned.py` parses those snapshots and checks the JSON and
Python transcriptions set-for-set and tuple-for-tuple. The optional command

```bash
python3 scripts/verify_upstream_online.py
```

compares the bundled snapshots with the same files at the pinned upstream
commit. This is the only verification command that requires network access.

## Continuous integration

The workflow in `.github/workflows/verify.yml` runs the arithmetic,
certificate derivation, live-derived arithmetic, parameter table, and full test
suite on every push and pull request.

## Release manifest

To generate `RELEASE_MANIFEST.json` containing SHA-256 hashes of all repository
files and the verification commands, run:

```bash
python3 scripts/make_manifest.py
```

Run this after the repository contents are final for a release.

## Scope

This repository provides exact finite verification and exact-integer
arithmetic for the constructions in the paper. It is not a Lean formalization
of the new heterogeneous theorem or of the complete new result. The bundled
BPZ source snapshots originate from a Lean-verified upstream development; the
new constructions here are checked by the Python verification described above.

## License

Apache License 2.0. See `LICENSE` and `NOTICE`. The redistributed BPZ-derived
base data and rule tables retain their required attribution.

## Citation

Please cite the accompanying paper when using this repository. The underlying
base constructions and recursive frameworks build on work of Polak and
Schrijver, Itty et al., Gao, and Buys, Polak, and Zuiddam.
