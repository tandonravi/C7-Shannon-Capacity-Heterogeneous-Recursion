"""Derive the certificates from the base data and check them against the ledger.

These are the slow tests (about a minute).  Run just this file with

    python3 -m pytest tests/test_certificates.py -q

The recomputation is exact and self-contained: everything comes from
data/base_c7.json, which is the BPZ base gadget transcribed into JSON.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shannon.certificates import compute_all, load_base

ROOT = Path(__file__).resolve().parent.parent
LEDGER = json.load(open(ROOT / "data" / "certificates.json"))


@pytest.fixture(scope="module")
def derived():
    return compute_all(verbose=False)


def test_base_gadget_satisfies_every_axiom(derived):
    for name, ok in derived["base_axioms"].items():
        assert ok, name


def test_base_gadget_profile(derived):
    assert derived["base_profile"] == tuple(LEDGER["base_gadget"]["profile"])


def test_transversal_split_matches_the_paper():
    """The paper's J_0 = {0, 5, 6} are the pairs whose centre lies in P^H."""
    base, _, _, _ = load_base()
    pairs = json.load(open(ROOT / "data" / "base_c7.json"))["private_pairs"]
    PH = set(json.load(open(ROOT / "data" / "base_c7.json"))["PH"])
    assert sorted(j for j, (r, _) in enumerate(pairs) if r in PH) == [0, 5, 6]
    assert base.profile[4] == 26 and base.profile[5] == 19


def test_exchange_construction_of_J_plus():
    """Appendix B.1: build() asserts each exchange preserves independence."""
    from shannon.certificates import build
    base, T, exchanges, _ = load_base()
    assert len(exchanges) == 8
    build(base, T, exchanges, verbose=False)  # raises if any check fails


@pytest.mark.parametrize("name", ["C1", "C2", "C3", "C4", "q15"])
def test_confirmed_certificates(derived, name):
    assert derived[name] == LEDGER["certificates"][name]["value"]


def test_C4_component_decomposition(derived):
    """Every entry of Table 5, both columns, plus the two totals."""
    recorded = {k: v for k, v in
                LEDGER["certificates"]["C4"]["component_decomposition"].items()
                if not k.startswith("_")}
    assert derived["C4_components"] == recorded
    assert sum(derived["C4_components"].values()) == derived["C4"]
    assert derived["I30pp"] == 2455726444728097
    assert derived["X30_0L"] == 1446768903083453
    sizes = json.load(open(ROOT / "data" / "paper_claims.json"))["table5_components"]
    assert derived["component_sizes"] == {k: v for k, v in sizes.items()
                                          if not k.startswith("_")}
