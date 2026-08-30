"""Re-derive the BPZ framework rather than trusting the transcription.

Admissibility (Definition 6) is checked from scratch for all ten combining
rules; pairwise separation is checked for all three terminal codes; and
Proposition 1 is verified at the level of actual sets on a small example.
"""

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shannon import bpz
from shannon.graph import is_independent, separated
from test_gadget_setlevel import base_gadget, set_product


def test_perp_is_symmetric_and_irreflexive():
    for x in bpz.LABELS:
        assert not bpz.perp(x, x), f"{x} _|_ {x} would break admissibility"
        for y in bpz.LABELS:
            assert bpz.perp(x, y) == bpz.perp(y, x)


def test_perp_has_eleven_pairs():
    pairs = {frozenset((x, y)) for x in bpz.LABELS for y in bpz.LABELS
             if bpz.perp(x, y)}
    assert len(pairs) == 11


def test_sigma_is_an_automorphism_of_the_separation_relation():
    """A <-> D, H <-> V preserves _|_.  This is why S3h is admissible."""
    swap = {"A": "D", "D": "A", "H": "V", "V": "H", "B": "B", "N": "N", "O": "O"}
    for x in bpz.LABELS:
        for y in bpz.LABELS:
            assert bpz.perp(x, y) == bpz.perp(swap[x], swap[y])


def test_all_combining_rules_are_admissible():
    for name, rule in bpz.RULES.items():
        problems = bpz.rule_problems(rule)
        assert not problems, f"{name}: {problems[:3]}"


def test_all_terminal_codes_are_pairwise_separated():
    for name, code in bpz.CODES.items():
        assert len(code) == len(set(code)), f"{name}: duplicate words"
        problems = bpz.code_problems(code)
        assert not problems, f"{name}: {problems[:3]}"


def test_word_counts_match_the_paper():
    import json
    claims = json.load(open(Path(__file__).resolve().parent.parent
                            / "data" / "paper_claims.json"))
    expected = claims["rule_word_counts"]
    for name, rule in bpz.RULES.items():
        assert bpz.word_count(rule) == expected[name], name
    for name, code in bpz.CODES.items():
        assert len(code) == expected[name], name


def test_s2a_reproduces_gao_product_as_sets():
    """S2a is not merely numerically equal to Lemma 1 -- it gives the same sets."""
    g = base_gadget()
    F = g.seven_family()

    combined = {}
    for lam in bpz.LABELS:
        acc = set()
        for word in bpz.S2a[lam]:
            for left, right in itertools.product(F[word[0]], F[word[1]]):
                acc.add(left + right)
        combined[lam] = acc

    product_families = set_product(g, g).seven_family()
    assert combined == product_families

    assert bpz.apply_rule(bpz.S2a, [tuple(len(F[c]) for c in bpz.LABELS)] * 2) \
        == tuple(len(product_families[c]) for c in bpz.LABELS)


def test_proposition1_output_is_a_seven_family():
    g = base_gadget()
    F = g.seven_family()
    out = {}
    for lam in bpz.LABELS:
        acc = set()
        for word in bpz.S2a[lam]:
            for a, b in itertools.product(F[word[0]], F[word[1]]):
                acc.add(a + b)
        out[lam] = acc
    for lam in bpz.LABELS:
        assert is_independent(out[lam]), f"F_{lam}^out not independent"
        for mu in bpz.LABELS:
            if bpz.perp(lam, mu):
                assert separated(out[lam], out[mu])


def test_terminal_code_yields_an_independent_set():
    """K3a applied to three copies of the base family, checked as a set."""
    g = base_gadget()
    F = g.seven_family()
    acc = set()
    for word in bpz.K3a:
        for combo in itertools.product(*(F[c] for c in word)):
            acc.add(tuple(x for part in combo for x in part))
    assert is_independent(acc)
    sizes = tuple(len(F[c]) for c in bpz.LABELS)
    assert len(acc) == bpz.apply_code(bpz.K3a, [sizes] * 3)
