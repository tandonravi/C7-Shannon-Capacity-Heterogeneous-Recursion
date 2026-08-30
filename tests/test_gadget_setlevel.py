"""Check the propagation formulas against honest set-level computation.

This is the test that matters most.  shannon/gadget.py implements Lemma 1 and
Theorem 1 as arithmetic on six-tuples; nothing in that file knows what a graph
is.  Here we build *actual* Gao gadgets as sets of vertices in small strong
powers of C_7, form the products as sets, verify every axiom of Definition 2 on
the output, and compare the resulting profile against the formula.

If Theorem 1's propagation rules were wrong, this would catch it.
"""

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shannon.graph import (closed_neighborhood, confusable, is_independent,
                           separated, vertices)
from shannon.gadget import (NeutralCodebook, Profile, SideCodebook, gao,
                            hetgao, Phi, sigma)


# --------------------------------------------------------------------------
# reference implementation: Gao gadgets as sets
# --------------------------------------------------------------------------

class SetGadget:
    """A Gao gadget (Definition 2) represented by its actual vertex sets."""

    def __init__(self, d, I, pairs, PH, PV, X):
        self.d = d
        self.I, self.pairs = set(I), list(pairs)
        self.PH, self.PV, self.X = set(PH), set(PV), set(X)
        self.R = {r for r, _ in self.pairs}
        self.Q = {q for _, q in self.pairs}
        self.B = self.I - self.R
        self.NH = closed_neighborhood(self.PH, d)
        self.NV = closed_neighborhood(self.PV, d)
        self.X0 = self.X - (self.NH | self.NV)
        self.XH = self.X & self.NH
        self.XV = self.X & self.NV

    def profile(self, name=""):
        return Profile(len(self.I), len(self.pairs), len(self.X),
                       len(self.X0), len(self.XH), len(self.XV), name=name)

    def axioms(self):
        """Every clause of Definition 2, as (name, bool) pairs."""
        pairs_ok = all(
            r in self.I and q not in self.I
            and {w for w in self.I if confusable(w, q)} == {r}
            for r, q in self.pairs)
        endpoints = [e for p in self.pairs for e in p]
        transversal_ok = (
            len(self.PH) == len(self.PV) == len(self.pairs)
            and all((r in self.PH and q in self.PV)
                    or (q in self.PH and r in self.PV)
                    for r, q in self.pairs))
        return [
            ("I independent", is_independent(self.I)),
            ("private pairs", pairs_ok),
            ("endpoint-disjoint", len(endpoints) == len(set(endpoints))),
            ("P^H independent", is_independent(self.PH)),
            ("P^V independent", is_independent(self.PV)),
            ("complementary transversals", transversal_ok),
            ("X independent", is_independent(self.X)),
            ("X cap N(P^H) cap N(P^V) empty",
             not (self.X & self.NH & self.NV)),
        ]

    def seven_family(self):
        """Table 2(a): the BPZ families of this gadget."""
        return {"B": set(self.B), "N": set(self.X0), "A": set(self.XV),
                "D": set(self.XH), "O": set(self.R), "H": set(self.PH),
                "V": set(self.PV)}


def _left_route(x, right, left):
    if x in right.X0:
        return left.R
    return left.PH if x in right.XH else left.PV


def _right_route(y, left, right):
    if y in left.X0:
        return right.R
    return right.PV if y in left.XH else right.PH


def set_product(left, right, X_new=None):
    """Lemma 1 at the level of sets; X_new overrides X_L x X_R (Theorem 1)."""
    AH = {p + x for x in right.X for p in _left_route(x, right, left)}
    AV = {y + q for y in left.X for q in _right_route(y, left, right)}
    core = {b1 + b2 for b1 in left.B for b2 in right.B}
    pairs = ([(r + x, q + x) for (r, q) in left.pairs for x in right.X0]
             + [(y + r, y + q) for y in left.X0 for (r, q) in right.pairs])
    PH = ({p + x for p in left.PH for x in right.X0}
          | {y + q for y in left.X0 for q in right.PV})
    PV = ({p + x for p in left.PV for x in right.X0}
          | {y + q for y in left.X0 for q in right.PH})
    X = X_new if X_new is not None else {a + b for a in left.X for b in right.X}
    return SetGadget(left.d + right.d, core | AH | AV, pairs, PH, PV, X)


def het_set_product(left, right, J0, JH, JV):
    X = ({u + x for u in J0 for x in right.X0}
         | {u + x for u in JH for x in right.XH}
         | {u + x for u in JV for x in right.XV})
    return set_product(left, right, X_new=X)


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------

def base_gadget():
    """Example 1 of the paper: I = {0,2,4} in C_7, private pair (0,6)."""
    return SetGadget(1, [(0,), (2,), (4,)], [((0,), (6,))],
                     [(0,)], [(6,)], [(1,), (3,), (5,)])


def random_maximal_independent(d, avoid=frozenset(), rng=random):
    pool = [v for v in vertices(d) if v not in avoid]
    rng.shuffle(pool)
    chosen = []
    for v in pool:
        if all(not confusable(v, c) for c in chosen):
            chosen.append(v)
    return set(chosen)


def codebook_decomposition(J, gadget, source=""):
    return NeutralCodebook(len(J), len(J - (gadget.NH | gadget.NV)),
                           len(J & gadget.NH), len(J & gadget.NV),
                           source=source)


def side_codebook(J, gadget, source=""):
    q = len(J - closed_neighborhood(gadget.X0, gadget.d))
    return SideCodebook(len(J), q, source=source)


# --------------------------------------------------------------------------
# tests
# --------------------------------------------------------------------------

def test_base_gadget_satisfies_definition_2():
    g = base_gadget()
    assert all(ok for _, ok in g.axioms()), g.axioms()
    assert g.profile().as_tuple() == (3, 1, 3, 1, 1, 1)


def test_lemma1_reproduces_example_1():
    g = base_gadget()
    prod = set_product(g, g)
    assert all(ok for _, ok in prod.axioms())
    assert prod.profile().as_tuple() == (10, 2, 9, 5, 2, 2)
    assert len(prod.I) == 10, "should recover an optimal 10-word code in C_7^2"
    assert is_independent(prod.I)


def test_lemma1_formula_matches_sets():
    g1 = base_gadget()
    g2 = set_product(g1, g1)
    for left, right in [(g1, g1), (g1, g2), (g2, g1), (g2, g2)]:
        prod = set_product(left, right)
        assert all(ok for _, ok in prod.axioms())
        assert prod.profile().as_tuple() == \
            gao(left.profile(), right.profile()).as_tuple()


def test_theorem1_formula_matches_sets():
    """Random heterogeneous codebooks, checked against brute force."""
    rng = random.Random(20260828)
    g1 = base_gadget()
    g2 = set_product(g1, g1)

    for left, right in [(g1, g1), (g1, g2), (g2, g1), (g2, g2)]:
        for _ in range(5):
            # J_0 must avoid N(P^H) cap N(P^V); the others need only be independent
            J0 = random_maximal_independent(left.d, left.NH & left.NV, rng)
            JH = random_maximal_independent(left.d, rng=rng)
            JV = random_maximal_independent(left.d, rng=rng)
            assert is_independent(J0) and not (J0 & left.NH & left.NV)

            prod = het_set_product(left, right, J0, JH, JV)
            assert all(ok for _, ok in prod.axioms()), \
                f"Theorem 1 output violated Definition 2: {prod.axioms()}"

            predicted = hetgao(
                left.profile(), right.profile(),
                J0=codebook_decomposition(J0, left),
                JH=side_codebook(JH, left),
                JV=side_codebook(JV, left))
            assert prod.profile().as_tuple() == predicted.as_tuple()


def test_theorem1_reduces_to_lemma1():
    """J_0 = J_H = J_V = X_L must recover Gao's product exactly."""
    g1 = base_gadget()
    g2 = set_product(g1, g1)
    for left, right in [(g1, g2), (g2, g2)]:
        prod = het_set_product(left, right, left.X, left.X, left.X)
        assert prod.profile().as_tuple() == \
            gao(left.profile(), right.profile()).as_tuple()


def test_trivial_q_lemma():
    """X independent and X^0 a subset of X implies |X \\ N(X^0)| = s - o."""
    g1 = base_gadget()
    for g in [g1, set_product(g1, g1), set_product(g1, set_product(g1, g1))]:
        q = len(g.X - closed_neighborhood(g.X0, g.d))
        assert q == len(g.X) - len(g.X0)


def test_sigma_swaps_orientation_only():
    """sigma exchanges the transversals, leaving I, the private pairs and X^0 alone.

    The private pairs must NOT be reversed: the centre of a private pair has to
    stay in I, so (q, r) is not a private pair.  An earlier version of this test
    reversed them and so built a gadget that violates Definition 2 -- which the
    test did not notice, because it only inspected X.
    """
    g1 = base_gadget()
    g = set_product(g1, g1)
    flipped = SetGadget(g.d, g.I, g.pairs, g.PV, g.PH, g.X)
    assert all(ok for _, ok in flipped.axioms()), flipped.axioms()
    assert flipped.X0 == g.X0
    assert flipped.XH == g.XV and flipped.XV == g.XH
    assert flipped.profile().as_tuple() == sigma(g.profile()).as_tuple()


def test_table2_conversion_is_a_seven_family():
    """Every one of the eleven separations of Definition 4 holds."""
    from shannon.bpz import LABELS, perp
    g1 = base_gadget()
    for g in [g1, set_product(g1, g1)]:
        F = g.seven_family()
        for lam in LABELS:
            assert is_independent(F[lam]), f"F_{lam} not independent"
        for lam in LABELS:
            for mu in LABELS:
                if perp(lam, mu):
                    assert separated(F[lam], F[mu]), f"F_{lam} not _|_ F_{mu}"
        assert tuple(len(F[c]) for c in LABELS) == Phi(g.profile())
