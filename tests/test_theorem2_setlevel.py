"""Rebuild the Section 4.2 wiring at base dimension 1 and check every profile.

test_c4_bruteforce.py validates the C4 branch; this test validates the rest of
the Theorem 2 pipeline.  Every gadget of Step 1 and the reachable products of
Steps 2-4 are built twice, with the exact codebook wiring of Table 6 translated
to the k = 1 analogue: once as honest vertex sets with the reference
implementation, and once through the hetgao/gao formula engine with codebook
parameters MEASURED from the actual sets.  The two must agree, and every
set-built gadget of dimension <= 6 must satisfy all of Definition 2.

Scope: G40^(8) lives in dimension 8 at toy scale, where the pairwise
independence re-check is too slow, so for it only the profile identity is
asserted (its auxiliary machinery is the same hetgao path already axiom-checked
at d <= 6).  The final G55 product is not rebuilt here.  Its one novel
ingredient -- using the auxiliary set of the heterogeneous sibling Ghat_30 as
J0 for the ordinary sibling G_30^L -- is validated explicitly below by
test_final_sibling_J0_transfer; test_c4_bruteforce.py
separately validates the G_30^L neutral reference and G_30^{++} source counts
against enumeration.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from shannon.gadget import gao, hetgao, sigma
from test_gadget_setlevel import (SetGadget, base_gadget,
                                  codebook_decomposition, het_set_product,
                                  set_product, side_codebook)


def T1(c):
    return (2 - c) % 7


def Tw(word):
    return tuple(T1(x) for x in word)


def sigma_set(g):
    return SetGadget(g.d, g.I, g.pairs, g.PV, g.PH, g.X)


def _measured_hetgao(left_set, right_set, J0, JH, JV, name):
    """Theorem 1's formula, with every codebook parameter read off the sets."""
    return hetgao(left_set.profile(), right_set.profile(),
                  J0=codebook_decomposition(J0, left_set),
                  JH=side_codebook(JH, left_set),
                  JV=side_codebook(JV, left_set), name=name)


_CACHE = {}


def build_pipeline():
    if "P" in _CACHE:
        return _CACHE["P"]
    g1 = base_gadget()
    sg1 = sigma_set(g1)

    # ---- Step 1: d = 2 and d = 3 ingredients, wired as in eq. (108)/(111) ----
    g2 = set_product(g1, g1)                    # G_10
    g2A = set_product(sg1, g1)                  # G_10^A
    g2D = set_product(g1, sg1)                  # G_10^D
    Jp = {Tw(w) for w in g2.I}                  # J^+ (no exchanges at toy scale)

    g3X = set_product(g2, g1)
    g3het = het_set_product(g2, g1, g2.X, Jp, Jp)
    g3AX = set_product(g2A, g1)
    g3Ahet = het_set_product(g2A, g1, g2A.X, Jp, g2.X)
    g3DX = set_product(g2D, sg1)
    g3Dhet = het_set_product(g2D, sg1, g2D.X, g2.X, Jp)

    # ---- codebooks at the left dimension of Steps 2-4 ----
    J3 = {Tw(w) for w in g3X.I}                 # J_15 = T^{x3}(I_15)

    # ---- Steps 2-4, wired exactly as Table 6 ----
    g6_6 = het_set_product(g3X, g3het, g3het.X, J3, J3)            # G_30^(6)
    g5_8 = het_set_product(g3AX, g2D, g3Ahet.X, J3, J3)            # G_25^(8)
    g8_8 = het_set_product(g3DX, g5_8, g3Dhet.X, g3het.X, J3)      # G_40^(8)
    g6L = set_product(g3AX, sigma_set(g3Ahet))                     # G_30^L
    ghat6 = het_set_product(g3AX, sigma_set(g3Ahet),
                            g3Ahet.X, J3, g3AX.X)                  # Ghat_30
    ghat5 = het_set_product(g3AX, g2D, g3Ahet.X, J3, g3het.X)      # Ghat_25
    g6pp = set_product(g3het, g3het)                               # G_30^{++}

    _CACHE["P"] = dict(g1=g1, sg1=sg1, g2=g2, g2A=g2A, g2D=g2D, Jp=Jp, J3=J3,
                g3X=g3X, g3het=g3het, g3AX=g3AX, g3Ahet=g3Ahet,
                g3DX=g3DX, g3Dhet=g3Dhet, g6_6=g6_6, g5_8=g5_8, g8_8=g8_8,
                g6L=g6L, ghat6=ghat6, ghat5=ghat5, g6pp=g6pp)
    return _CACHE["P"]


def test_every_low_dimensional_gadget_satisfies_definition_2():
    P = build_pipeline()
    for name in ("g2", "g2A", "g2D", "g3X", "g3het", "g3AX", "g3Ahet",
                 "g3DX", "g3Dhet", "g5_8", "ghat5", "g6_6", "g6L",
                 "ghat6", "g6pp"):
        g = P[name]
        assert all(ok for _, ok in g.axioms()), (name, g.axioms())


def test_step1_profiles_match_the_formula_engine():
    P = build_pipeline()
    p1, sp1 = P["g1"].profile(), P["sg1"].profile()
    assert P["g2"].profile().as_tuple() == gao(p1, p1).as_tuple()
    assert P["g2A"].profile().as_tuple() == gao(sp1, p1).as_tuple()
    assert P["g2D"].profile().as_tuple() == gao(p1, sp1).as_tuple()
    assert P["g3X"].profile().as_tuple() == \
        gao(P["g2"].profile(), p1).as_tuple()
    for het_name, left, right, J0, JH, JV in [
            ("g3het", P["g2"], P["g1"], P["g2"].X, P["Jp"], P["Jp"]),
            ("g3Ahet", P["g2A"], P["g1"], P["g2A"].X, P["Jp"], P["g2"].X),
            ("g3Dhet", P["g2D"], P["sg1"], P["g2D"].X, P["g2"].X, P["Jp"])]:
        predicted = _measured_hetgao(left, right, J0, JH, JV, het_name)
        assert P[het_name].profile().as_tuple() == predicted.as_tuple(), \
            het_name


def test_orientation_identities():
    """G_15^{D,het} = sigma(G_15^{A,het}) and the shared physical sets."""
    P = build_pipeline()
    assert P["g3Dhet"].profile().as_tuple() == \
        sigma(P["g3Ahet"].profile()).as_tuple()
    assert P["g3Ahet"].X == P["g3Dhet"].X, \
        "X_15^{A,het} and X_15^{D,het} must be the same physical set"
    assert P["g3X"].X == P["g3AX"].X == P["g3DX"].X
    assert P["g3X"].X0 == P["g3AX"].X0 == P["g3DX"].X0, \
        "X_15^{0,X} must be common to the three ordinary gadgets"


def test_steps_2_to_4_profiles_match_the_formula_engine():
    P = build_pipeline()
    cases = [
        ("g6_6", P["g3X"], P["g3het"], P["g3het"].X, P["J3"], P["J3"]),
        ("g5_8", P["g3AX"], P["g2D"], P["g3Ahet"].X, P["J3"], P["J3"]),
        ("g8_8", P["g3DX"], P["g5_8"], P["g3Dhet"].X, P["g3het"].X, P["J3"]),
        ("ghat6", P["g3AX"], sigma_set(P["g3Ahet"]),
         P["g3Ahet"].X, P["J3"], P["g3AX"].X),
        ("ghat5", P["g3AX"], P["g2D"], P["g3Ahet"].X, P["J3"], P["g3het"].X),
    ]
    for name, left, right, J0, JH, JV in cases:
        predicted = _measured_hetgao(left, right, J0, JH, JV, name)
        assert P[name].profile().as_tuple() == predicted.as_tuple(), name
    assert P["g6L"].profile().as_tuple() == \
        gao(P["g3AX"].profile(), sigma(P["g3Ahet"].profile())).as_tuple()
    assert P["g6pp"].profile().as_tuple() == \
        gao(P["g3het"].profile(), P["g3het"].profile()).as_tuple()


def test_sibling_J0_admissibility_conditions():
    """The sibling-auxiliary transfers of Table 6 require shared transversals
    and the eq. 33 separation; check both as sets, per application."""
    P = build_pipeline()
    for J0, host in [(P["g3het"].X, P["g3X"]),      # G_30^(6)
                     (P["g3Ahet"].X, P["g3AX"]),    # G_25^(8), Ghat_30, Ghat_25
                     (P["g3Dhet"].X, P["g3DX"])]:   # G_40^(8)
        assert not (J0 & host.NH & host.NV), "eq. 33 fails for a J0 transfer"
    assert P["g3het"].PH == P["g3X"].PH and P["g3het"].PV == P["g3X"].PV
    assert P["g3Ahet"].PH == P["g3AX"].PH and P["g3Ahet"].PV == P["g3AX"].PV
    assert P["g3Dhet"].PH == P["g3DX"].PH and P["g3Dhet"].PV == P["g3DX"].PV


def test_final_sibling_J0_transfer():
    """The G55 step takes J0 = Xhat_30 (auxiliary of the heterogeneous sibling)
    inside the ordinary sibling G_30^L. Assert the transfer directly --
    identical main code, private pairs, and transversals, the
    eq. (33) separation, and the declared codebook decomposition."""
    P = build_pipeline()
    ghat, gL = P["ghat6"], P["g6L"]
    assert ghat.I == gL.I
    assert {tuple(p) for p in ghat.pairs} == {tuple(p) for p in gL.pairs}
    assert ghat.PH == gL.PH and ghat.PV == gL.PV
    assert not (ghat.X & gL.NH & gL.NV), "eq. 33 fails for Xhat_30 in G_30^L"
    j0 = len(ghat.X)
    o0 = len(ghat.X - (gL.NH | gL.NV))
    h0 = len(ghat.X & gL.NH)
    v0 = len(ghat.X & gL.NV)
    assert o0 + h0 + v0 == j0
    # measured against G_30^L equals Ghat_30's own declared decomposition,
    # because the two siblings share their transversals
    assert (o0, h0, v0) == (len(ghat.X0), len(ghat.XH), len(ghat.XV))
