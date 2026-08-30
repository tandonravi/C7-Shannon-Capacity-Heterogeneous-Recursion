"""Run the whole C4 pipeline at base dimension k = 1 and check it by brute force.

This is the test that matters for the certificates.  `shannon/certificates.py`
counts |J_30^{++} \\ N(X_30^{0,L})| over a set of 2.4e15 words by decomposing it
into products of five-dimensional atoms.  Nothing about that decomposition is
specific to k = 5, so the identical code is run here with the one-dimensional
gadget of Example 1.  The whole pipeline then lives in C_7^{box 6}, which has
117649 vertices, so every object can be built explicitly with the set-level
reference implementation and the answer obtained by direct enumeration.

An earlier version of this repository reported C4 = 852176977598432 instead of
841760069965664.  The cause was a neighbourhood predicate that had been paired
by hand with the wrong set: the atom was X_10^{A,H} (the 12236 auxiliary words
confusable with the A-oriented H-transversal) while the test evaluated
membership in N(P_10^{A,H}), the neighbourhood of the 5152-word transversal
itself.  The two agree on X_10 and differ off it, so the trivial-q self-test --
which only probes inside X_30^L -- could not see it.  This test can: it compares
the final number against enumeration, with no shared code path.

Atoms now derive their own neighbourhoods, so a predicate can no longer be
paired with the wrong set at all.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from shannon.certificates import (Base, Ring, apply_T_blocks, build,
                                  components, count_outside, reference, size)
from shannon.graph import closed_neighborhood, is_independent
from test_gadget_setlevel import (SetGadget, base_gadget, het_set_product,
                                  set_product)


def T1(c):
    """w -> 2 - w, an automorphism of C_7."""
    return (2 - c) % 7


def sigma(g):
    """Swap the two transversals.  Private pairs are NOT reversed: the centre
    of a private pair must stay in I."""
    return SetGadget(g.d, g.I, g.pairs, g.PV, g.PH, g.X)


def _explicit_pipeline():
    """The same construction as Section 4.2, built as honest vertex sets."""
    g1 = base_gadget()
    g2 = set_product(g1, g1)                                  # G_10
    g2A = set_product(sigma(g1), g1)                          # G_10^A
    Jp = {(T1(a), T1(b)) for a, b in g2.I}                    # J^+
    g3het = het_set_product(g2, g1, g2.X, Jp, Jp)             # G_15^het
    g3AX = set_product(g2A, g1)                               # G_15^{A,X}
    g3Ahet = het_set_product(g2A, g1, g2A.X, Jp, g2.X)        # G_15^{A,het}
    g6L = set_product(g3AX, sigma(g3Ahet))                    # G_30^L
    g6pp = set_product(g3het, g3het)                          # G_30^{++}
    return dict(g1=g1, g2=g2, g2A=g2A, Jp=Jp, g3het=g3het, g3AX=g3AX,
                g3Ahet=g3Ahet, g6L=g6L, g6pp=g6pp)


def _block_pipeline():
    ring = Ring(1)
    base = Base(ring, {0, 2, 4}, [(0, 6)], {0}, {6}, {1, 3, 5})
    G, _, _ = build(base, T1, (), verbose=False)
    return G


def _expand(blocks):
    """Materialise a block set as explicit words."""
    out = set()
    for blk in blocks:
        stack = [()]
        for atom in blk:
            items = (atom.elements if hasattr(atom, "elements") else atom.pairs)
            nxt = []
            for prefix in stack:
                for e in items:
                    nxt.append(prefix + (e if isinstance(e, tuple) else (e,)))
            stack = nxt
        out |= set(stack)
    return out


def test_every_toy_gadget_satisfies_definition_2():
    for name, g in _explicit_pipeline().items():
        if isinstance(g, SetGadget):
            assert all(ok for _, ok in g.axioms()), (name, g.axioms())


def test_J_plus_is_independent():
    assert is_independent(_explicit_pipeline()["Jp"])


def test_block_decompositions_equal_the_explicit_sets():
    """Not just the right size -- the right set, coordinate by coordinate."""
    E, G = _explicit_pipeline(), _block_pipeline()
    cases = [
        ("B15", G["B15"], E["g3het"].B), ("R15", G["R15"], E["g3het"].R),
        ("P15^H", G["PH15"], E["g3het"].PH), ("P15^V", G["PV15"], E["g3het"].PV),
        ("X15^0", G["X15_0"], E["g3het"].X0),
        ("X15^H", G["X15_H"], E["g3het"].XH),
        ("X15^V", G["X15_V"], E["g3het"].XV),
        ("X15^{0,X}", G["X15X_0"], E["g3AX"].X0),
        ("X15^{X,rest}", G["X15X_r"], E["g3AX"].XH | E["g3AX"].XV),
        ("X15^{A,het,0}", G["X15A_0"], E["g3Ahet"].X0),
        ("X15^{A,het,rest}", G["X15A_r"], E["g3Ahet"].XH | E["g3Ahet"].XV),
    ]
    for name, blocks, explicit in cases:
        assert _expand(blocks) == set(explicit), name
        assert size(blocks) == len(explicit), f"{name}: blocks not disjoint"


def test_reference_set_is_correct():
    """X_30^{0,L} must be exactly the neutral part of G_30^L."""
    E, G = _explicit_pipeline(), _block_pipeline()
    ref = reference(G)
    assert _expand(ref) == set(E["g6L"].X0)
    assert size(ref) == len(E["g6L"].X0)


def test_components_reconstruct_I30pp():
    E, G = _explicit_pipeline(), _block_pipeline()
    comps = components(G)
    assert sum(size(b) for _, b in comps) == len(E["g6pp"].I)
    union = set()
    for _, blocks in comps:
        union |= _expand(blocks)
    assert union == set(E["g6pp"].I)


def test_C4_matches_brute_force():
    """The headline check: the block count equals direct enumeration."""
    E, G = _explicit_pipeline(), _block_pipeline()
    J6 = {tuple(T1(x) for x in w) for w in E["g6pp"].I}
    assert len(J6) == len(E["g6pp"].I) and is_independent(J6)

    brute = len(J6 - closed_neighborhood(E["g6L"].X0, 6))
    machinery = sum(count_outside(apply_T_blocks(blocks, T1), reference(G))
                    for _, blocks in components(G))
    assert machinery == brute, f"machinery {machinery} != brute force {brute}"
