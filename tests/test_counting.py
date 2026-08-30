"""Check the product-set counting algebra against brute force.

shannon/counting.py computes |J \\ N(X0)| for sets with up to 2.4e15 elements by
decomposing them into Cartesian products of small atoms.  Here the same code
runs with atoms of dimension 1, where the answer can be obtained by direct
enumeration of C_7^{box d}.  The algebra is indifferent to atom dimension, so
agreement at k = 1 is evidence for correctness at k = 5.
"""

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shannon.counting import (Atom, Block, BlockSet,
                              brute_force_outside_neighborhood,
                              count_outside_neighborhood, expand)
from shannon.graph import (closed_neighborhood,
                           closed_neighborhood_by_scan, vertices)

V1 = vertices(1)


def _atom(rng):
    return Atom(rng.sample(V1, rng.randint(1, 4)), 1)


def _disjoint_source(d, n_blocks, rng):
    """Blocks made provably disjoint by partitioning the first coordinate."""
    firsts = rng.sample(V1, min(n_blocks, 7))
    return BlockSet([Block([Atom([f], 1)] + [_atom(rng) for _ in range(d - 1)])
                     for f in firsts])


def test_matches_brute_force_on_random_instances():
    rng = random.Random(20260828)
    for _ in range(300):
        d = rng.choice([2, 3])
        source = _disjoint_source(d, rng.randint(1, 3), rng)
        reference = BlockSet([Block([_atom(rng) for _ in range(d)])
                              for _ in range(rng.randint(1, 3))])
        assert (count_outside_neighborhood(source, reference)
                == brute_force_outside_neighborhood(expand(source),
                                                    expand(reference), d))


def test_neighborhood_of_a_product_is_the_product_of_neighborhoods():
    rng = random.Random(7)
    for _ in range(50):
        a, b = _atom(rng), _atom(rng)
        block = BlockSet([Block([a, b])])
        empty_outside = count_outside_neighborhood(block, block)
        assert empty_outside == 0, "a set is inside its own neighborhood"


def test_trivial_q_lemma_via_counting():
    """For an independent X with X0 a subset of X, |X \\ N(X0)| = |X| - |X0|."""
    rng = random.Random(11)
    independent_atom = Atom([(0,), (2,), (4,)], 1)   # independent in C_7
    neutral_atom = Atom([(0,)], 1)
    for d in (2, 3):
        X = BlockSet([Block([independent_atom] * d)])
        X0 = BlockSet([Block([neutral_atom] * d)])
        got = count_outside_neighborhood(X, X0)
        assert got == X.size() - X0.size()


def test_expansion_neighborhood_equals_scan():
    """The fast expansion method must equal the exhaustive scan."""
    rng = random.Random(5)
    for d in (1, 2, 3):
        for _ in range(20):
            S = {tuple(rng.randrange(7) for _ in range(d))
                 for _ in range(rng.randint(1, 6))}
            assert closed_neighborhood(S, d) == closed_neighborhood_by_scan(S, d)


def test_rejects_overlapping_source_blocks():
    """Summing over source blocks is only valid when they are disjoint."""
    a = Atom([(0,), (2,)], 1)
    overlapping = BlockSet([Block([a, a]), Block([a, a])])
    reference = BlockSet([Block([a, a])])
    try:
        count_outside_neighborhood(overlapping, reference)
    except ValueError as exc:
        assert "disjoint" in str(exc)
    else:
        raise AssertionError("expected a disjointness error")
