"""Confusability primitives for strong powers of the odd cycle C_n (n = 7).

Everything here is brute force over Z_n^d and is intended for small d only
(d <= 4 is instant, d = 5 is fine, beyond that use shannon.counting).
It exists so that the propagation formulas in shannon.gadget can be checked
against honest set-level computation in the tests.
"""

from itertools import product

N_CYCLE = 7


def confusable_symbol(a, b, n=N_CYCLE):
    """a ~ b in C_n, where '~' includes equality (closed relation)."""
    return (a - b) % n in (0, 1, n - 1)


def confusable(x, y, n=N_CYCLE):
    """x ~ y in C_n^{box d}: confusable in *every* coordinate."""
    if len(x) != len(y):
        raise ValueError(f"words of unequal length: {len(x)} vs {len(y)}")
    return all(confusable_symbol(a, b, n) for a, b in zip(x, y))


def vertices(d, n=N_CYCLE):
    return list(product(range(n), repeat=d))


def is_independent(S, n=N_CYCLE):
    S = list(S)
    for i, x in enumerate(S):
        for y in S[i + 1:]:
            if confusable(x, y, n):
                return False
    return True


def closed_neighborhood(S, d, n=N_CYCLE):
    """N(S) = {u : u ~ s for some s in S}, as a set of vertices of C_n^{box d}.

    Computed by expanding each word of S to its 3^d confusable words.  Exact,
    and linear in |S| * 3^d rather than n^d * |S|; the two methods are compared
    against each other in tests/test_counting.py.
    """
    out = set()
    for s in S:
        if len(s) != d:
            raise ValueError(f"word of length {len(s)} in dimension-{d} call")
        stack = [()]
        for c in s:
            stack = [p + ((c + delta) % n,) for p in stack for delta in (-1, 0, 1)]
        out.update(stack)
    return out


def closed_neighborhood_by_scan(S, d, n=N_CYCLE):
    """Reference implementation by full scan of the ambient graph.  Tests only."""
    S = list(S)
    out = set()
    for u in vertices(d, n):
        for s in S:
            if confusable(u, s, n):
                out.add(u)
                break
    return out


def separated(A, B, n=N_CYCLE):
    """A _|_ B: no vertex of A is confusable with a vertex of B."""
    return not any(confusable(a, b, n) for a in A for b in B)


def independent_of_neighborhood_count(J, X0, d, n=N_CYCLE):
    """|J \\ N(X0)| -- the shape every certificate C1..C4 and q15 takes."""
    return len(set(J) - closed_neighborhood(X0, d, n))
