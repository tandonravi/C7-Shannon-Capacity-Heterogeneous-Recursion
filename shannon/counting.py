"""Counting |J \\ N(X0)| for sets far too large to enumerate.

The certificates live in C_7^{box 10}, ^{box 15} and ^{box 30}, where the sets
have up to 2.4e15 elements.  Brute force is impossible.  But every set involved
is a *union of Cartesian products of atoms*, where an atom is a subset of
C_7^{box k} for a small k (k = 5 in the paper).  That structure is enough.

Two facts do all the work:

  1.  N(A_1 x ... x A_d) = N(A_1) x ... x N(A_d)
      The closed neighborhood of a product is the product of the neighborhoods,
      because confusability in the strong product is coordinatewise.

  2.  For a union, inclusion-exclusion:
          |S \\ N(union_j T_j)| = sum_{J} (-1)^{|J|} |S cap intersect_{j in J} N(T_j)|
      and when S and every T_j is a single product block, each term factors
      coordinatewise into atom-level intersections.

So a count in dimension 30 reduces to a handful of set operations on subsets of
C_7^{box 5}, which are tiny.  The number of inclusion-exclusion terms is
2^(number of blocks in the reference set), and the reference sets here have two
or three blocks, so this is cheap.

Everything is exercised in tests/test_counting.py at atom dimension 1, where
brute force is available for comparison; the algebra does not care whether an
atom lives in C_7^{box 1} or C_7^{box 5}.
"""

from itertools import combinations

from .graph import closed_neighborhood, vertices


class Atom:
    """A subset of C_7^{box k}, with its closed neighborhood cached."""

    __slots__ = ("elements", "k", "_nbhd", "name")

    def __init__(self, elements, k, name=""):
        self.elements = frozenset(elements)
        self.k = k
        self.name = name
        self._nbhd = None

    @property
    def neighborhood(self):
        if self._nbhd is None:
            self._nbhd = frozenset(closed_neighborhood(self.elements, self.k))
        return self._nbhd

    def __len__(self):
        return len(self.elements)

    def __repr__(self):
        return f"Atom({self.name or len(self.elements)}, k={self.k})"


class Block:
    """A Cartesian product Atom_1 x ... x Atom_d, of dimension sum(k_i)."""

    __slots__ = ("atoms",)

    def __init__(self, atoms):
        self.atoms = tuple(atoms)

    @property
    def dimension(self):
        return sum(a.k for a in self.atoms)

    def size(self):
        n = 1
        for a in self.atoms:
            n *= len(a)
        return n

    def provably_disjoint(self, other) -> bool:
        """True if some coordinate has disjoint atoms, which forces disjointness."""
        return any(not a.elements & b.elements
                   for a, b in zip(self.atoms, other.atoms))

    def __repr__(self):
        return " x ".join(a.name or str(len(a)) for a in self.atoms)


class BlockSet:
    """A union of Blocks, assumed pairwise disjoint.

    Every set appearing in the constructions is of this form, and the blocks are
    pairwise separated (hence disjoint) by the gadget axioms.  `size()` checks
    nothing; `assert_disjoint_shape` is provided for the caller who wants it.
    """

    def __init__(self, blocks, name=""):
        self.blocks = list(blocks)
        self.name = name
        dims = {b.dimension for b in self.blocks}
        if len(dims) > 1:
            raise ValueError(f"BlockSet({name}): blocks disagree on dimension")
        self.dimension = dims.pop() if dims else 0
        widths = {len(b.atoms) for b in self.blocks}
        if len(widths) > 1:
            raise ValueError(f"BlockSet({name}): blocks disagree on atom count")

    def size(self):
        return sum(b.size() for b in self.blocks)

    def __repr__(self):
        return (f"BlockSet({self.name}, {len(self.blocks)} blocks, "
                f"|.|={self.size()})")


def _block_minus_neighborhood(source: Block, reference: BlockSet) -> int:
    """|source \\ N(reference)| by inclusion-exclusion over reference blocks."""
    total = 0
    blocks = reference.blocks
    for r in range(len(blocks) + 1):
        for chosen in combinations(blocks, r):
            term = 1
            for i, atom in enumerate(source.atoms):
                live = atom.elements
                for ref_block in chosen:
                    live = live & ref_block.atoms[i].neighborhood
                    if not live:
                        break
                term *= len(live)
                if term == 0:
                    break
            total += (-1) ** r * term
    return total


def count_outside_neighborhood(source: BlockSet, reference: BlockSet) -> int:
    """|source \\ N(reference)|.

    Both arguments must be unions of Cartesian products over the *same* atom
    decomposition (same number of coordinates, matching atom dimensions).
    """
    # Total dimension and atom count alone do not determine the ordered
    # coordinate decomposition: (1, 2) and (2, 1) must not be zipped together.
    # Inspect EVERY block, including nonleading source/reference blocks.
    shapes = {tuple(atom.k for atom in block.atoms)
              for block in source.blocks + reference.blocks}
    if len(shapes) > 1:
        raise ValueError(f"atom decompositions differ: {sorted(shapes)}")
    if not source.blocks:
        return 0
    # The per-block counts are summed, so the source blocks must be disjoint.
    # (The reference blocks may overlap freely; inclusion-exclusion handles it.)
    for b1, b2 in combinations(source.blocks, 2):
        if not b1.provably_disjoint(b2):
            raise ValueError(
                f"source blocks '{b1}' and '{b2}' are not provably disjoint; "
                "counts would be double-counted")
    if not reference.blocks:
        return source.size()  # N(empty) is empty; disjointness was checked above.
    return sum(_block_minus_neighborhood(b, reference) for b in source.blocks)


def brute_force_outside_neighborhood(source_elements, reference_elements, d):
    """Reference implementation, for tests only.  Enumerates C_7^{box d}."""
    return len(set(source_elements) - closed_neighborhood(reference_elements, d))


def expand(bs: BlockSet):
    """Materialize a BlockSet as explicit vertices.  Tests and small d only."""
    out = set()
    for block in bs.blocks:
        parts = [list(a.elements) for a in block.atoms]
        stack = [()]
        for part in parts:
            stack = [prefix + piece for prefix in stack for piece in part]
        out.update(stack)
    return out


def atoms_from_partition(sets, k, names=None):
    """Convenience: turn a list of vertex sets into Atoms of dimension k."""
    names = names or [""] * len(sets)
    return [Atom(s, k, name=n) for s, n in zip(sets, names)]


def all_vertices_atom(k):
    return Atom(vertices(k), k, name="V")
