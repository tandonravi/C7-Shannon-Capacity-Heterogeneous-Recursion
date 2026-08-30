"""Derivation of the certificates C1-C4 and q15 from the five-dimensional base data.

Written to be validated end-to-end: the SAME code runs at base dimension k = 1,
where the whole d = 6k pipeline is small enough to enumerate by brute force
(tests/test_c4_bruteforce.py), and at k = 5, which is the real construction.
Nothing here is memoized by id(); every atom carries a stable uid, and every
neighbourhood predicate is derived from the atom itself, so a predicate cannot
be paired with the wrong set.

Representation.  Every set is a list of blocks; a block is a tuple of atoms
whose dimensions sum to the block dimension.  An atom is an explicit set of
words of dimension k (A1) or 2k (A2).  The only facts used are

    N(A x B) = N(A) x N(B)          and          N(union) = union of N,

so a word (w_1, ..., w_m) lies in N(Ref) iff some block j of Ref has
w_i in N(Ref[j][i]) for every i.  Counting the flag vectors that occur in each
coordinate separately and combining is therefore exact.
"""

import json
from collections import Counter
from itertools import product
from pathlib import Path


class Ring:
    """Confusability structure of C_7^{box k}."""

    def __init__(self, k):
        self.k = k
        self.n = 7 ** k
        deltas = list(product((-1, 0, 1), repeat=k))
        self.NB = []
        for c in range(self.n):
            w = [(c // 7 ** j) % 7 for j in range(k)]
            self.NB.append(frozenset(
                sum(((w[j] + d[j]) % 7) * 7 ** j for j in range(k)) for d in deltas))
        self.BITS = []
        for c in range(self.n):
            m = 0
            for x in self.NB[c]:
                m |= 1 << x
            self.BITS.append(m)

    def decode(self, c):
        return tuple((c // 7 ** j) % 7 for j in range(self.k))

    def nbhd(self, S):
        out = set()
        for c in S:
            out |= self.NB[c]
        return out


_UID = [0]


def _uid():
    _UID[0] += 1
    return _UID[0]


class A1:
    """An explicit subset of C_7^{box k}."""
    width = 1

    def __init__(self, ring, elements, name=""):
        self.ring, self.name, self.uid = ring, name, _uid()
        self.elements = frozenset(elements)
        self._nb = None

    def __len__(self):
        return len(self.elements)

    def in_nbhd(self, c):
        if self._nb is None:
            self._nb = self.ring.nbhd(self.elements)
        return c in self._nb

    def apply(self, T):
        return A1(self.ring, {T(c) for c in self.elements}, "T" + self.name)


class A2:
    """An explicit subset of C_7^{box 2k}, stored as pairs of k-dim codes."""
    width = 2

    def __init__(self, ring, pairs, name=""):
        self.ring, self.name, self.uid = ring, name, _uid()
        self.pairs = frozenset(pairs)
        self._by_first = None
        self._dilated = None

    def __len__(self):
        return len(self.pairs)

    def in_nbhd(self, w):
        a, b = w
        if self._dilated is None:
            by_first = {}
            for x, y in self.pairs:
                by_first[x] = by_first.get(x, 0) | (1 << y)
            dil = {}
            for c in range(self.ring.n):
                m = 0
                for c2 in self.ring.NB[c]:
                    v = by_first.get(c2)
                    if v:
                        m |= v
                if m:
                    dil[c] = m
            self._dilated = dil
        m = self._dilated.get(a)
        return m is not None and (self.ring.BITS[b] & m) != 0

    def apply(self, T):
        return A2(self.ring, {(T(x), T(y)) for x, y in self.pairs}, "T" + self.name)


def prod2(ring, A, B, name=""):
    """A1 x A1 -> A2."""
    return A2(ring, {(a, b) for a in A.elements for b in B.elements}, name)


# --------------------------------------------------------------------------
# block sets
# --------------------------------------------------------------------------

def size(blocks):
    total = 0
    for blk in blocks:
        n = 1
        for atom in blk:
            n *= len(atom)
        total += n
    return total


def cross(left, right):
    """Blocks of a Cartesian product of two block sets."""
    return [l + r for l in left for r in right]


def apply_T_blocks(blocks, T):
    """Apply T to every k-coordinate block; atoms are rebuilt, not cached."""
    cache = {}
    out = []
    for blk in blocks:
        new = []
        for atom in blk:
            if atom.uid not in cache:
                cache[atom.uid] = atom.apply(T)
            new.append(cache[atom.uid])
        out.append(tuple(new))
    return out


def count_outside(source, reference):
    """|union(source) \\ N(union(reference))|.

    Source blocks must be pairwise disjoint (checked by the caller against an
    explicit size where possible).  Reference blocks may overlap freely.
    """
    shapes = {tuple(a.width for a in b) for b in source} | \
             {tuple(a.width for a in b) for b in reference}
    if len(shapes) != 1:
        raise ValueError(f"inconsistent atom shapes: {shapes}")
    # enforce, not merely assume, pairwise disjointness of the source blocks:
    # two product blocks are disjoint iff some coordinate has disjoint atoms.
    dmemo = {}
    def _atoms_disjoint(x, y):
        if x.uid == y.uid:
            return False
        key = (x.uid, y.uid) if x.uid < y.uid else (y.uid, x.uid)
        if key not in dmemo:
            ex = x.elements if isinstance(x, A1) else x.pairs
            ey = y.elements if isinstance(y, A1) else y.pairs
            dmemo[key] = not (ex & ey) if len(ex) <= len(ey) else not (ey & ex)
        return dmemo[key]
    for i, b1 in enumerate(source):
        for b2 in source[i + 1:]:
            if not any(_atoms_disjoint(a1, a2) for a1, a2 in zip(b1, b2)):
                raise ArithmeticError(
                    "source blocks are not provably disjoint in any coordinate: "
                    f"{[a.name for a in b1]} vs {[a.name for a in b2]}")
    m = len(reference[0])
    memo = {}
    total = 0
    for blk in source:
        counters = []
        for i in range(m):
            atom = blk[i]
            ref_atoms = tuple(r[i] for r in reference)
            key = (atom.uid, tuple(r.uid for r in ref_atoms))
            if key not in memo:
                items = atom.elements if isinstance(atom, A1) else atom.pairs
                memo[key] = Counter(tuple(r.in_nbhd(e) for r in ref_atoms)
                                    for e in items)
            counters.append(memo[key])
        acc = Counter({(True,) * len(reference): 1})
        for c in counters:
            nxt = Counter()
            for k1, n1 in acc.items():
                for k2, n2 in c.items():
                    nxt[tuple(a and b for a, b in zip(k1, k2))] += n1 * n2
            acc = nxt
        for flags, n in acc.items():
            if not any(flags):
                total += n
    return total


# --------------------------------------------------------------------------
# the construction, expressed once, parameterised by the base gadget
# --------------------------------------------------------------------------

class Base:
    """A base gadget at dimension k, as explicit A1 atoms."""

    def __init__(self, ring, I, pairs, PH, PV, X):
        r = ring
        self.ring = r
        R = {a for a, _ in pairs}
        NH, NV = r.nbhd(PH), r.nbhd(PV)
        X0 = set(X) - (NH | NV)
        XH = set(X) & NH
        XV = set(X) & NV
        self.I = A1(r, I, "I1")
        self.B = A1(r, set(I) - R, "B1")
        self.R = A1(r, R, "R1")
        self.PH = A1(r, PH, "PH1")
        self.PV = A1(r, PV, "PV1")
        self.X = A1(r, X, "X1")
        self.X0 = A1(r, X0, "X0_1")
        self.XH = A1(r, XH, "XH1")
        self.XV = A1(r, XV, "XV1")
        self.Xr = A1(r, XH | XV, "Xr1")
        self.profile = (len(I), len(pairs), len(X), len(X0), len(XH), len(XV))


def build(base, T, exchanges=(), verbose=True):
    """Build every object needed for C4, as block sets of shape (2k, k, 2k, k)."""
    r = base.ring
    g = base

    def say(*a):
        if verbose:
            print(*a)

    # ---- dimension 2k: G_10 = Gao(G5, G5), plus the two oriented variants ----
    P = lambda A, B, nm="": prod2(r, A, B, nm)
    R2 = A2(r, P(g.R, g.X0).pairs | P(g.X0, g.R).pairs, "R2")
    PH2 = A2(r, P(g.PH, g.X0).pairs | P(g.X0, g.PV).pairs, "PH2")
    PV2 = A2(r, P(g.PV, g.X0).pairs | P(g.X0, g.PH).pairs, "PV2")
    X2 = P(g.X, g.X, "X2")
    X2_0 = A2(r, P(g.X0, g.X0).pairs | P(g.Xr, g.Xr).pairs, "X2_0")
    X2_H = A2(r, P(g.XH, g.X0).pairs | P(g.X0, g.XV).pairs, "X2_H")
    X2_V = A2(r, P(g.XV, g.X0).pairs | P(g.X0, g.XH).pairs, "X2_V")
    X2_r = A2(r, X2.pairs - X2_0.pairs, "X2_r")
    I2 = A2(r, (P(g.B, g.B).pairs | P(g.R, g.X0).pairs | P(g.PH, g.XH).pairs
                | P(g.PV, g.XV).pairs | P(g.X0, g.R).pairs | P(g.XH, g.PV).pairs
                | P(g.XV, g.PH).pairs), "I2")
    B2 = A2(r, I2.pairs - R2.pairs, "B2")

    # G_10^A = Gao(sigma G5, G5): transversals (PV,X0)u(X0,PV) and (PH,X0)u(X0,PH)
    PA_H = A2(r, P(g.PV, g.X0).pairs | P(g.X0, g.PV).pairs, "PA_H")
    PA_V = A2(r, P(g.PH, g.X0).pairs | P(g.X0, g.PH).pairs, "PA_V")
    X2A_0 = A2(r, {p for p in X2.pairs
                   if not PA_H.in_nbhd(p) and not PA_V.in_nbhd(p)}, "X2A_0")
    X2A_H = A2(r, {p for p in X2.pairs if PA_H.in_nbhd(p)}, "X2A_H")
    X2A_V = A2(r, {p for p in X2.pairs if PA_V.in_nbhd(p)}, "X2A_V")
    say(f"  d=2k: |I|={len(I2)} |X^0|={len(X2_0)} |X^H|={len(X2_H)} |X^V|={len(X2_V)}")
    say(f"        oriented A: |X^0|={len(X2A_0)} |X^H|={len(X2A_H)} |X^V|={len(X2A_V)}")
    if X2A_0.pairs != X2_0.pairs:
        raise ArithmeticError("X_10^0 must be orientation-independent")

    # ---- J^+ and its split against N(X_2^0) ----
    raw = {(T(a), T(b)) for a, b in I2.pairs}
    removed = {tuple(x) for x, _ in exchanges}
    inserted = {tuple(y) for _, y in exchanges}
    if exchanges:
        if not removed <= raw:
            raise ArithmeticError("a removed word is not in (TxT)(I_10)")
        if inserted & raw:
            raise ArithmeticError("an inserted word is already in (TxT)(I_10)")
        index = {}
        for a, b in raw:
            index.setdefault(a, set()).add(b)
        for (rem, ins) in [(tuple(x), tuple(y)) for x, y in exchanges]:
            nb = {(a2, b2) for a2 in r.NB[ins[0]] for b2 in index.get(a2, ())
                  if b2 in r.NB[ins[1]]}
            if nb != {rem}:
                raise ArithmeticError(
                    f"inserted {ins}: J-neighbours {nb} != {{{rem}}}")
            if not X2_0.in_nbhd(rem):
                raise ArithmeticError(
                    f"removed {rem} is not inside N(X_10^0); exchange gains nothing")
            if X2_0.in_nbhd(ins):
                raise ArithmeticError(
                    f"inserted {ins} lies inside N(X_10^0); exchange gains nothing")
        ins_l = sorted(inserted)
        for i, p in enumerate(ins_l):
            for q in ins_l[i + 1:]:
                if p[0] in r.NB[q[0]] and p[1] in r.NB[q[1]]:
                    raise ArithmeticError(
                        "inserted words are confusable with each other")
    Jp = A2(r, (raw - removed) | inserted, "Jp")
    Jq = A2(r, {p for p in Jp.pairs if not X2_0.in_nbhd(p)}, "Jq")
    JpN = A2(r, {p for p in Jp.pairs if X2_0.in_nbhd(p)}, "JpN")
    pre_exchange = sum(1 for p in raw if not X2_0.in_nbhd(p))
    if len(Jq) != pre_exchange + len(exchanges):
        raise ArithmeticError(
            f"|J^+ \\ N(X_10^0)| = {len(Jq)} != pre-exchange {pre_exchange} "
            f"+ {len(exchanges)} exchanges")
    say(f"        |J^+|={len(Jp)}  |Jq|={len(Jq)}  |J^+ cap N(X^0)|={len(JpN)}"
        f"  (pre-exchange {pre_exchange})")

    S = lambda *bl: list(bl)
    # ---- dimension 3k blocks, shape (2k, k) ----
    G = {}
    G["B15"] = S((B2, g.B), (PH2, g.XH), (PV2, g.XV), (X2_H, g.PV), (X2_V, g.PH))
    G["R15"] = S((R2, g.X0), (X2_0, g.R))
    G["PH15"] = S((PH2, g.X0), (X2_0, g.PV))
    G["PV15"] = S((PV2, g.X0), (X2_0, g.PH))
    G["X15_0"] = S((X2_0, g.X0), (Jq, g.XH), (Jq, g.XV))
    G["X15_H"] = S((X2_H, g.X0), (JpN, g.XV))
    G["X15_V"] = S((X2_V, g.X0), (JpN, g.XH))
    # ordinary 15-dim auxiliary set and its neutral part / rest
    G["_atoms"] = S((X2, g.X), (Jp, g.X))   # keeps X2 and Jp reachable by name
    G["X15X_0"] = S((X2_0, g.X0), (X2_r, g.Xr))
    G["X15X_r"] = S((X2_r, g.X0), (X2_0, g.Xr))
    # auxiliary of G_15^{A,het} = HetGao(G_10^A, G5)[(X_10, J^+, X_10); X_5]
    X2A_r = A2(r, X2A_H.pairs | X2A_V.pairs, "X2A_r")
    G["X15A_0"] = S((X2A_0, g.X0), (Jq, g.XH), (X2_r, g.XV))
    G["X15A_r"] = S((X2A_H, g.X0), (X2_0, g.XV), (X2A_V, g.X0), (JpN, g.XH))
    return G, T, {"C1_pre_exchange": pre_exchange}


def _nbhd_pairs(ring, atom2):
    """N(atom2) as an explicit set of pairs, for small instances only."""
    out = set()
    for a, b in atom2.pairs:
        for a2 in ring.NB[a]:
            for b2 in ring.NB[b]:
                out.add((a2, b2))
    return out


def components(G):
    """The seven Gao-product components of I_30^{++} = Gao(G_15^het, G_15^het)."""
    return [
        ("B15 x B15", cross(G["B15"], G["B15"])),
        ("R15 x X15^0", cross(G["R15"], G["X15_0"])),
        ("P15^H x X15^H", cross(G["PH15"], G["X15_H"])),
        ("P15^V x X15^V", cross(G["PV15"], G["X15_V"])),
        ("X15^0 x R15", cross(G["X15_0"], G["R15"])),
        ("X15^H x P15^V", cross(G["X15_H"], G["PV15"])),
        ("X15^V x P15^H", cross(G["X15_V"], G["PH15"])),
    ]


def reference(G):
    """X_30^{0,L}, the neutral part of G_30^L = Gao(G_15^{A,X}, sigma G_15^{A,het})."""
    return cross(G["X15X_0"], G["X15A_0"]) + cross(G["X15X_r"], G["X15A_r"])


def compute_C4(base, T, exchanges=(), verbose=True):
    G, _, _ = build(base, T, exchanges, verbose)
    comps = components(G)
    ref = reference(G)
    if verbose:
        print(f"  |I_30^++| = {sum(size(b) for _, b in comps)}")
        print(f"  |X_30^{{0,L}}| = {size(ref)}")
    parts = []
    for name, blocks in comps:
        tb = apply_T_blocks(blocks, T)
        v = count_outside(tb, ref)
        parts.append((name, size(blocks), v))
        if verbose:
            print(f"    {name:<16} size {size(blocks):>19}  contribution {v:>18}")
    return parts, G, ref


# --------------------------------------------------------------------------
# the five certificates
# --------------------------------------------------------------------------

DATA = Path(__file__).resolve().parent.parent / "data"


def load_base(path=None):
    """Rebuild the d = 5 gadget from data/base_c7.json and check Definition 2."""
    raw = json.load(open(path or DATA / "base_c7.json"))
    ring = Ring(5)
    I = set(raw["I0"]); X = set(raw["X"])
    pairs = [tuple(p) for p in raw["private_pairs"]]
    PH, PV = set(raw["PH"]), set(raw["PV"])
    NB = ring.NB
    axioms = {
        "I independent": all(not (NB[u] & (I - {u})) for u in I),
        "X independent": all(not (NB[u] & (X - {u})) for u in X),
        "private pairs valid": all(r in I and q not in I
                                   and {w for w in I if w in NB[q]} == {r}
                                   for r, q in pairs),
        "endpoint-disjoint": len({e for p in pairs for e in p}) == 2 * len(pairs),
        "P^H independent": all(not (NB[u] & (PH - {u})) for u in PH),
        "P^V independent": all(not (NB[u] & (PV - {u})) for u in PV),
        "complementary transversals": (
            len(PH) == len(PV) == len(pairs)
            and all((r in PH and q in PV) or (q in PH and r in PV)
                    for r, q in pairs)
            and PH | PV == {e for p in pairs for e in p}
            and not (PH & PV)),
        "X cap N(P^H) cap N(P^V) empty":
            not (X & ring.nbhd(PH) & ring.nbhd(PV)),
    }
    base = Base(ring, I, pairs, PH, PV, X)
    exchanges = [(tuple(a), tuple(b)) for a, b in raw["exchanges"]["pairs"]]

    def T(n):
        w = ring.decode(n)
        return sum(c * 7 ** j for j, c in
                   enumerate(((2 - w[1]) % 7, w[3], w[0], (2 - w[2]) % 7, w[4])))

    # T must be a bijective automorphism of C_7^{box 5}.  Checked exhaustively:
    # bijectivity on all 16807 codes, and T(N(c)) = N(T(c)) for every code c,
    # which is exactly preservation of the closed confusability relation.
    images = [T(c) for c in range(ring.n)]
    axioms["T is a bijection"] = len(set(images)) == ring.n
    axioms["T preserves confusability"] = all(
        {T(x) for x in ring.NB[c]} == ring.NB[images[c]] for c in range(ring.n))
    return base, T, exchanges, axioms


def compute_all(path=None, verbose=True):
    def say(*a):
        if verbose:
            print(*a)

    base, T, exchanges, axioms = load_base(path)
    say("Base gadget (d = 5)")
    for name, ok in axioms.items():
        say(f"  [{'ok' if ok else 'FAIL'}] {name}")
        if not ok:
            raise ArithmeticError(f"base gadget axiom failed: {name}")
    say(f"  profile {base.profile}")
    out = {"base_profile": base.profile, "base_axioms": axioms}

    G, _, extras = build(base, T, exchanges, verbose=verbose)
    out.update(extras)
    g = base
    A = {a.name: a for key in G for blk in G[key] for a in blk}

    # ---- d = 10 ----
    out["C1"] = len(A["Jq"])
    if out["C1"] != out["C1_pre_exchange"] + len(exchanges):
        raise ArithmeticError("C1 != pre-exchange count + number of exchanges")
    raw_J = A["Jp"]
    say(f"\n  C1 = |J^+ \\ N(X_10^0)| = {out['C1']}"
        f"  ({out['C1_pre_exchange']} before the exchanges)")

    # ---- d = 15, shape (10, 5) ----
    I15 = [(A["B2"], g.B), (A["R2"], g.X0), (A["PH2"], g.XH), (A["PV2"], g.XV),
           (A["X2_0"], g.R), (A["X2_H"], g.PV), (A["X2_V"], g.PH)]
    X15het = [(A["X2"], g.X0), (A["Jp"], g.XH), (A["Jp"], g.XV)]
    J15 = apply_T_blocks(I15, T)
    say(f"\n  |I_15| = {size(I15)}   |X_15| = {size(X15het)}")
    out["C3"] = count_outside(X15het, G["X15X_0"])
    out["C2"] = count_outside(J15, G["X15X_0"])
    out["q15"] = count_outside(J15, G["X15_0"])
    say(f"  C3  = |X_15 \\ N(X_15^{{0,X}})| = {out['C3']}")
    say(f"  C2  = |J_15 \\ N(X_15^{{0,X}})| = {out['C2']}")
    say(f"  q15 = |J_15 \\ N(X_15^0)|      = {out['q15']}")

    # ---- d = 30 ----
    comps = components(G)
    ref = reference(G)
    say(f"\n  |I_30^{{++}}| = {sum(size(b) for _, b in comps)}")
    say(f"  |X_30^{{0,L}}| = {size(ref)}")
    parts = {}
    for name, blocks in comps:
        v = count_outside(apply_T_blocks(blocks, T), ref)
        parts[name] = v
        say(f"    {name:<16} size {size(blocks):>19}  contribution {v:>18}")
    out["C4_components"] = parts
    out["C4"] = sum(parts.values())
    out["component_sizes"] = {n: size(b) for n, b in comps}
    out["I30pp"] = sum(size(b) for _, b in comps)
    out["X30_0L"] = size(ref)
    say(f"  C4 = {out['C4']}")
    return out
