"""Gao gadget profiles and their propagation rules.

Two operations live here:

  gao(L, R)       Lemma 1  -- Gao's binary product.
  hetgao(L, R, ...)  Theorem 1 -- the heterogeneous refinement.

`hetgao` deliberately *requires* every codebook parameter to be passed by
keyword as an explicit Codebook object.  Filling those in is the whole content
of an application of Theorem 1, and leaving any of them implicit is how
verification errors hide.  Every call site therefore reads as a complete,
auditable statement of which sets were used.
"""

from dataclasses import dataclass
from typing import Optional


# --------------------------------------------------------------------------
# profiles
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Profile:
    """pi(G) = (a, t, s, o, h, v), the six-parameter Gao profile.

    a = |I|          size of the main independent set
    t                number of selected private pairs
    s = |X|          size of the auxiliary independent set
    o = |X^0|        auxiliary points confusable with neither transversal
    h = |X^H|        auxiliary points confusable with P^H only
    v = |X^V|        auxiliary points confusable with P^V only
    """
    a: int
    t: int
    s: int
    o: int
    h: int
    v: int
    name: str = ""

    def __post_init__(self):
        if any(x < 0 for x in (self.a, self.t, self.s, self.o, self.h, self.v)):
            raise ValueError(f"{self.name or 'profile'}: negative entry")
        if self.t > self.a:
            raise ValueError(f"{self.name or 'profile'}: t > a, but R is a subset of I")
        if self.o + self.h + self.v != self.s:
            raise ValueError(f"{self.name or 'profile'}: o + h + v != s")

    def as_tuple(self):
        return (self.a, self.t, self.s, self.o, self.h, self.v)

    def renamed(self, name):
        return Profile(*self.as_tuple(), name=name)

    def __str__(self):
        return f"{self.name or 'pi'} = {self.as_tuple()}"


def sigma(p: Profile, name: Optional[str] = None) -> Profile:
    """Reverse orientation: swap the two transversals P^H <-> P^V.

    The main code, private pairs and auxiliary set are unchanged; X^0 is
    unchanged because it depends only on the union of the two neighborhoods,
    while X^H and X^V are exchanged.  (Paper, eq. 79.)
    """
    a, t, s, o, h, v = p.as_tuple()
    return Profile(a, t, s, o, v, h, name=name or f"sigma({p.name})")


# --------------------------------------------------------------------------
# codebooks for Theorem 1
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class NeutralCodebook:
    """J_0: used over X_R^0.  Must be independent AND satisfy

        J_0 cap N(P_L^H) cap N(P_L^V) = empty            (paper, eq. 33)

    so it decomposes as J_0 = J_0^0 + J_0^H + J_0^V relative to the *left*
    gadget's transversals.  Sizes: (j0, o0, h0, v0).
    """
    j0: int
    o0: int
    h0: int
    v0: int
    source: str = ""

    def __post_init__(self):
        if min(self.j0, self.o0, self.h0, self.v0) < 0:
            raise ValueError(f"NeutralCodebook({self.source}): negative entry")
        if self.o0 + self.h0 + self.v0 != self.j0:
            raise ValueError(f"NeutralCodebook({self.source}): o0+h0+v0 != j0")


@dataclass(frozen=True)
class SideCodebook:
    """J_H or J_V: used over X_R^H / X_R^V.  Only needs to be independent.

    The single derived quantity is

        q = |J \\ N(X_L^0)|                                (paper, eq. 35)

    which splits J into the part contributing to the new neutral class and the
    part contributing to a one-sided class.
    """
    j: int
    q: int
    source: str = ""
    certificate: str = ""      # ledger key, or "" if derived structurally

    def __post_init__(self):
        if not 0 <= self.q <= self.j:
            raise ValueError(f"SideCodebook({self.source}): need 0 <= q <= j")


def own_auxiliary(p: Profile, source: str = "") -> SideCodebook:
    """The auxiliary set X of a gadget, used as a side codebook over itself.

    Trivial-q lemma.  X is independent and X^0 is a subset of X, so for x in X
    we have x ~ y for some y in X^0 only if x = y.  Hence X cap N(X^0) = X^0
    exactly, and

        q = |X \\ N(X^0)| = s - o

    with no computation required.  This is why several codebook choices in the
    paper need no entry in the certificate ledger.
    """
    return SideCodebook(p.s, p.s - p.o, source=source or f"X({p.name})",
                        certificate="")


def neutral_from(p: Profile, source: str = "") -> NeutralCodebook:
    """Use the auxiliary set of a gadget G' as J_0 for a product whose left
    input G shares G''s product code, private pairs and transversals.

    Legitimate exactly when G and G' arise from the same ordered input pair,
    because Theorem 1 leaves code, pairs and transversals unchanged; then
    X(G') satisfies eq. 33 for G's transversals and decomposes with G''s
    own (o, h, v).
    """
    return NeutralCodebook(p.s, p.o, p.h, p.v, source=source or f"X({p.name})")


# --------------------------------------------------------------------------
# products
# --------------------------------------------------------------------------

def gao(L: Profile, R: Profile, name: str = "") -> Profile:
    """Lemma 1 (Gao's product lemma), eq. 22."""
    a1, t1, s1, o1, h1, v1 = L.as_tuple()
    a2, t2, s2, o2, h2, v2 = R.as_tuple()
    return Profile(
        a=(a1 - t1) * (a2 - t2) + t1 * s2 + s1 * t2,
        t=t1 * o2 + o1 * t2,
        s=s1 * s2,
        o=o1 * o2 + (h1 + v1) * (h2 + v2),
        h=h1 * o2 + o1 * v2,
        v=v1 * o2 + o1 * h2,
        name=name or f"Gao({L.name},{R.name})",
    )


def hetgao(L: Profile, R: Profile, *,
           J0: NeutralCodebook,
           JH: SideCodebook,
           JV: SideCodebook,
           name: str = "",
           trace: Optional[list] = None) -> Profile:
    """Theorem 1 (heterogeneous refinement of Gao's product), eqs. 38-39.

    The main code size and private-pair count are those of gao(L, R); only the
    auxiliary parameters (s, o, h, v) change.
    """
    aL, tL, sL, oL, hL, vL = L.as_tuple()
    aR, tR, sR, oR, hR, vR = R.as_tuple()

    out = Profile(
        a=(aL - tL) * (aR - tR) + tL * sR + sL * tR,
        t=tL * oR + oL * tR,
        s=J0.j0 * oR + JH.j * hR + JV.j * vR,
        o=J0.o0 * oR + JH.q * hR + JV.q * vR,
        h=J0.h0 * oR + (JV.j - JV.q) * vR,
        v=J0.v0 * oR + (JH.j - JH.q) * hR,
        name=name or f"HetGao({L.name},{R.name})",
    )

    # Theorem 1 reduces to Lemma 1 when J0 = JH = JV = X_L.  Check it whenever
    # the caller happens to have made that choice: a free regression test.
    # (An explicit exception, not `assert`: python -O must not silence it.)
    if (J0 == neutral_from(L, J0.source) and JH.j == sL and JV.j == sL
            and JH.q == sL - oL and JV.q == sL - oL):
        if out.as_tuple() != gao(L, R).as_tuple():
            raise ArithmeticError("Theorem 1 did not reduce to Lemma 1 "
                                  f"for {out.name}")

    if trace is not None:
        trace.append(_TraceRow(out.name, L.name, R.name, J0, JH, JV, out))
    return out


@dataclass(frozen=True)
class _TraceRow:
    output: str
    left: str
    right: str
    J0: NeutralCodebook
    JH: SideCodebook
    JV: SideCodebook
    profile: Profile


# --------------------------------------------------------------------------
# Gao gadget -> BPZ seven-family cardinality vector
# --------------------------------------------------------------------------

def Phi(p: Profile):
    """Gao-to-BPZ cardinality map, eq. 77.

    Returned in the fixed order (B, N, A, D, O, H, V).  Note the crossed
    assignment: F_A = X^V and F_D = X^H.  It is forced -- the alternative
    would require X^H _|_ P^H, which is false.
    """
    a, t, s, o, h, v = p.as_tuple()
    return (a - t, o, v, h, t, t, t)


def format_codebook_table(trace):
    """Emit the per-application codebook table (Table 6 of the manuscript).

    Every application of Theorem 1 needs (j0,o0,h0,v0), (jH,qH), (jV,qV) stated
    explicitly; this generates the table straight from the code, and
    tests/test_paper_claims.py asserts it against the transcription of Table 6
    in data/paper_claims.json, so code and manuscript cannot drift apart.
    """
    lines = ["output                left            right           "
             "(j0, o0, h0, v0)                                   "
             "(jH, qH)                    (jV, qV)                    cert"]
    for r in trace:
        j0 = f"({r.J0.j0}, {r.J0.o0}, {r.J0.h0}, {r.J0.v0})"
        jh = f"({r.JH.j}, {r.JH.q})"
        jv = f"({r.JV.j}, {r.JV.q})"
        cert = ",".join(x for x in (r.JH.certificate, r.JV.certificate) if x)
        lines.append(f"{r.output:<21} {r.left:<15} {r.right:<15} "
                     f"{j0:<50} {jh:<27} {jv:<27} {cert or '-'}")
    return "\n".join(lines)
