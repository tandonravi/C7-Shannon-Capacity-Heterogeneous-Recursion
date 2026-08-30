"""The constructions of the paper, built from the certificate ledger.

  warmup()        Section 3.1, the binary heterogeneous recursion -> d = 200
  theorem2()      Section 4.2, the heterogeneous BPZ recursion    -> d = 500
  baselines()     Gao, BPZ v1, BPZ v2 -- reproduced independently

Every function returns a dict of named intermediate objects so the tests can
assert against the paper equation by equation.
"""

import json
from pathlib import Path

from .gadget import (Profile, gao, hetgao, sigma, Phi,
                     NeutralCodebook, SideCodebook, own_auxiliary, neutral_from)
from . import bpz

DATA = Path(__file__).resolve().parent.parent / "data"


def load_certificates(path=None, recomputed=False):
    """Load the certificate ledger.

    With `recomputed=True` the certificate layer is DERIVED LIVE from
    data/base_c7.json by shannon.certificates.compute_all (about half a
    minute) and validated in full against the ledger: the five counts, the
    base profile, every base-gadget axiom, the pre-exchange C1 value, all
    seven C4 component contributions, and the |I_30^{++}| / |X_30^{0,L}|
    totals.  Any disagreement raises.
    """
    with open(path or DATA / "certificates.json") as f:
        led = json.load(f)
    if recomputed:
        from .certificates import compute_all
        derived = compute_all(verbose=False)
        for name in ("C1", "C2", "C3", "C4", "q15"):
            entry = led["certificates"][name]
            if derived[name] != entry["value"]:
                raise ArithmeticError(
                    f"{name}: derived {derived[name]} != ledger {entry['value']}")
            entry["value"] = derived[name]
            entry["live_derived"] = True
        if not all(derived["base_axioms"].values()):
            raise ArithmeticError("a base gadget axiom failed on live derivation")
        if tuple(derived["base_profile"]) != tuple(led["base_gadget"]["profile"]):
            raise ArithmeticError(
                f"base profile {derived['base_profile']} != ledger "
                f"{tuple(led['base_gadget']['profile'])}")
        c1 = led["certificates"]["C1"]
        if derived["C1_pre_exchange"] != c1["pre_exchange_value"]:
            raise ArithmeticError(
                f"pre-exchange count {derived['C1_pre_exchange']} != ledger "
                f"{c1['pre_exchange_value']}")
        c4 = led["certificates"]["C4"]
        recorded = {k: v for k, v in c4["component_decomposition"].items()
                    if not k.startswith("_")}
        if derived["C4_components"] != recorded:
            raise ArithmeticError("C4 component contributions disagree with ledger")
        if derived["I30pp"] != c4["size_I30pp"]:
            raise ArithmeticError(
                f"|I_30^++| {derived['I30pp']} != ledger {c4['size_I30pp']}")
        if derived["X30_0L"] != c4["size_X30_0L"]:
            raise ArithmeticError(
                f"|X_30^0L| {derived['X30_0L']} != ledger {c4['size_X30_0L']}")
    return led


def _cert(led, key):
    return led["certificates"][key]["value"]


def base_profiles(led):
    G5 = Profile(*led["base_gadget"]["profile"], name="G5")
    return G5, sigma(G5, name="sG5")


# --------------------------------------------------------------------------
# baselines
# --------------------------------------------------------------------------

def baselines(led=None):
    """Reproduce the three published bounds this paper builds on."""
    led = led or load_certificates()
    out = {}

    def binary_chain(p):
        g2 = gao(p, p)
        g3 = gao(p, g2)
        g5 = gao(g2, g3)
        g10 = gao(g5, g5)
        g20 = gao(g10, g10)
        return gao(g20, g20)

    gao_base = Profile(367, 8, 367, 321, 26, 20, name="Gao G5")
    out["gao_a200"] = binary_chain(gao_base).a
    out["bpz1_a200"] = binary_chain(base_profiles(led)[0]).a

    # BPZ v2 multi-gadget recursion, eq. 103
    G5, sG5 = base_profiles(led)
    w, sw = Phi(G5), Phi(sG5)
    S2a, S3a, S3b = bpz.S2a, bpz.S3a, bpz.S3b
    q2 = bpz.apply_rule(S2a, [w, sw])
    n3 = bpz.apply_rule(S3a, [w, w, w])
    n4 = bpz.apply_rule(S3a, [q2, w, w])
    n5 = bpz.apply_rule(S3a, [n3, w, w])
    n6 = bpz.apply_rule(S3a, [n3, w, q2])
    n8 = bpz.apply_rule(S3a, [n4, q2, q2])
    n11 = bpz.apply_rule(S3a, [n5, n3, n3])
    n25 = bpz.apply_rule(S3b, [n6, n11, n8])
    out["bpz2_M"] = bpz.apply_code(bpz.K4a, [n25] * 4)
    out["bpz2_vectors"] = dict(w=w, sw=sw, q2=q2, n3=n3, n4=n4, n5=n5,
                               n6=n6, n8=n8, n11=n11, n25=n25)
    return out


# --------------------------------------------------------------------------
# Section 3.1: the warmup binary recursion
# --------------------------------------------------------------------------

def warmup(led=None):
    """Section 3.1.  Uses C1 and q15.  Yields an independent set in C_7^{box 200}."""
    led = led or load_certificates()
    G5, _ = base_profiles(led)
    q10, q15 = _cert(led, "C1"), _cert(led, "q15")
    trace = []
    R = {}

    # Step 1
    G10 = gao(G5, G5, name="G10")

    # J^+ : automorphic image of I_10 plus eight exchanges (Appendix B.1)
    Jplus = SideCodebook(G10.a, q10, source="J^+", certificate="C1")

    # Step 2
    G15h = hetgao(G10, G5, J0=neutral_from(G10, "X_10"), JH=Jplus, JV=Jplus,
                  name="G15^het", trace=trace)

    # J_15 : pure automorphic image of I_15, measured against X_15^0
    J15 = SideCodebook(G15h.a, q15, source="J_15", certificate="q15")
    X15 = neutral_from(G15h, "X_15")

    # Steps 3-6
    G25h = hetgao(G15h, G10, J0=X15, JH=J15, JV=J15, name="G25^het", trace=trace)
    G40h = hetgao(G15h, G25h, J0=X15, JH=J15, JV=J15, name="G40^het", trace=trace)
    G30h = hetgao(G15h, G15h, J0=X15, JH=J15, JV=J15, name="G30^het", trace=trace)
    G60h = gao(G30h, G30h, name="G60^het")

    # Step 7: HetGao(G60^het, G40^het)[(X_60, I_60, I_60); X_40].  Only
    # (a, t, s) of the result are needed by the final Gao product, and none of
    # them depends on an overlap count, so no certificate is required -- and no
    # full Profile is fabricated here: the o/h/v split of the 100-dimensional
    # auxiliary set is simply not determined by the available data.
    a100 = (G60h.a - G60h.t) * (G40h.a - G40h.t) + G60h.t * G40h.s + G60h.s * G40h.t
    t100 = G60h.t * G40h.o + G60h.o * G40h.t
    s100 = G60h.s * G40h.o + G60h.a * (G40h.h + G40h.v)   # j0*o_R + j_H*h_R + j_V*v_R

    # Step 8: a200 for Gao(G100^het, G100^het) uses only (a, t, s) of Step 7.
    a200 = (a100 - t100) ** 2 + 2 * t100 * s100

    R.update(G10=G10, G15h=G15h, G25h=G25h, G40h=G40h, G30h=G30h, G60h=G60h,
             a100=a100, t100=t100, s100=s100, trace=trace, a200=a200)
    return R


# --------------------------------------------------------------------------
# Section 4.2: Theorem 2
# --------------------------------------------------------------------------

def theorem2(led=None):
    """Theorem 2.  Uses C1-C4.  Yields an independent set in C_7^{box 500}."""
    led = led or load_certificates()
    G5, sG5 = base_profiles(led)
    C1, C2, C3, C4 = (_cert(led, k) for k in ("C1", "C2", "C3", "C4"))
    trace = []

    # ---- Step 1: ten- and fifteen-dimensional ingredients -----------------
    G10 = gao(G5, G5, name="G10")
    G10A = gao(sG5, G5, name="G10^A")
    G10D = gao(G5, sG5, name="G10^D")
    # All three share the physical auxiliary set X_10 = X_5 x X_5 and, because
    # X^0 depends only on the union N(P^H) u N(P^V), the same neutral part.

    Jplus = SideCodebook(G10.a, C1, source="J^+", certificate="C1")
    X10 = own_auxiliary(G10, "X_10")            # trivial q = s - o = 28980

    G15X = gao(G10, G5, name="G15^X")
    G15h = hetgao(G10, G5, J0=neutral_from(G10, "X_10"),
                  JH=Jplus, JV=Jplus, name="G15^het", trace=trace)
    G15AX = gao(G10A, G5, name="G15^{A,X}")
    G15Ah = hetgao(G10A, G5, J0=neutral_from(G10A, "X_10"),
                   JH=Jplus, JV=X10, name="G15^{A,het}", trace=trace)
    G15DX = gao(G10D, sG5, name="G15^{D,X}")
    G15Dh = hetgao(G10D, sG5, J0=neutral_from(G10D, "X_10"),
                   JH=X10, JV=Jplus, name="G15^{D,het}", trace=trace)

    if sigma(G15Ah).as_tuple() != G15Dh.as_tuple():
        raise ArithmeticError("G15^{D,het} != sigma(G15^{A,het})")

    # codebooks living in C_7^{box 15}, all measured against X_15^{0,X}
    J15 = SideCodebook(G15X.a, C2, source="J_15", certificate="C2")
    X15 = SideCodebook(G15h.s, C3, source="X_15", certificate="C3")
    X15X = own_auxiliary(G15X, "X_15^X")        # trivial q = 14088465

    # ---- Step 2: n6 -------------------------------------------------------
    G30_6 = hetgao(G15X, G15h, J0=neutral_from(G15h, "X_15"),
                   JH=J15, JV=J15, name="G30^(6)", trace=trace)
    n6 = Phi(G30_6)

    # ---- Step 3: n8 -------------------------------------------------------
    G25_8 = hetgao(G15AX, G10D, J0=neutral_from(G15Ah, "X_15^{A,het}"),
                   JH=J15, JV=J15, name="G25^(8)", trace=trace)
    G40_8 = hetgao(G15DX, G25_8, J0=neutral_from(G15Dh, "X_15^{D,het}"),
                   JH=X15, JV=J15, name="G40^(8)", trace=trace)
    n8 = Phi(G40_8)

    # ---- Step 4: n11 ------------------------------------------------------
    G30L = gao(G15AX, sigma(G15Ah), name="G30^L")
    Ghat30 = hetgao(G15AX, sigma(G15Ah), J0=neutral_from(G15Ah, "X_15^{A,het}"),
                    JH=J15, JV=X15X, name="Ghat30", trace=trace)
    Ghat25 = hetgao(G15AX, G10D, J0=neutral_from(G15Ah, "X_15^{A,het}"),
                    JH=J15, JV=X15, name="Ghat25", trace=trace)
    G25R = sigma(Ghat25, name="G25^R")

    Gpp30 = gao(G15h, G15h, name="G30^{++}")
    J30pp = SideCodebook(Gpp30.a, C4, source="J_30^{++}", certificate="C4")

    G55 = hetgao(G30L, G25R, J0=neutral_from(Ghat30, "Xhat_30"),
                 JH=J30pp, JV=J30pp, name="G55", trace=trace)
    n11 = Phi(sigma(G55))

    # ---- Step 5: unchanged BPZ top level ----------------------------------
    n25 = bpz.apply_rule(bpz.S3b, [n6, n11, n8])
    M = bpz.apply_code(bpz.K4a, [n25] * 4)

    return dict(G10=G10, G10A=G10A, G10D=G10D, G15X=G15X, G15h=G15h,
                G15AX=G15AX, G15Ah=G15Ah, G15DX=G15DX, G15Dh=G15Dh,
                G30_6=G30_6, G25_8=G25_8, G40_8=G40_8, G30L=G30L,
                Ghat30=Ghat30, Ghat25=Ghat25, G25R=G25R, Gpp30=Gpp30,
                G55=G55, n6=n6, n8=n8, n11=n11, n25=n25, M=M, trace=trace)


# --------------------------------------------------------------------------
# exact roots
# --------------------------------------------------------------------------

def nth_root_digits(M, n, digits=60):
    """Truncate the nonnegative nth root to ``digits`` decimal places.

    M, n, and digits must be integers, with M >= 0, n >= 1, digits >= 0.
    Integer bisection certifies the returned lower endpoint; no floating
    point or rounding-to-nearest is used, including at zero precision.
    """
    from operator import index

    try:
        M, n, digits = index(M), index(n), index(digits)
    except TypeError as exc:
        raise TypeError("M, n, and digits must be integers") from exc
    if M < 0:
        raise ValueError("M must be nonnegative")
    if n <= 0:
        raise ValueError("root degree n must be positive")
    if digits < 0:
        raise ValueError("digits must be nonnegative")

    scale = 10 ** digits
    target = M * scale ** n
    # This is an upper bound for the root, also when target == 0.
    lo, hi = 0, 1 << ((target.bit_length() + n - 1) // n)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** n <= target:
            lo = mid
        else:
            hi = mid - 1
    if digits == 0:
        return str(lo)
    whole, fractional = divmod(lo, scale)
    return f"{whole}.{fractional:0{digits}d}"
