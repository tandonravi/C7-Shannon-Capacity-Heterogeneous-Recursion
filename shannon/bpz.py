"""The BPZ multi-gadget framework: seven-family representations.

Rules and terminal codes are transcribed from
    ShannonBounds/Substitutions.lean, ShannonBounds/TerminalCodes.lean
at commit aa21eeb12b75b0413d3fa9fb4208b5d0bf2c4d65 of
    github.com/spectra-research/shannon-capacity-lean

Admissibility is re-checked here from Definition 6 rather than assumed; see
tests/test_bpz_rules.py.
"""

from functools import reduce
from itertools import combinations

# fixed coordinate order for every cardinality vector in this repo
LABELS = ("B", "N", "A", "D", "O", "H", "V")
INDEX = {c: i for i, c in enumerate(LABELS)}

# Definition 4: for each label, the labels it is required to be separated from.
PERP = {
    "B": {"O", "H", "V"},
    "N": {"A", "D", "O", "H", "V"},
    "A": {"N", "D", "H"},
    "D": {"N", "A", "V"},
    "O": {"B", "N"},
    "H": {"B", "N", "A"},
    "V": {"B", "N", "D"},
}


def perp(x, y):
    return y in PERP[x]


def _p(s):
    """Parse '(B,B),(H,D)' into [('B','B'), ('H','D')]."""
    s = s.replace(" ", "").strip()
    if not s:
        return []
    return [tuple(tok.strip("()").split(",")) for tok in s.split("),(")]


# --------------------------------------------------------------------------
# binary combining rules
# --------------------------------------------------------------------------

S2a = {
    "B": _p("(B,B),(H,D),(V,A),(D,V),(A,H)"),
    "N": _p("(N,N),(A,A),(A,D),(D,A),(D,D)"),
    "A": _p("(A,N),(N,D)"),
    "D": _p("(D,N),(N,A)"),
    "O": _p("(O,N),(N,O)"),
    "H": _p("(H,N),(N,V)"),
    "V": _p("(V,N),(N,H)"),
}

S2b = {
    "B": _p("(A,V),(B,B),(D,H)"),
    "N": _p("(A,A),(N,N)"),
    "A": _p("(A,N),(B,D),(D,V),(V,H)"),
    "D": _p("(D,N),(N,A)"),
    "O": [],
    "H": _p("(H,B),(N,V)"),
    "V": _p("(H,A),(N,H),(V,D),(V,N)"),
}

# --------------------------------------------------------------------------
# ternary combining rules
# --------------------------------------------------------------------------

S3a = {
    "B": _p("(A,H,N),(A,N,V),(B,A,H),(B,B,B),(B,D,V),(B,H,D),(B,V,A),"
            "(D,N,H),(D,V,N),(H,D,N),(H,N,A),(V,A,N),(V,N,D)"),
    "N": _p("(A,A,N),(A,D,N),(A,N,A),(A,N,D),(D,A,N),(D,D,N),(D,N,A),"
            "(D,N,D),(N,A,A),(N,A,D),(N,D,A),(N,D,D),(N,N,N)"),
    "A": _p("(A,A,A),(A,A,D),(A,D,A),(A,D,D),(A,N,N),(N,A,N),(N,N,A)"),
    "D": _p("(D,A,A),(D,A,D),(D,D,A),(D,D,D),(D,N,N),(N,D,N),(N,N,D)"),
    "O": _p("(O,A,A),(O,A,D),(O,D,A),(O,D,D)"),
    "H": _p("(H,A,A),(H,A,D),(H,D,A),(H,D,D),(H,N,N),(N,H,N),(N,N,H)"),
    "V": _p("(N,N,V),(N,V,N),(V,A,A),(V,A,D),(V,D,A),(V,D,D),(V,N,N)"),
}

S3b = {
    "B": _p("(B,B,B),(B,H,A),(N,A,V),(N,D,H),(N,V,D),(A,B,H),(A,H,N),"
            "(A,H,D),(D,B,V),(D,V,N),(D,V,D),(H,N,D),(H,D,N),(H,D,D),"
            "(V,N,A),(V,A,B),(V,D,A)"),
    "N": _p("(B,V,D),(N,N,N),(N,V,A),(A,N,H),(A,D,H),(A,H,N),(D,N,V),"
            "(D,D,H),(D,H,N),(H,V,D)"),
    "A": _p("(N,B,A),(N,A,N),(A,B,B),(A,A,V),(A,H,A)"),
    "D": _p("(N,B,D),(N,D,N),(D,B,B),(D,A,V),(D,H,A)"),
    "O": [],
    "H": _p("(N,N,H),(N,H,N),(H,N,N),(H,V,V)"),
    "V": _p("(N,N,V),(N,V,N),(V,N,N),(V,V,V)"),
}

S3c = {
    "B": _p("(A,B,V),(A,H,H),(A,V,A),(B,B,B),(D,A,H),(D,H,A),(D,N,H),"
            "(H,A,A),(H,A,N),(H,H,H),(N,H,D),(N,V,A),(V,A,D),(V,D,B)"),
    "N": _p("(A,A,D),(A,D,A),(B,H,D),(D,A,B),(D,B,H),(H,H,A),(H,V,H),"
            "(N,N,N),(V,A,A),(V,D,H)"),
    "A": _p("(A,B,N),(A,N,D),(B,N,A),(B,V,V),(H,D,N),(N,A,B),(N,D,A),"
            "(V,N,V),(V,V,N)"),
    "D": _p("(D,N,N),(N,D,N),(N,N,D)"),
    "O": [],
    "H": _p("(B,H,N),(H,N,B),(N,B,H)"),
    "V": _p("(A,V,N),(D,H,N),(H,N,A),(N,A,H),(N,D,V),(N,N,V),(N,V,N),"
            "(V,N,D),(V,N,N)"),
}

S3d = {
    "B": _p("(A,B,V),(B,B,B),(B,H,D),(B,V,A),(D,A,H),(D,D,V),(D,N,H),"
            "(H,A,A),(H,A,N),(V,A,D),(V,D,N),(V,H,A),(V,H,H)"),
    "N": _p("(A,A,D),(B,H,D),(D,A,B),(D,B,H),(H,V,H),(N,N,N),(V,A,A),"
            "(V,D,H),(V,H,A)"),
    "A": _p("(A,B,N),(A,N,D),(B,D,A),(B,N,A),(B,V,V),(H,D,N),(N,A,B),"
            "(V,N,V),(V,V,N)"),
    "D": _p("(D,N,N),(N,D,N),(N,N,D)"),
    "O": [],
    "H": _p("(B,H,N),(H,N,B),(N,B,H)"),
    "V": _p("(A,V,N),(D,H,N),(H,N,A),(N,A,H),(N,D,V),(N,N,V),(N,V,N),"
            "(V,N,D),(V,N,N)"),
}

S3e = {
    "B": _p("(A,H,N),(B,A,H),(B,B,B),(B,D,V),(B,H,D),(B,V,A),(D,V,N),"
            "(H,D,N),(H,N,A),(V,A,N),(V,N,D)"),
    "N": _p("(A,A,N),(A,H,N),(D,A,N),(D,D,N),(H,N,A),(N,A,A),(N,A,D),"
            "(N,D,A),(N,D,D),(N,N,N)"),
    "A": _p("(A,A,A),(A,A,D),(A,D,A),(A,D,D),(A,N,N),(B,N,A),(H,N,D),"
            "(N,A,N)"),
    "D": _p("(A,N,D),(D,A,A),(D,A,D),(D,D,A),(D,D,D),(D,N,N),(N,D,N),"
            "(N,N,D)"),
    "O": [],
    "H": _p("(B,N,H),(H,A,A),(H,A,D),(H,D,D),(H,H,A),(H,N,N),(N,H,N)"),
    "V": _p("(A,N,V),(D,N,H),(N,N,V),(N,V,N),(V,A,A),(V,A,D),(V,D,A),"
            "(V,D,D),(V,N,N)"),
}

S3f = {
    "B": _p("(A,H,D),(A,H,N),(B,A,H),(B,B,B),(B,D,V),(B,V,A),(H,D,B),"
            "(H,H,H),(H,N,A),(N,H,D),(V,N,D)"),
    "N": _p("(A,H,B),(B,H,H),(D,H,A),(H,D,D),(H,N,A),(N,A,A),(N,A,H),"
            "(N,H,A),(N,N,N),(V,H,V)"),
    "A": _p("(A,A,A),(A,A,H),(A,B,N),(B,N,A),(H,H,N),(H,N,H),(N,A,N)"),
    "D": _p("(A,N,D),(D,A,A),(D,A,D),(D,B,N),(N,D,N),(N,N,D)"),
    "O": [],
    "H": _p("(B,N,H),(H,A,B),(H,A,V),(H,N,N),(N,H,N)"),
    "V": _p("(A,N,V),(D,N,H),(D,V,N),(N,N,V),(N,V,N),(V,A,B),(V,A,V),"
            "(V,N,N)"),
}

S3g = {
    "B": _p("(A,H,N),(B,A,H),(B,B,B),(B,D,V),(B,H,D),(B,V,A),(H,D,N),"
            "(H,N,A),(V,H,A),(V,H,H),(V,N,D)"),
    "N": _p("(A,H,B),(B,H,H),(D,H,A),(H,D,D),(H,N,A),(N,A,A),(N,A,H),"
            "(N,H,A),(N,N,N),(V,H,V)"),
    "A": _p("(A,A,A),(A,A,H),(A,B,N),(B,N,A),(H,H,N),(H,N,H),(N,A,N)"),
    "D": _p("(A,N,D),(D,A,A),(D,A,H),(D,B,N),(N,D,N),(N,N,D)"),
    "O": [],
    "H": _p("(B,N,H),(H,A,B),(H,A,V),(H,N,N),(N,H,N)"),
    "V": _p("(A,N,V),(D,N,H),(D,V,N),(N,N,V),(N,V,N),(V,A,B),(V,A,V),"
            "(V,N,N)"),
}

S3h = {
    "B": S3b["B"],
    "N": S3b["N"],
    "A": _p("(N,B,D),(N,D,N),(D,B,B),(D,A,V),(D,H,A)"),
    "D": _p("(N,B,A),(N,A,N),(A,B,B),(A,A,V),(A,H,A)"),
    "O": [],
    "H": _p("(N,N,V),(N,V,N),(V,N,N),(V,V,V)"),
    "V": _p("(N,N,H),(N,H,N),(H,N,N),(H,V,V)"),
}

RULES = {"S2a": S2a, "S2b": S2b, "S3a": S3a, "S3b": S3b, "S3c": S3c,
         "S3d": S3d, "S3e": S3e, "S3f": S3f, "S3g": S3g, "S3h": S3h}

# --------------------------------------------------------------------------
# terminal codes
# --------------------------------------------------------------------------

K3a = _p("(A,B,V),(A,V,H),(A,H,N),(B,A,B),(B,H,A),(B,D,D),(B,N,B),(D,D,N),"
         "(D,A,H),(D,N,H),(H,H,H),(H,B,A),(H,V,N),(N,H,N),(N,B,V),(N,V,H),"
         "(V,B,N),(V,V,A),(V,B,D)")

K4a = _p("(B,B,B,B),(B,N,V,B),(B,A,H,B),(B,D,V,B),(B,H,N,B),(B,H,D,B),"
         "(B,V,A,B),(N,B,B,V),(N,N,V,H),(N,A,V,H),(N,D,H,H),(N,H,N,H),"
         "(N,H,A,H),(N,V,D,V),(A,B,B,V),(A,N,H,H),(A,A,V,H),(A,D,H,H),"
         "(A,H,A,H),(A,V,N,H),(A,V,D,H),(D,B,B,H),(D,N,V,V),(D,A,V,V),"
         "(D,D,H,V),(D,H,N,V),(D,H,A,V),(D,V,D,V),(H,A,H,D),(H,B,B,A),"
         "(H,N,V,D),(H,D,V,D),(H,A,V,N),(H,D,H,N),(H,H,A,N),(H,H,N,N),"
         "(H,H,N,D),(H,H,D,D),(H,V,D,N),(H,V,A,D),(V,N,V,N),(V,B,B,N),"
         "(V,B,B,D),(V,D,V,A),(V,N,V,A),(V,A,H,A),(V,V,A,A),(V,H,D,A),"
         "(V,V,N,A)")

K4b = _p("(A,A,A,V),(A,A,N,V),(A,B,H,A),(A,H,A,A),(A,H,A,N),(A,H,N,A),"
         "(A,H,N,N),(A,N,A,V),(A,N,D,H),(A,N,H,N),(A,N,N,H),(A,V,D,B),"
         "(B,A,D,H),(B,A,V,N),(B,B,B,B),(B,B,V,D),(B,D,B,H),(B,D,H,N),"
         "(B,H,H,H),(B,V,A,D),(B,V,N,D),(D,A,V,A),(D,H,H,A),(D,N,B,V),"
         "(D,N,V,A),(D,N,V,N),(D,V,B,A),(D,V,B,N),(H,B,V,H),(H,A,B,H),"
         "(H,B,D,B),(H,N,N,D),(H,D,A,B),(H,D,N,B),(N,A,A,V),(H,V,H,D),"
         "(H,V,H,N),(N,H,A,N),(N,A,N,V),(N,B,H,A),(N,N,D,H),(N,H,B,A),"
         "(N,N,A,V),(N,V,D,N),(N,N,N,V),(N,N,V,N),(V,A,A,B),(N,V,H,A),"
         "(N,V,N,N),(V,H,A,H),(V,A,N,B),(V,D,D,H),(V,N,N,A),(V,H,N,H),"
         "(V,N,A,B),(V,N,N,N),(V,V,H,V)")

CODES = {"K3a": K3a, "K4a": K4a, "K4b": K4b}


# --------------------------------------------------------------------------
# operations on cardinality vectors
# --------------------------------------------------------------------------

def sigma_vector(vec):
    """sigma on cardinality vectors, eq. 81: (B,N,A,D,O,H,V) -> (B,N,D,A,O,V,H)."""
    B, N, A, D, O, H, V = vec
    return (B, N, D, A, O, V, H)


def _word_value(word, vectors):
    return reduce(lambda acc, kv: acc * vectors[kv[0]][INDEX[kv[1]]],
                  enumerate(word), 1)


def apply_rule(rule, vectors):
    """Proposition 1, eq. 85."""
    arity = {len(w) for tuples in rule.values() for w in tuples}
    if arity and arity != {len(vectors)}:
        raise ValueError("rule arity does not match number of inputs")
    return tuple(sum(_word_value(w, vectors) for w in rule[lam])
                 for lam in LABELS)


def apply_code(code, vectors):
    """Cardinality of the independent set produced by a terminal code."""
    if {len(w) for w in code} != {len(vectors)}:
        raise ValueError("terminal code arity does not match number of inputs")
    return sum(_word_value(w, vectors) for w in code)


# --------------------------------------------------------------------------
# admissibility (Definition 6)
# --------------------------------------------------------------------------

def _label_separated(w1, w2):
    return any(perp(a, b) for a, b in zip(w1, w2))


def rule_problems(rule):
    """Return [] if the rule is admissible, else a list of violations.

    Also validates the schema first: every output label present, all word
    letters legal, and a single uniform arity across the whole rule.
    """
    bad = []
    if set(rule) != set(LABELS):
        bad.append(f"output labels {sorted(rule)} != {sorted(LABELS)}")
    arities = {len(w) for words in rule.values() for w in words}
    if len(arities) > 1:
        bad.append(f"mixed arities {sorted(arities)}")
    for lam, words in rule.items():
        for w in words:
            for c in w:
                if c not in INDEX:
                    bad.append(f"T_{lam}: illegal letter {c!r} in {w}")
    for lam, words in rule.items():
        if len(words) != len(set(words)):
            bad.append(f"T_{lam}: duplicate words")
        for w1, w2 in combinations(words, 2):
            if not _label_separated(w1, w2):
                bad.append(f"cond (i) T_{lam}: {w1} vs {w2}")
    for lam in LABELS:
        for mu in LABELS:
            if perp(lam, mu):
                for w1 in rule[lam]:
                    for w2 in rule[mu]:
                        if not _label_separated(w1, w2):
                            bad.append(f"cond (ii) {lam}|{mu}: {w1} vs {w2}")
    return bad


def code_problems(code):
    """Terminal code words must be pairwise separated in some coordinate."""
    return [f"{w1} vs {w2}" for w1, w2 in combinations(code, 2)
            if not _label_separated(w1, w2)]


def word_count(rule):
    return sum(len(v) for v in rule.values())
