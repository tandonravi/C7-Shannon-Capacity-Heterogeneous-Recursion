"""Assert that the code reproduces every numbered quantity in the manuscript.

data/paper_claims.json is a transcription of what the paper prints.  If any
assertion here fails, either the code or the manuscript is wrong, and the two
cannot silently drift apart across revisions.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shannon.construction import (baselines, load_certificates,
                                  nth_root_digits, theorem2, warmup)

ROOT = Path(__file__).resolve().parent.parent
CLAIMS = json.load(open(ROOT / "data" / "paper_claims.json"))


def _profiles():
    t2, wu = theorem2(), warmup()
    return {
        "eq49_G10": t2["G10"], "eq55_G15het": wu["G15h"],
        "eq60_G25het": wu["G25h"], "eq62_G40het": wu["G40h"],
        "eq64_G30het": wu["G30h"], "eq66_G60het": wu["G60h"],
        "eq109_G10A": t2["G10A"], "eq109_G10D": t2["G10D"],
        "eq112_G15X": t2["G15X"], "eq112_G15AX": t2["G15AX"],
        "eq112_G15Ahet": t2["G15Ah"], "eq112_G15DX": t2["G15DX"],
        "eq112_G15Dhet": t2["G15Dh"], "eq118_G30_6": t2["G30_6"],
        "eq122_G25_8": t2["G25_8"], "eq126_G40_8": t2["G40_8"],
        "eq129_G30L": t2["G30L"], "eq131_Ghat30": t2["Ghat30"],
        "eq134_Ghat25": t2["Ghat25"], "eq135_G25R": t2["G25R"],
        "eq141_G55": t2["G55"],
    }


def test_every_printed_profile():
    computed = _profiles()
    for key, expected in CLAIMS["profiles"].items():
        assert computed[key].as_tuple() == tuple(expected), key


def test_step7_partial_profile():
    """Step 7 yields only (a, t, s); no full profile exists to compare."""
    wu = warmup()
    exp = CLAIMS["step7_partial"]
    assert (wu["a100"], wu["t100"], wu["s100"]) \
        == (exp["a100"], exp["t100"], exp["s100"])


def test_cardinality_vectors():
    t2 = theorem2()
    base = baselines()
    exp = CLAIMS["cardinality_vectors"]
    assert base["bpz2_vectors"]["w"] == tuple(exp["eq101_w"])
    assert base["bpz2_vectors"]["sw"] == tuple(exp["eq102_sw"])
    assert t2["n6"] == tuple(exp["eq119_n6"])
    assert t2["n8"] == tuple(exp["eq127_n8"])
    assert t2["n11"] == tuple(exp["eq142_n11"])
    assert t2["n25"] == tuple(exp["eq163_n25"])


def test_n25_has_zero_O_coordinate():
    """S3b has an empty T_O, so the O-entry of n25 must vanish."""
    from shannon import bpz
    assert bpz.S3b["O"] == []
    assert theorem2()["n25"][bpz.INDEX["O"]] == 0


def test_code_sizes():
    exp = CLAIMS["code_sizes"]
    assert warmup()["a200"] == exp["eq72_a200_warmup"]
    t2 = theorem2()
    assert t2["M"] == exp["eq106_Mstar"]
    assert t2["Gpp30"].a == exp["eq137_I30pp"]


def test_bounds_to_every_printed_digit():
    exp = CLAIMS["bounds"]
    base, t2, wu = baselines(), theorem2(), warmup()
    cases = [
        (base["gao_a200"], 200, exp["eq26_gao"]),
        (base["bpz1_a200"], 200, exp["bpz_v1"]),
        (base["bpz2_M"], 500, exp["eq104_bpz_v2"]),
        (wu["a200"], 200, exp["eq73_warmup"]),
        (t2["M"], 500, exp["eq165_theorem2"]),
    ]
    for value, n, printed in cases:
        digits = len(printed.split(".")[1])
        assert nth_root_digits(value, n, digits) == printed, printed


def test_improvement_is_strictly_positive():
    """Theorem 2 must beat BPZ v2, and the warmup must sit between v1 and v2."""
    t2, wu, base = theorem2(), warmup(), baselines()
    d = 80
    thm2 = nth_root_digits(t2["M"], 500, d)
    bpz2 = nth_root_digits(base["bpz2_M"], 500, d)
    bpz1 = nth_root_digits(base["bpz1_a200"], 200, d)
    gao = nth_root_digits(base["gao_a200"], 200, d)
    warm = nth_root_digits(wu["a200"], 200, d)
    assert gao < bpz1 < warm < bpz2 < thm2, "ordering of published bounds"


def test_table5_components():
    """The seven Gao-product components of I_30^{++} and their total."""
    wu = warmup()
    a, t, s, o, h, v = wu["G15h"].as_tuple()
    computed = {
        "B15 x B15": (a - t) ** 2,
        "R15 x X15^0": t * o, "P15^H x X15^H": t * h, "P15^V x X15^V": t * v,
        "X15^0 x R15": o * t, "X15^H x P15^V": h * t, "X15^V x P15^H": v * t,
    }
    expected = {k: v for k, v in CLAIMS["table5_components"].items()
                if not k.startswith("_")}
    assert computed == expected
    assert sum(computed.values()) == CLAIMS["code_sizes"]["eq137_I30pp"]


def test_certificate_contributions_are_consistent():
    """Each Table 5 contribution to C4 must not exceed its component size."""
    led = load_certificates()
    c4 = led["certificates"]["C4"]
    parts = {k: v for k, v in c4["component_decomposition"].items()
             if not k.startswith("_")}
    assert sum(parts.values()) == c4["value"]
    sizes = {k: v for k, v in CLAIMS["table5_components"].items()
             if not k.startswith("_")}
    for name, contribution in parts.items():
        assert contribution <= sizes[name], name


def test_table6_trace_matches_the_manuscript():
    """Every row of Table 6: input pair, J0 decomposition, both q-pairs."""
    trace = theorem2()["trace"]
    table = {k: v for k, v in CLAIMS["table6"].items() if not k.startswith("_")}
    assert len(trace) == len(table) == 9
    for row in trace:
        left, right, j0, jh, jv = table[row.output]
        assert (row.left, row.right) == (left, right), row.output
        assert [row.J0.j0, row.J0.o0, row.J0.h0, row.J0.v0] == j0, row.output
        assert [row.JH.j, row.JH.q] == jh and [row.JV.j, row.JV.q] == jv, \
            row.output



def test_Mstar_appendix_copy_agrees():
    """eq. (164) repeats M* in Appendix B.4; both copies must equal the
    computed value."""
    led = load_certificates()
    t2 = theorem2(led)
    cs = CLAIMS["code_sizes"]
    assert cs["eq106_Mstar"] == cs["eq164_Mstar_appendix"]
    assert t2["M"] == cs["eq164_Mstar_appendix"]


def test_exchange_list_matches_base_data():
    """The eight exchanges as printed in Appendix B.1 (mirrored in
    paper_claims.json) must equal the machine input in base_c7.json."""
    base = json.load(open(ROOT / "data" / "base_c7.json"))
    assert CLAIMS["appendixB1_exchanges"]["pairs"] == base["exchanges"]["pairs"]
    assert len(base["exchanges"]["pairs"]) == 8


def test_certificates_within_range():
    led = load_certificates()
    for name, entry in led["certificates"].items():
        assert 0 <= entry["value"] <= entry["upper_bound"], name


def test_c3_forces_the_documented_cross_check():
    """|X_15 cap N(X_15^{0,X})| = 34038298 + 45 m, so C3 pins m = 29992."""
    led = load_certificates()
    C3 = led["certificates"]["C3"]["value"]
    s15 = warmup()["G15h"].s
    m, remainder = divmod(s15 - C3 - 105709 * 322, 45)
    assert remainder == 0 and m == 29992
