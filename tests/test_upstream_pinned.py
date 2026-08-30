"""Parse the pinned upstream Lean files and assert the transcriptions match.

This converts the transcription trust root into a machine check.  The
combining rules and terminal codes in shannon/bpz.py, and the base data in
data/base_c7.json, were transcribed from the BPZ Lean repository at commit
aa21eeb12b75b0413d3fa9fb4208b5d0bf2c4d65; verbatim snapshots of the three
source files live in upstream/ with their SHA-256 pins.  Here they are
re-parsed with an independent parser and compared set-for-set and
word-for-word.  What remains on trust is only that the snapshots themselves
equal the upstream commit -- which anyone can check against SHA256SUMS.
"""

import hashlib
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shannon import bpz

ROOT = Path(__file__).resolve().parent.parent
UP = ROOT / "upstream"


def _words(block):
    return [tuple(x.strip().lstrip(".") for x in w.split(","))
            for w in re.findall(r"!\[([^\]]*)\]", block)]


def _parse_rule(text, table_name):
    m = re.search(rf"def {table_name} : Letter . Finset \(Fin \d+ . Letter\)"
                  r"(.*?)(?=\n/--|\ndef |\Z)", text, re.S)
    out = {}
    for lam, block in re.findall(
            r"\|\s*\.?([BNADOHV])\s*=>\s*(.*?)(?=\n\s*\|\s*\.?[BNADOHV]\s*=>|\Z)",
            m.group(1), re.S):
        out[lam] = [] if "∅" in block else _words(block)
    return out


def _parse_code(text, name):
    m = re.search(rf"def {name} : Finset \(Fin \d+ . Letter\) :="
                  r"(.*?)(?=\n/--|\ndef |\Z)", text, re.S)
    return _words(m.group(1))


def _ints(text, name):
    m = re.search(rf"def {name} : List Code :=\s*\[(.*?)\]", text, re.S)
    return [int(x) for x in re.findall(r"\d+", m.group(1))]


def test_snapshot_hashes_match_the_pin_file():
    pins = dict(line.split()[::-1] for line in
                (UP / "SHA256SUMS").read_text().strip().split("\n"))
    for rel, expected in pins.items():
        actual = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
        assert actual == expected, rel


def test_all_rules_match_upstream_word_for_word():
    text = (UP / "Substitutions.lean").read_text()
    for name in ("S2a", "S2b", "S3a", "S3b", "S3c", "S3d",
                 "S3e", "S3f", "S3g", "S3h"):
        upstream = _parse_rule(text, "T" + name[1:])
        ours = bpz.RULES[name]
        for lam in bpz.LABELS:
            assert set(map(tuple, upstream.get(lam, []))) == set(ours[lam]), \
                (name, lam)


def test_all_terminal_codes_match_upstream():
    text = (UP / "TerminalCodes.lean").read_text()
    for ours_name, lean_name in (("K3a", "C3a"), ("K4a", "C4a"), ("K4b", "C4b")):
        assert set(map(tuple, _parse_code(text, lean_name))) \
            == set(bpz.CODES[ours_name]), ours_name


def test_base_data_matches_upstream():
    text = (UP / "BaseC7Data.lean").read_text()
    data = json.load(open(ROOT / "data" / "base_c7.json"))
    assert sorted(_ints(text, "Ilist")) == sorted(data["I0"])
    assert sorted(_ints(text, "Xlist")) == sorted(data["X"])

    lean_pairs = [(int(a), int(b), side == "true") for a, b, side in
                  re.findall(r"\((\d+),\s*(\d+),\s*(true|false)\)", text)]
    assert [[r, q] for r, q, _ in lean_pairs] == data["private_pairs"]

    # transversal naming: P^H is the `true` side of BPZ's `bep`, i.e. it takes
    # the centre r when side is true and the alternative q when side is false.
    PH = {r if side else q for r, q, side in lean_pairs}
    PV = {q if side else r for r, q, side in lean_pairs}
    assert PH == set(data["PH"]) and PV == set(data["PV"])
    # ...which puts exactly the pairs {0, 5, 6} centre-first in P^H (eq. 23).
    assert sorted(i for i, (r, _, s) in enumerate(lean_pairs) if s) == [0, 5, 6]
