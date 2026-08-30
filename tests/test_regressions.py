"""Regression tests for root formatting and product-counting edge cases."""
import subprocess
import sys
from pathlib import Path

import pytest

from shannon.construction import nth_root_digits
from shannon.counting import Atom, Block, BlockSet, count_outside_neighborhood, expand
from shannon.graph import closed_neighborhood

ROOT = Path(__file__).resolve().parent.parent


def test_zero_root():
    assert nth_root_digits(0, 2, 3) == "0.000"
    assert nth_root_digits(0, 500, 0) == "0"
    assert nth_root_digits(0, 500, 68) == "0." + "0" * 68


def test_zero_precision():
    assert nth_root_digits(10, 2, 0) == "3"
    assert nth_root_digits(1, 500, 0) == "1"
    assert nth_root_digits(16, 2, 0) == "4"


def test_small_root_grid_by_integer_inequalities():
    # 12,642 independently checked intervals, including exact powers and zero.
    for M in range(301):
        for degree in range(1, 8):
            for digits in range(6):
                s = nth_root_digits(M, degree, digits)
                assert ("." in s) == (digits > 0)
                if digits:
                    assert len(s.split(".")[1]) == digits
                k = int(s.replace(".", ""))
                target = M * 10 ** (digits * degree)
                assert k ** degree <= target < (k + 1) ** degree


@pytest.mark.parametrize("args", [(-1, 2, 3), (10, 0, 3), (10, -1, 3), (10, 2, -1)])
def test_invalid_root_ranges_rejected(args):
    with pytest.raises(ValueError):
        nth_root_digits(*args)


@pytest.mark.parametrize("args", [(10.0, 2, 3), (10, 2.0, 3), (10, 2, 3.0), ("10", 2, 3)])
def test_noninteger_root_arguments_rejected(args):
    with pytest.raises(TypeError):
        nth_root_digits(*args)


def test_cli_negative_precision_rejected_before_recomputation():
    run = subprocess.run([sys.executable, str(ROOT / "verify_arithmetic.py"),
                          "--digits", "-1", "--recomputed"],
                         capture_output=True, text=True, timeout=15)
    assert run.returncode == 2
    assert "--digits must be a nonnegative integer" in run.stderr
    assert "This paper" not in run.stdout


def test_cli_zero_precision():
    run = subprocess.run([sys.executable, str(ROOT / "verify_arithmetic.py"),
                          "--digits", "0"], capture_output=True, text=True, timeout=15)
    assert run.returncode == 0, run.stdout + run.stderr
    lines = [line for line in run.stdout.splitlines()
             if "Section 3.1 (d = 200)" in line or "Theorem 2   (d = 500)" in line]
    assert len(lines) == 2
    assert all(line.split()[-1] == "3" for line in lines)


def _block(shape, symbol=0):
    return Block([Atom({(symbol,) * k}, k) for k in shape])


@pytest.mark.parametrize("source_shapes,reference_shapes", [
    ([(1, 2)], [(2, 1)]),
    ([(1, 2), (2, 1)], [(1, 2)]),
    ([(1, 2)], [(1, 2), (2, 1)]),
    ([(1, 2)], [(1, 1, 1)]),
])
def test_every_block_shape_is_checked(source_shapes, reference_shapes):
    source = BlockSet([_block(s, i) for i, s in enumerate(source_shapes)])
    reference = BlockSet([_block(s, i) for i, s in enumerate(reference_shapes)])
    with pytest.raises(ValueError, match="atom decompositions differ"):
        count_outside_neighborhood(source, reference)


def test_compatible_mixed_atom_dimensions_match_brute_force():
    source = BlockSet([_block((1, 2), 0), _block((1, 2), 3)])
    reference = BlockSet([_block((1, 2), 0)])
    expected = len(expand(source) - closed_neighborhood(expand(reference), 3))
    assert expected == 1
    assert count_outside_neighborhood(source, reference) == expected


def test_empty_source_and_reference():
    empty = BlockSet([])
    one = BlockSet([_block((1, 2), 0)])
    assert count_outside_neighborhood(empty, one) == 0
    assert count_outside_neighborhood(one, empty) == 1
    assert count_outside_neighborhood(empty, empty) == 0


def test_empty_reference_does_not_bypass_source_disjointness():
    block = _block((1, 2), 0)
    with pytest.raises(ValueError, match="not provably disjoint"):
        count_outside_neighborhood(BlockSet([block, block]), BlockSet([]))
