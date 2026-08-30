# Pinned upstream snapshots

Verbatim copies of three files from the BPZ Lean repository,

    github.com/spectra-research/shannon-capacity-lean
    commit aa21eeb12b75b0413d3fa9fb4208b5d0bf2c4d65

redistributed under Apache-2.0 (see ../NOTICE and ../LICENSE).

These make the transcription trust root machine-checkable:
`tests/test_upstream_pinned.py` re-parses them and asserts equality with
`shannon/bpz.py` (all ten combining rules, all three terminal codes, word for
word) and with `data/base_c7.json` (the 367-word sets, the eight private
pairs, and the transversal split). `SHA256SUMS` pins the exact bytes.
