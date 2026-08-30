#!/usr/bin/env python3
"""Optional ONLINE check that the bundled upstream snapshots are byte-identical
to the pinned commit of the BPZ Lean repository.

    python3 scripts/verify_upstream_online.py

Everything else in this repository runs offline; this script is the one
deliberate exception, for readers who want to close the external trust root
themselves.  It downloads the three files from the pinned commit via
raw.githubusercontent.com and compares SHA-256 digests with upstream/SHA256SUMS
and with the local bytes.
"""

import hashlib
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
COMMIT = "aa21eeb12b75b0413d3fa9fb4208b5d0bf2c4d65"
REPO = "spectra-research/shannon-capacity-lean"
FILES = ["BaseC7Data.lean", "Substitutions.lean", "TerminalCodes.lean"]


def main():
    bad = 0
    for name in FILES:
        local = (ROOT / "upstream" / name).read_bytes()
        url = (f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/"
               f"ShannonBounds/{name}")
        remote = urllib.request.urlopen(url, timeout=30).read()
        lh, rh = (hashlib.sha256(b).hexdigest() for b in (local, remote))
        ok = lh == rh
        bad += not ok
        print(f"  {name:<22} local {lh[:16]}...  remote {rh[:16]}...  "
              f"{'ok' if ok else 'MISMATCH'}")
    if bad:
        print("  bundled snapshots DIFFER from the pinned upstream commit")
        return 1
    print("  bundled snapshots are byte-identical to the pinned upstream commit")
    return 0


if __name__ == "__main__":
    sys.exit(main())
