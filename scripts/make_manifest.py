#!/usr/bin/env python3
"""Generate RELEASE_MANIFEST.json for a repository release.

The manifest records SHA-256 hashes for all release files, the Python and
platform versions used to create it, and the recommended verification commands.
"""

import hashlib
import json
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    files = sorted(
        path
        for path in ROOT.rglob("*")
        if path.is_file()
        and ".git" not in path.parts
        and "__pycache__" not in path.parts
        and ".pytest_cache" not in path.parts
        and path.suffix != ".pyc"
        and path.name != "RELEASE_MANIFEST.json"
    )
    manifest = {
        "package": "c7-heterogeneous",
        "python": sys.version,
        "platform": platform.platform(),
        "verification_commands": [
            "python3 verify_arithmetic.py",
            "python3 verify_certificates.py",
            "python3 verify_arithmetic.py --recomputed",
            "python3 verify_arithmetic.py --table",
            "python3 -m pytest tests -q",
            "python3 scripts/verify_upstream_online.py  # optional, online",
        ],
        "files": {str(path.relative_to(ROOT)): sha256(path) for path in files},
    }
    output = ROOT / "RELEASE_MANIFEST.json"
    output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {output} ({len(files)} files)")


if __name__ == "__main__":
    main()
