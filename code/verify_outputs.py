"""Check regenerated result tables and the analysis plan against recorded hashes.

Usage (from the repository root):
    python code/verify_outputs.py
"""

import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def digest(path: Path) -> str:
    data = path.read_bytes().replace(bytes([13, 10]), bytes([10]))
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    failures = 0
    for line in (ROOT / "EXPECTED_SHA256.txt").read_text().splitlines():
        expected, rel = line.split("  ", 1)
        path = ROOT / rel
        ok = path.exists() and digest(path) == expected
        failures += not ok
        print(f"{'OK  ' if ok else 'FAIL'}  {rel}")
    frozen = re.search(r"SHA-256: ([0-9a-f]{64})",
                       (ROOT / "governance" / "PROTOCOL_FREEZE.txt").read_text()).group(1)
    ok = hashlib.sha256((ROOT / "PROTOCOL.md").read_bytes()).hexdigest() == frozen
    failures += not ok
    print(f"{'OK  ' if ok else 'FAIL'}  PROTOCOL.md matches the freeze record")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
