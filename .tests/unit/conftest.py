"""Make workflow/scripts importable as plain modules.

The scripts are standalone CLIs, not a package -- they are invoked as
`python3 workflow/scripts/foo.py` by the rules. They each do a
`sys.path.insert(0, dirname(__file__))` so they can import their siblings
(gz_io), which means adding that directory here is enough to import them.
"""
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
