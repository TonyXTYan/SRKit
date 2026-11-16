#!/usr/bin/env python3
"""Build SRKit and join it with the Navigator's adapter into src/srkit-backend.js.

    python3 -I navigator/src/build_srkit.py

Runs tools/build_js.py (the pinned toolchain and locked dependencies compile the Rust engine to WebAssembly and
embed it in SRKit's JavaScript loader), then appends src/srkit-navigator.js, which maps the page's calls onto SRKit.
The result is one script, generated and not tracked; build_nav.py embeds it in the page.
"""
import subprocess
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
repo = here.parents[1]
subprocess.run([sys.executable, '-I', str(repo / 'tools/build_js.py')], check=True)
engine = (repo / 'dist/srkit.js').read_text()
adapter = (here / 'srkit-navigator.js').read_text()
(here / 'srkit-backend.js').write_text(engine + '\n' + adapter)
print(f'Navigator backend: {len(engine) + len(adapter) + 1} bytes (SRKit {len(engine)}, adapter {len(adapter)})')
