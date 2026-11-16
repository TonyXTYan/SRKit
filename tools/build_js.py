#!/usr/bin/env python3
"""Build SRKit for JavaScript: compile srkit-wasm, embed it in the loader, write dist/.

    python3 -I srkit/tools/build_js.py [--copy-to PATH ...]

Writes srkit/dist/srkit.js (one classic script with the module embedded), srkit/dist/srkit.wasm and
srkit/dist/srkit.d.ts. Each --copy-to PATH also writes the script there, for an application that carries its own
copy. The toolchain is pinned by rust-toolchain.toml and the dependencies by Cargo.lock; source paths are remapped,
so the same commit gives the same bytes on any machine. The SHA-256 of the module is recorded in SRKit.info.
"""
import base64
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
env = os.environ.copy()
cargo_home = Path(env.get('CARGO_HOME', Path.home() / '.cargo'))
remap = [f'--remap-path-prefix={cargo_home}=/cargo', f'--remap-path-prefix={root}=/srkit']
rustup_home = Path(env.get('RUSTUP_HOME', Path.home() / '.rustup'))
remap.append(f'--remap-path-prefix={rustup_home}=/rustup')
env['RUSTFLAGS'] = ' '.join(filter(None, [env.get('RUSTFLAGS', ''), *remap]))
target = root / 'target' / 'js'
subprocess.run(['cargo', 'build', '-p', 'srkit-wasm', '--release', '--target', 'wasm32-unknown-unknown',
                '--target-dir', str(target), '--locked'], cwd=root, env=env, check=True)
wasm = (target / 'wasm32-unknown-unknown' / 'release' / 'srkit_wasm.wasm').read_bytes()
digest = hashlib.sha256(wasm).hexdigest()
version = re.search(r'^version\s*=\s*"([^"]+)"', (root / 'Cargo.toml').read_text(), re.M).group(1)

src = (root / 'js' / 'srkit.src.js').read_text()
for mark in ('/*WASM_BASE64*/', '/*VERSION*/', '/*WASM_SHA256*/', '/*WASM_BYTES*/0'):
    assert src.count(mark) == 1, mark
js = (src.replace('/*WASM_BASE64*/', base64.b64encode(wasm).decode('ascii'))
         .replace('/*VERSION*/', version).replace('/*WASM_SHA256*/', digest).replace('/*WASM_BYTES*/0', str(len(wasm))))
js = js.replace("This file is the source; tools/build_js.py embeds the compiled module and writes dist/srkit.js.",
                f"Built by srkit/tools/build_js.py from srkit {version}; do not edit. Module SHA-256 {digest}.")
dist = root / 'dist'
dist.mkdir(exist_ok=True)
(dist / 'srkit.js').write_text(js)
(dist / 'srkit.wasm').write_bytes(wasm)
shutil.copyfile(root / 'js' / 'srkit.d.ts', dist / 'srkit.d.ts')
print(f'srkit {version}: module {len(wasm)} bytes, sha256 {digest}; dist/srkit.js {len(js)} bytes')
args = sys.argv[1:]
while args:
    if args[0] != '--copy-to' or len(args) < 2:
        sys.exit('usage: build_js.py [--copy-to PATH ...]')
    out = Path(args[1])
    out.write_text(js)
    print(f'  copied to {out}')
    args = args[2:]
