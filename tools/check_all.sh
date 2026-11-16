#!/usr/bin/env bash
# Everything that must pass before a change is committed.
#     tools/check_all.sh                the engine: format, lint, tests, the JavaScript build and its tests
#     tools/check_all.sh --navigator    also the Navigator: its build, the committed page, parity and browser checks
# Run from anywhere. Needs the pinned Rust toolchain, Python 3 and Node; --navigator also needs Playwright.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
step() { printf '\n== %s\n' "$*"; }

cd "$root"
step "format";            cargo fmt --all -- --check
step "lint";              cargo clippy --workspace --all-targets --locked -- -D warnings
step "tests";             cargo test --workspace --locked 2>&1 | grep -E '^test result|FAILED|panicked|error' | awk '/^test result/ { p += $4; f += $6 } !/^test result: ok/ { print } END { print p " passed, " f " failed"; if (f > 0) exit 1 }'
step "JavaScript build";  python3 -I tools/build_js.py
step "binding tests";     node --expose-gc tools/test_js.js
step "native and wasm32 agree to the bit"; node tools/bits.js

if [[ "${1:-}" == "--navigator" ]]; then
  check="${TMPDIR:-/tmp}/navigator-check.html"
  step "Navigator: build";  python3 -I navigator/src/build_srkit.py
  python3 -I navigator/src/build_nav.py navigator/src/navigator.src.html navigator/src/stars.json "$check" --standalone
  step "Navigator: the committed page is this build"; cmp "$check" navigator/relativistic-navigator.html
  step "Navigator: parity"; node navigator/test/backend-parity.js; node navigator/test/reach-traveller.js
  step "Navigator: browser"; node navigator/test/browser-smoke.js
fi
printf '\nall checks passed\n'
