#!/usr/bin/env bash
# Everything that must pass before a change under srkit/ is committed (engine rule 7).
#     srkit/tools/check_all.sh            the engine: format, lint, tests, the JavaScript build and its tests
#     srkit/tools/check_all.sh --apps     also the applications that carry a built copy (Navigator, alpha 3.3e3)
# Run from anywhere. Needs the pinned Rust toolchain, Python 3 and Node; --apps also needs Playwright for the browser checks.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
step() { printf '\n== %s\n' "$*"; }

cd "$root/srkit"
step "format";            cargo fmt --all -- --check
step "lint";              cargo clippy --workspace --all-targets --locked -- -D warnings
step "tests";             cargo test --workspace --locked 2>&1 | grep -E '^test result|FAILED|panicked|error' | awk '/^test result/ { p += $4; f += $6 } !/^test result: ok/ { print } END { print p " passed, " f " failed"; if (f > 0) exit 1 }'
cd "$root"
step "JavaScript build";  python3 -I srkit/tools/build_js.py
step "binding tests";     node --expose-gc srkit/tools/test_js.js
step "native and wasm32 agree to the bit"; node srkit/tools/bits.js

if [[ "${1:-}" == "--apps" ]]; then
  step "the alpha's copy is this build"; cmp srkit/dist/srkit.js prototypes/alpha-3-3e3/shared/srkit.js
  step "Navigator: build";  python3 -I docs/design/mockups/src/build_srkit.py
  python3 -I docs/design/mockups/src/build_nav.py docs/design/mockups/src/navigator.src.html docs/design/mockups/src/stars.json "${TMPDIR:-/tmp}/navigator-check.html" --standalone
  cmp "${TMPDIR:-/tmp}/navigator-check.html" docs/design/mockups/relativistic-navigator.html
  step "Navigator: parity"; node docs/design/mockups/test/backend-parity.js; node docs/design/mockups/test/reach-traveller.js
  step "Navigator: browser"; node docs/design/mockups/test/browser-smoke.js
  step "alpha 3.3e3: parity"; node prototypes/alpha-3-3e3/shared/test/srkit-parity.js | tail -1
  cd prototypes/alpha-3-3e3
  for t in sim-selftest sim-22-test sim-23-test sim-30-test sim-31-test sim-32-test ships-test content-check; do
    step "alpha 3.3e3: $t"; node "shared/test/$t.js" --no-browser | tail -1
  done
fi
printf '\nall checks passed\n'
