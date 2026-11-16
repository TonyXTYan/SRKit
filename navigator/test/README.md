# Navigator backend checks

Build the page first (the commands are in the Navigator's [README](../README.md)). Then, from the repository root:

```sh
node navigator/test/backend-parity.js
node navigator/test/reach-traveller.js
node navigator/test/browser-smoke.js
```

The browser test needs Playwright and uses its own Chromium; set `SRG_CHROMIUM` to use another browser, and `NODE_PATH` if Playwright is installed outside the checkout.

The three files were written for the first port of the Navigator, which used a Navigator-only Rust module. The backend they load now is SRKit's general engine behind `src/srkit-navigator.js`, so they are the acceptance test for the port: 7,145 and 926 assertions against the original equations, and 61 browser checks.

`backend-parity.js` loads the generated embedded WebAssembly adapter in a clean
Node VM. It compares trip plans, ship clocks before/during/after a trip,
light-cone inverses and intersections, orbits, optics and Lorentz boosts with
an independent copy of the original Navigator formulas. It verifies that the
compiled backend is called and that unavailable WebAssembly is an explicit
failure. The legacy fixture is test-only and never loaded by the application.

`reach-traveller.js` checks the conservative reach contour against the original
construction and its rest-frame flip-and-burn limit. Test-traveller emission
times are checked with independent monotone bisection in three dimensions,
including their waiting period before launch.

`browser-smoke.js` exercises actual controls, clock states and chart switches
on desktop and a 390 px viewport. It also opens the single-file build with
network requests blocked to verify that its embedded physics backend works
without a separate Wasm download. If a managed browser blocks the file scheme
with `ERR_BLOCKED_BY_ADMINISTRATOR`, the runner explicitly reports that direct
file navigation is unavailable and instead boots the same standalone HTML with
all HTTP requests blocked. Other navigation failures remain failures. It checks
page errors and visible numerical readouts. Device sensors, headset hardware and the optional external three.js
renderer are outside this backend-port smoke test.

The comparison tolerances accommodate binary64 differences between native
JavaScript transcendental functions and Rust `libm`, especially null-coordinate
inversion at high gamma. They are regression tolerances rather than certified
physical error bounds or a frozen bitwise determinism contract.
