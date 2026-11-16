# Relativistic Navigator

A single web page for exploring relativistic travel: plan a ride at constant proper acceleration to a real star, the Galactic centre or Andromeda, and see the clocks, the sky, the maps and the spacetime charts as the traveller would, each drawn from what light has actually reached the ship. It doubles as a shareable explainer, so its text is written for general readers.

The page is [`relativistic-navigator.html`](relativistic-navigator.html): one file, with the engine and the star table embedded. It is published from this repository with GitHub Pages (`.github/workflows/pages.yml`) at <https://tonyxtyan.github.io/SRKit/>.

## On SRKit

All of the Navigator's physics runs on SRKit, compiled to WebAssembly and embedded in the page: ride planning and state, orbit clocks and motion, light-cone intersections, the null-coordinate inverses, Doppler and aberration, boosts, reach contours and test travellers. [`src/srkit-navigator.js`](src/srkit-navigator.js) is the whole of the Navigator's side: fifty lines that convert its units (drives in g) and field names. JavaScript keeps the controls, the charts, the sampling and the drawing; the GPU sky shader is a rendering path.

The page needs no separate download and has no JavaScript physics to fall back on: without WebAssembly it says so and stops. The Debug panel names the engine's version and the first digits of its module's SHA-256.

The page began as a design mockup for *The Special Relativity Game* (`docs/design/maps-and-navigation.md` in the game's repository, where its design rationale stays) with its physics in JavaScript. That JavaScript is kept in `test/fixtures/` as an independent reference, and the tests compare the engine with it.

| File | What it is |
| --- | --- |
| `relativistic-navigator.html` | The built page, standalone. Committed so that it can be opened or served directly; CI checks that it is exactly what the sources build |
| `src/navigator.src.html` | Its source: page, styles and script, with placeholders for the star list and the engine |
| `src/srkit-navigator.js` | The adapter: the Navigator's calls mapped onto SRKit's general API |
| `src/build_srkit.py` | Builds SRKit and joins it with the adapter into the generated `src/srkit-backend.js` (not tracked) |
| `src/build_nav.py` | Injects the star list and the backend; `--standalone` adds the HTML skeleton |
| `src/build_stars.py` | Builds `src/stars.json` from SIMBAD and the IAU Catalog of Star Names |
| `src/stars.json` | 9,715 stars (V ≤ 6.5, plus everything SIMBAD has within 50 ly) in galactic coordinates. Its licence is not the code's: [`STAR-DATA.md`](STAR-DATA.md) |
| `test/` | Comparisons with the original equations, and browser checks: [`test/README.md`](test/README.md) |

## Building and checking

From the repository root, with Python 3, Node and the pinned Rust toolchain:

```sh
python3 -I navigator/src/build_srkit.py
python3 -I navigator/src/build_nav.py navigator/src/navigator.src.html navigator/src/stars.json \
  navigator/relativistic-navigator.html --standalone
node navigator/test/backend-parity.js
node navigator/test/reach-traveller.js
node navigator/test/browser-smoke.js        # needs Playwright
```

`tools/check_all.sh --navigator` does the same after the engine's own checks, and fails if the committed page is not what the sources build. Options of `build_nav.py`: without `--standalone` it writes the page body alone, for a viewer that adds its own HTML skeleton; `--stamp VERSION` writes a build stamp for the Debug panel (version, source commit and branch, build time, source hash; commit the source first so that the stamp names the commit).

To rebuild the star table from fresh downloads (about 4 MB; the second SIMBAD address is the mirror the script falls back to):

```sh
python3 -I navigator/src/build_stars.py fetch /path/to/empty/dir
python3 -I navigator/src/build_stars.py build /path/to/empty/dir navigator/src/stars.json
```

## Opening it

- **Desktop:** open `relativistic-navigator.html` in a browser. Everything works except the phone gyro and WebXR. The 360° view loads three.js r128 from cdnjs, so it needs a network connection.
- **Phone gyro and WebXR headsets:** these need the page served over HTTPS (or `localhost` on the same device) [K: browsers expose device orientation and WebXR only in secure contexts]. Options:
  - serve the folder (`python3 -m http.server 8000`) and put an HTTPS tunnel in front of it;
  - host the file on any HTTPS static host.

  Then:
  - **Phone:** open the page, press "Open 360° view", then "Phone gyro". iOS asks for motion permission. Your heading when you switch it on becomes "ahead"; "Look ahead" re-centres.
  - **Headset:** open the page in the headset's browser (Quest Browser, Safari on Vision Pro) and press "Enter VR". Use the floating panel: touch a button with a fingertip, or look at it and pinch (Vision Pro), or point and pull the trigger. The button lights up and acts when you let go; pinches just off a button snap to it. Pinch on the sky to pick a star. The panel's last line shows where your last pinch landed.
- **Inside the claude.ai viewer,** motion sensors and WebXR are blocked: the page explains this when you press those buttons. Drag-to-look and full screen still work there (full window on phones).

## Status

The table under Status records what was checked while the page was a mockup in the game's repository, mostly before the port to SRKit; `navigation_maps.py` and the design document it mentions are there. Headset and phone-sensor hardware have not been tried since the port.

| Item | Status |
| --- | --- |
| 360° view (WebGL sky sphere) | Checked in headless Chromium: dead ahead at the screen centre, not mirrored (looking right brings Hadar, Vela, the LMC and Canopus into view where they belong) |
| 360° sky on the GPU | Checked in headless Chromium (software GL): same sky at rest and squeezed ahead in flight, exposure −5.4 stops at γ ≈ 82 and about −31 at the flip of the Galactic-centre ride. JavaScript about 1–2 ms per frame, worst scan slice 3.4 ms; the old CPU texture took 390–600 ms per update. The Debug panel reports both. Frame rate on real GPUs not yet measured |
| Headset control panel | Shown in the 360° view through the Debug panel, operated with the mouse, simulated rays and simulated fingertips. A press picks its button where it starts; drifting 5 px into the next button keeps it, moving to the middle of another switches. Pinches in a gap or 50 px off the panel snap to the nearest button (2.2° away); a pinch 40 px below the panel does nothing rather than pick a star. Fingertip: hover at 5 cm lights the button, touching and lifting presses it, a hand coming through from behind does nothing, and a fingertip drags the bar (30% → 70%). visionOS gives no gaze before a pinch, so there is no hover highlight for gaze. **Not yet tried on a real headset since the first report.** |
| Playback in a headset | With the page's own frame loop stopped (as in immersive VR), the 360° loop alone advanced ship time 2.35 → 4.67 yr in 2 s |
| Spacetime in 3D / VR | Galactic, momentary-rest and light-cone frames render in the 3D view and match the 2D panel; Enter VR without a headset leaves you on the page with the reason; with simulated rays, dragging 0.8 rad round the ring turned the model 0.80 rad and dragging the plate moved it, and it stayed put; "Highlight" (default on) and "Travellers" (default off) stay in sync between the page checkboxes, the 3D view's buttons and the headset panel; turning travellers off while viewing from one falls back to the galactic frame |
| Keep playing / round trip | "⟳ Loop" beside Play on the page and in the 360° view (and on the headset panel), always in sync. Off: a ride stops on arrival. On: a transit becomes a repeating round trip (19.8 → 39.5 yr to the Galactic centre) and wraps, also with page frames stopped as in a headset. The loop stops after one lap unless it is on. Return mirrors the outbound leg to 4 × 10⁻¹⁵ |
| Sharp point stars and level of detail | Checked in headless Chromium: point stars stay sharp when zoomed; arriving inside the bulge, nearby cells resolve into individual seeded stars over the glow |
| Spacetime explorer 2+1D | Checked in headless Chromium: galactic, ship, test-traveller and picked-star frames; γ shown matches the ride (51,616 at the flip, 2.294 at 0.9c); at rest every frame coincides. Own frame (light-cone coordinates): before departure and on arrival every seen dot lands on the same pixel as in the galactic frame (1 g and 0.01 g, Proxima and the 100,000 ly crossing); in flight its current cone agrees with an independent 3D boost to 10⁻¹⁵ |
| Spacetime 1+1D (distance from you, along an axis) | Checked in headless Chromium and against `navigation_maps.py`. Distance from you (default): every body placed exactly; off the line, reception and ping times match brute-force bisection and the light-cone chart a 3D boost (~10⁻¹⁵); on the line it equals the signed chart folded. Along an axis: toward the destination it equals the signed chart (radar coordinates vs bisection and light-cone chart vs boost, ~10⁻¹⁵); toward another body the Maps panel draws the axis and bodies fade with their distance from it. At the flip of the 1 g crossing, home's radar distance 5.256 ly (SR-100) and light age in your frame 0.4843 yr (SR-97). Chart switches morph point by point in both modes; grids drawn as contours stay right at γ ≈ 13,000; 5–15 ms a draw while playing. Phone width: shorter labels, no overflow. Circular orbits: every body placed exactly (reception and ping by Newton in brackets that always hold; ~10⁻¹⁵ against brute force); the centre's light age γr and radar distance r/γ to rounding; grids and "now"s along your heading, ahead and dotted astern. Hover read-outs for every line (grid lines with value and side, the three "now"s, cones, flashes, worldlines with the time under the pointer) |
| Spacetime: scales, aspect, 2D plane | Checked in headless Chromium and against `navigation_maps.py`. Scales: log axes, log with light at 45°, Penrose (arctan, tanh) in all four charts and both 1D modes; light rays keep ΔY′ ∓ X′ to 6 × 10⁻¹⁴; the Penrose diamond with i⁺, i⁻, i⁰ and ℐ±. The time stretch (light off 45°) by slider, Shift+wheel and Alt+wheel. 2D plane: the four charts on a transit and on both orbits; radar rings along 300 random directions sit on their ping and echo cones (3 × 10⁻¹⁵ ride, 2 × 10⁻¹⁴ orbit); drag turns, Shift-drag pans, wheel zooms about the pointer, double-click fits, chart switches morph; 198 of 493 hover points over the plot name a line or body; no horizontal scroll at 390 px. A transit draw takes 4–21 ms settled and 3–12 ms playing; on an orbit up to 41 ms settled and 32 ms playing (radar chart). The orbit "floating point" bands (axis mode) are gone: bodies are traced by reception time, with directions taken where their light arrives; after arriving at Sirius (light-cone chart, line of flight) home stays on the home side, and on the line every body matches the signed chart to 4 × 10⁻¹⁵ on a round trip too. A zoomed Penrose or log view keeps you at the same place on screen while playing |
| Separate scales, orbit lines, flash densities, Things to try | Every scale and every pairing of separate space and time functions (21), in all four charts and three ways across, on the Proxima ride and the 100 ly orbit: 504 draws, no errors. In each, you are on the plot, the hover's inverse lands within 0.05 px, and every seen body sits on the drawn past cone. On orbits under a scale, redrawing with 25 times the budget changes at most 0.3 % of pixels (the orbit radar chart along the axis: 53 pixels, before 35,387). Flash densities per source up to 10 per grid step: about 25 ms settled and 15 ms playing. All 15 "Try it" presets apply their ride, time and panel settings and scroll to their panel, with no errors and no horizontal scroll at 390 px |
| Ship time steps and beyond the ride | Steps of 0.1 yr to 1 Myr either way, Play/Pause, Launch, Flip and Arrival in the right-hand column. Ship time before launch (waiting at home) and after arrival (at rest there), from −1 Myr to +1 Myr, on a one-way ride, a round trip and an orbit: every panel draws, with no errors and no NaN in the read-outs. Leaving the 2D plane lands exactly on the 1D chart (0.0 px after panning and zooming) |
| Animated changes of the way across, build stamp | Filmed frame by frame in headless Chromium (frozen at 0, ¼, ½, ¾ and 1): distance ↔ axis morph point by point; axis, distance → 2D plane in stages (the lines the 1D chart has move into the plane, then it tilts, then the 3D-only lines fade in; the start frame puts every body on the 1D chart's pixels, 0.00 px); 2D plane → axis or distance the same backwards (lands on the 1D chart to 0.00 px, also after panning and zooming); on a transit and an orbit. Each ends in the mode asked for, also in real time and when the choice changes mid-tilt; no errors. The Debug panel shows the stamp (or "unstamped" for local builds) |
| Layout and theme | The ride strip sits below the sky sphere. Two columns at 1,100 px and up, with the ride, drive, common destinations, theme, trip card and panel switches in a sticky right-hand column (stays at the top while the page scrolls); one column below, with no horizontal scroll at 390 px or 1,000 px. Auto, Dark and Light themes; the maps, sky and spacetime canvases stay on black, and so do the map's zoom buttons and the 360° view's controls |
| Chart zoom and pan, seat-time bar | Checked in headless Chromium: wheel zoom about the cursor, drag pan (a drag never picks a star), a synthetic two-finger pinch, the + and − buttons, "Reset view"; the bar drives every panel in transit and loop mode; no horizontal scroll at 390 px |
| Maps: ship frame, now | Four frames on the Maps panel (renamed from Chart). The "now" positions match an independent Lorentz boost: t′ = 0 on the slice, relative position error ≤ 2 × 10⁻⁸ at γ = 13,421; at rest the view equals the galactic frame. Mid-ride with heading up the galaxy squashes to a line through you; the four frame buttons wrap cleanly at 390 px |
| Maps: ring labels | Light-years on every map except the galactic frame centred on you, which keeps "seen N yr ago". Ship frame, as seen: rings and labels read "looks N ly away" (Galactic centre at 4% of the ride: relay 26,000 yr ago, looks 58,786 ly away = 26,000 × D 2.26); the hover card adds the same line |
| Maps: ship frame, outside view | Galactic positions with the ship's light-age ellipses. At rest the contours are the galactic frame's circles. In flight they flatten ahead and stream astern; on the log scale contours start where their front tip clears the 0.1 ly centre, and labels stagger along the down-right diagonal like the circles'. The ellipse geometry (semi-major γA, semi-minor A, t′ = −A on every point) is checked in `navigation_maps.py`. No horizontal scroll at 390 px |
| Maps: seen from above | Over one lap of the 100 ly orbit round the Sun (300 ly zoom, galactic up), six named stars 20–94 ly off the plane drifted up to 89 px relative to Earth with the old projection and 0 px now, in the galactic, outside-view and orbit-centre frames. At rest all five frames still put M31, the Galactic centre and Omega Centauri at the same pixels. Brightness → Height off the plane dims off-plane stars and fades or drops their labels |
| Maps: future light cone | Checkbox shown only on the outside view, off by default. Orange dashed ellipses mirrored from the past ones, labelled "+N yr" up and to the right, never on the abeam spot of the past label of the same age. Ring labels now stay clear of the zoom buttons |
| Map plane and reach contours | The map plane holds the line of flight, with least tilt from the galactic plane: 54.5° for υ PsA, 86.5° for HIP 3056, 0° for the Galactic centre and for circular orbits. It swings to a new destination over 0.9 s. Ride to the Galactic centre: identical to the previous build to 0 px in all five frames. Six off-plane destinations, all five frames, four stages of the ride: 0 crossings between consecutive reach contours in 467 pairs (46 in the build that sliced by the galactic plane). Plan isochrone through the destination at every stage (gap ~10⁻¹³ px). Circles at launch (radius constant to 0.01 px) |
| Maps: free 3D camera | "3D" button: orthographic camera starting 34° above the map plane; drag turns and tilts, wheel and pinch zoom, "Reset view" levels it. Rings become ellipses; stalks drop from labelled stars to the plane |
| Maps: ring transitions | Frozen at 0, 25, 50, 75 and 100% of the transition: zooming from 120 kly to 3 kly, the old rings expand away and fade while the new ones grow from the centre; switching from the galactic frame to the outside view, circles morph into the ship's ellipses; switching linear to log, shared values glide to their new radii |
| Landmarks and waypoints | Real objects at rounded literature positions and distances; Andromeda's waypoints use the same disc frame as its seeded points (NGC 206 lands beside the 10 kpc waypoint at the south-west end) |
| Reach from here | Checked in headless Chromium and against `navigation_maps.py`: at rest the contours equal the seat-years rings to machine precision; the plan isochrone passes through the destination (relative miss about 10⁻¹⁶); after the flip the brake-now point is the destination. Dynamic density keeps at least three contours astern of the ship at every tested stage, frame and destination |
| Clock table | Checked in headless Chromium: per ship-second, day, month and year for home, destination and the launch-site calendar; at the flip of the Galactic-centre ride one ship-day shows 3.22 s of home and 73.5 yr of the centre |
| Panel toggles | Hiding the chart and sky views cuts a redraw from about 230 ms to 5 ms |
| Lazy drawing, Debug panel | Checked in headless Chromium: off-screen panels skip seat-time redraws and catch up when scrolled into view; a seat-time change at the top costs about 46 ms instead of 427 ms; the slider keeps the same width at all 61 tested seat times (809.6 px at 1280 px wide, 306 px at 390 px) |
| Phone gyro, WebXR | Code paths written to the standard APIs, **not yet tried on a real phone or headset** |
| Physics | The numbers come from `navigation_maps.py` and the navigation rows of the oracle (`oracle/oracle.csv`: SR-74 to SR-79 and SR-88 to SR-110; the status notes above quote them) |

## Credits

- **Star data:** the SIMBAD database, operated at CDS, Strasbourg, France (Wenger et al. 2000, A&AS 143, 9), under the Open Database License 1.0. Star names: the IAU Catalog of Star Names (IAU Working Group on Star Names), CC BY. Details and terms: [`STAR-DATA.md`](STAR-DATA.md). Masses in the table are a placeholder main-sequence fit from luminosity, used only as label priority.
- **Seeded galaxies:** the Milky Way and Andromeda are procedural component models (described in the game's design document), with landmarks at rounded literature positions.
- **3D library:** three.js r128 (MIT), loaded from cdnjs at run time for the 360° view.
- **Code:** MIT, as the rest of this repository; written with Claude Code.
