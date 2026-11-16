"""Third-party check of the Kerr orbit references against KerrGeoPy (Black Hole Perturbation Toolkit).

tests/crosscheck.rs shows that the crate reproduces tests/data/crosscheck.csv to 2e-15. This script closes the
loop for the Kerr orbit rows of that file by recomputing them with an independent package:

  * ISCO radius                kerrgeopy.constants.separatrix(a, 0, +-1)     (bracketed root finding)
  * specific energy            kerrgeopy.constants.energy(a, p, 0, +-1)
  * specific angular momentum  kerrgeopy.constants.angular_momentum(a, p, 0, +-1)
  * dt/dtau and angular rate   kerrgeopy.frequencies.mino_frequencies: Gamma = dt/dlambda with
                               dtau/dlambda = r^2 on the equator, and Omega = Upsilon_phi / Gamma

KerrGeoPy works in binary64 with root finders and elliptic integrals, so agreement is expected to about 1e-9,
not to the last bit. It does not support the extremal hole and describes only stable orbits, so rows at a = 0
or 1 and rows inside the ISCO are skipped.

Run by hand (it is not part of `cargo test`):

    pip install kerrgeopy mpmath
    python srkit/crates/srkit-gr/tools/crosscheck_kerrgeopy.py
"""
import csv
import sys
from pathlib import Path

from kerrgeopy.constants import angular_momentum, energy, separatrix
from kerrgeopy.frequencies import mino_frequencies

C = 299_792_458.0
G = 6.674_30e-11
TOLERANCE = 1e-9
# Within 1e-10 of the extremal hole KerrGeoPy's own double-precision root finder loses digits: the ISCO sits at
# 1 + (4 (1 - a))^(1/3), a cube-root singularity. At a = 1 - 1e-12 and 1 - 1e-14 an 80-digit root of the ISCO
# condition r^2 - 6r + 8a sqrt(r) - 3a^2 = 0 agrees with the crate to 2e-16 and with KerrGeoPy to 1e-8 and 3e-7.
# Those rows are therefore held to KerrGeoPy's accuracy, not the crate's.
NEAR_EXTREMAL = 1e-10
TOLERANCE_NEAR = 1e-6

path = Path(__file__).resolve().parent.parent / 'tests' / 'data' / 'crosscheck.csv'
rows = [r for r in csv.reader(l for l in open(path, encoding='utf-8') if not l.startswith('#'))][1:]
worst = {}
worst_near = {}
skipped = 0


def compare(name, ours, theirs, a, sense):
    error = abs(ours - theirs) / abs(theirs)
    near = 1.0 - a <= NEAR_EXTREMAL
    book = worst_near if near else worst
    book[name] = max(book.get(name, 0.0), error)
    if error > (TOLERANCE_NEAR if near else TOLERANCE):
        print(f'DISAGREE {name} at a = {a}, sense {sense:+.0f}: reference {ours!r}, KerrGeoPy {theirs!r}')


isco_x = {}
for function, *args, expected, _ in rows:
    if function == 'kerr::isco_radius':
        mass, a, sense = (float(v) for v in args[:3])
        isco_x[(a, sense)] = float(expected) / (G * mass / C ** 2)

for function, *args, expected, _ in rows:
    if not function.startswith('kerr::'):
        continue
    x = [float(v) for v in args if v]
    ours = float(expected)
    if function == 'kerr::isco_radius':
        mass, a, sense = x
        if not 0.0 < a < 1.0:
            skipped += 1
            continue
        compare(function, ours / (G * mass / C ** 2), separatrix(a, 0.0, sense), a, sense)
    elif function in ('kerr::isco_specific_energy', 'kerr::isco_specific_angular_momentum'):
        mass = x[0] if len(x) == 3 else None
        a, sense = x[-2:]
        if not 0.0 < a < 1.0:
            skipped += 1
            continue
        # KerrGeoPy refuses an orbit on the separatrix itself. The energy and the angular momentum are stationary in
        # the radius at the ISCO, so a step of 1e-7 outwards moves them by 1e-14.
        p = isco_x[(a, sense)] * (1 + 1e-7)
        if function.endswith('energy'):
            compare(function, ours, energy(a, p, 0.0, sense), a, sense)
        else:
            compare(function, ours / (G * mass / C), abs(angular_momentum(a, p, 0.0, sense)), a, sense)
    elif function.startswith('kerr::circular_orbit_'):
        mass, a, sense, r = x
        m = G * mass / C ** 2
        p = r / m
        if not 0.0 < a < 1.0:
            skipped += 1
            continue
        if p <= isco_x[(a, sense)] * (1 + 1e-6):
            skipped += 1   # KerrGeoPy only describes stable orbits; the crate's values inside the ISCO are checked by mpmath
        elif function.endswith('specific_energy'):
            compare(function, ours, energy(a, p, 0.0, sense), a, sense)
        elif function.endswith('specific_angular_momentum'):
            compare(function, ours / (m * C), abs(angular_momentum(a, p, 0.0, sense)), a, sense)
        elif True:
            _, _, upsilon_phi, gamma = mino_frequencies(a, p, 0.0, sense)
            if function.endswith('dt_dtau'):
                compare(function, ours, gamma / p ** 2, a, sense)
            else:
                compare(function, ours * m / C, abs(upsilon_phi / gamma), a, sense)
        else:
            skipped += 1

for name, error in sorted(worst.items()):
    print(f'{name:<52} worst relative difference {error:.1e}')
for name, error in sorted(worst_near.items()):
    print(f'{name:<52} worst relative difference {error:.1e}   (1 - a <= {NEAR_EXTREMAL:g}: KerrGeoPy is the less accurate here)')
print(f'{skipped} rows skipped (a = 0 or 1, or an orbit inside the ISCO)')
if not worst or max(worst.values()) > TOLERANCE or max(worst_near.values(), default=0.0) > TOLERANCE_NEAR:
    sys.exit(1)
print(f'all compared rows agree with KerrGeoPy to {TOLERANCE:g} ({TOLERANCE_NEAR:g} within {NEAR_EXTREMAL:g} of extremal)')
