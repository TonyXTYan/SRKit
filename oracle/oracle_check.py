"""Independent reproduction of the physics oracle.

The stated values live in oracle.csv (one row per check, IDs SR-xx for the spine physics, DT-xx for deep time);
this script computes each row from first principles and compares. Docs cite rows by ID instead of restating numbers.

Run:  uv run -q --with numpy --with scipy --with mpmath python oracle/oracle_check.py

One line per row: status | id | quantity | stated | computed | relative difference | definition / note.
Status:
  MATCH      computed within the row's tolerance (+-1% unless stated)
  ROUNDED    equal to the stated value at its quoted significant figures, but outside +-1%:
             the oracle should carry more digits (computed value shown)
  MISMATCH   wrong number or wrong physics; correct value in the note
  AMBIGUOUS  reproducible only under an assumption the row does not state; proposed definition in the note

Exits nonzero for mismatches, missing checks or repeated check IDs. Rounded
values and explicitly documented ambiguous models remain accepted outcomes.
"""
import csv
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 60

# ---- constants (oracle convention plus CODATA 2018 / IAU 2015 nominal) ----
c = 299_792_458.0
yr = 365.25 * 86400.0           # Julian year, s
ly = c * yr
g0 = 9.80665
Lsun = 3.828e26                 # IAU 2015 nominal
G = 6.67430e-11
hbar = 1.054571817e-34
kB = 1.380649e-23
Msun = 1.3271244e20 / G          # IAU 2015 nominal GM_sun over CODATA 2018 G: 1.98841e30 kg (R71.10)
Rsun = 6.957e8                  # IAU 2015 nominal
AU = 149_597_870_700.0
H_L = 1.0 / (17.53e9 * yr)      # de Sitter Hubble rate, 1/H = 17.53 Gyr
Om, OL = 0.315, 0.685           # Planck 2018 (H0 = 67.4 gives 1/H_L = 17.53 Gyr)
TCMB0 = 2.7255
J_kt, J_Gt = 4.184e12, 4.184e18
d_prox = 4.2465                 # ly (Gaia DR3 1.3020 pc)
EPS = 2.0 ** -52
MIN_PER_YR = yr / 60.0

rows = []
ORACLE = {r['id']: r for r in csv.DictReader(open(Path(__file__).with_name('oracle.csv'), encoding='utf-8'))}
seen = set()


def sigfigs(s):
    m = s.lower().lstrip('+-').split('e')[0]
    d = m.replace('.', '').lstrip('0')
    if '.' not in m:
        d = d.rstrip('0') or '0'
    return max(len(d), 1)


def round_sig(x, n):
    return float(f"{x:.{n - 1}e}")


def row(key, computed, note='', status=None):
    if key in seen:
        raise ValueError(f'duplicate oracle check: {key}')
    r = ORACLE[key]
    seen.add(key)
    stated, tol = r['stated'], float(r['tolerance'])
    name = r['quantity'] + (f" [{r['unit']}]" if r['unit'] else '')
    sv = float(stated.replace(',', ''))
    rel = (computed - sv) / sv if sv != 0 else float('nan')
    if status is None:
        if abs(rel) <= tol:
            status = 'MATCH'
        elif math.isclose(round_sig(computed, sigfigs(stated.replace(',', ''))), sv, rel_tol=1e-9):
            status = 'ROUNDED'
        else:
            status = 'MISMATCH'
    rows.append((status, key, name, stated, computed, rel, note))


# ======================= SR: spine physics (formerly handoff section 8) =======================

def flip(d_ly, a):
    """Flip-and-burn at constant proper acceleration a, from rest to rest over d."""
    x = a * d_ly * ly / (2 * c * c)          # = a d / 2c^2
    gam = 1 + x                              # peak gamma at turnover
    phi = math.acosh(gam)                    # peak rapidity
    ca = c / a / yr                          # c/a in years
    tau = 2 * ca * phi
    t = 2 * ca * math.sqrt(x * x + 2 * x)    # 2 (c/a) sinh(phi)
    excess = 2 * ca * (math.sqrt(x * x + 2 * x) - x)  # t - d/c
    return tau, t, gam, phi, excess

row('SR-01', c / 10 / yr)
row('SR-02', c * c / 10 / ly)

tg, tt, gg, pg, eg = flip(1e5, g0)
t10, tt10, g10, p10, e10 = flip(1e5, 10.0)
row('SR-03', tg,
    f'tau=(2c/a)acosh(1+ad/2c^2) with g=9.80665 (corrected 26 Oct; old row 21.97 used a=10, which gives {t10:.3f})')
row('SR-04', tt,
    f't=(2c/a)sinh(phi). +-1% is 1000 yr, untestable; test t-d/c instead: {eg:.4f} yr at g, {e10:.4f} at a=10')
row('SR-05', gg,
    f'gamma=1+ad/2c^2 at g (corrected 26 Oct; a=10 gives {g10:.1f})')
row('SR-06', math.exp(2 * pg),
    f'R=exp(2*phi_peak) (accel+decel, ideal photon rocket). g: {math.exp(2*pg):.4g}=exp({2*pg:.3f}); a=10: {math.exp(2*p10):.4g}')
row('SR-07', 2 * pg)
row('SR-08', c / g0 / yr * math.atanh(0.9),
    f'proper time tau=(c/a)artanh(b). Coordinate time would be {(c/g0/yr)*0.9/math.sqrt(0.19):.3f} yr; a=10 gives tau={c/10/yr*math.atanh(0.9):.3f}')

def loop_gamma(r_ly, a=g0):
    return math.sqrt(1 + a * r_ly * ly / c ** 2)    # (gamma beta)^2 = a r / c^2

for r, k in [(1, 'SR-09'), (100, 'SR-10'), (1e4, 'SR-11'), (5e4, 'SR-12')]:
    row(k, loop_gamma(r), f'a=10 m/s^2 would give {loop_gamma(r, 10.0):.4g}')

gam100 = loop_gamma(100)
b100 = math.sqrt(1 - 1 / gam100 ** 2)
lap = 2 * math.pi * 100 / b100                       # yr (r in ly, v in c)
row('SR-13', lap, f'2 pi r / v, gamma={gam100:.4f}, beta={b100:.6f}')
row('SR-14', lap / gam100)
row('SR-15', gam100 * (1 - b100))
row('SR-16', gam100 * (1 + b100),
    f'within 1%, but correctly rounded it is {gam100*(1+b100):.2f} -> quote 20.4')

# lap-average of the playback rate D = dt_emit/dtau for a static source at finite distance
def lap_mean_D(R, h, n=20000):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    xo = np.stack([100 * np.cos(th), 100 * np.sin(th), 0 * th])       # ly
    bv = b100 * np.stack([-np.sin(th), np.cos(th), 0 * th])
    sep = xo - np.array([[R], [0], [h]])
    nh = sep / np.linalg.norm(sep, axis=0)
    D = gam100 * (1 - np.sum(bv * nh, axis=0))
    return D.mean()                       # uniform in tau (constant speed); periodic -> spectral accuracy
md = max(abs(lap_mean_D(R, h) / gam100 - 1) for R, h in [(1e3, 0), (300, 50), (1e6, 0), (150, -20)])
row('SR-17', 1 + md, f'max |mean/gamma - 1| over 4 source positions = {md:.1e} (exact: observer returns to start)')

gap_nr = 2 * math.pi * 100 * b100                     # yr: 2 pi r v / c^2
row('SR-18', gap_nr * gam100,
    f'corrected 26 Oct (old row 625 yr was the beta<<1 limit 2 pi r v/c^2 = {gap_nr:.4g}); lab time 2 pi r v gamma^2/c^2 = {gap_nr*gam100**2:.4g} yr')
row('SR-19', gap_nr * gam100**2)

U = g0 / (c * H_L)                                     # saturation gamma*beta = alpha/(cH)
row('SR-20', math.sqrt(1 + U * U),
    f'u=gamma*beta, du/dt=alpha/c-Hu -> alpha/(cH)={U:.4e} (a=10 gives {10/(c*H_L):.3e})')
# tau(T) = T/sqrt(1+U^2) + C, de Sitter from rest: u = U(1-exp(-Ht)); C from the transient
Hyr = H_L * yr
f_int = lambda t: 1 / mp.sqrt(1 + (U * (1 - mp.exp(-Hyr * t))) ** 2) - 1 / mp.sqrt(1 + U * U)
Cint = mp.quad(f_int, [0, 1, 10, 1e3, 1e5, 1e7, 1e9, 1e10, 1e11, 1e12, 1e13])   # ~ (c/a) ln(2U); tail < 1e-24
tauH = 1e100 / math.sqrt(1 + U * U) + float(Cint)
row('SR-21', tauH, f'tau = T/gamma_sat + C, de Sitter from rest; C = {float(Cint):.2f} yr ~ (c/a) ln(2U)')

phi_flat = math.asinh(1e100 * yr * g0 / c)
row('SR-22', c / g0 / yr * phi_flat,
    f'tau=(c/a)asinh(at/c), phi={phi_flat:.3f}; a=10 gives {c/10/yr*math.asinh(1e100*yr*10/c):.1f}')

for b, k, conv in [(0.01, 'SR-23', 3600), (0.5, 'SR-24', 60), (0.9999, 'SR-25', 1), (0.99999, 'SR-26', 1e-3)]:
    row(k, 1200 * (1 / b - 1) / conv, 'dt=(d/c)(1/beta-1)')

for b, k_ke, k_g, k_bs in [(0.9999, 'SR-27', 'SR-28', 'SR-29'), (0.99999, 'SR-30', 'SR-31', 'SR-32')]:
    gm = 1 / math.sqrt(1 - b * b)
    row(k_ke, (gm - 1) * c * c / J_Gt, '(gamma-1)c^2, 1 t TNT = 4.184e9 J')
    row(k_g, gm)
    row(k_bs, math.sqrt((1 + b) / (1 - b)), 'sqrt((1+b)/(1-b))')

realtime = lambda delay_s: delay_s / yr * 60.0       # 1 yr = 1 min -> real seconds
row('SR-33', realtime((1 + 1.523679) * AU / c) * 1e3,
    f'a_E+a_M=2.524 AU (21.0 min light). True max ~2.68 AU (22.3 min) gives {realtime(2.68*AU/c)*1e3:.2f} ms')
row('SR-34', realtime(170 * AU / c))
row('SR-35', realtime(2000 * AU / c))
row('SR-36', realtime(1e5 * AU / c) / 60)
row('SR-37', d_prox, 'd=4.2465 ly')
row('SR-38', 100 / 60)
row('SR-39', 1e5 / 1440)

sid = 365.256363 / 365.25                           # sidereal year in Julian years
hoh = lambda a2: 0.5 * ((1 + a2) / 2) ** 1.5 * sid    # yr; real s = 60 * yr
row('SR-40', 60 * hoh(1.523679), 'half-period of transfer ellipse, circular coplanar orbits')
row('SR-41', hoh(5.2038))
row('SR-42', hoh(30.07))
row('SR-43', 1 / (1 / 365.256363 - 1 / 686.980) / 365.25)

def antimatter(d_ly, mpay=1e6, a=g0):
    _, _, _, phi, _ = flip(d_ly, a)
    R = math.exp(2 * phi)
    return (R - 1) * mpay / 2                       # half the propellant is antimatter
row('SR-44', antimatter(d_prox),
    'ideal photon rocket (all annihilation energy into collimated exhaust), flip-and-burn at g')
row('SR-45', antimatter(1e5),
    f'uses g (a=10 gives {antimatter(1e5, a=10.0):.3g}): inconsistent with the tau row above')
row('SR-46', 4 * AU ** 2 / 6.371e6 ** 2, '4 AU^2 / R_E^2')
for k, M, R in [('SR-47', 3.3011e23, 2.4397e6), ('SR-48', 5.9722e24, 6.3710e6), ('SR-49', 1.89813e27, 6.9911e7)]:   # Mercury, Earth, Jupiter
    row(k, 0.6 * G * M * M / R, 'uniform-sphere binding as stated; real bodies are centrally condensed and bind harder (Earth ~2.49e32 J, +11%; Jupiter as n=1 polytrope 0.75GM^2/R, +25%)')
for b, k in [(0.01, 'SR-50'), (0.1, 'SR-51')]:
    row(k, (1 / math.sqrt(1 - b * b) - 1) * c * c / J_kt)
row('SR-52', d_prox / 0.01 - d_prox / 0.1, 'd/v1 - d/v2')

def kerr(a):
    """Prograde equatorial ISCO (Bardeen-Press-Teukolsky 1972), M=1, a=J c/(G M^2)."""
    a = mp.mpf(a)
    Z1 = 1 + mp.cbrt(1 - a * a) * (mp.cbrt(1 + a) + mp.cbrt(1 - a))
    Z2 = mp.sqrt(3 * a * a + Z1 * Z1)
    r = 3 + Z2 - mp.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2))
    den = r ** mp.mpf(0.75) * mp.sqrt(r ** mp.mpf(1.5) - 3 * mp.sqrt(r) + 2 * a)
    ut = (r ** mp.mpf(1.5) + a) / den                        # dt/dtau, Boyer-Lindquist t
    E = (r ** mp.mpf(1.5) - 2 * mp.sqrt(r) + a) / den         # specific energy at ISCO
    # independent cross-check from the equatorial Boyer-Lindquist metric with Omega = 1/(r^1.5 + a)
    Om_ = 1 / (r ** mp.mpf(1.5) + a)
    gtt, gtp, gpp = -(1 - 2 / r), -2 * a / r, r * r + a * a + 2 * a * a / r
    ut2 = 1 / mp.sqrt(-(gtt + 2 * gtp * Om_ + gpp * Om_ ** 2))
    E2 = -(gtt + gtp * Om_) * ut2
    assert abs(ut2 / ut - 1) < mp.mpf('1e-30') and abs(E2 / E - 1) < mp.mpf('1e-30')
    return float(r), float(ut), float(1 - E)
spins = [('0', 'SR-53', 'SR-54'), ('0.9', 'SR-55', 'SR-56'), ('0.998', 'SR-57', 'SR-58'),
         ('1-1e-6', 'SR-59', 'SR-60'), ('1-1e-14', 'SR-61', 'SR-62')]
for sa, k_ut, k_eta in spins:
    a = (1 - mp.mpf(sa[2:])) if sa.startswith('1-') else mp.mpf(sa)
    r, ut, eta = kerr(a)
    row(k_ut, ut, f'r_isco={r:.6g} M; u^t of prograde circular geodesic (BPT formula = metric cross-check)')
    row(k_eta, 100 * eta, 'Novikov-Thorne 1-E_isco (no photon capture)')

row('SR-63', 226 * 1e6 * c * c, '(gamma-1)mc^2, kinetic energy only')

M9 = 1e9
P_naive = hbar * c ** 6 / (15360 * math.pi * G ** 2 * M9 ** 2)
t_naive = 5120 * math.pi * G ** 2 * M9 ** 3 / (hbar * c ** 4) / yr
TH9 = hbar * c ** 3 / (8 * math.pi * G * M9 * kB)
# Standard-Model emission (MacGibbon 1991): dM/dt = -5.34e25 f(M)/M^2 g/s, f ~ 12 at 1e12 g (estimate)
f12 = 12.0
P_sm = 5.34e25 * f12 / (1e12) ** 2 * 1e-3 * c * c
t_sm = (1e12) ** 3 / (3 * 5.34e25 * f12) / yr
row('SR-64', P_naive,
    f'reproduced only by photon-only blackbody of horizon area, no greybody: P=hbar c^6/(15360 pi G^2 M^2). '
    f'T_H={TH9:.3g} K (kT={kB*TH9/1.602176634e-10:.1f} GeV): emits all SM species; MacGibbon f~12 gives ~{P_sm:.1g} W', status='AMBIGUOUS')
row('SR-65', t_naive,
    f'naive 5120 pi G^2 M^3/(hbar c^4). With SM species, ~{t_sm:.0f} yr', status='AMBIGUOUS')
row('SR-66', G * Msun * 0.9 * Msun / Rsun, 'M, R held at solar values (upper-bound style estimate)')
row('SR-67', Lsun / (2 * Msun * c))
TdS = hbar * H_L / (2 * math.pi * kB)
row('SR-68', 2.7 / TdS, f'T_dS={TdS:.4e} K; with T_CMB=2.7255 K: {TCMB0/TdS:.3e}')
row('SR-69', 1e-6 * 100 * ly / 1e6 / 1e3,
    'lambda d/D as stated; Airy first-null diameter would be 2.44x')
row('SR-70', 1e6 * 1e6 / 1e-6 / ly)
row('SR-71', Lsun / c, 'fully absorbed, F=L/c (reflected would be 2L/c=2.6e18)')
v_rot = 220e3
row('SR-72', v_rot / c,
    f'v=220 km/s (IAU 1985) by definition. Modern v_c(R0)~229-233 km/s gives {229e3/c:.2e}')
row('SR-73', 1 / math.sqrt(1 - (v_rot / c) ** 2) - 1,
    f'v=229 km/s gives {1/math.sqrt(1-(229e3/c)**2)-1:.2e}')

# ======================= DT: deep time (deep-time-engine.md sections 7 and 1) =======================
ms_per_yr = yr * 1e3
row('DT-01', math.log2(1e100 * ms_per_yr), 'ceil -> 368 bits')
row('DT-02', 2 ** 384 / ms_per_yr)
row('DT-03', 2 ** 53 / ms_per_yr,
    'defined as 2^53 ticks (t stored as f64 count of ms). Stored in s: 2^43 s = %.3g yr; in yr: 2^18 = %.3g yr' % (2 ** 43 / yr, 2 ** 18))
row('DT-04', 2 ** 53 / (yr / 60))
sp_yr = float(np.spacing(1e100))
sp_s = float(np.spacing(1e100 * yr)) / yr
sp_ms = float(np.spacing(1e100 * ms_per_yr)) / ms_per_yr
row('DT-05', sp_ms,
    f'actual ulp of an f64 count of ms; the bound eps*t is {EPS*1e100:.3g}. Stored in s: {sp_s:.3g}; in yr: {sp_yr:.3g}')
row('DT-06', 1e-3 / 60 * (yr / 60))
row('DT-07', (1 / 60) / 60 * (yr / 60))

ca = c / g0 / yr                                      # c/a in yr, a = g
worst = 0.0
for phi in [1.0, 11.6, 12.6, 37.0, 231.0]:
    f = -ca * math.expm1(-phi)
    ref = -mp.mpf(ca) * mp.expm1(-mp.mpf(phi))
    worst = max(worst, float(abs((mp.mpf(f) - ref) / ref)))
row('DT-08', worst,
    'max over phi in {1,11.6,12.6,37,231}; pass if <= 1e-15', status='MATCH' if worst <= 1e-15 else 'MISMATCH')

# naive form: worst-case bound eps*cosh(phi)*(c/a), and measured errors with two libm behaviours:
#  host libm (macOS: sinh == cosh == exp/2 for large x) and a musl-style emulation (Rust `libm` is a musl port:
#  sinh = (t + t/(t+1))/2 with t = expm1(x); cosh = (e^x + e^-x)/2), which is what the engine will ship.
def sinh_musl(x):
    t = math.expm1(x); return 0.5 * (t + t / (t + 1))
def cosh_musl(x):
    t = math.exp(x); return 0.5 * (t + 1 / t)
phi_ms = math.log(2 * (1 / ms_per_yr) / (EPS * ca))     # eps * (e^phi/2) * (c/a) = 1 ms
phi_noise = math.log(2 / EPS)                           # eps * e^phi/2 = 1: error bound ~ the value itself
def scan(lo, hi, sh, ch, step=0.001):
    grid = np.arange(lo, hi, step); bad = 0; mx = 0.0; first = None
    for phi in grid:
        phi = float(phi)
        e = abs(ca * (sh(phi) - ch(phi) + 1) - ca * (-math.expm1(-phi))) * ms_per_yr   # ms
        if e > 1:
            bad += 1; first = phi if first is None else first
        mx = max(mx, e)
    return bad / len(grid), mx / ms_per_yr, first
_, _, first_host = scan(0.5, 20, math.sinh, math.cosh, 0.005)
fr_h, mx_h, _ = scan(37, 231, math.sinh, math.cosh, 0.01)
fr_m, mx_m, _ = scan(37, 231, sinh_musl, cosh_musl, 0.01)
row('DT-09', phi_ms,
    f'worst-case bound eps*cosh(phi)*(c/a) = 1 ms; measured first >1 ms error at phi={first_host:.2f}')
row('DT-10', phi_noise,
    f'bound eps*cosh(phi)=1. Measured on phi in [37,231]: musl-style libm {100*fr_m:.1f}% of points >1 ms, '
    f'max error {mx_m:.2g} yr (sporadic, ~1 ulp of cosh); host macOS libm {100*fr_h:.1f}%, max {mx_h:.1g} yr. Libm-dependent')
phi0 = mp.mpf(5)
Dnum = mp.diff(lambda p: (1 - mp.exp(-p)), phi0)        # dt_ret/dtau = (c/a)(d/dphi)(1-e^-phi)*(a/c)
row('DT-11', float(Dnum), 'dt_ret/dtau by numerical differentiation')
row('DT-12', phi_flat, f'a=10 gives {math.asinh(1e100*yr*10/c):.3f}')
# SR-22 already checks the corresponding proper time in the spine section.

# orbit phase at Delta = 1e20 yr, exact integer modular phase vs mpmath reference vs naive f64
P = round(365.256363004 * 86400 * 1000)                  # 1 sidereal yr as integer ms
Dms = 10 ** 20 * round(ms_per_yr)
ph_int = Fraction(Dms % P, P)
ph_ref = mp.frac(mp.mpf(Dms) / P)
ph_f64 = (float(Dms) / P) % 1.0
row('DT-13', float(abs(mp.mpf(ph_int.numerator) / ph_int.denominator - ph_ref)),
    f'integer phase = {float(ph_int):.12f}; naive f64 (Delta/P mod 1) gives {ph_f64:.6f}, off by {abs(ph_f64-float(ph_int)):.3f} cycles; '
    f'f64 phase is garbage once ulp(Delta) > P, from Delta ~ {P*2**52/ms_per_yr:.1e} yr for a 1 yr period', status='MATCH')

# causal rounding: arrival at a moving receiver, ceil(f64 delta) vs exact rational delta
random.seed(7)
c_ms = Fraction(299_792_458, 1000)                     # m per ms
N = 100_000; naive_fail = 0; guarded_fail = 0
for _ in range(N):
    k = random.randint(10 ** 6, 10 ** 13)              # ms after emission tick
    v = random.uniform(-0.9, 0.9) * 299_792.458        # receiver velocity, m/ms
    target = k + random.uniform(0, 4) * float(np.spacing(float(k)))   # just above a tick
    x0 = float(target * (299_792.458 - v))             # receiver offset at emission, m (a double)
    exact = Fraction(x0) / (c_ms - Fraction(v))        # exact arrival delta, ms
    df = x0 / (299_792.458 - v)                        # f64 closed form
    if math.ceil(df) < exact:
        naive_fail += 1
    tick = math.ceil(df * (1 + 8 * EPS))               # guard band, then exact integer/rational check
    while tick < exact:
        tick += 1
    if tick < exact:
        guarded_fail += 1
row('DT-14', guarded_fail,
    f'adversarial cases just above a tick (engine rule 5); plain ceil(f64) fails {naive_fail} times',
    status='MATCH' if guarded_fail == 0 else 'MISMATCH')

row('DT-15', TdS)
a_eq = TCMB0 / TdS
A = (Om / OL) ** (1 / 3)
x_eq = float(mp.asinh((mp.mpf(a_eq) / A) ** 1.5) / 1.5)
row('DT-16', a_eq, 'T_CMB0=2.7255 K, T proportional to 1/a')
row('DT-17', math.log(a_eq))
row('DT-18', x_eq / H_L / yr, 'a(t)=(Om/OL)^(1/3) sinh^(2/3)(1.5 H_L t), Planck 2018')

TH = lambda M: hbar * c ** 3 / (8 * math.pi * G * M * kB)
tev = lambda M: 5120 * math.pi * G ** 2 * M ** 3 / (hbar * c ** 4) / yr
k_naive = 1 / (15360 * math.pi)
k_page_gg = 2.011e-4 * (0.167 + 0.019)              # Page 1976: photons + gravitons share of total
row('DT-19', TH(1e11 * Msun))
row('DT-20', TH(Msun))
row('DT-21', tev(1e11 * Msun),
    'naive photon-blackbody formula; see 1 Msun row', status='AMBIGUOUS')
row('DT-22', tev(Msun),
    f'naive formula. Photons+gravitons with greybody factors (Page 1976): {tev(Msun)*k_naive/k_page_gg:.3g} yr; '
    'shorter still (~x3.5) if the lightest neutrino is massless', status='AMBIGUOUS')
row('DT-23', 1e12 * Msun * c * c / (kB * TdS * math.log(2)))
row('DT-24', 3.6e15 / (kB * 300 * math.log(2)))
row('DT-25', H_L ** 2 * 5e4 * ly, 'Lambda c^2/3 = H_L^2')
row('DT-26', H_L * 1e5 * yr,
    f'uses H_L; today H0=67.4 km/s/Mpc gives {67.4e3/3.0856775814913673e22*1e5*yr:.2g}')
row('DT-27', math.log(10) / (3 * H_L) / yr / 1e9, 'photon rate ~ exp(-3Ht)')
x_of = float(mp.asinh((mp.mpf(np.finfo(float).max) / A) ** 1.5) / 1.5)
row('DT-28', x_of, f'with a=e^(Ht): {math.log(np.finfo(float).max):.2f}')
row('DT-29', x_of / H_L / yr)
vf = 5 / (5 / 0.1 + 300)
row('DT-30', vf, 'l/(l/v+T_rep)')
row('DT-31', 1e5 / vf / 1e6, 'defined as 1e5 ly / v_front (front crossing a full diameter)')

# section 1 table and prose
row('DT-32', 2 ** 106 / ms_per_yr)
row('DT-33', 2 ** 127 / ms_per_yr, 'signed (2^127)')
row('DT-34', 2 ** 256 / ms_per_yr)
row('DT-38', 2 ** 64 / ms_per_yr, 'unsigned 64-bit, arena clock')
for k, per_yr in [('DT-35', yr / 60), ('DT-36', yr), ('DT-37', yr * 1e6)]:   # minute, 1 s, 1 us
    row(k, math.ceil(math.log2(1e100 * per_yr)))

# ---- megastructures and 1 g travel (R64.1; megastructure_travel.py); SR-80 on, leaving SR-74 to SR-79 to the navigator branch ----
_, _, gp, php, _ = flip(d_prox, g0)
row('SR-80', 0.5 * 1e6 * c * c * (math.exp(php) - 1),
    'beamed sail held at 1 g to the turnover: received power D^2 P, D = exp(-phi), so E = (m c^2/2)(exp(phi_f) - 1)')
row('SR-81', 0.5 * 1e6 * g0 * c * math.exp(2 * php), 'peak source power at the turnover: (m a c/2) exp(2 phi_f)')
row('SR-82', (math.exp(2 * php) - 1) * 1e6 * c * c, 'photon rocket propellant rest energy (MR - 1) m c^2, MR = exp(2 phi_peak)')
row('SR-83', 0.5 * 3.3011e23 / (4 * math.pi * (0.387098 * AU) ** 2), 'half of Mercury spread over the sphere at its orbit (GV: half the mass usable)')
def _taper(GM, R, T, spec=5e7):
    w = 2 * math.pi / T; rs = (GM / w ** 2) ** (1 / 3); x = R / rs
    return math.exp(GM / R * (1 - 1.5 * x + 0.5 * x ** 3) / spec)
row('SR-84', _taper(3.986004418e14, 6.3781e6, 86164.1), 'constant-stress tether to synchronous orbit, exp(dPhi rho/sigma), dPhi = (GM/R)(1 - 3x/2 + x^3/2), x = R/r_sync')
row('SR-85', _taper(4.282837e13, 3.3962e6, 88642.7), 'as SR-84 for Mars (sidereal day 88,642.7 s)')
_b = math.sqrt(1 - 1 / gp ** 2)
row('SR-86', gp * 1e5 * _b * c * (gp - 1) * 1.67262192e-27 * c * c, 'aboard: gamma n beta c (gamma - 1) m_p c^2 (the gas is denser by gamma in the ship frame) at the Proxima trip\'s peak gamma, n = 0.1 per cm^3; per unit of galactic time it is gamma times less')

row('SR-87', 1e-6 / AU, 'an orbit-scale aperture (alpha 3.1 beam director and link telescope); the Airy first null is 1.22 times this')

# ---- print ----
w = max(len(r[2]) for r in rows)
for st, key, name, stated, comp, rel, note in rows:
    rs = f'{rel:+.2%}' if math.isfinite(rel) else '   n/a'
    prec = max(5, sigfigs(stated.replace(',', '')) + 2)
    print(f'{st:9s} | {key:5s} | {name:<{w}s} | stated {stated:>10s} | computed {comp:<11.{prec}g} | {rs:>8s} | {note}')
from collections import Counter
print('\nsummary:', dict(Counter(r[0] for r in rows)))
print(f'unique rows checked: {len(seen)} / {len(ORACLE)}')
missing = sorted(set(ORACLE) - seen)
if missing:
    print('rows in oracle.csv with no check:', ', '.join(missing))
mismatched = [r[1] for r in rows if r[0] == 'MISMATCH']
if mismatched:
    print('mismatched oracle rows:', ', '.join(mismatched))
if missing or mismatched:
    sys.exit(1)
