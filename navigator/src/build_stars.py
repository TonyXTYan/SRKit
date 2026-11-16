"""Build the Navigator's star list from SIMBAD and the IAU Catalog of Star Names.

    python3 -I build_stars.py fetch <dir>             download the three source files into <dir>
    python3 -I build_stars.py build <dir> <out.json>  build the star list from them

Selection: every star SIMBAD lists with V <= 6.5 and a parallax (the naked-eye sky), and everything with a parallax
of 65.2 mas or more (within 50 light-years) that has a magnitude in V, G or R.

Sources and licences. The output is a derived database and carries the licences of its sources, not the code's:
  * SIMBAD (CDS, Strasbourg; Wenger et al. 2000, A&AS 143, 9), under the Open Database License (ODbL) 1.0:
    positions, parallaxes, magnitudes, spectral types and designations. "Use of SIMBAD service is free without any
    registration under ODbL licence" (simbad.cds.unistra.fr). The star list is therefore under ODbL 1.0.
  * IAU Catalog of Star Names (IAU WGSN), Creative Commons Attribution: proper names.
The earlier list came from the HYG database (CC BY-SA 4.0), which is compiled from Hipparcos; ESA distributes
Hipparcos and Gaia data under CC BY-NC 3.0 IGO. SIMBAD is the publisher whose licence this list relies on.

Masses are a placeholder main-sequence fit from luminosity, used only to decide which labels show first.
Code AI-generated (Claude Code)."""
import csv
import json
import math
import re
import subprocess
import sys

# SIMBAD's query service at CDS, then its mirror at Harvard, tried in turn
TAP = ('https://simbad.cds.unistra.fr/simbad/sim-tap/sync', 'https://simbad.harvard.edu/simbad/sim-tap/sync')
COLUMNS = ('b.oid, b.main_id, b.ra, b.dec, b.plx_value, b.plx_err, b.plx_bibcode, b.sp_type, b.otype, f.V, f.B, f.G, f.R, i.ids '
           'FROM basic AS b {join} allfluxes AS f ON f.oidref = b.oid JOIN ids AS i ON i.oidref = b.oid')
QUERIES = {'bright.csv': 'SELECT ' + COLUMNS.format(join='JOIN') + ' WHERE f.V <= 6.5 AND b.plx_value > 0',
           'near.csv': 'SELECT ' + COLUMNS.format(join='LEFT JOIN') + ' WHERE b.plx_value >= 65.2'}
IAU_CSN = 'https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt'

GREEK = {'alf': 'α', 'bet': 'β', 'gam': 'γ', 'del': 'δ', 'eps': 'ε', 'zet': 'ζ', 'eta': 'η', 'tet': 'θ', 'iot': 'ι',
         'kap': 'κ', 'lam': 'λ', 'mu.': 'μ', 'nu.': 'ν', 'ksi': 'ξ', 'omi': 'ο', 'pi.': 'π', 'rho': 'ρ', 'sig': 'σ',
         'tau': 'τ', 'ups': 'υ', 'phi': 'φ', 'chi': 'χ', 'psi': 'ψ', 'ome': 'ω'}
# not single stars or systems of stars: clusters, planets, and SIMBAD's marker for a doubtful entry
NOT_STARS = {'OpC', 'GlC', 'Cl*', 'As*', 'Pl', 'Pl?', 'err'}
# J2000 equatorial -> galactic rotation (Hipparcos convention)
R = [[-0.0548755604, -0.8734370902, -0.4838350155],
     [+0.4941094279, -0.4448296300, +0.7469822445],
     [-0.8676661490, -0.1980763734, +0.4559837762]]
PC_LY = 3.261563777
# B - V where SIMBAD has no B magnitude, by the first letter of the spectral type
CI_BY_CLASS = {'O': -0.3, 'B': -0.15, 'A': 0.1, 'F': 0.45, 'G': 0.65, 'K': 1.0, 'M': 1.5, 'L': 2.0, 'T': 2.0, 'Y': 2.0, 'D': 0.2, 'C': 1.6, 'S': 1.6}


def galactic(ra_deg, dec_deg):
    a, d = math.radians(ra_deg), math.radians(dec_deg)
    v = (math.cos(d) * math.cos(a), math.cos(d) * math.sin(a), math.sin(d))
    g = [sum(R[i][j] * v[j] for j in range(3)) for i in range(3)]
    return math.degrees(math.atan2(g[1], g[0])) % 360, math.degrees(math.asin(max(-1, min(1, g[2]))))


def mass_from_lum(L):
    """Placeholder priority: main-sequence mass-luminosity fit (overestimates giants)."""
    if L <= 0:
        return 0.1
    m = L ** 0.25 if L < 16 else (L / 1.4) ** (1 / 3.5)
    return max(0.08, min(60.0, m))


def tidy(s):
    return re.sub(r'\s+', ' ', s).strip()


def designation(ids):
    """A Bayer or Flamsteed designation from SIMBAD's '* alf CMa', '* tet01 Ori C' or '*  61 Cyg B'."""
    bayer = flamsteed = ''
    for i in ids:
        m = re.match(r'^\*\s+([a-z]{2}[a-z.])(\d\d)?\s+(\S+)(\s+\S+)?$', i)
        if m and m.group(1) in GREEK and not bayer:
            bayer = GREEK[m.group(1)] + (str(int(m.group(2))) if m.group(2) else '') + ' ' + m.group(3) + (m.group(4) or '')
        m = re.match(r'^\*\s+(\d+)\s+(\S+)(\s+\S+)?$', i)
        if m and not flamsteed:
            flamsteed = m.group(1) + ' ' + m.group(2) + (m.group(3) or '')
    return tidy(bayer or flamsteed)


def ident(ids, main_id):
    for prefix in ('HIP', 'HD', 'GJ', 'HR'):
        for i in ids:
            if re.match(rf'^{prefix}\s+\S', i):
                return tidy(i)
    return tidy(main_id.lstrip('*').replace('NAME ', ''))


def fetch(folder):
    for name, query in QUERIES.items():
        for attempt, tap in enumerate(TAP):
            done = subprocess.run(['curl', '-sS', '--fail', '-m', '900', '-G', tap, '--data-urlencode', 'request=doQuery',
                                   '--data-urlencode', 'lang=adql', '--data-urlencode', 'format=csv',
                                   '--data-urlencode', 'query=' + query, '-o', f'{folder}/{name}'])
            # an unavailable service can answer 200 with an HTML page: accept only the CSV header
            if done.returncode == 0 and open(f'{folder}/{name}', 'rb').read(3) == b'oid':
                break
            if attempt == len(TAP) - 1:
                sys.exit(f'{name}: no SIMBAD service answered')
    subprocess.run(['curl', '-sS', '--fail', '-m', '120', '-o', f'{folder}/iau-csn.txt', IAU_CSN], check=True)


def build(folder, out):
    # proper names, by Hipparcos and Henry Draper number (columns 11 and 12 of the IAU list)
    names = {}
    for line in open(f'{folder}/iau-csn.txt', encoding='utf-8'):
        if line[:1] in '#$' or not line.strip():
            continue
        # fixed-width names, then columns that are counted back from the date (the designation can hold a space)
        name, f = line[18:36].strip(), line[36:].split()
        k = next(i for i, tok in enumerate(f) if re.fullmatch(r'\d{4}-\d{2}-\d{2}', tok))
        hip, hd = f[k - 4], f[k - 3]
        if hip != '_':
            names['HIP ' + hip] = name
        if hd != '_':
            names.setdefault('HD ' + hd, name)

    objects = {}
    for name in ('bright.csv', 'near.csv'):
        with open(f'{folder}/{name}', newline='', encoding='utf-8') as fh:
            for r in csv.DictReader(fh):
                objects[r['oid']] = r
    stars, skipped = [], {'not a star': 0, 'no magnitude': 0}
    for r in objects.values():
        if r['otype'] in NOT_STARS:
            skipped['not a star'] += 1
            continue
        mag = next((float(r[k]) for k in ('V', 'G', 'R') if r[k]), None)
        if mag is None:
            skipped['no magnitude'] += 1
            continue
        plx = float(r['plx_value'])
        d_pc = 1000.0 / plx
        if not (mag <= 6.5 or d_pc * PC_LY <= 50):
            continue
        ids = [tidy(i) for i in r['ids'].split('|')]
        sp = tidy(r['sp_type'])
        ci = float(r['B']) - float(r['V']) if r['B'] and r['V'] else CI_BY_CLASS.get(sp[:1].upper(), 0.65)
        lum = 10 ** ((4.83 - (mag - 5 * math.log10(d_pc / 10))) / 2.5)
        proper = next((names[i] for i in ids if i in names), '')
        l, b = galactic(float(r['ra']), float(r['dec']))
        stars.append({'row': [round(l, 2), round(b, 2), float(f'{d_pc * PC_LY:.4g}'), round(mag, 2), round(max(-0.4, min(2.5, ci)), 2),
                              float(f'{mass_from_lum(lum):.2g}'), proper or designation(ids), ident(ids, r['main_id']), sp[:12]],
                      'system': r['otype'] == '**', 'ra': float(r['ra']), 'dec': float(r['dec']), 'plx': plx, 'main': tidy(r['main_id'])})

    # SIMBAD lists a multiple system and its components as separate objects. Where a system's entry has a selected
    # component within three arcminutes (positions of fast nearby stars differ by epoch) at the same distance, keep the
    # components and drop the system.
    cells = {}
    for s in stars:
        cells.setdefault((round(s['ra'] * 20), round(s['dec'] * 20)), []).append(s)

    def has_component(s):
        cx, cy = round(s['ra'] * 20), round(s['dec'] * 20)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for o in cells.get((cx + dx, cy + dy), ()):
                    if o is s or o['system']:
                        continue
                    sep = math.hypot((o['ra'] - s['ra']) * math.cos(math.radians(s['dec'])), o['dec'] - s['dec'])
                    if sep < 3 / 60 and abs(o['plx'] / s['plx'] - 1) < 0.2:
                        return True
        return False

    # The same for a system that SIMBAD types otherwise (alpha Centauri is a spectroscopic binary there): if its
    # primary is listed as "<name> A", the system's own entry goes. Sirius stays: only its companion is listed apart.
    mains = {s['main'] for s in stars}
    kept = [s for s in stars if not (s['system'] and has_component(s)) and s['main'] + ' A' not in mains]
    rows = sorted((s['row'] for s in kept), key=lambda x: (-x[5], x[7]))
    json.dump({'fields': ['l', 'b', 'd_ly', 'mag', 'ci', 'mass_est', 'name', 'id', 'spect'], 'stars': rows,
               'licence': 'Open Database License (ODbL) 1.0. Contains information from the SIMBAD database, operated at CDS, Strasbourg, France '
                          '(Wenger et al. 2000), made available under the ODbL. Proper names: IAU Catalog of Star Names (IAU WGSN), CC BY.'},
              open(out, 'w', encoding='utf-8'), ensure_ascii=False, separators=(',', ':'))
    named = sum(1 for x in rows if x[6])
    print(f'{len(rows)} stars ({len(stars) - len(kept)} system entries dropped for their components; skipped: {skipped}); {named} named; '
          f'nearest {min(x[2] for x in rows)} ly; farthest {max(x[2] for x in rows)} ly')
    for x in rows:
        if x[6] in ('Proxima Centauri', 'Sirius', 'Vega', 'Rigil Kentaurus', 'Toliman', 'Deneb', 'Betelgeuse', "Barnard's Star", 'τ Cet', 'Arcturus', '61 Cyg A'):
            print('  ', x)


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == 'fetch':
        fetch(sys.argv[2])
    elif len(sys.argv) == 4 and sys.argv[1] == 'build':
        build(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
