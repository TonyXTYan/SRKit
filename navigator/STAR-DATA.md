# The Navigator's star table: sources and licence

`src/stars.json`, and the copy of it embedded in `relativistic-navigator.html`, is **not** under the repository's MIT licence. It is a database derived from other databases and carries their terms.

**Licence: Open Database License (ODbL) 1.0**, <https://opendatacommons.org/licenses/odbl/1-0/>.

> Contains information from the SIMBAD database, operated at CDS, Strasbourg, France, which is made available under the ODbL.

| Source | What comes from it | Terms |
|---|---|---|
| [SIMBAD](https://simbad.cds.unistra.fr/simbad/) (CDS, Strasbourg; Wenger et al. 2000, A&AS 143, 9) | Positions, parallaxes, magnitudes, colours, spectral types, designations | ODbL 1.0. From the SIMBAD site: "Use of SIMBAD service is free without any registration under ODbL licence" |
| [IAU Catalog of Star Names](https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt) (IAU Working Group on Star Names) | Proper names | Creative Commons Attribution |

What the ODbL asks of anyone who passes the table on, in short (the licence text governs):

- **Attribution.** Keep the notice above with the table and with anything made from it that is shown publicly.
- **Share-alike.** A changed or extended version of the table that is used publicly is offered under the ODbL too.
- **Keep it open.** If you publish the table with technical restrictions, also offer a copy without them.

The ODbL permits commercial use. A work made *with* the table (a rendered sky, a chart, an application that reads it) is a "produced work": it carries the attribution notice and need not itself be under the ODbL.

## How it is made

`src/build_stars.py` holds the whole method, which meets the ODbL's condition that a derived database is offered with its alterations or the means to make them. It selects every star SIMBAD lists with V ≤ 6.5 and a parallax, and everything with a parallax of 65.2 mas or more (within 50 light-years); converts to galactic coordinates and light-years; takes B − V from SIMBAD where it has both magnitudes and from the spectral class otherwise; removes clusters, planets and system entries that duplicate their components; and names stars from the IAU catalogue, then by Bayer or Flamsteed designation. `mass_est` is not data: it is a rough main-sequence fit from luminosity, used only to decide which labels show first.

Downloaded 16 November 2026. 9,715 stars.

## Why SIMBAD

The first table came from the HYG database (CC BY-SA 4.0), which is compiled from the Hipparcos, Yale and Gliese catalogues. ESA distributes Hipparcos and Gaia data under CC BY-NC 3.0 IGO, with commercial use by prior authorisation, so a table that descends from them was not clearly free to carry into an MIT-licensed project or a commercial product. SIMBAD publishes its own compilation under the ODbL, which is the licence this table relies on. Individual measurements in SIMBAD come from the literature, the Hipparcos and Gaia catalogues among it; the terms above are SIMBAD's for its database, and this note is not legal advice.
