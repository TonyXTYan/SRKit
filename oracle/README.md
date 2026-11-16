# The physics oracle

`oracle.csv` is the table of reference numbers for SRKit and for *The Special Relativity Game*: one row per quantity, with the value as stated in the game's design documents, its unit and a tolerance. Rows are named `SR-nn` (special relativity, optics, navigation, gravity and the engineering built on them) and `DT-nn` (deep time: cosmology, exact clocks, black holes). This folder is its home: a new row is added here, with its check.

Two independent things reproduce it:

- `oracle_check.py` computes every row from its formula with mpmath, NumPy and SciPy, and reports each as MATCH, ROUNDED, AMBIGUOUS or MISMATCH. It shares no code with the crates, and fails if a row has no check.

  ```sh
  python oracle/oracle_check.py        # needs numpy, scipy and mpmath
  ```

- The crates' `tests/oracle.rs` read the same file at build time and recompute their rows through the public API: `srkit-core` 47 rows, `srkit-flrw` the cosmology rows, `srkit-gr` the gravity rows. The Python binding's tests do the same for the core's rows. Rows that are game design (megastructures, power budgets) and the navigation rows (SR-74 to SR-79 and SR-88 to SR-110, which the Navigator's page quotes) belong to no crate yet and are left to `oracle_check.py`.

`scripts/` holds the four reference scripts that crate tests cite by name and reproduce to the digits the scripts print (`srkit-flrw/tests/scripts.rs`, `srkit-gr/tests/oracle.rs`). The game's repository keeps the scripts behind its design numbers (`tools/physics-checks/` there).

## Rules

- **Columns:** `id`, `group`, `quantity`, `stated`, `unit`, `tolerance` (relative; 0.01 means ±1%).
- A stated value carries at least three significant figures, so that it can meet a 1% tolerance.
- Constants are the official ones: exact `c`, the Julian year, CODATA `G`, the IAU nominal `GM☉`. The cosmology is Planck 2018 (`H0` = 67.4 km/s/Mpc, `Ω_m` = 0.315, flat).
- A stated value is never edited to make a test pass. A disagreement is settled by a derivation, and the row and both checks change together.
- IDs are never reused. The game's documents cite rows by ID.
