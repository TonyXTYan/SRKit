# The physics oracle

`oracle.csv` is a table of reference numbers: one row per quantity, with the value as stated in the design documents of *The Special Relativity Game*, its unit and a tolerance. Rows are named `SR-nn` (special relativity, optics, gravity and the engineering built on them) and `DT-nn` (deep time: cosmology, exact clocks, black holes).

Two independent things reproduce it:

- `oracle_check.py` computes every row from its formula with mpmath, NumPy and SciPy, and reports each as MATCH, ROUNDED, AMBIGUOUS or MISMATCH. It shares no code with the crates.

  ```sh
  python oracle/oracle_check.py        # needs numpy, scipy and mpmath
  ```

- The crates' `tests/oracle.rs` read the same file at build time and recompute their rows through the public API: `srkit-core` 47 rows, `srkit-flrw` the cosmology rows, `srkit-gr` the gravity rows. The Python binding's tests do the same for the core's rows. Rows that are game design (megastructures, power budgets) belong to no crate and are left to `oracle_check.py`.

`scripts/` holds the four reference scripts that crate tests cite by name and reproduce to the digits the scripts print (`srkit-flrw/tests/scripts.rs`, `srkit-gr/tests/oracle.rs`).

## Where rows come from

The game's repository keeps the working copy of the oracle (`tools/physics-checks/` there) and adds rows as its design grows; this folder is the copy SRKit is tested against. A row that exists in both must be identical in both, and the game's repository checks that. A stated value is never edited to make a test pass: a disagreement is settled by a derivation, and the row and both checks change together.
