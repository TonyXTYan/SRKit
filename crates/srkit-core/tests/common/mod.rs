//! Shared by the oracle tests: the physics oracle table and its comparison rule.
#![allow(dead_code)]
use std::collections::BTreeMap;

pub const ORACLE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../../tools/physics-checks/oracle.csv"
));

/// One oracle row: its stated value as written (the number of significant figures matters) and its tolerance.
pub struct Row {
    pub quantity: String,
    pub stated: String,
    pub tolerance: f64,
}

/// Split one CSV record. Quantities are quoted where they contain commas; doubled quotes are accepted.
fn fields(record: &str) -> Vec<String> {
    let (mut values, mut value, mut quoted) = (Vec::new(), String::new(), false);
    let mut chars = record.chars().peekable();
    while let Some(ch) = chars.next() {
        match ch {
            '"' if quoted && chars.peek() == Some(&'"') => {
                value.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            ',' if !quoted => values.push(std::mem::take(&mut value)),
            _ => value.push(ch),
        }
    }
    assert!(!quoted, "unterminated CSV field in {record}");
    values.push(value);
    values
}

pub fn rows() -> BTreeMap<String, Row> {
    let mut out = BTreeMap::new();
    for record in ORACLE.lines().skip(1).filter(|l| !l.is_empty()) {
        let c = fields(record);
        assert_eq!(c.len(), 6, "invalid oracle record: {record}");
        let row = Row {
            quantity: c[2].clone(),
            stated: c[3].clone(),
            tolerance: c[5].parse().unwrap(),
        };
        assert!(
            out.insert(c[0].clone(), row).is_none(),
            "duplicate oracle id {}",
            c[0]
        );
    }
    out
}

fn sig_figs(stated: &str) -> usize {
    let mantissa = stated.to_lowercase();
    let mantissa = mantissa
        .trim_start_matches(['+', '-'])
        .split('e')
        .next()
        .unwrap()
        .to_string();
    let mut digits = mantissa
        .replace('.', "")
        .trim_start_matches('0')
        .to_string();
    if !mantissa.contains('.') {
        digits = digits.trim_end_matches('0').to_string();
    }
    digits.len().max(1)
}

/// The oracle's rule (as `oracle_check.py`): within the row's tolerance, or equal to the stated value once
/// rounded to the significant figures it is stated to.
pub fn agrees(row: &Row, computed: f64) -> bool {
    let stated: f64 = row.stated.replace(',', "").parse().unwrap();
    if stated == 0.0 {
        return computed == 0.0;
    }
    if ((computed - stated) / stated).abs() <= row.tolerance {
        return true;
    }
    let rounded: f64 = format!(
        "{:.*e}",
        sig_figs(&row.stated.replace(',', "")) - 1,
        computed
    )
    .parse()
    .unwrap();
    ((rounded - stated) / stated).abs() <= 1e-9
}

/// Check every `(id, value)` against the table and return how many rows the table holds in all.
pub fn check(computed: &[(&str, f64)]) -> usize {
    let table = rows();
    let mut seen = std::collections::BTreeSet::new();
    let mut failures = Vec::new();
    for (id, value) in computed {
        assert!(seen.insert(*id), "oracle row {id} checked twice");
        let row = table
            .get(*id)
            .unwrap_or_else(|| panic!("no oracle row {id}"));
        if !value.is_finite() || !agrees(row, *value) {
            failures.push(format!(
                "{id} ({}): stated {}, computed {value:e}",
                row.quantity, row.stated
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "oracle mismatches:\n{}",
        failures.join("\n")
    );
    table.len()
}
