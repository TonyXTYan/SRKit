//! Helpers shared by the integration tests.
#![allow(dead_code)]

/// Split one CSV record into fields. Fields may be quoted; a quoted field may hold commas and doubled quotes.
/// Records in this repository never contain line breaks.
pub fn csv_fields(record: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = record.chars().peekable();
    while let Some(ch) = chars.next() {
        match ch {
            '"' if quoted && chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            ',' if !quoted => fields.push(std::mem::take(&mut field)),
            _ => field.push(ch),
        }
    }
    assert!(!quoted, "unterminated quoted field in: {record}");
    fields.push(field);
    fields
}

/// Relative difference `|got − expected| / |expected|`, or the absolute difference when `expected` is zero.
pub fn relative_error(got: f64, expected: f64) -> f64 {
    if got == expected {
        0.0
    } else if expected == 0.0 {
        got.abs()
    } else {
        ((got - expected) / expected).abs()
    }
}

#[test]
fn csv_fields_handle_quotes() {
    assert_eq!(csv_fields("a,b,,d"), ["a", "b", "", "d"]);
    assert_eq!(
        csv_fields(r#"SR-06,spine,"photon-rocket mass ratio, that trip",1.066e10,,0.01"#),
        [
            "SR-06",
            "spine",
            "photon-rocket mass ratio, that trip",
            "1.066e10",
            "",
            "0.01"
        ]
    );
    assert_eq!(csv_fields(r#""say ""hi""",2"#), [r#"say "hi""#, "2"]);
}
