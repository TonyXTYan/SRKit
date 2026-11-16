//! Helpers shared by the integration tests.
#![allow(dead_code)]

/// Split one CSV record into its fields.
///
/// Fields may be quoted, quoted fields may contain commas, and a doubled quote inside a quoted field is a
/// literal quote (RFC 4180). The oracle's quantities contain commas, so splitting on commas is not enough.
/// Records with embedded newlines are not supported; neither file has any.
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

/// `|got − want| / |want|`, or the absolute difference when `want` is zero.
pub fn relative_error(got: f64, want: f64) -> f64 {
    if want == 0.0 {
        (got - want).abs()
    } else {
        ((got - want) / want).abs()
    }
}

/// Assert that `got` is within the relative tolerance `tol` of `want`.
#[track_caller]
pub fn assert_close(got: f64, want: f64, tol: f64, what: &str) {
    let error = relative_error(got, want);
    assert!(
        error <= tol,
        "{what}: got {got:e}, want {want:e}, relative error {error:e} > {tol:e}"
    );
}
