//! A small SVG writer: escaped text, numbers with a fixed number of decimals, and a width estimate for text.

use std::fmt::Write as _;

/// Decimals written for every pixel coordinate and length.
pub(crate) const DECIMALS: usize = 2;

/// `value` with exactly `decimals` digits after the point. A value that rounds to zero is written without a
/// sign, so `-0.001` and `0.001` give the same text.
pub(crate) fn fixed(value: f64, decimals: usize) -> String {
    let text = format!("{value:.decimals$}");
    match text.strip_prefix('-') {
        Some(rest) if rest.bytes().all(|b| b == b'0' || b == b'.') => rest.to_string(),
        _ => text,
    }
}

/// Append `text` with the five XML special characters escaped. Tabs and line breaks become spaces, and
/// characters XML 1.0 does not allow become U+FFFD.
pub(crate) fn escape_into(out: &mut String, text: &str) {
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            '\t' | '\n' | '\r' => out.push(' '),
            c if (c as u32) < 0x20 || c == '\u{fffe}' || c == '\u{ffff}' => out.push('\u{fffd}'),
            c => out.push(c),
        }
    }
}

/// An SVG document under construction.
pub(crate) struct Svg {
    out: String,
}

impl Svg {
    pub fn new() -> Self {
        Self {
            out: String::with_capacity(16 * 1024),
        }
    }
    pub fn finish(self) -> String {
        self.out
    }
    /// Start a tag: `<name`. Follow with attributes, then [`Svg::children`], [`Svg::content`] or [`Svg::empty`].
    pub fn open(&mut self, name: &str) -> &mut Self {
        self.out.push('<');
        self.out.push_str(name);
        self
    }
    /// An attribute with an escaped text value.
    pub fn attr(&mut self, name: &str, value: &str) -> &mut Self {
        self.out.push(' ');
        self.out.push_str(name);
        self.out.push_str("=\"");
        escape_into(&mut self.out, value);
        self.out.push('"');
        self
    }
    /// An attribute with a number, written with [`DECIMALS`] decimals.
    pub fn num(&mut self, name: &str, value: f64) -> &mut Self {
        let _ = write!(self.out, " {name}=\"{}\"", fixed(value, DECIMALS));
        self
    }
    /// A `points` attribute from pixel positions.
    pub fn points(&mut self, points: impl IntoIterator<Item = (f64, f64)>) -> &mut Self {
        self.out.push_str(" points=\"");
        for (i, (x, y)) in points.into_iter().enumerate() {
            if i > 0 {
                self.out.push(' ');
            }
            let _ = write!(self.out, "{},{}", fixed(x, DECIMALS), fixed(y, DECIMALS));
        }
        self.out.push('"');
        self
    }
    /// End the start tag of an element that holds other elements.
    pub fn children(&mut self) {
        self.out.push_str(">\n");
    }
    /// End the start tag of an element that holds text.
    pub fn content(&mut self) -> &mut Self {
        self.out.push('>');
        self
    }
    /// End an element that holds nothing.
    pub fn empty(&mut self) {
        self.out.push_str("/>\n");
    }
    /// Escaped text.
    pub fn text(&mut self, text: &str) -> &mut Self {
        escape_into(&mut self.out, text);
        self
    }
    /// The end tag `</name>`.
    pub fn close(&mut self, name: &str) {
        self.out.push_str("</");
        self.out.push_str(name);
        self.out.push_str(">\n");
    }
    /// The end tag of an element inside a line of text: no line break after it.
    pub fn close_inline(&mut self, name: &str) -> &mut Self {
        self.out.push_str("</");
        self.out.push_str(name);
        self.out.push('>');
        self
    }
}

/// An estimate of the width of `text` in pixels at a font size, for laying out margins and the legend.
///
/// No font is measured (the output names only system fonts), so this leans wide: the widths are those of a
/// typical sans-serif, rounded up.
pub(crate) fn text_width(text: &str, size: f64) -> f64 {
    let ems: f64 = text
        .chars()
        .map(|c| match c {
            'i' | 'j' | 'l' | 'I' | '!' | '|' | '\'' | '.' | ',' | ':' | ';' => 0.30,
            ' ' | 'f' | 't' | 'r' | '(' | ')' | '[' | ']' | '-' | '/' | '"' => 0.36,
            'm' | 'w' | 'M' | 'W' | '@' | '%' => 0.90,
            'A'..='Z' => 0.70,
            '0'..='9' => 0.58,
            'a'..='z' => 0.57,
            c if (c as u32) < 0x2e80 => 0.62,
            _ => 1.0,
        })
        .sum();
    ems * size
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_decimals_and_no_negative_zero() {
        assert_eq!(fixed(1.0, 2), "1.00");
        assert_eq!(fixed(-1.005, 2), "-1.00");
        assert_eq!(fixed(-0.001, 2), "0.00");
        assert_eq!(fixed(-0.0, 2), "0.00");
        assert_eq!(fixed(0.0, 0), "0");
        assert_eq!(fixed(-0.4, 0), "0");
        assert_eq!(fixed(12345.678, 2), "12345.68");
        assert_eq!(fixed(-0.006, 2), "-0.01");
    }

    #[test]
    fn escapes_markup_and_control_characters() {
        let mut out = String::new();
        escape_into(&mut out, "a<b & \"c\" > 'd'\n\u{1}");
        assert_eq!(
            out,
            "a&lt;b &amp; &quot;c&quot; &gt; &apos;d&apos; \u{fffd}"
        );
    }

    #[test]
    fn writes_elements() {
        let mut svg = Svg::new();
        svg.open("g").attr("class", "a\"b").children();
        svg.open("line").num("x1", 1.0).num("y1", -0.0).empty();
        svg.open("polyline")
            .points([(0.0, 1.5), (2.25, 3.0)])
            .empty();
        svg.open("text").content().text("x<y").close("text");
        svg.close("g");
        assert_eq!(
            svg.finish(),
            "<g class=\"a&quot;b\">\n<line x1=\"1.00\" y1=\"0.00\"/>\n\
             <polyline points=\"0.00,1.50 2.25,3.00\"/>\n<text>x&lt;y</text>\n</g>\n"
        );
    }
}
