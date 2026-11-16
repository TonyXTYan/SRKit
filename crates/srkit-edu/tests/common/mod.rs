//! A strict little XML reader for the tests: it checks that a document is well formed and hands back its
//! elements, so that the tests can read coordinates out of the SVG the crate wrote. No XML crate is used.
#![allow(dead_code)]

use srkit_core::{segment::Segment, worldline::WorldlineBuf, Vec3};

/// One element of a parsed document.
#[derive(Clone, Debug)]
pub struct Element {
    pub name: String,
    /// Attributes in document order, values unescaped.
    pub attrs: Vec<(String, String)>,
    /// The `class` of every enclosing element, outermost first, then this element's own.
    pub classes: Vec<String>,
    /// Names of the enclosing elements, outermost first.
    pub ancestors: Vec<String>,
    /// The text directly inside the element, unescaped.
    pub text: String,
}

impl Element {
    pub fn attr(&self, name: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v.as_str())
    }
    pub fn num(&self, name: &str) -> f64 {
        self.attr(name)
            .unwrap_or_else(|| panic!("<{}> has no {name}", self.name))
            .parse()
            .unwrap()
    }
    /// The `points` of a polyline or polygon, in pixels.
    pub fn points(&self) -> Vec<(f64, f64)> {
        self.attr("points")
            .unwrap()
            .split(' ')
            .map(|pair| {
                let (x, y) = pair.split_once(',').unwrap();
                (x.parse().unwrap(), y.parse().unwrap())
            })
            .collect()
    }
    /// The two ends of a `<line>`.
    pub fn ends(&self) -> ((f64, f64), (f64, f64)) {
        (
            (self.num("x1"), self.num("y1")),
            (self.num("x2"), self.num("y2")),
        )
    }
    pub fn centre(&self) -> (f64, f64) {
        (self.num("cx"), self.num("cy"))
    }
}

/// The elements called `name` inside a group of class `class`.
pub fn find<'a>(elements: &'a [Element], class: &str, name: &str) -> Vec<&'a Element> {
    elements
        .iter()
        .filter(|e| e.name == name && e.classes.iter().any(|c| c == class))
        .collect()
}

fn is_name(name: &str) -> bool {
    let mut chars = name.chars();
    chars
        .next()
        .is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
        && chars.all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.' | ':'))
}

/// Replace the entity references in `text`. A bare `&`, an unknown entity or a control character is an error.
fn unescape(text: &str) -> Result<String, String> {
    let mut out = String::new();
    let mut rest = text;
    while let Some(amp) = rest.find('&') {
        out.push_str(&rest[..amp]);
        let end = rest[amp..]
            .find(';')
            .ok_or_else(|| format!("bare & in {text:?}"))?;
        let entity = &rest[amp + 1..amp + end];
        out.push(match entity {
            "amp" => '&',
            "lt" => '<',
            "gt" => '>',
            "quot" => '"',
            "apos" => '\'',
            _ => {
                let code = entity
                    .strip_prefix("#x")
                    .map(|hex| u32::from_str_radix(hex, 16))
                    .or_else(|| {
                        entity
                            .strip_prefix('#')
                            .map(|decimal| decimal.parse::<u32>())
                    })
                    .ok_or_else(|| format!("unknown entity &{entity};"))?
                    .map_err(|_| format!("bad character reference &{entity};"))?;
                char::from_u32(code).ok_or_else(|| format!("bad character reference &{entity};"))?
            }
        });
        rest = &rest[amp + end + 1..];
    }
    out.push_str(rest);
    if let Some(c) = out
        .chars()
        .find(|&c| (c as u32) < 0x20 && !"\t\n\r".contains(c))
    {
        return Err(format!("control character {:#x} in {text:?}", c as u32));
    }
    Ok(out)
}

fn attributes(mut rest: &str) -> Result<Vec<(String, String)>, String> {
    let mut out: Vec<(String, String)> = Vec::new();
    while !rest.trim_start().is_empty() {
        if !rest.starts_with(char::is_whitespace) {
            return Err(format!("no space before an attribute at {rest:?}"));
        }
        rest = rest.trim_start();
        let eq = rest
            .find('=')
            .ok_or_else(|| format!("attribute without a value at {rest:?}"))?;
        let name = &rest[..eq];
        if !is_name(name) {
            return Err(format!("bad attribute name {name:?}"));
        }
        let value = rest[eq + 1..]
            .strip_prefix('"')
            .ok_or_else(|| format!("{name} is not quoted"))?;
        let close = value
            .find('"')
            .ok_or_else(|| format!("{name} has no closing quote"))?;
        if out.iter().any(|(n, _)| n == name) {
            return Err(format!("{name} is given twice"));
        }
        out.push((name.to_string(), unescape(&value[..close])?));
        rest = &value[close + 1..];
    }
    Ok(out)
}

/// Parse a document, checking that it is well formed: one root, every tag closed in order, attributes quoted
/// and given once, and every `&` and `<` escaped.
pub fn parse(document: &str) -> Result<Vec<Element>, String> {
    let mut out: Vec<Element> = Vec::new();
    let mut open: Vec<usize> = Vec::new();
    let mut roots = 0;
    let mut rest = document;
    loop {
        let lt = rest.find('<');
        let text = &rest[..lt.unwrap_or(rest.len())];
        match open.last() {
            Some(&top) => out[top].text.push_str(&unescape(text)?),
            None if text.trim().is_empty() => {}
            None => return Err(format!("text outside the root: {text:?}")),
        }
        let Some(lt) = lt else { break };
        rest = &rest[lt + 1..];
        let gt = rest.find('>').ok_or("a tag is never closed")?;
        let tag = &rest[..gt];
        rest = &rest[gt + 1..];
        if tag.contains('<') {
            return Err(format!("< inside a tag: {tag:?}"));
        }
        if let Some(name) = tag.strip_prefix('/') {
            let top = open
                .pop()
                .ok_or_else(|| format!("</{name}> closes nothing"))?;
            if out[top].name != name {
                return Err(format!("</{name}> closes <{}>", out[top].name));
            }
            continue;
        }
        let (body, empty) = match tag.strip_suffix('/') {
            Some(body) => (body, true),
            None => (tag, false),
        };
        let name_end = body.find(char::is_whitespace).unwrap_or(body.len());
        let name = &body[..name_end];
        if !is_name(name) {
            return Err(format!("bad element name {name:?}"));
        }
        let attrs = attributes(&body[name_end..])?;
        if open.is_empty() {
            roots += 1;
        }
        let (mut classes, mut ancestors) = match open.last() {
            Some(&top) => (out[top].classes.clone(), out[top].ancestors.clone()),
            None => (Vec::new(), Vec::new()),
        };
        if let Some(&top) = open.last() {
            ancestors.push(out[top].name.clone());
        }
        if let Some((_, class)) = attrs.iter().find(|(n, _)| n == "class") {
            classes.push(class.clone());
        }
        out.push(Element {
            name: name.to_string(),
            attrs,
            classes,
            ancestors,
            text: String::new(),
        });
        if !empty {
            open.push(out.len() - 1);
        }
    }
    if let Some(&top) = open.last() {
        return Err(format!("<{}> is never closed", out[top].name));
    }
    if roots != 1 {
        return Err(format!("{roots} root elements"));
    }
    Ok(out)
}

/// Distance from a point to a segment.
pub fn distance_to_segment(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (bx, by) = (b.0 - a.0, b.1 - a.1);
    let (px, py) = (p.0 - a.0, p.1 - a.1);
    let length2 = bx * bx + by * by;
    let u = if length2 > 0.0 {
        ((px * bx + py * by) / length2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (px - u * bx).hypot(py - u * by)
}

/// Distance from a point to the nearest of several polylines.
pub fn distance_to_polylines(p: (f64, f64), polylines: &[Vec<(f64, f64)>]) -> f64 {
    polylines
        .iter()
        .flat_map(|line| line.windows(2))
        .map(|pair| distance_to_segment(p, pair[0], pair[1]))
        .fold(f64::INFINITY, f64::min)
}

/// A body at rest at `x` on the axis from `t0` to `t1`.
pub fn at_rest(x: f64, t0: f64, t1: f64) -> WorldlineBuf {
    WorldlineBuf::new([Segment::rest(t0, t1, Vec3::new(x, 0.0, 0.0))], 0.0, 0.0).unwrap()
}

/// A body moving at `v` along the axis from `t0` to `t1`, passing the origin at `t = 0` with its clock at zero.
pub fn moving(v: f64, t0: f64, t1: f64) -> WorldlineBuf {
    WorldlineBuf::new(
        [Segment::inertial(
            t0,
            t1,
            Vec3::new(v * t0, 0.0, 0.0),
            Vec3::new(v, 0.0, 0.0),
        )],
        0.0,
        0.0,
    )
    .unwrap()
}

/// A single segment as a worldline with its clock at zero at the segment's start.
pub fn single(segment: Segment) -> WorldlineBuf {
    WorldlineBuf::new([segment], segment.t0(), 0.0).unwrap()
}
