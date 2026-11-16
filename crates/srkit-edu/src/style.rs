//! Colours, dash patterns and the two colour schemes.

/// The colour scheme of a diagram.
///
/// The palette of each theme is fixed and was checked by computation, not by eye. The four series colours and
/// the colour for light have a WCAG contrast of at least 3:1 against the surface. Every pair of them is at
/// least 9 apart (distance in OKLab, times 100) under simulated protanopia and deuteranopia, and at least 15
/// apart for full colour vision. Text is at least 7:1.
///
/// Colour is never the only cue: [`Style::series`] also gives each series its own dash pattern, light is
/// drawn dashed, and every labelled line is named in the legend.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Theme {
    /// Dark ink on a near-white surface.
    #[default]
    Light,
    /// Light ink on a near-black surface.
    Dark,
}

/// A stroke colour, named by its job so that it follows the [`Theme`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Colour {
    /// One of the four series colours, for worldlines: 0 blue, 1 magenta, 2 green, 3 purple. Larger indices
    /// wrap round, so tell a fifth series apart by its dash pattern.
    Series(u8),
    /// The colour for light: cones, signals and lines of sight. A dark gold.
    Light,
    /// The colour of the text.
    Ink,
    /// A mid grey, for guide lines that should stay in the background.
    Muted,
    /// A colour of the caller's choosing, the same in both themes. Its contrast is the caller's concern.
    Rgb(u8, u8, u8),
}

/// A dash pattern. The lengths are in pixels and do not grow with the stroke width.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Dash {
    /// An unbroken line.
    #[default]
    Solid,
    /// Long dashes: 12 on, 4 off.
    LongDash,
    /// Dashes: 6 on, 4 off.
    Dashed,
    /// A dash and a dot: 9 on, 3 off, 2 on, 3 off.
    DashDot,
    /// A dash and two dots: 9 on, 3 off, 2 on, 3 off, 2 on, 3 off.
    DashDotDot,
    /// Dots: 2 on, 3 off.
    Dotted,
}

impl Dash {
    /// The SVG `stroke-dasharray` value, or `None` for a solid line.
    pub const fn pattern(self) -> Option<&'static str> {
        match self {
            Self::Solid => None,
            Self::LongDash => Some("12 4"),
            Self::Dashed => Some("6 4"),
            Self::DashDot => Some("9 3 2 3"),
            Self::DashDotDot => Some("9 3 2 3 2 3"),
            Self::Dotted => Some("2 3"),
        }
    }
}

/// How a line or a mark is drawn: stroke colour, stroke width in pixels and dash pattern.
///
/// Dots (events and ticks) are filled with the colour; the width and the dash do not apply to them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Style {
    /// Stroke colour.
    pub colour: Colour,
    /// Stroke width in pixels: positive, finite and at most 64.
    pub width: f64,
    /// Dash pattern.
    pub dash: Dash,
}

/// The dash patterns [`Style::series`] hands out, in order.
const SERIES_DASHES: [Dash; 5] = [
    Dash::Solid,
    Dash::LongDash,
    Dash::DashDot,
    Dash::Dashed,
    Dash::DashDotDot,
];

impl Style {
    /// A style from its three parts.
    pub const fn new(colour: Colour, width: f64, dash: Dash) -> Self {
        Self {
            colour,
            width,
            dash,
        }
    }
    /// The style of the `index`-th worldline: 2 pixels wide, in series colour `index` (of four) and dash pattern
    /// `index` (of five), so twenty series differ in colour or in dash and neighbours differ in both.
    pub const fn series(index: u8) -> Self {
        Self::new(
            Colour::Series(index % 4),
            2.0,
            SERIES_DASHES[(index % 5) as usize],
        )
    }
    /// The style for light: thin dashes in the colour for light.
    pub const fn light() -> Self {
        Self::new(Colour::Light, 1.25, Dash::Dashed)
    }
    /// The style for guide lines, such as another frame's grid: thin grey dots.
    pub const fn guide() -> Self {
        Self::new(Colour::Muted, 1.0, Dash::Dotted)
    }
    /// The style for a horizon: a dash-dot line in the colour of the text.
    pub const fn horizon() -> Self {
        Self::new(Colour::Ink, 1.5, Dash::DashDot)
    }
    /// The same style in another colour.
    pub const fn with_colour(mut self, colour: Colour) -> Self {
        self.colour = colour;
        self
    }
    /// The same style at another stroke width, in pixels.
    pub const fn with_width(mut self, width: f64) -> Self {
        self.width = width;
        self
    }
    /// The same style with another dash pattern.
    pub const fn with_dash(mut self, dash: Dash) -> Self {
        self.dash = dash;
        self
    }
}

/// The fixed colours of a theme, as `[red, green, blue]`.
pub(crate) struct Palette {
    pub surface: [u8; 3],
    pub ink: [u8; 3],
    /// Axis titles, tick labels and the numbers beside proper-time ticks.
    pub secondary: [u8; 3],
    pub muted: [u8; 3],
    pub grid: [u8; 3],
    pub light: [u8; 3],
    pub series: [[u8; 3]; 4],
}

const LIGHT: Palette = Palette {
    surface: [0xfc, 0xfc, 0xfb],
    ink: [0x0b, 0x0b, 0x0b],
    secondary: [0x52, 0x51, 0x4e],
    muted: [0x89, 0x87, 0x81],
    grid: [0xe1, 0xe0, 0xd9],
    light: [0xc7, 0x83, 0x00],
    series: [
        [0x2a, 0x78, 0xd6],
        [0xd5, 0x51, 0x81],
        [0x00, 0x6d, 0x00],
        [0x6a, 0x28, 0xb0],
    ],
};

const DARK: Palette = Palette {
    surface: [0x1a, 0x1a, 0x19],
    ink: [0xf0, 0xef, 0xec],
    secondary: [0xc3, 0xc2, 0xb7],
    muted: [0x89, 0x87, 0x81],
    grid: [0x2c, 0x2c, 0x2a],
    light: [0xc9, 0x85, 0x00],
    series: [
        [0x39, 0x87, 0xe5],
        [0xd5, 0x51, 0x81],
        [0x42, 0x76, 0x37],
        [0x89, 0x55, 0xac],
    ],
};

impl Theme {
    pub(crate) const fn palette(self) -> &'static Palette {
        match self {
            Self::Light => &LIGHT,
            Self::Dark => &DARK,
        }
    }
    /// The `[red, green, blue]` this theme gives a colour.
    pub const fn rgb(self, colour: Colour) -> [u8; 3] {
        let p = self.palette();
        match colour {
            Colour::Series(i) => p.series[(i % 4) as usize],
            Colour::Light => p.light,
            Colour::Ink => p.ink,
            Colour::Muted => p.muted,
            Colour::Rgb(r, g, b) => [r, g, b],
        }
    }
    /// The `[red, green, blue]` of the background.
    pub const fn surface(self) -> [u8; 3] {
        self.palette().surface
    }
}

/// `#rrggbb`.
pub(crate) fn hex(rgb: [u8; 3]) -> String {
    format!("#{:02x}{:02x}{:02x}", rgb[0], rgb[1], rgb[2])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// WCAG 2 relative luminance.
    fn luminance(rgb: [u8; 3]) -> f64 {
        let lin = |c: u8| {
            let c = f64::from(c) / 255.0;
            if c <= 0.04045 {
                c / 12.92
            } else {
                ((c + 0.055) / 1.055).powf(2.4)
            }
        };
        0.2126 * lin(rgb[0]) + 0.7152 * lin(rgb[1]) + 0.0722 * lin(rgb[2])
    }

    fn contrast(a: [u8; 3], b: [u8; 3]) -> f64 {
        let (la, lb) = (luminance(a), luminance(b));
        (la.max(lb) + 0.05) / (la.min(lb) + 0.05)
    }

    #[test]
    fn every_mark_colour_clears_three_to_one_and_text_seven_to_one() {
        for theme in [Theme::Light, Theme::Dark] {
            let p = theme.palette();
            let marks = [
                Colour::Series(0),
                Colour::Series(1),
                Colour::Series(2),
                Colour::Series(3),
                Colour::Light,
                Colour::Muted,
                Colour::Ink,
            ];
            for colour in marks {
                let ratio = contrast(theme.rgb(colour), p.surface);
                assert!(ratio >= 3.0, "{theme:?} {colour:?}: {ratio:.2}");
            }
            for text in [p.ink, p.secondary] {
                let ratio = contrast(text, p.surface);
                assert!(ratio >= 7.0, "{theme:?} text {text:?}: {ratio:.2}");
            }
        }
    }

    /// OKLab coordinates of a linear-light colour.
    fn oklab(rgb: [f64; 3]) -> [f64; 3] {
        let [r, g, b] = rgb;
        let l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b).cbrt();
        let m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b).cbrt();
        let s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b).cbrt();
        [
            0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
            1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s,
        ]
    }

    /// Distance in OKLab, times 100, between two colours as seen through `vision`, a matrix on linear RGB.
    fn distance(a: [u8; 3], b: [u8; 3], vision: &[[f64; 3]; 3]) -> f64 {
        let seen =
            |rgb: [u8; 3]| {
                let lin = rgb.map(|c| {
                    let c = f64::from(c) / 255.0;
                    if c <= 0.04045 {
                        c / 12.92
                    } else {
                        ((c + 0.055) / 1.055).powf(2.4)
                    }
                });
                oklab(vision.map(|row| {
                    (row[0] * lin[0] + row[1] * lin[1] + row[2] * lin[2]).clamp(0.0, 1.0)
                }))
            };
        let (p, q) = (seen(a), seen(b));
        100.0 * ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
    }

    /// The claim in the documentation of [`Theme`], checked for every pair of the five line colours. Colour
    /// blindness is simulated with the matrices of Machado, Oliveira and Fernandes (2009) at full severity.
    #[test]
    fn line_colours_stay_apart_with_and_without_colour_blindness() {
        const NORMAL: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        const PROTAN: [[f64; 3]; 3] = [
            [0.152286, 1.052583, -0.204868],
            [0.114503, 0.786281, 0.099216],
            [-0.003882, -0.048116, 1.051998],
        ];
        const DEUTAN: [[f64; 3]; 3] = [
            [0.367322, 0.860646, -0.227968],
            [0.280085, 0.672501, 0.047413],
            [-0.011820, 0.042940, 0.968881],
        ];
        for theme in [Theme::Light, Theme::Dark] {
            let colours = [
                Colour::Series(0),
                Colour::Series(1),
                Colour::Series(2),
                Colour::Series(3),
                Colour::Light,
            ]
            .map(|c| theme.rgb(c));
            for (i, &a) in colours.iter().enumerate() {
                for &b in &colours[..i] {
                    let blind = distance(a, b, &PROTAN).min(distance(a, b, &DEUTAN));
                    let normal = distance(a, b, &NORMAL);
                    assert!(blind >= 9.0, "{theme:?} {a:?} {b:?}: {blind:.1}");
                    assert!(normal >= 15.0, "{theme:?} {a:?} {b:?}: {normal:.1}");
                }
            }
        }
    }

    #[test]
    fn neighbouring_series_differ_in_colour_and_in_dash() {
        for i in 0..20u8 {
            for j in 0..i {
                let (a, b) = (Style::series(i), Style::series(j));
                assert!(a.colour != b.colour || a.dash != b.dash, "{i} and {j}");
            }
        }
        for i in 0..19u8 {
            let (a, b) = (Style::series(i), Style::series(i + 1));
            assert!(a.colour != b.colour && a.dash != b.dash);
        }
    }

    #[test]
    fn colours_wrap_and_print_as_hex() {
        assert_eq!(
            Theme::Light.rgb(Colour::Series(5)),
            Theme::Light.rgb(Colour::Series(1))
        );
        assert_eq!(hex(Theme::Light.rgb(Colour::Rgb(1, 171, 255))), "#01abff");
        assert_eq!(hex(Theme::Dark.surface()), "#1a1a19");
    }
}
