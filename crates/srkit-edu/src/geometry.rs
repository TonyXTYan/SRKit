//! Points and clipping in the coordinates of the drawing frame.

/// A point of the diagram: position `x` along the axis and time `t`, in the frame the diagram is drawn in.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Pt {
    pub x: f64,
    pub t: f64,
}

impl Pt {
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.t.is_finite()
    }
}

/// The view window: `x0..=x1` across and `t0..=t1` up.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Rect {
    pub x0: f64,
    pub x1: f64,
    pub t0: f64,
    pub t1: f64,
}

impl Rect {
    pub fn contains(&self, p: Pt) -> bool {
        p.x >= self.x0 && p.x <= self.x1 && p.t >= self.t0 && p.t <= self.t1
    }
    /// The point nearest `p` inside the window. Clipped ends are passed through this so that rounding cannot
    /// leave one a hair outside.
    fn clamp(&self, p: Pt) -> Pt {
        Pt {
            x: p.x.max(self.x0).min(self.x1),
            t: p.t.max(self.t0).min(self.t1),
        }
    }
    /// The part of `p + u d`, for `u` in `u0..=u1`, that lies inside the window (Liang and Barsky's method).
    /// The bounds may be infinite.
    fn clip(&self, p: Pt, d: Pt, mut u0: f64, mut u1: f64) -> Option<(f64, f64)> {
        if !p.is_finite() || !d.is_finite() {
            return None;
        }
        if d.x == 0.0 && d.t == 0.0 {
            return (self.contains(p) && u0 <= u1).then_some((u0, u1));
        }
        for (pos, dir, lo, hi) in [(p.x, d.x, self.x0, self.x1), (p.t, d.t, self.t0, self.t1)] {
            if dir == 0.0 {
                if !(pos >= lo && pos <= hi) {
                    return None;
                }
            } else {
                let (a, b) = ((lo - pos) / dir, (hi - pos) / dir);
                u0 = u0.max(a.min(b));
                u1 = u1.min(a.max(b));
            }
        }
        (u0 <= u1).then_some((u0, u1))
    }
    fn at(&self, p: Pt, d: Pt, u: f64) -> Pt {
        self.clamp(Pt {
            x: p.x + u * d.x,
            t: p.t + u * d.t,
        })
    }
    /// The part of the segment from `a` to `b` inside the window.
    pub fn clip_segment(&self, a: Pt, b: Pt) -> Option<(Pt, Pt)> {
        let d = Pt {
            x: b.x - a.x,
            t: b.t - a.t,
        };
        let (u0, u1) = self.clip(a, d, 0.0, 1.0)?;
        Some((self.end(a, b, d, u0), self.end(a, b, d, u1)))
    }
    /// The point at parameter `u` of the segment from `a` to `b`, exact at both ends.
    fn end(&self, a: Pt, b: Pt, d: Pt, u: f64) -> Pt {
        if u <= 0.0 {
            a
        } else if u >= 1.0 {
            b
        } else {
            self.at(a, d, u)
        }
    }
    /// The part of the whole line through `p` along `d` inside the window. `None` if it misses the window or
    /// `d` is zero.
    pub fn clip_line(&self, p: Pt, d: Pt) -> Option<(Pt, Pt)> {
        self.clip_between(p, d, f64::NEG_INFINITY, f64::INFINITY)
    }
    /// The part of the ray from `p` along `d` inside the window.
    pub fn clip_ray(&self, p: Pt, d: Pt) -> Option<(Pt, Pt)> {
        self.clip_between(p, d, 0.0, f64::INFINITY)
    }
    fn clip_between(&self, p: Pt, d: Pt, u0: f64, u1: f64) -> Option<(Pt, Pt)> {
        if d.x == 0.0 && d.t == 0.0 {
            return None;
        }
        let (u0, u1) = self.clip(p, d, u0, u1)?;
        (u0.is_finite() && u1.is_finite()).then(|| (self.at(p, d, u0), self.at(p, d, u1)))
    }
    /// The parts of a polyline inside the window, as separate polylines of two or more points.
    pub fn clip_polyline(&self, points: &[Pt]) -> Vec<Vec<Pt>> {
        let mut out = Vec::new();
        let mut run: Vec<Pt> = Vec::new();
        let mut flush = |run: &mut Vec<Pt>| {
            if run.len() >= 2 {
                out.push(std::mem::take(run));
            } else {
                run.clear();
            }
        };
        for pair in points.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            let d = Pt {
                x: b.x - a.x,
                t: b.t - a.t,
            };
            match self.clip(a, d, 0.0, 1.0) {
                Some((u0, u1)) => {
                    if run.is_empty() || u0 > 0.0 {
                        flush(&mut run);
                        run.push(self.end(a, b, d, u0));
                    }
                    run.push(self.end(a, b, d, u1));
                    if u1 < 1.0 {
                        flush(&mut run);
                    }
                }
                None => flush(&mut run),
            }
        }
        flush(&mut run);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const R: Rect = Rect {
        x0: 0.0,
        x1: 10.0,
        t0: 0.0,
        t1: 5.0,
    };

    fn pt(x: f64, t: f64) -> Pt {
        Pt { x, t }
    }

    #[test]
    fn segments_are_cut_at_the_window() {
        assert_eq!(
            R.clip_segment(pt(-5.0, 2.0), pt(15.0, 2.0)),
            Some((pt(0.0, 2.0), pt(10.0, 2.0)))
        );
        // ends inside are returned untouched
        assert_eq!(
            R.clip_segment(pt(1.0, 1.0), pt(2.0, 4.0)),
            Some((pt(1.0, 1.0), pt(2.0, 4.0)))
        );
        assert_eq!(R.clip_segment(pt(-5.0, 6.0), pt(15.0, 6.0)), None);
        assert_eq!(R.clip_segment(pt(11.0, 1.0), pt(12.0, 2.0)), None);
        // a single point
        assert_eq!(
            R.clip_segment(pt(3.0, 3.0), pt(3.0, 3.0)),
            Some((pt(3.0, 3.0), pt(3.0, 3.0)))
        );
        assert_eq!(R.clip_segment(pt(30.0, 3.0), pt(30.0, 3.0)), None);
        assert_eq!(R.clip_segment(pt(f64::NAN, 3.0), pt(3.0, 3.0)), None);
    }

    #[test]
    fn lines_and_rays_run_to_the_edges() {
        // a 45 degree line through (2, 0)
        assert_eq!(
            R.clip_line(pt(2.0, 0.0), pt(1.0, 1.0)),
            Some((pt(2.0, 0.0), pt(7.0, 5.0)))
        );
        // the same line from a point far outside
        let (a, b) = R.clip_line(pt(-98.0, -100.0), pt(1.0, 1.0)).unwrap();
        assert!((a.x - 2.0).abs() < 1e-12 && a.t == 0.0 && (b.x - 7.0).abs() < 1e-12 && b.t == 5.0);
        assert_eq!(
            R.clip_ray(pt(4.0, 1.0), pt(-1.0, 1.0)),
            Some((pt(4.0, 1.0), pt(0.0, 5.0)))
        );
        assert_eq!(R.clip_ray(pt(4.0, 6.0), pt(-1.0, 1.0)), None);
        assert_eq!(R.clip_line(pt(4.0, 1.0), pt(0.0, 0.0)), None);
        // vertical and horizontal lines
        assert_eq!(
            R.clip_line(pt(3.0, 99.0), pt(0.0, 1.0)),
            Some((pt(3.0, 0.0), pt(3.0, 5.0)))
        );
        assert_eq!(R.clip_line(pt(13.0, 99.0), pt(0.0, 1.0)), None);
    }

    #[test]
    fn polylines_split_where_they_leave() {
        let line = [
            pt(1.0, 1.0),
            pt(2.0, 2.0),
            pt(12.0, 2.0),
            pt(12.0, 4.0),
            pt(8.0, 4.0),
            pt(8.0, 4.5),
        ];
        let parts = R.clip_polyline(&line);
        assert_eq!(
            parts,
            vec![
                vec![pt(1.0, 1.0), pt(2.0, 2.0), pt(10.0, 2.0)],
                vec![pt(10.0, 4.0), pt(8.0, 4.0), pt(8.0, 4.5)],
            ]
        );
        assert!(R.clip_polyline(&[pt(1.0, 1.0)]).is_empty());
        assert!(R.clip_polyline(&[]).is_empty());
        assert!(R.clip_polyline(&[pt(-1.0, 1.0), pt(-2.0, 2.0)]).is_empty());
    }
}
