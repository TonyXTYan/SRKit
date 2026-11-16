//! Why a diagram refused its input.

/// Why a diagram, or something added to it, could not be built.
///
/// Every function that takes numbers from the caller checks them and returns one of these. Nothing in the crate
/// panics on input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    /// A view range is empty, reversed or too narrow to resolve next to its own offset, or the plot is smaller
    /// than one pixel.
    EmptyView,
    /// An input is NaN or infinite, or it leads to a position that is.
    NonFinite,
    /// A view bound beyond ±1e150, or a plot larger than [`crate::MAX_PIXELS`] on a side.
    TooLarge,
    /// The axis direction has no length.
    ZeroDirection,
    /// A stroke width that is not a positive, finite number of at most 64 pixels.
    InvalidStyle,
    /// A parameter outside the range the function accepts: a step, a distance or an acceleration that is not
    /// positive, or a speed at or above `c`.
    OutOfRange,
    /// The two events of a signal are not joined by light: the time between them differs from their distance.
    NotLightLike,
    /// More than [`crate::MAX_MARKS`] ticks or grid lines would be drawn. Choose a larger step.
    TooManyMarks,
    /// A worldline needs more than [`crate::MAX_POINTS`] points to stay within the curve tolerance, as a circle
    /// with thousands of turns in view does.
    TooDetailed,
    /// The observer sees nothing of the source: no light from the drawn part of its worldline reaches the event.
    NothingSeen,
    /// The segment has no Rindler horizon: its proper acceleration is not constant.
    NoHorizon,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::EmptyView => "the view is empty or too narrow to draw",
            Self::NonFinite => "an input is not a finite number",
            Self::TooLarge => "the view or the plot is too large",
            Self::ZeroDirection => "the axis direction has no length",
            Self::InvalidStyle => "the stroke width is not a positive number of pixels",
            Self::OutOfRange => "a parameter is outside the accepted range",
            Self::NotLightLike => "the two events are not joined by light",
            Self::TooManyMarks => "too many ticks or grid lines; choose a larger step",
            Self::TooDetailed => "the worldline needs too many points to draw",
            Self::NothingSeen => "no light from the source reaches the observer",
            Self::NoHorizon => "the segment has no Rindler horizon",
        })
    }
}

impl std::error::Error for Error {}
