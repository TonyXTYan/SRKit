//! Print the Rindler-chase diagram as SVG.
//!
//! ```sh
//! cargo run -p srkit-edu --example rindler_chase > chase.svg
//! cargo run -p srkit-edu --example rindler_chase -- 0.5 dark > chase.svg
//! ```
//!
//! The arguments are the proper acceleration in light-years per year squared (default 1.0323, one gravity)
//! and `dark` for the dark theme.

use std::process::ExitCode;

use srkit_edu::{rindler_chase_diagram, Theme};

fn main() -> ExitCode {
    let mut acceleration = 1.0323;
    let mut theme = Theme::Light;
    for argument in std::env::args().skip(1) {
        if argument == "dark" {
            theme = Theme::Dark;
        } else if let Ok(value) = argument.parse::<f64>() {
            acceleration = value;
        } else {
            eprintln!("usage: rindler_chase [acceleration] [dark]");
            return ExitCode::FAILURE;
        }
    }
    match rindler_chase_diagram(acceleration) {
        Ok(mut diagram) => {
            print!("{}", diagram.theme(theme).render());
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("rindler_chase: {error}");
            ExitCode::FAILURE
        }
    }
}
