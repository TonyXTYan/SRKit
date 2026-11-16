//! Print the twin-paradox diagram as SVG.
//!
//! ```sh
//! cargo run -p srkit-edu --example twin_paradox > twins.svg
//! cargo run -p srkit-edu --example twin_paradox -- 10 0.5 dark > twins.svg
//! ```
//!
//! The arguments are the one-way distance in light-years (default 4.2465, Proxima Centauri), the proper
//! acceleration in light-years per year squared (default 1.0323, one gravity) and `dark` for the dark theme.

use std::process::ExitCode;

use srkit_edu::{twin_paradox_diagram, Theme};

fn main() -> ExitCode {
    let mut numbers = [4.2465, 1.0323];
    let mut theme = Theme::Light;
    let mut next = 0;
    for argument in std::env::args().skip(1) {
        if argument == "dark" {
            theme = Theme::Dark;
        } else if let (Ok(value), true) = (argument.parse::<f64>(), next < numbers.len()) {
            numbers[next] = value;
            next += 1;
        } else {
            eprintln!("usage: twin_paradox [distance] [acceleration] [dark]");
            return ExitCode::FAILURE;
        }
    }
    match twin_paradox_diagram(numbers[0], numbers[1]) {
        Ok(mut diagram) => {
            print!("{}", diagram.theme(theme).render());
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("twin_paradox: {error}");
            ExitCode::FAILURE
        }
    }
}
