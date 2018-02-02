//! The primary function of this library is `integrate`, witch uses the double exponential algorithm.
//! It is a port of the [Fast Numerical Integration](https://www.codeproject.com/kb/recipes/fastnumericalintegration.aspx) from c++ to rust. The original code is by John D. Cook, and is licensed under the BSD.
//!
//! Other Algorithms are provided in modules.
#![no_std]

extern crate num_traits;

pub mod double_exponential;
pub mod clenshaw_curtis;
mod traits;

#[derive(Clone, Copy, Debug)]
pub struct Output {
    pub num_function_evaluations: u32,
    pub error_estimate: f64,
    pub integral: f64,
}

impl Output {
    fn scale(self, c: f64) -> Self {
        Output {
            num_function_evaluations: self.num_function_evaluations,
            error_estimate: c * self.error_estimate,
            integral: c * self.integral,
        }
    }
}

pub use double_exponential::integrate;
