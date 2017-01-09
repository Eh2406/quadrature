pub mod double_exponential;
pub mod clenshaw_curtis;

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
