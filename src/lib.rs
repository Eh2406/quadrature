mod double_exponential;

#[derive(Clone, Copy, Debug)]
pub struct Output {
    pub num_function_evaluations: u32,
    pub error_estimate: f64,
    pub integral: f64,
}

pub use double_exponential::integrate;
