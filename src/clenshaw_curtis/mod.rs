//! The `clenshaw_curtis` module provides a `integrate` function with the same signature as `quadrature::integrate`.
//! The implemented variant of clenshaw curtis quadrature is adaptive, however the weights change for each adaptation. This unfortunately means that the sum needs to be recalculated for each layer of adaptation.
//! It also does not allocate on the heap, however it does use a `[f64; 129]` to store the function values. It has a hard coded maximum of approximately 257 function evaluations. This guarantees that the algorithm will return.
//! The clenshaw curtis algorithm exactly integrates polynomials of order N. This implementation starts with an N of approximately 5 and increases up to an N of approximately 257. In general the error in the algorithm decreases exponentially in the number of function evaluations. In summery clenshaw curtis will in general use **more stack space** and **run slower** than the double exponential algorithm, unless clenshaw curtis can get the exact solution.

use core::f64;

mod constants;

use self::constants::*;

use super::Output;

/// Integrate an analytic function over a finite interval.
/// f is the function to be integrated.
/// a is left limit of integration.
/// b is right limit of integration
/// target_absolute_error is the desired bound on error
///
/// # Examples
///
/// ```
/// use quadrature::clenshaw_curtis::integrate;
/// fn integrand(x: f64) -> f64 {
///     (-x / 5.0).exp() * x.powf(-1.0 / 3.0)
/// }
/// let o = integrate(integrand , 0.0, 10.0, 1e-6);
/// assert!((o.integral - 3.6798142583691758).abs() <= o.error_estimate)
/// ```
pub fn integrate<F>(f: F, a: f64, b: f64, target_absolute_error: f64) -> Output
    where F: Fn(f64) -> f64
{
    // Apply the linear change of variables x = ct + d
    // $$\int_a^b f(x) dx = c \int_{-1}^1 f( ct + d ) dt$$
    // c = (b-a)/2, d = (a+b)/2
    let c = 0.5 * (b - a);
    let d = 0.5 * (a + b);
    integrate_core(|x| {
        let out = f(c * x + d);
        if out.is_finite() { out } else { 0.0 }
    },
                   0.25 * target_absolute_error / c)
        .scale(c)
}

/// Integrate f(x) from [-1.0, 1.0]
fn integrate_core<F>(f: F, target_absolute_error: f64) -> Output
    where F: Fn(f64) -> f64
{
    let mut f_value = [::core::f64::NAN; 129];
    debug_assert_eq!(f_value.len(), ABCISSAS.len());
    let mut max_x_idx = 1;
    f_value[0] = f(0.0);

    let mut error_estimate = ::core::f64::MAX;
    let mut integral = ::core::f64::MAX;
    for &w in WEIGHTS.iter() {
        for (v, &x) in f_value[max_x_idx..w.len()]
            .iter_mut()
            .zip(&ABCISSAS[max_x_idx..w.len()]) {
            *v = f(x) + f(-x);
        }
        max_x_idx = w.len();

        let last_integral = integral;
        integral = f_value.iter().zip(w.iter()).fold(0.0, |sum, x| sum + (x.0 * x.1));
        error_estimate = (last_integral - integral).abs();

        if error_estimate < target_absolute_error {
            break;
        }
    }

    Output {
        num_function_evaluations: (max_x_idx * 2 - 1) as u32,
        error_estimate: error_estimate.abs(),
        integral: integral,
    }
}

