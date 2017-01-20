//! The double exponential algorithm is naturally adaptive, it stops calling the integrand when the error is reduced to below the desired threshold.
//! It also does not allocate. No box, no vec, etc.
//! It has a hard coded maximum of approximately 350 function evaluations. This guarantees that the algorithm will return.
//! The error in the algorithm decreases exponentially in the number of function evaluations, specifically O(exp(-cN/log(N))). So if 350 function evaluations is not giving the desired accuracy than the programmer probably needs to give some guidance by splitting up the range at singularities or [other preparation techniques](http://www.johndcook.com/blog/2012/02/21/care-and-treatment-of-singularities/).
//!
//! This is a port of the [Fast Numerical Integration](https://www.codeproject.com/kb/recipes/fastnumericalintegration.aspx) from c++ to rust. The original code is by John D. Cook, and is licensed under the [BSD](https://opensource.org/licenses/bsd-license.php).

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
/// use quadrature::integrate;
/// fn integrand(x: f64) -> f64 {
///     (-x / 5.0).exp() * x.powf(-1.0 / 3.0)
/// }
/// let o = integrate(integrand , 0.0, 10.0, 1e-6);
/// assert!((o.integral - 3.6798142583691758).abs() <= 1e-6);
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
    let mut error_estimate = ::std::f64::MAX;
    let mut num_function_evaluations = 1;
    let mut current_delta = ::std::f64::MAX;

    let mut integral = 2.0 * ::std::f64::consts::FRAC_PI_2 * f(0.0);

    for &weight in &WEIGHTS {
        let new_contribution = weight.iter()
            .map(|&(w, x)| w * (f(x) + f(-x)))
            .fold(0.0, |sum, x| sum + x);
        num_function_evaluations += 2 * weight.len();

        // difference in consecutive integral estimates
        let previous_delta_ln = current_delta.ln();
        current_delta = (0.5 * integral - new_contribution).abs();
        integral = 0.5 * integral + new_contribution;

        // Once convergence kicks in, error is approximately squared at each step.
        // Determine whether we're in the convergent region by looking at the trend in the error.
        if num_function_evaluations <= 13 {
            // level <= 1
            continue; // previousDelta meaningless, so cannot check convergence.
        }

        // Exact comparison with zero is harmless here.  Could possibly be replaced with
        // a small positive upper limit on the size of currentDelta, but determining
        // that upper limit would be difficult.  At worse, the loop is executed more
        // times than necessary.  But no infinite loop can result since there is
        // an upper bound on the loop variable.
        if current_delta == 0.0 {
            error_estimate = 0.0;
            break;
        }
        // previousDelta != 0 or would have been kicked out previously
        let r = current_delta.ln() / previous_delta_ln;

        if r > 1.9 && r < 2.1 {
            // If convergence theory applied perfectly, r would be 2 in the convergence region.
            // r close to 2 is good enough. We expect the difference between this integral estimate
            // and the next one to be roughly delta^2.
            error_estimate = current_delta * current_delta;
        } else {
            // Not in the convergence region.  Assume only that error is decreasing.
            error_estimate = current_delta;
        }

        if error_estimate < target_absolute_error {
            break;
        }
    }

    Output {
        num_function_evaluations: num_function_evaluations as u32,
        error_estimate: error_estimate,
        integral: integral,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    unit_test!(trivial_function_works = |_| 0.5; -1.0..1.0; 1e-14);
    unit_test!(demo_function1_works = |x| (-x / 5.0).exp() * x.powf(-1.0 / 3.0); 0.0..10.0; 1e-6 => 3.6798142583691758);
    unit_test!(demo_function2_works = |x| (1.0 - x).powf(5.0) * x.powf(-1.0 / 3.0); 0.0..1.0; 1e-6 => 0.41768525592055004);
    unit_test!(demo_function3_works = |x| (-x / 5000.0).exp() * (x / 1000.0).powf(-1.0 / 3.0); 0.0..10000.0; 1e-6);
    unit_test!(demo_bad_function1_works = |x| (1.0 - x).powf(0.99); 0.0..1.0; 1e-6 => 0.50251256281407035);
    unit_test!(demo_bad_function2_works = |x| x.abs(); -1.0..1.0; 1e-6 => 1.0; 385);
    unit_test!(demo_bad_function3_works = |x| (0.5 - x.abs()).abs(); -1.0..1.0; 1e-6 => 0.5; 385);
}
