//! The double exponential rule is based on the observation that the trapezoid rule converges
//! very rapidly for functions on the entire real line that go to zero like exp( - exp(t) ).
//! The change of variables x = tanh( pi sinh(t) /2) transforms an integral over [-1, 1]
//! into an integral with integrand suited to the double exponential rule.
//!
//! The transformed integral is infinite, but we truncate the domain of integration to [-3, 3].
//! The limit '3' was chosen for two reasons: for t = 3, the transformed x values
//! are nearly equal to 1 for 12 or more significant figures.  Also, for t = 3, the
//! smallest weights are 12 orders of magnitude smaller than the largest weights; setting
//! the cutoff larger than 3 would not have a significant impact on the integral value
//! unless there is a strong singularity at one of the end points.
//!
//! The change of variables x(t) is an odd function, i.e. x(-t) = -x(t), and so we need only
//! store the positive x values.  Also, the derivative w(t) = x'(t) is even, i.e. w(-t) = w(t),
//! and so we need only store the weights corresponding to positive values of x.
//!
//! The integration first applies the trapezoid rule to [-3, 3] in steps of size 1.
//! Then it subsequently cuts the step size in half each time, comparing the results.
//! Integration stops when subsequent iterations are close enough together or the maximum
//! integration points have been used.
//! By cutting h in half, the previous integral can be reused; we only need evaluate the
//! integrand at the newly added points.
//!
//! Finally, note that we're not strictly using the trapezoid rule: we don't treat the
//! end points differently.  This is because we assume the values at the ends of the interval
//! hardly matter due to the rapid decay of the integrand.

include!(concat!(env!("OUT_DIR"), "/double_exponential_constants.rs"));

use super::Output;

/// Integrate an analytic function over a finite interval.
/// f is the function to be integrated.
/// a is left limit of integration.
/// b is right limit of integration
/// target_absolute_error is the desired bound on error
pub fn integrate<F>(f: F, a: f64, b: f64, target_absolute_error: f64) -> Output
    where F: Fn(f64) -> f64
{
    // Apply the linear change of variables x = ct + d
    // $$\int_a^b f(x) dx = c \int_{-1}^1 f( ct + d ) dt$$
    // c = (b-a)/2, d = (a+b)/2
    let c = 0.5 * (b - a);
    let d = 0.5 * (a + b);
    integrate_core(|x| f(c * x + d), 0.25 * target_absolute_error / c).scale(c)
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
    #[test]
    fn trivial_function_works() {
        let o = integrate(|_| 0.5, -1.0, 1.0, 1e-14);
        assert!(o.error_estimate <= 1e-14,
                "error_estimate larger then asked");
    }

    #[test]
    fn demo_function1_works() {
        let o = integrate(|x| (-x / 5.0).exp() * x.powf(-1.0 / 3.0), 0.0, 10.0, 1e-6);
        assert!(o.error_estimate <= 1e-6, "error_estimate larger then asked");
        assert!((o.integral - 3.6798142583691758).abs() <= 1e-6,
                "error larger then error_estimate");
    }

    #[test]
    fn demo_function2_works() {
        let o = integrate(|x| (1.0 - x).powf(5.0) * x.powf(-1.0 / 3.0), 0.0, 1.0, 1e-6);
        assert!(o.error_estimate <= 1e-6, "error_estimate larger then asked");
        assert!((o.integral - 0.41768525592055004).abs() <= o.error_estimate,
                "error larger then error_estimate");
    }

    #[test]
    fn demo_function3_works() {
        let o = integrate(|x| (-x / 5000.0).exp() * (x / 1000.0).powf(-1.0 / 3.0),
                          0.0,
                          10000.0,
                          1e-6);
        assert!(o.error_estimate <= 1e-6, "error_estimate larger then asked");
    }

    #[test]
    fn demo_bad_function1_works() {
        let o = integrate(|x| (1.0 - x).powf(0.99), 0.0, 1.0, 1e-6);
        assert!(o.error_estimate <= 1e-6, "error_estimate larger then asked");
        assert!((o.integral - 0.50251256281407035).abs() <= o.error_estimate,
                "error larger then error_estimate");
    }

    #[test]
    fn demo_bad_function2_works() {
        let o = integrate(|x| x.abs(), -1.0, 1.0, 1e-6);
        assert!((o.integral - 1.0).abs() <= o.error_estimate,
                "error larger then error_estimate");
    }

    #[test]
    fn demo_bad_function3_works() {
        let o = integrate(|x| (0.5 - x.abs()).abs(), -1.0, 1.0, 1e-6);
        assert!((o.integral - 0.5).abs() <= o.error_estimate,
                "error larger then error_estimate");
    }
}
