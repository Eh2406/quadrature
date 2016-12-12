mod double_exponential_constants;
use double_exponential_constants::*;

#[derive(Clone, Copy, Debug)]
pub struct Output {
    pub num_function_evaluations: u32,
    pub error_estimate: f64,
    pub integral: f64,
}

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
    integrate_core(f, c, 0.5 * (a + b), target_absolute_error)
}

/// Integrate f(cx + d) with the given integration constants
fn integrate_core<F>(f: F, c: f64, d: f64, target_absolute_error: f64) -> Output
    where F: Fn(f64) -> f64
{
    let target_absolute_error = 0.25 * target_absolute_error / c;
    let mut error_estimate = ::std::f64::MAX;
    let mut h = 1.0;
    let mut num_function_evaluations = 0;
    let mut previous_delta;
    let mut current_delta = ::std::f64::MAX;

    let mut integral = WEIGHTS[0] * f(c * ABCISSAS[0] + d);
    num_function_evaluations += 1;

    let func = |i| WEIGHTS[i] * (f(c * ABCISSAS[i] + d) + f(-c * ABCISSAS[i] + d));

    for i in OFFSETS[0]..OFFSETS[1] {
        integral += func(i);
    }
    num_function_evaluations += OFFSETS[1] - OFFSETS[0];

    for level in 1..(OFFSETS.len() - 1) {
        h *= 0.5;
        let mut new_contribution = 0.0;
        for i in OFFSETS[level]..OFFSETS[level + 1] {
            new_contribution += func(i);
        }
        num_function_evaluations += OFFSETS[level + 1] - OFFSETS[level];
        new_contribution *= h;

        // difference in consecutive integral estimates
        previous_delta = current_delta;
        current_delta = (0.5 * integral - new_contribution).abs();
        integral = 0.5 * integral + new_contribution;

        // Once convergence kicks in, error is approximately squared at each step.
        // Determine whether we're in the convergent region by looking at the trend in the error.
        if level <= 1 {
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
        let r = current_delta.ln() / previous_delta.ln();

        if r > 1.9 && r < 2.1 {
            // If convergence theory applied perfectly, r would be 2 in the convergence region.
            // r close to 2 is good enough. We expect the difference between this integral estimate
            // and the next one to be roughly delta^2.
            error_estimate = previous_delta * previous_delta;
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
        error_estimate: c * error_estimate,
        integral: c * integral,
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
        assert!((o.integral - 3.6798142583691758).abs() <= o.error_estimate,
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
