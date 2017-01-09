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
    let mut f_value = [::std::f64::NAN; 129];
    debug_assert_eq!(f_value.len(), ABCISSAS.len());
    let mut max_x_idx = 1;
    f_value[0] = f(0.0);

    let mut error_estimate = ::std::f64::MAX;
    let mut integral = ::std::f64::MAX;
    for &w in WEIGHTS.iter() {
        for (v, &x) in f_value[max_x_idx..w.len()]
            .iter_mut()
            .zip(&ABCISSAS[max_x_idx..w.len()]) {
            *v = f(x) + f(-x);
        }
        max_x_idx = w.len();

        let last_integral = integral;
        debug_assert_eq!(f_value[..w.len()].len(), w.len());
        integral = f_value[..w.len()].iter().zip(w.iter()).fold(0.0, |sum, x| sum + (x.0 * x.1));
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn trivial_function_works() {
        let o = integrate(|_| 0.5, -1.0, 1.0, 1e-14);
        assert!(o.error_estimate <= 1e-14,
                "error_estimate larger then asked. estimate: {:#?}, asked: {:#?}",
                o.error_estimate,
                1e-14);
    }

    #[test]
    fn demo_function1_works() {
        let o = integrate(|x| (-x / 5.0).exp() * x.powf(-1.0 / 3.0), 0.0, 10.0, 1e-6);
        assert!((o.integral - 3.6798142583691758).abs() <= o.error_estimate,
                "error larger then error_estimate");
    }

    #[test]
    fn demo_function2_works() {
        let o = integrate(|x| (1.0 - x).powf(5.0) * x.powf(-1.0 / 3.0), 0.0, 1.0, 1e-6);
        assert!((o.integral - 0.41768525592055004).abs() <= o.error_estimate,
                "error larger then error_estimate");
    }

    #[test]
    fn demo_function3_works() {
        let o = integrate(|x| (-x / 5000.0).exp() * (x / 1000.0).powf(-1.0 / 3.0),
                          0.0,
                          10000.0,
                          1e-6);
        assert!((o.integral - 3679.81425836918).abs() <= o.error_estimate,
                "error larger then error_estimate");
    }

    #[test]
    fn demo_bad_function1_works() {
        let o = integrate(|x| (1.0 - x).powf(0.99), 0.0, 1.0, 1e-6);
        assert!(o.error_estimate <= 1e-6,
                "error_estimate larger then asked. estimate: {:#?}, asked: {:#?}",
                o.error_estimate,
                1e-6);
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
