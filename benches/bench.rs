#![feature(test)]
extern crate test;
extern crate quadrature;

mod double_exponential {
    use test::Bencher;
    use quadrature::double_exponential::integrate;
    use std::f64::consts::PI;

    #[bench]
    fn trivial_function(b: &mut Bencher) {
        b.iter(|| integrate(|_| 0.5, -1.0, 1.0, 1e-14));
    }

    #[bench]
    fn fn_4_iii(b: &mut Bencher) {
        b.iter(|| integrate(|x| ((PI * x).cos()) * ((1.0 - x).sqrt()), -1.0, 1.0, 1e-6));
    }

    #[bench]
    fn demo_function1(b: &mut Bencher) {
        b.iter(|| integrate(|x| (-x / 5.0).exp() * x.powf(-1.0 / 3.0), 0.0, 10.0, 1e-6));
    }

    #[bench]
    fn demo_function2(b: &mut Bencher) {
        b.iter(|| integrate(|x| (1.0 - x).powf(5.0) * x.powf(-1.0 / 3.0), 0.0, 1.0, 1e-6));
    }

    #[bench]
    fn demo_circle(b: &mut Bencher) {
        b.iter(|| integrate(|x| ((1.0 - (x.powi(2))).sqrt()).abs(), 0.0, 1.0, 1e-6));
    }

    #[bench]
    fn demo_bad_circle(b: &mut Bencher) {
        b.iter(|| integrate(|x| ((1.0 - (x.powi(2))).sqrt() - 0.7).abs(), 0.0, 1.0, 1e-6));
    }
}

mod clenshaw_curtis {
    use test::Bencher;
    use quadrature::clenshaw_curtis::integrate;
    use std::f64::consts::PI;

    #[bench]
    fn trivial_function(b: &mut Bencher) {
        b.iter(|| integrate(|_| 0.5, -1.0, 1.0, 1e-14));
    }

    #[bench]
    fn fn_4_iii(b: &mut Bencher) {
        b.iter(|| integrate(|x| ((PI * x).cos()) * ((1.0 - x).sqrt()), -1.0, 1.0, 1e-6));
    }

    #[bench]
    fn demo_function1(b: &mut Bencher) {
        b.iter(|| integrate(|x| (-x / 5.0).exp() * x.powf(-1.0 / 3.0), 0.0, 10.0, 1e-6));
    }

    #[bench]
    fn demo_function2(b: &mut Bencher) {
        b.iter(|| integrate(|x| (1.0 - x).powf(5.0) * x.powf(-1.0 / 3.0), 0.0, 1.0, 1e-6));
    }

    #[bench]
    fn demo_circle(b: &mut Bencher) {
        b.iter(|| integrate(|x| ((1.0 - (x.powi(2))).sqrt()).abs(), 0.0, 1.0, 1e-6));
    }

    #[bench]
    fn demo_bad_circle(b: &mut Bencher) {
        b.iter(|| integrate(|x| ((1.0 - (x.powi(2))).sqrt() - 0.7).abs(), 0.0, 1.0, 1e-6));
    }
}