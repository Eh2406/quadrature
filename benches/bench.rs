#![feature(test)]
extern crate test;
extern crate quadrature;


macro_rules! unit_bench {
    ($name:ident = $inta:expr ; $lim:expr; $eps:expr) => {
        #[bench]
        fn $name(b: &mut Bencher) {
            b.iter(|| integrate($inta, $lim.start, $lim.end, $eps));
        }
    }
}
macro_rules! make_bench {
    ($($name:ident = $inta:expr ; $lim:expr; $eps:expr);*;) => {
        mod double_exponential {
            use test::Bencher;
            use quadrature::double_exponential::integrate;
            use std::f64::consts::PI;
            $(unit_bench!{$name = $inta; $lim; $eps})*
        }

        mod clenshaw_curtis {
            use test::Bencher;
            use quadrature::clenshaw_curtis::integrate;
            use std::f64::consts::PI;
            $(unit_bench!{$name = $inta; $lim; $eps})*
        }
    }
}

make_bench! {
    trivial_function = |_| 0.5; - 1.0..1.0; 1e-14;
    fn_4_iii = |x| ((PI * x).cos()) * ((1.0 - x).sqrt()); - 1.0..1.0; 1e-14;
    demo_function1 = |x| (-x / 5.0).exp() * x.powf(-1.0 / 3.0); 0.0..10.0; 1e-6;
    demo_function2 = |x| (1.0 - x).powf(5.0) * x.powf(-1.0 / 3.0); 0.0..1.0; 1e-6;
    demo_circle = |x| ((1.0 - (x.powi(2))).sqrt()).abs(); 0.0..1.0; 1e-6;
    demo_bad_circle = |x| ((1.0 - (x.powi(2))).sqrt() - 0.7).abs(); 0.0..1.0; 1e-6;
}
