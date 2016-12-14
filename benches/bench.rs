#![feature(test)]
extern crate test;
use test::Bencher;

extern crate quadrature;
use quadrature::integrate;

use std::f64::consts::PI;

#[bench]
fn bench_trivial_function(b: &mut Bencher) {
    b.iter(|| integrate(|_| 0.5, -1.0, 1.0, 1e-14));
}

#[bench]
fn bench_4_iii(b: &mut Bencher) {
    b.iter(|| integrate(|x| ((PI * x).cos()) * ((1.0 - x).sqrt()), -1.0, 1.0, 1e-6));
}

#[bench]
fn bench_demo_function1(b: &mut Bencher) {
    b.iter(|| integrate(|x| (-x / 5.0).exp() * x.powf(-1.0 / 3.0), 0.0, 10.0, 1e-6));
}

#[bench]
fn bench_demo_function2(b: &mut Bencher) {
    b.iter(|| integrate(|x| (1.0 - x).powf(5.0) * x.powf(-1.0 / 3.0), 0.0, 1.0, 1e-6));
}
