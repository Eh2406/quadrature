#![feature(test)]
extern crate test;
use test::Bencher;

extern crate quadrature;
use quadrature::{double_exponential, clenshaw_curtis};

use std::f64::consts::PI;

#[bench]
fn double_exponential_trivial_function(b: &mut Bencher) {
    b.iter(|| double_exponential::integrate(|_| 0.5, -1.0, 1.0, 1e-14));
}

#[bench]
fn double_exponential_4_iii(b: &mut Bencher) {
    b.iter(|| double_exponential::integrate(|x| ((PI * x).cos()) * ((1.0 - x).sqrt()), -1.0, 1.0, 1e-6));
}

#[bench]
fn double_exponential_demo_function1(b: &mut Bencher) {
    b.iter(|| double_exponential::integrate(|x| (-x / 5.0).exp() * x.powf(-1.0 / 3.0), 0.0, 10.0, 1e-6));
}

#[bench]
fn double_exponential_demo_function2(b: &mut Bencher) {
    b.iter(|| double_exponential::integrate(|x| (1.0 - x).powf(5.0) * x.powf(-1.0 / 3.0), 0.0, 1.0, 1e-6));
}

#[bench]
fn clenshaw_curtis_trivial_function(b: &mut Bencher) {
    b.iter(|| clenshaw_curtis::integrate(|_| 0.5, -1.0, 1.0, 1e-14));
}

#[bench]
fn clenshaw_curtis_4_iii(b: &mut Bencher) {
    b.iter(|| clenshaw_curtis::integrate(|x| ((PI * x).cos()) * ((1.0 - x).sqrt()), -1.0, 1.0, 1e-6));
}

#[bench]
fn clenshaw_curtis_demo_function1(b: &mut Bencher) {
    b.iter(|| clenshaw_curtis::integrate(|x| (-x / 5.0).exp() * x.powf(-1.0 / 3.0), 0.0, 10.0, 1e-6));
}

#[bench]
fn clenshaw_curtis_demo_function2(b: &mut Bencher) {
    b.iter(|| clenshaw_curtis::integrate(|x| (1.0 - x).powf(5.0) * x.powf(-1.0 / 3.0), 0.0, 1.0, 1e-6));
}