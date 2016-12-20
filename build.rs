// build.rs

use std::env;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::f64::consts::FRAC_PI_2;

fn abcissas (x: f64) -> f64 {(FRAC_PI_2 * x.sinh()).tanh()}

fn weights (x: f64, h:f64) -> f64 {FRAC_PI_2 * h *(x.cosh()/ ((FRAC_PI_2 * x.sinh()).cosh()).powi(2))}

fn make_first_row(h: f64, d: f64) -> Vec<(f64, f64)> {
    let mut row = Vec::with_capacity((d / h) as usize);
    let mut i = h;
    while i <= d {
        row.push((weights(i, h), abcissas(i)));
        i += h;
    }
    row
}

fn make_row(h: f64, d: f64) -> Vec<(f64, f64)> {
    let mut row = Vec::with_capacity((d / h) as usize);
    let mut i = h;
    while i <= d {
        row.push((weights(i, h), abcissas(i)));
        i += 2.0 * h;
    }
    row
}

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("double_exponential_constants.rs");
    let mut f = File::create(&dest_path).unwrap();

    f.write_all(b"const WEIGHTS: [&'static [(f64, f64)]; 7] =\n[").unwrap();

    let inf = 3.0;

    f.write_all(&format!("&{:#?},", &make_first_row(1.0, inf)).into_bytes()).unwrap();
    f.write_all(&format!("&{:#?},", &make_row(0.5, inf)).into_bytes()).unwrap();
    f.write_all(&format!("&{:#?},", &make_row(0.25, inf)).into_bytes()).unwrap();
    f.write_all(&format!("&{:#?},", &make_row(0.125, inf)).into_bytes()).unwrap();
    f.write_all(&format!("&{:#?},", &make_row(0.0625, inf)).into_bytes()).unwrap();
    f.write_all(&format!("&{:#?},", &make_row(0.03125, inf)).into_bytes()).unwrap();
    f.write_all(&format!("&{:#?},", &make_row(0.015625, inf)).into_bytes()).unwrap();

    f.write_all(b"];").unwrap();
}