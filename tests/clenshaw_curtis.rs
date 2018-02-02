extern crate quadrature;

use quadrature::clenshaw_curtis::integrate;

#[macro_use]
mod test_macro;

unit_test!(trivial_function_works = |_| 0.5; -1.0..1.0; 1e-14);
unit_test!(demo_function1_works = |x| (-x / 5.0).exp() * x.powf(-1.0 / 3.0); 0.0..10.0; 1e-6 => 3.6798142583691758; 257);
unit_test!(demo_function2_works = |x| (1.0 - x).powf(5.0) * x.powf(-1.0 / 3.0); 0.0..1.0; 1e-6 => 0.41768525592055004; 257);
unit_test!(demo_function3_works = |x| (-x / 5000.0).exp() * (x / 1000.0).powf(-1.0 / 3.0); 0.0..10000.0; 1e-6 => 3679.81425836918; 257);
unit_test!(demo_bad_function1_works = |x| (1.0 - x).powf(0.99); 0.0..1.0; 1e-6 => 0.50251256281407035);
unit_test!(demo_bad_function2_works = |x| x.abs(); -1.0..1.0; 1e-6 => 1.0; 257);
unit_test!(demo_bad_function3_works = |x| (0.5 - x.abs()).abs(); -1.0..1.0; 1e-6 => 0.5; 257);
unit_test!(demo_circle = |x| ((1.0-(x.powi(2))).sqrt()).abs(); -1.0..1.0; 1e-6 => 1.5707963267949);
unit_test!(demo_bad_circle = |x| ((1.0-(x.powi(2))).sqrt()-0.7).abs(); -1.0..1.0; 1e-4 => 0.420201353577392);
