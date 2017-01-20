macro_rules! unit_test {
($name:ident = $inta:expr ; $lim:expr; $eps:expr => $out:expr; $max:expr) => (
    #[test]
    fn $name() {
        let o = integrate($inta, $lim.start, $lim.end, $eps);
        println!("\n{:#?}", o);
        assert!(o.num_function_evaluations == $max,
                "num_function_evaluations is not maxed out. evaluations: {:#?}, max: {:#?}",
                o.num_function_evaluations,
                $max);
        assert!((o.integral - $out).abs() <= o.error_estimate,
                "error larger then error_estimate. error: {:#?}, estimate: {:#?}",
                (o.integral - $out).abs(),
                o.error_estimate);
    }
);
($name:ident = $inta:expr ; $lim:expr; $eps:expr => $out:expr) => (
    #[test]
    fn $name() {
        let o = integrate($inta, $lim.start, $lim.end, $eps);
        println!("\n{:#?}", o);
        assert!(o.error_estimate <= $eps,
                "error_estimate larger then asked. estimate: {:#?}, asked: {:#?}",
                o.error_estimate,
                $eps);
        assert!((o.integral - $out).abs() <= $eps,
                "error larger then asked. error: {:#?}, estimate: {:#?}",
                (o.integral - $out).abs(),
                $eps);
    }
);
($name:ident = $inta:expr ; $lim:expr; $eps:expr) => (
    #[test]
    fn $name() {
        let o = integrate($inta, $lim.start, $lim.end, $eps);
        println!("\n{:#?}", o);
        assert!(o.error_estimate <= $eps,
                "error_estimate larger then asked. estimate: {:#?}, asked: {:#?}",
                o.error_estimate,
                $eps);
    }
)
}
