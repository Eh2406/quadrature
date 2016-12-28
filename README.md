This is a port of the [Fast Numerical Integration](https://www.codeproject.com/kb/recipes/fastnumericalintegration.aspx) from c++ to rust. The original code is by John D. Cook, and is licensed under the [BSD](https://opensource.org/licenses/bsd-license.php) Other quadrature techniques may be added in the future.

This library is pure safe rust, and cross-platform. The double exponential algorithm is naturally adaptive, it stops calling the integrand when the error is reduced to below the desired threshold. It also does not allocate. No box, no vec, etc. It has a hard coded maximum of approximately 350 function evaluations. This guarantees that the algorithm will return. The error in the algorithm decreases exponentially in the number of function evaluations, specifically O(exp(-cN/log(N))). So if 350 function evaluations is not giving the desired accuracy than the programmer probably needs to give some guidance by splitting up the range at singularities or [other preparation techniques](http://www.johndcook.com/blog/2012/02/21/care-and-treatment-of-singularities/).

I have been testing against `rust 1.0` with:
```cmd
> rustup override set 1.0.0
> cargo test
> rustup override set nightly
> cargo test
> cargo bench
```
It is also tested on travis [![Build Status](https://travis-ci.org/Eh2406/quadrature.svg?branch=master)](https://travis-ci.org/Eh2406/quadrature).
I think bumping rust version is a breaking change, and will be respected in semver.
