This is a port of the [Fast Numerical Integration](https://www.codeproject.com/kb/recipes/fastnumericalintegration.aspx) from c++ to rust. The original code is by John D. Cook, and is licensed under the [BSD](https://opensource.org/licenses/bsd-license.php) Other quadrature techniques may be added in the future.


I have been testing against `rust 1.0` with:
```cmd
> rustup override set 1.0.0
> cargo test
```
I think bumping rust version is a breaking change, and will be respected in semver.
