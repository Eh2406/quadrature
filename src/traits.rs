use core::ops::Neg;

use num_traits::Num;

/// Basic floating point operations that work with `core`.
pub trait BasicFloat: Num + Neg<Output = Self> + PartialOrd + Copy {
    /// Returns positive infinity.
    #[inline]
    fn infinity() -> Self;

    /// Returns negative infinity.
    #[inline]
    fn neg_infinity() -> Self;

    /// Returns NaN.
    #[inline]
    fn nan() -> Self;

    /// Returns `true` if the number is NaN.
    #[inline]
    fn is_nan(self) -> bool {
        self != self
    }

    /// Returns `true` if the number is infinite.
    #[inline]
    fn is_infinite(self) -> bool {
        self == Self::infinity() || self == Self::neg_infinity()
    }

    /// Returns `true` if the number is neither infinite or NaN.
    #[inline]
    fn is_finite(self) -> bool {
        !(self.is_nan() || self.is_infinite())
    }

    /// Computes the absolute value of `self`. Returns `BasicFloat::nan()` if the
    /// number is `BasicFloat::nan()`.
    #[inline]
    fn abs(self) -> Self {
        if self.is_sign_positive() {
            return self;
        }
        if self.is_sign_negative() {
            return -self;
        }
        Self::nan()
    }

    /// Returns a number that represents the sign of `self`.
    ///
    /// - `1.0` if the number is positive, `+0.0` or `BasicFloat::infinity()`
    /// - `-1.0` if the number is negative, `-0.0` or `BasicFloat::neg_infinity()`
    /// - `BasicFloat::nan()` if the number is `BasicFloat::nan()`
    #[inline]
    fn signum(self) -> Self {
        if self.is_sign_positive() {
            return Self::one();
        }
        if self.is_sign_negative() {
            return -Self::one();
        }
        Self::nan()
    }

    /// Returns `true` if `self` is positive, including `+0.0` and
    /// `BasicFloat::infinity()`.
    #[inline]
    fn is_sign_positive(self) -> bool {
        self > Self::zero() || (Self::one() / self) == Self::infinity()
    }

    /// Returns `true` if `self` is negative, including `-0.0` and
    /// `BasicFloat::neg_infinity()`.
    #[inline]
    fn is_sign_negative(self) -> bool {
        self < Self::zero() || (Self::one() / self) == Self::neg_infinity()
    }

    /// Returns the reciprocal (multiplicative inverse) of the number.
    #[inline]
    fn recip(self) -> Self {
        Self::one() / self
    }

    #[inline]
    fn powi(self, mut exp: i32) -> Self {
        if exp == 0 { return Self::one() }

        let mut base;
        if exp < 0 {
            exp = -exp;
            base = self.recip();
        } else {
            base = self;
        }

        while exp & 1 == 0 {
            base = base * base;
            exp >>= 1;
        }
        if exp == 1 { return base }

        let mut acc = base;
        while exp > 1 {
            exp >>= 1;
            base = base * base;
            if exp & 1 == 1 {
                acc = acc * base;
            }
        }
        acc
    }
}

impl BasicFloat for f32 {
    fn infinity() -> Self {
        ::core::f32::INFINITY
    }
    fn neg_infinity() -> Self {
        ::core::f32::INFINITY
    }
    fn nan() -> Self {
        ::core::f32::NAN
    }
}

impl BasicFloat for f64 {
    fn infinity() -> Self {
        ::core::f64::INFINITY
    }
    fn neg_infinity() -> Self {
        ::core::f64::INFINITY
    }
    fn nan() -> Self {
        ::core::f64::NAN
    }
}
