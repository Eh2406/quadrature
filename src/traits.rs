use core::mem;
use core::num::FpCategory;
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

    /// Returns `true` if the number is neither zero, infinite, subnormal or NaN.
    #[inline]
    fn is_normal(self) -> bool {
        self.classify() == FpCategory::Normal
    }

    /// Returns the floating point category of the number. If only one property
    /// is going to be tested, it is generally faster to use the specific
    /// predicate instead.
    #[inline]
    fn classify(self) -> FpCategory;

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

    /// Converts to degrees, assuming the number is in radians.
    #[inline]
    fn to_degrees(self) -> Self;

    /// Converts to radians, assuming the number is in degrees.
    #[inline]
    fn to_radians(self) -> Self;
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
    fn classify(self) -> FpCategory {
        const EXP_MASK: u32 = 0x7f800000;
        const MAN_MASK: u32 = 0x007fffff;

        let bits: u32 = unsafe { mem::transmute(self) };
        match (bits & MAN_MASK, bits & EXP_MASK) {
            (0, 0) => FpCategory::Zero,
            (_, 0) => FpCategory::Subnormal,
            (0, EXP_MASK) => FpCategory::Infinite,
            (_, EXP_MASK) => FpCategory::Nan,
            _ => FpCategory::Normal,
        }
    }
    fn to_degrees(self) -> Self {
        self * (180.0 / ::core::f32::consts::PI)
    }
    fn to_radians(self) -> Self {
        self * (::core::f32::consts::PI / 180.0)
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
    fn classify(self) -> FpCategory {
        const EXP_MASK: u64 = 0x7ff0000000000000;
        const MAN_MASK: u64 = 0x000fffffffffffff;

        let bits: u64 = unsafe { mem::transmute(self) };
        match (bits & MAN_MASK, bits & EXP_MASK) {
            (0, 0) => FpCategory::Zero,
            (_, 0) => FpCategory::Subnormal,
            (0, EXP_MASK) => FpCategory::Infinite,
            (_, EXP_MASK) => FpCategory::Nan,
            _ => FpCategory::Normal,
        }
    }
    fn to_degrees(self) -> Self {
        self * (180.0 / ::core::f64::consts::PI)
    }
    fn to_radians(self) -> Self {
        self * (::core::f64::consts::PI / 180.0)
    }
}
