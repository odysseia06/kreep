use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::algebra::field::Field;
use crate::algebra::ring::Ring;
use crate::utils::is_prime;

/// Prime field GF(p) where `p` is a `u64`-sized modulus.
#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct Fp<const P: u64> {
    value: u64,
}

impl<const P: u64> Fp<P> {
    /// Create a new field element, automatically reduced modulo `P`.
    ///
    /// In debug builds, this asserts that `P` is prime.
    pub fn new(value: u64) -> Self {
        debug_assert!(is_prime(P), "Fp modulus P={} is not prime", P);
        Self { value: value % P }
    }

    /// Get the representative in `[0, P-1]`.
    pub const fn value(self) -> u64 {
        self.value
    }

    /// The modulus `p`.
    pub const fn modulus() -> u64 {
        P
    }

    /// Validate that the modulus `P` is prime.
    ///
    /// Returns `Ok(())` if `P` is prime, or an error message otherwise.
    /// Call this at application startup for early failure on misconfiguration.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::Fp;
    ///
    /// type F17 = Fp<17>;
    /// assert!(F17::validate_prime().is_ok());
    ///
    /// type F15 = Fp<15>;
    /// assert!(F15::validate_prime().is_err());
    /// ```
    pub const fn validate_prime() -> Result<(), &'static str> {
        if !is_prime(P) {
            return Err("modulus P is not prime");
        }
        Ok(())
    }
}

impl<const P: u64> fmt::Debug for Fp<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Fp<{}>({})", P, self.value)
    }
}

/* ---- standard arithmetic operators ---- */

impl<const P: u64> Add for Fp<P> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        let mut sum = self.value + rhs.value;
        if sum >= P {
            sum -= P;
        }
        Self { value: sum }
    }
}

impl<const P: u64> Sub for Fp<P> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        // self - rhs ≡ self + (P - rhs)  (mod P)
        if self.value >= rhs.value {
            Self {
                value: self.value - rhs.value,
            }
        } else {
            Self {
                value: self.value + P - rhs.value,
            }
        }
    }
}

impl<const P: u64> Mul for Fp<P> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        // Use u128 to avoid overflow for moderately sized P.
        let prod = (self.value as u128) * (rhs.value as u128);
        let m = P as u128;
        Self {
            value: (prod % m) as u64,
        }
    }
}

impl<const P: u64> Neg for Fp<P> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        if self.value == 0 {
            self
        } else {
            Self {
                value: P - self.value,
            }
        }
    }
}

/// Division implemented via multiplicative inverse.
impl<const P: u64> Div for Fp<P> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inverse().expect("division by zero in Fp")
    }
}

/* ---- implement Ring ---- */

impl<const P: u64> Ring for Fp<P> {
    const ZERO: Self = Self { value: 0 };
    const ONE: Self = Self { value: 1 };
}

/* ---- implement Field ---- */

impl<const P: u64> Field for Fp<P> {
    fn inverse(self) -> Option<Self> {
        if self.value == 0 {
            return None;
        }

        // Extended Euclidean algorithm to find x such that a*x ≡ 1 (mod P)
        let a = self.value as i128;
        let m = P as i128;

        let (g, x, _) = egcd(a, m);
        if g != 1 {
            // Not invertible: gcd(a, P) != 1 (this can happen if P is not prime)
            return None;
        }

        // x may be negative, so bring it into [0, P-1]
        let mut x = x % m;
        if x < 0 {
            x += m;
        }
        Some(Self::new(x as u64))
    }
}

/* ---- internal helper: extended Euclidean algorithm ---- */

/// Returns `(g, x, y)` such that `g = gcd(a, b)` and `a*x + b*y = g`.
fn egcd(a: i128, b: i128) -> (i128, i128, i128) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x1, y1) = egcd(b, a % b);
        (g, y1, x1 - (a / b) * y1)
    }
}

/* ---- basic tests ---- */

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::field::Field;
    use crate::algebra::ring::Ring;

    type F17 = Fp<17>;

    #[test]
    fn add_basic() {
        let a = F17::new(5);
        let b = F17::new(13);
        assert_eq!((a + b).value(), 1); // 5 + 13 = 18 ≡ 1 (mod 17)
    }

    #[test]
    fn mul_basic() {
        let a = F17::new(3);
        let b = F17::new(7);
        assert_eq!((a * b).value(), (21 % 17) as u64);
    }

    #[test]
    fn neg_basic() {
        let a = F17::new(3);
        let minus_a = -a;
        assert_eq!((a + minus_a).value(), F17::ZERO.value());
    }

    #[test]
    fn inverse_exists_for_nonzero() {
        for x in 1u64..17 {
            let a = F17::new(x);
            let inv = a.inverse().expect("should be invertible in prime field");
            let prod = a * inv;
            assert_eq!(prod, F17::ONE);
        }
    }

    #[test]
    fn inverse_none_for_zero() {
        let zero = F17::ZERO;
        assert!(zero.inverse().is_none());
    }

    #[test]
    fn validate_prime_ok() {
        assert!(Fp::<2>::validate_prime().is_ok());
        assert!(Fp::<3>::validate_prime().is_ok());
        assert!(Fp::<17>::validate_prime().is_ok());
        assert!(Fp::<101>::validate_prime().is_ok());
    }

    #[test]
    fn validate_prime_err() {
        assert!(Fp::<1>::validate_prime().is_err());
        assert!(Fp::<4>::validate_prime().is_err());
        assert!(Fp::<15>::validate_prime().is_err());
        assert!(Fp::<100>::validate_prime().is_err());
    }
}
