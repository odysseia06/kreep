use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};

use crate::algebra::ring::Ring;
use crate::structures::fp::Fp;

/// Polynomial over a prime field Fp.
///
/// Coefficients are stored in ascending order of degree:
/// `coeffs[i]` is the coefficient of `x^i`.
///
/// The zero polynomial is represented as an empty coefficient vector.
#[derive(Clone, PartialEq, Eq)]
pub struct Poly<const P: u64> {
    coeffs: Vec<Fp<P>>,
}

impl<const P: u64> Poly<P> {
    /// Create a polynomial from coefficients in ascending order.
    ///
    /// `coeffs[i]` is the coefficient of `x^i`.
    /// Trailing zeros are automatically removed.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // 3 + 2x + x^2
    /// let p = Poly::new(vec![F17::new(3), F17::new(2), F17::new(1)]);
    /// assert_eq!(p.degree(), Some(2));
    /// ```
    pub fn new(coeffs: Vec<Fp<P>>) -> Self {
        let mut poly = Self { coeffs };
        poly.normalize();
        poly
    }

    /// Create the zero polynomial.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let zero = Poly::<17>::zero();
    /// assert!(zero.is_zero());
    /// assert_eq!(zero.degree(), None);
    /// ```
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    /// Create a constant polynomial.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let c = Poly::constant(F17::new(5));
    /// assert_eq!(c.degree(), Some(0));
    /// assert_eq!(c.eval(F17::new(10)), F17::new(5));
    /// ```
    pub fn constant(c: Fp<P>) -> Self {
        if c == Fp::ZERO {
            Self::zero()
        } else {
            Self { coeffs: vec![c] }
        }
    }

    /// Create the polynomial `x`.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let x = Poly::<17>::x();
    /// assert_eq!(x.degree(), Some(1));
    /// assert_eq!(x.eval(F17::new(5)), F17::new(5));
    /// ```
    pub fn x() -> Self {
        Self {
            coeffs: vec![Fp::ZERO, Fp::ONE],
        }
    }

    /// Create a monomial `c * x^n`.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let m = Poly::monomial(F17::new(3), 2);  // 3x^2
    /// assert_eq!(m.degree(), Some(2));
    /// assert_eq!(m.eval(F17::new(2)), F17::new(12));  // 3 * 4 = 12
    /// ```
    pub fn monomial(c: Fp<P>, n: usize) -> Self {
        if c == Fp::ZERO {
            return Self::zero();
        }
        let mut coeffs = vec![Fp::ZERO; n + 1];
        coeffs[n] = c;
        Self { coeffs }
    }

    /// Check if this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Get the degree of the polynomial.
    ///
    /// Returns `None` for the zero polynomial, `Some(n)` otherwise
    /// where `n` is the highest power with a non-zero coefficient.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let p = Poly::new(vec![F17::new(1), F17::new(2), F17::new(3)]);
    /// assert_eq!(p.degree(), Some(2));
    ///
    /// let zero = Poly::<17>::zero();
    /// assert_eq!(zero.degree(), None);
    /// ```
    pub fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    /// Get the leading coefficient.
    ///
    /// Returns `None` for the zero polynomial.
    pub fn leading_coeff(&self) -> Option<Fp<P>> {
        self.coeffs.last().copied()
    }

    /// Get the coefficient of `x^i`.
    ///
    /// Returns zero if `i` is beyond the polynomial's degree.
    pub fn coeff(&self, i: usize) -> Fp<P> {
        self.coeffs.get(i).copied().unwrap_or(Fp::ZERO)
    }

    /// Get a slice of all coefficients.
    pub fn coefficients(&self) -> &[Fp<P>] {
        &self.coeffs
    }

    /// Evaluate the polynomial at a point using Horner's method.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // p(x) = 1 + 2x + 3x^2
    /// let p = Poly::new(vec![F17::new(1), F17::new(2), F17::new(3)]);
    ///
    /// // p(2) = 1 + 4 + 12 = 17 ≡ 0 (mod 17)
    /// assert_eq!(p.eval(F17::new(2)), F17::new(0));
    ///
    /// // p(1) = 1 + 2 + 3 = 6
    /// assert_eq!(p.eval(F17::new(1)), F17::new(6));
    /// ```
    pub fn eval(&self, x: Fp<P>) -> Fp<P> {
        if self.coeffs.is_empty() {
            return Fp::ZERO;
        }

        // Horner's method: p(x) = a_0 + x(a_1 + x(a_2 + ... + x*a_n))
        let mut result = Fp::ZERO;
        for &coeff in self.coeffs.iter().rev() {
            result = result * x + coeff;
        }
        result
    }

    /// Remove trailing zero coefficients.
    fn normalize(&mut self) {
        while self.coeffs.last() == Some(&Fp::ZERO) {
            self.coeffs.pop();
        }
    }

    /// Create a polynomial from its roots: `(x - r1)(x - r2)...(x - rn)`.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // (x - 2)(x - 5) = x^2 - 7x + 10
    /// let p = Poly::from_roots(&[F17::new(2), F17::new(5)]);
    ///
    /// // Roots should evaluate to zero
    /// assert_eq!(p.eval(F17::new(2)), F17::new(0));
    /// assert_eq!(p.eval(F17::new(5)), F17::new(0));
    /// ```
    pub fn from_roots(roots: &[Fp<P>]) -> Self {
        if roots.is_empty() {
            return Self::constant(Fp::ONE);
        }

        let mut result = Self::new(vec![-roots[0], Fp::ONE]); // (x - r0)

        for &root in &roots[1..] {
            let factor = Self::new(vec![-root, Fp::ONE]); // (x - ri)
            result = result * factor;
        }

        result
    }

    /// Lagrange interpolation: find the unique polynomial of degree < n
    /// passing through the given points.
    ///
    /// Returns `None` if any x-coordinates are duplicated.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // Find polynomial passing through (0, 1), (1, 3), (2, 7)
    /// let points = [
    ///     (F17::new(0), F17::new(1)),
    ///     (F17::new(1), F17::new(3)),
    ///     (F17::new(2), F17::new(7)),
    /// ];
    /// let p = Poly::interpolate(&points).unwrap();
    ///
    /// // Verify it passes through all points
    /// for (x, y) in &points {
    ///     assert_eq!(p.eval(*x), *y);
    /// }
    /// ```
    pub fn interpolate(points: &[(Fp<P>, Fp<P>)]) -> Option<Self> {
        use crate::algebra::field::Field;

        let n = points.len();
        if n == 0 {
            return Some(Self::zero());
        }

        // Check for duplicate x-coordinates
        for i in 0..n {
            for j in (i + 1)..n {
                if points[i].0 == points[j].0 {
                    return None;
                }
            }
        }

        let mut result = Self::zero();

        for i in 0..n {
            let (xi, yi) = points[i];

            // Build the Lagrange basis polynomial L_i(x)
            // L_i(x) = product over j != i of (x - x_j) / (x_i - x_j)
            let mut basis = Self::constant(Fp::ONE);
            let mut denom = Fp::ONE;

            for j in 0..n {
                if i != j {
                    let xj = points[j].0;
                    // Multiply by (x - x_j)
                    basis = basis * Self::new(vec![-xj, Fp::ONE]);
                    // Accumulate denominator (x_i - x_j)
                    denom = denom * (xi - xj);
                }
            }

            // Divide by denominator and multiply by y_i
            let coeff = yi * denom.inverse()?;
            result = result + basis * coeff;
        }

        Some(result)
    }
}

/* ---- Arithmetic operators ---- */

impl<const P: u64> Add for Poly<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        let mut coeffs = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let a = self.coeff(i);
            let b = rhs.coeff(i);
            coeffs.push(a + b);
        }

        Self::new(coeffs)
    }
}

impl<const P: u64> Add<&Poly<P>> for Poly<P> {
    type Output = Self;

    fn add(self, rhs: &Poly<P>) -> Self::Output {
        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        let mut coeffs = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let a = self.coeff(i);
            let b = rhs.coeff(i);
            coeffs.push(a + b);
        }

        Self::new(coeffs)
    }
}

impl<const P: u64> Neg for Poly<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let coeffs = self.coeffs.into_iter().map(|c| -c).collect();
        Self { coeffs }
    }
}

impl<const P: u64> Sub for Poly<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<const P: u64> Sub<&Poly<P>> for Poly<P> {
    type Output = Self;

    fn sub(self, rhs: &Poly<P>) -> Self::Output {
        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        let mut coeffs = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let a = self.coeff(i);
            let b = rhs.coeff(i);
            coeffs.push(a - b);
        }

        Self::new(coeffs)
    }
}

impl<const P: u64> Mul for Poly<P> {
    type Output = Self;

    /// Polynomial multiplication using naive O(n*m) convolution.
    fn mul(self, rhs: Self) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }

        let n = self.coeffs.len();
        let m = rhs.coeffs.len();
        let mut coeffs = vec![Fp::ZERO; n + m - 1];

        for i in 0..n {
            for j in 0..m {
                coeffs[i + j] = coeffs[i + j] + self.coeffs[i] * rhs.coeffs[j];
            }
        }

        Self::new(coeffs)
    }
}

impl<const P: u64> Mul<&Poly<P>> for Poly<P> {
    type Output = Self;

    fn mul(self, rhs: &Poly<P>) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }

        let n = self.coeffs.len();
        let m = rhs.coeffs.len();
        let mut coeffs = vec![Fp::ZERO; n + m - 1];

        for i in 0..n {
            for j in 0..m {
                coeffs[i + j] = coeffs[i + j] + self.coeffs[i] * rhs.coeffs[j];
            }
        }

        Self::new(coeffs)
    }
}

/// Scalar multiplication: polynomial * field element
impl<const P: u64> Mul<Fp<P>> for Poly<P> {
    type Output = Self;

    fn mul(self, rhs: Fp<P>) -> Self::Output {
        if rhs == Fp::ZERO {
            return Self::zero();
        }
        let coeffs = self.coeffs.into_iter().map(|c| c * rhs).collect();
        Self::new(coeffs)
    }
}

impl<const P: u64> fmt::Debug for Poly<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        let mut first = true;
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            if coeff == Fp::ZERO {
                continue;
            }

            if !first {
                write!(f, " + ")?;
            }
            first = false;

            match i {
                0 => write!(f, "{}", coeff.value())?,
                1 if coeff == Fp::ONE => write!(f, "x")?,
                1 => write!(f, "{}*x", coeff.value())?,
                _ if coeff == Fp::ONE => write!(f, "x^{}", i)?,
                _ => write!(f, "{}*x^{}", coeff.value(), i)?,
            }
        }

        Ok(())
    }
}

impl<const P: u64> fmt::Display for Poly<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type F17 = Fp<17>;
    type P17 = Poly<17>;

    #[test]
    fn new_normalizes() {
        // Trailing zeros should be removed
        let p = P17::new(vec![F17::new(1), F17::new(2), F17::ZERO, F17::ZERO]);
        assert_eq!(p.degree(), Some(1));
        assert_eq!(p.coefficients().len(), 2);
    }

    #[test]
    fn new_all_zeros() {
        let p = P17::new(vec![F17::ZERO, F17::ZERO]);
        assert!(p.is_zero());
        assert_eq!(p.degree(), None);
    }

    #[test]
    fn zero() {
        let z = P17::zero();
        assert!(z.is_zero());
        assert_eq!(z.degree(), None);
        assert_eq!(z.leading_coeff(), None);
    }

    #[test]
    fn constant() {
        let c = P17::constant(F17::new(5));
        assert_eq!(c.degree(), Some(0));
        assert_eq!(c.leading_coeff(), Some(F17::new(5)));
        assert_eq!(c.eval(F17::new(100)), F17::new(5));
    }

    #[test]
    fn constant_zero() {
        let c = P17::constant(F17::ZERO);
        assert!(c.is_zero());
    }

    #[test]
    fn x_polynomial() {
        let x = P17::x();
        assert_eq!(x.degree(), Some(1));
        assert_eq!(x.eval(F17::new(5)), F17::new(5));
        assert_eq!(x.eval(F17::new(0)), F17::ZERO);
    }

    #[test]
    fn monomial() {
        let m = P17::monomial(F17::new(3), 2); // 3x^2
        assert_eq!(m.degree(), Some(2));
        assert_eq!(m.coeff(0), F17::ZERO);
        assert_eq!(m.coeff(1), F17::ZERO);
        assert_eq!(m.coeff(2), F17::new(3));
    }

    #[test]
    fn monomial_zero_coeff() {
        let m = P17::monomial(F17::ZERO, 5);
        assert!(m.is_zero());
    }

    #[test]
    fn coeff_out_of_range() {
        let p = P17::new(vec![F17::new(1), F17::new(2)]);
        assert_eq!(p.coeff(0), F17::new(1));
        assert_eq!(p.coeff(1), F17::new(2));
        assert_eq!(p.coeff(2), F17::ZERO);
        assert_eq!(p.coeff(100), F17::ZERO);
    }

    #[test]
    fn eval_zero_poly() {
        let z = P17::zero();
        assert_eq!(z.eval(F17::new(5)), F17::ZERO);
    }

    #[test]
    fn eval_constant() {
        let c = P17::constant(F17::new(7));
        assert_eq!(c.eval(F17::new(0)), F17::new(7));
        assert_eq!(c.eval(F17::new(100)), F17::new(7));
    }

    #[test]
    fn eval_linear() {
        // p(x) = 2 + 3x
        let p = P17::new(vec![F17::new(2), F17::new(3)]);
        assert_eq!(p.eval(F17::new(0)), F17::new(2));
        assert_eq!(p.eval(F17::new(1)), F17::new(5));
        assert_eq!(p.eval(F17::new(5)), F17::ZERO); // 2 + 15 = 17 ≡ 0
    }

    #[test]
    fn eval_quadratic() {
        // p(x) = 1 + 2x + 3x^2
        let p = P17::new(vec![F17::new(1), F17::new(2), F17::new(3)]);
        // p(0) = 1
        assert_eq!(p.eval(F17::new(0)), F17::new(1));
        // p(1) = 1 + 2 + 3 = 6
        assert_eq!(p.eval(F17::new(1)), F17::new(6));
        // p(2) = 1 + 4 + 12 = 17 ≡ 0
        assert_eq!(p.eval(F17::new(2)), F17::ZERO);
    }

    #[test]
    fn debug_format() {
        let p = P17::new(vec![F17::new(1), F17::new(2), F17::new(3)]);
        let s = format!("{:?}", p);
        assert_eq!(s, "1 + 2*x + 3*x^2");
    }

    #[test]
    fn debug_format_with_one_coeff() {
        let p = P17::new(vec![F17::ZERO, F17::ONE, F17::ONE]);
        let s = format!("{:?}", p);
        assert_eq!(s, "x + x^2");
    }

    #[test]
    fn debug_format_zero() {
        let z = P17::zero();
        assert_eq!(format!("{:?}", z), "0");
    }

    // ---- Arithmetic tests ----

    #[test]
    fn add_basic() {
        // (1 + 2x) + (3 + 4x) = 4 + 6x
        let p1 = P17::new(vec![F17::new(1), F17::new(2)]);
        let p2 = P17::new(vec![F17::new(3), F17::new(4)]);
        let sum = p1 + p2;
        assert_eq!(sum.coeff(0), F17::new(4));
        assert_eq!(sum.coeff(1), F17::new(6));
    }

    #[test]
    fn add_different_degrees() {
        // (1 + 2x) + (3x^2) = 1 + 2x + 3x^2
        let p1 = P17::new(vec![F17::new(1), F17::new(2)]);
        let p2 = P17::monomial(F17::new(3), 2);
        let sum = p1 + p2;
        assert_eq!(sum.degree(), Some(2));
        assert_eq!(sum.coeff(0), F17::new(1));
        assert_eq!(sum.coeff(1), F17::new(2));
        assert_eq!(sum.coeff(2), F17::new(3));
    }

    #[test]
    fn add_with_zero() {
        let p = P17::new(vec![F17::new(1), F17::new(2)]);
        let z = P17::zero();
        assert_eq!(p.clone() + z, p);
    }

    #[test]
    fn add_cancellation() {
        // (1 + 2x) + (16 + 15x) = 0 (mod 17)
        let p1 = P17::new(vec![F17::new(1), F17::new(2)]);
        let p2 = P17::new(vec![F17::new(16), F17::new(15)]);
        let sum = p1 + p2;
        assert!(sum.is_zero());
    }

    #[test]
    fn neg_basic() {
        let p = P17::new(vec![F17::new(1), F17::new(2)]);
        let neg_p = -p.clone();
        assert_eq!(neg_p.coeff(0), F17::new(16)); // -1 mod 17
        assert_eq!(neg_p.coeff(1), F17::new(15)); // -2 mod 17
    }

    #[test]
    fn neg_zero() {
        let z = P17::zero();
        assert!((-z).is_zero());
    }

    #[test]
    fn sub_basic() {
        // (5 + 3x) - (2 + x) = 3 + 2x
        let p1 = P17::new(vec![F17::new(5), F17::new(3)]);
        let p2 = P17::new(vec![F17::new(2), F17::new(1)]);
        let diff = p1 - p2;
        assert_eq!(diff.coeff(0), F17::new(3));
        assert_eq!(diff.coeff(1), F17::new(2));
    }

    #[test]
    fn sub_self_is_zero() {
        let p = P17::new(vec![F17::new(1), F17::new(2), F17::new(3)]);
        let diff = p.clone() - p;
        assert!(diff.is_zero());
    }

    #[test]
    fn mul_constants() {
        let p1 = P17::constant(F17::new(3));
        let p2 = P17::constant(F17::new(5));
        let prod = p1 * p2;
        assert_eq!(prod.degree(), Some(0));
        assert_eq!(prod.coeff(0), F17::new(15));
    }

    #[test]
    fn mul_by_x() {
        // (1 + 2x) * x = x + 2x^2
        let p = P17::new(vec![F17::new(1), F17::new(2)]);
        let x = P17::x();
        let prod = p * x;
        assert_eq!(prod.degree(), Some(2));
        assert_eq!(prod.coeff(0), F17::ZERO);
        assert_eq!(prod.coeff(1), F17::new(1));
        assert_eq!(prod.coeff(2), F17::new(2));
    }

    #[test]
    fn mul_linear() {
        // (1 + x) * (2 + 3x) = 2 + 5x + 3x^2
        let p1 = P17::new(vec![F17::new(1), F17::new(1)]);
        let p2 = P17::new(vec![F17::new(2), F17::new(3)]);
        let prod = p1 * p2;
        assert_eq!(prod.degree(), Some(2));
        assert_eq!(prod.coeff(0), F17::new(2));
        assert_eq!(prod.coeff(1), F17::new(5));
        assert_eq!(prod.coeff(2), F17::new(3));
    }

    #[test]
    fn mul_by_zero_poly() {
        let p = P17::new(vec![F17::new(1), F17::new(2)]);
        let z = P17::zero();
        assert!((p * z).is_zero());
    }

    #[test]
    fn mul_scalar() {
        // (1 + 2x) * 3 = 3 + 6x
        let p = P17::new(vec![F17::new(1), F17::new(2)]);
        let prod = p * F17::new(3);
        assert_eq!(prod.coeff(0), F17::new(3));
        assert_eq!(prod.coeff(1), F17::new(6));
    }

    #[test]
    fn mul_scalar_zero() {
        let p = P17::new(vec![F17::new(1), F17::new(2)]);
        let prod = p * F17::ZERO;
        assert!(prod.is_zero());
    }

    #[test]
    fn mul_degree_sum() {
        // degree(p * q) = degree(p) + degree(q)
        let p = P17::new(vec![F17::new(1), F17::new(2), F17::new(3)]); // degree 2
        let q = P17::new(vec![F17::new(4), F17::new(5)]); // degree 1
        let prod = p * q;
        assert_eq!(prod.degree(), Some(3));
    }

    #[test]
    fn arithmetic_eval_consistency() {
        // (p + q)(x) = p(x) + q(x)
        // (p * q)(x) = p(x) * q(x)
        let p = P17::new(vec![F17::new(1), F17::new(2)]);
        let q = P17::new(vec![F17::new(3), F17::new(4), F17::new(5)]);
        let x = F17::new(7);

        assert_eq!((p.clone() + q.clone()).eval(x), p.eval(x) + q.eval(x));
        assert_eq!((p.clone() * q.clone()).eval(x), p.eval(x) * q.eval(x));
    }

    // ---- from_roots tests ----

    #[test]
    fn from_roots_empty() {
        let p = P17::from_roots(&[]);
        assert_eq!(p.degree(), Some(0));
        assert_eq!(p.coeff(0), F17::ONE);
    }

    #[test]
    fn from_roots_single() {
        // (x - 3)
        let p = P17::from_roots(&[F17::new(3)]);
        assert_eq!(p.degree(), Some(1));
        assert_eq!(p.eval(F17::new(3)), F17::ZERO);
        assert_eq!(p.eval(F17::new(0)), F17::new(14)); // -3 mod 17
    }

    #[test]
    fn from_roots_two() {
        // (x - 2)(x - 5) should have roots at 2 and 5
        let p = P17::from_roots(&[F17::new(2), F17::new(5)]);
        assert_eq!(p.degree(), Some(2));
        assert_eq!(p.eval(F17::new(2)), F17::ZERO);
        assert_eq!(p.eval(F17::new(5)), F17::ZERO);
    }

    #[test]
    fn from_roots_multiple() {
        let roots = [F17::new(1), F17::new(3), F17::new(7), F17::new(11)];
        let p = P17::from_roots(&roots);
        assert_eq!(p.degree(), Some(4));

        for &r in &roots {
            assert_eq!(
                p.eval(r),
                F17::ZERO,
                "root {} should evaluate to zero",
                r.value()
            );
        }
    }

    #[test]
    fn from_roots_duplicate() {
        // (x - 2)^2 = x^2 - 4x + 4
        let p = P17::from_roots(&[F17::new(2), F17::new(2)]);
        assert_eq!(p.degree(), Some(2));
        assert_eq!(p.eval(F17::new(2)), F17::ZERO);
    }

    // ---- interpolate tests ----

    #[test]
    fn interpolate_empty() {
        let p = P17::interpolate(&[]).unwrap();
        assert!(p.is_zero());
    }

    #[test]
    fn interpolate_single_point() {
        // Constant polynomial through (3, 7)
        let p = P17::interpolate(&[(F17::new(3), F17::new(7))]).unwrap();
        assert_eq!(p.degree(), Some(0));
        assert_eq!(p.eval(F17::new(3)), F17::new(7));
        assert_eq!(p.eval(F17::new(0)), F17::new(7));
    }

    #[test]
    fn interpolate_two_points() {
        // Line through (0, 1) and (1, 3) -> y = 1 + 2x
        let p =
            P17::interpolate(&[(F17::new(0), F17::new(1)), (F17::new(1), F17::new(3))]).unwrap();

        assert_eq!(p.degree(), Some(1));
        assert_eq!(p.eval(F17::new(0)), F17::new(1));
        assert_eq!(p.eval(F17::new(1)), F17::new(3));
    }

    #[test]
    fn interpolate_three_points() {
        // Quadratic through (0, 1), (1, 2), (2, 5)
        // p(x) = 1 + 0x + x^2 = 1 + x^2
        let points = [
            (F17::new(0), F17::new(1)),
            (F17::new(1), F17::new(2)),
            (F17::new(2), F17::new(5)),
        ];
        let p = P17::interpolate(&points).unwrap();

        assert_eq!(p.degree(), Some(2));
        for (x, y) in &points {
            assert_eq!(p.eval(*x), *y);
        }
    }

    #[test]
    fn interpolate_all_zeros() {
        // Zero polynomial through (1, 0), (2, 0), (3, 0)
        let points = [
            (F17::new(1), F17::ZERO),
            (F17::new(2), F17::ZERO),
            (F17::new(3), F17::ZERO),
        ];
        let p = P17::interpolate(&points).unwrap();

        // Should be zero polynomial
        assert!(p.is_zero());
    }

    #[test]
    fn interpolate_duplicate_x_fails() {
        let points = [
            (F17::new(1), F17::new(2)),
            (F17::new(1), F17::new(3)), // duplicate x
        ];
        assert!(P17::interpolate(&points).is_none());
    }

    #[test]
    fn interpolate_roundtrip_with_from_roots() {
        // Create polynomial from roots, then interpolate points on it
        let roots = [F17::new(2), F17::new(5), F17::new(11)];
        let p = P17::from_roots(&roots);

        // Need n+1 points to recover a degree-n polynomial
        // p has degree 3, so we need 4 points
        let points: Vec<_> = (0..4)
            .map(|i| {
                let x = F17::new(i);
                (x, p.eval(x))
            })
            .collect();

        let q = P17::interpolate(&points).unwrap();

        // Both polynomials should agree on all field elements
        for i in 0..17 {
            let x = F17::new(i);
            assert_eq!(p.eval(x), q.eval(x), "mismatch at x={}", i);
        }
    }
}
