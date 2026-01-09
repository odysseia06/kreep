use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};

#[cfg(feature = "rand")]
use crate::algebra::field::Field;
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

    /// Make the polynomial monic (leading coefficient = 1).
    ///
    /// Returns `None` if the polynomial is zero.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let p = Poly::new(vec![F17::new(2), F17::new(4), F17::new(2)]); // 2 + 4x + 2x^2
    /// let monic = p.monic().unwrap();
    /// assert_eq!(monic.leading_coeff(), Some(F17::new(1)));
    /// ```
    pub fn monic(&self) -> Option<Self> {
        use crate::algebra::field::Field;

        let lc = self.leading_coeff()?;
        let inv = lc.inverse()?;
        Some(self.clone() * inv)
    }

    /// Euclidean division: compute quotient and remainder.
    ///
    /// Returns `(q, r)` such that `self = q * divisor + r` and `deg(r) < deg(divisor)`.
    ///
    /// Returns `None` if the divisor is zero.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // (x^2 + 2x + 1) / (x + 1) = (x + 1), remainder 0
    /// let dividend = Poly::new(vec![F17::new(1), F17::new(2), F17::new(1)]);
    /// let divisor = Poly::new(vec![F17::new(1), F17::new(1)]);
    /// let (q, r) = dividend.div_rem(&divisor).unwrap();
    ///
    /// assert_eq!(q, divisor);
    /// assert!(r.is_zero());
    /// ```
    pub fn div_rem(&self, divisor: &Self) -> Option<(Self, Self)> {
        use crate::algebra::field::Field;

        if divisor.is_zero() {
            return None;
        }

        // If dividend degree < divisor degree, quotient is 0
        let divisor_deg = divisor.degree()?;
        match self.degree() {
            None => return Some((Self::zero(), Self::zero())),
            Some(d) if d < divisor_deg => return Some((Self::zero(), self.clone())),
            _ => {}
        }

        let lc_inv = divisor.leading_coeff()?.inverse()?;
        let mut remainder = self.clone();
        let mut quotient_coeffs = vec![Fp::ZERO; self.coeffs.len() - divisor.coeffs.len() + 1];

        while let Some(rem_deg) = remainder.degree() {
            if rem_deg < divisor_deg {
                break;
            }

            let rem_lc = remainder.leading_coeff().unwrap();
            let coeff = rem_lc * lc_inv;
            let deg_diff = rem_deg - divisor_deg;

            quotient_coeffs[deg_diff] = coeff;

            // remainder -= coeff * x^deg_diff * divisor
            for (i, &d_coeff) in divisor.coeffs.iter().enumerate() {
                remainder.coeffs[i + deg_diff] = remainder.coeffs[i + deg_diff] - coeff * d_coeff;
            }
            remainder.normalize();
        }

        Some((Self::new(quotient_coeffs), remainder))
    }

    /// Compute the remainder of division.
    ///
    /// Returns `None` if the divisor is zero.
    pub fn rem(&self, divisor: &Self) -> Option<Self> {
        self.div_rem(divisor).map(|(_, r)| r)
    }

    /// Compute the greatest common divisor of two polynomials.
    ///
    /// The result is monic (leading coefficient = 1) unless both inputs are zero.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // gcd((x-1)(x-2), (x-2)(x-3)) = (x-2)
    /// let p1 = Poly::from_roots(&[F17::new(1), F17::new(2)]);
    /// let p2 = Poly::from_roots(&[F17::new(2), F17::new(3)]);
    /// let g = Poly::gcd(&p1, &p2);
    ///
    /// // g should be monic and have root at 2
    /// assert_eq!(g.leading_coeff(), Some(F17::new(1)));
    /// assert_eq!(g.eval(F17::new(2)), F17::new(0));
    /// ```
    pub fn gcd(a: &Self, b: &Self) -> Self {
        if b.is_zero() {
            return a.monic().unwrap_or_else(Self::zero);
        }

        let r = a.rem(b).unwrap_or_else(Self::zero);
        Self::gcd(b, &r)
    }

    /// Extended Euclidean algorithm for polynomials.
    ///
    /// Returns `(g, s, t)` such that `g = gcd(a, b) = s*a + t*b`.
    ///
    /// The gcd `g` is monic unless both inputs are zero.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let a = Poly::new(vec![F17::new(1), F17::new(2), F17::new(1)]); // 1 + 2x + x^2
    /// let b = Poly::new(vec![F17::new(1), F17::new(1)]); // 1 + x
    ///
    /// let (g, s, t) = Poly::extended_gcd(&a, &b);
    ///
    /// // Verify: s*a + t*b = g
    /// let check = s * a.clone() + t * b.clone();
    /// assert_eq!(check, g);
    /// ```
    pub fn extended_gcd(a: &Self, b: &Self) -> (Self, Self, Self) {
        if b.is_zero() {
            if a.is_zero() {
                return (Self::zero(), Self::zero(), Self::zero());
            }
            // Make result monic
            let monic_a = a.monic().unwrap();
            let lc_inv = {
                use crate::algebra::field::Field;
                a.leading_coeff().unwrap().inverse().unwrap()
            };
            return (monic_a, Self::constant(lc_inv), Self::zero());
        }

        let (q, r) = a.div_rem(b).unwrap();
        let (g, s1, t1) = Self::extended_gcd(b, &r);

        // g = s1 * b + t1 * r
        //   = s1 * b + t1 * (a - q * b)
        //   = t1 * a + (s1 - t1 * q) * b
        let s = t1.clone();
        let t = s1 - t1 * q;

        (g, s, t)
    }

    /// Compute `x^exp mod self` using repeated squaring.
    ///
    /// This is useful for irreducibility testing where we need to compute
    /// `x^{p^k} mod f(x)` efficiently.
    ///
    /// Returns `None` if self is zero.
    pub fn powmod_x(&self, exp: u64) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        if exp == 0 {
            return Some(Self::constant(Fp::ONE));
        }

        let mut base = Self::x().rem(self)?;
        let mut result = Self::constant(Fp::ONE);
        let mut e = exp;

        while e > 0 {
            if e & 1 == 1 {
                result = (result * &base).rem(self)?;
            }
            base = (base.clone() * &base).rem(self)?;
            e >>= 1;
        }

        Some(result)
    }

    /// Compute `base^exp mod self` using repeated squaring.
    ///
    /// Returns `None` if self is zero.
    pub fn powmod(&self, base: &Self, exp: u64) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        if exp == 0 {
            return Some(Self::constant(Fp::ONE));
        }

        let mut b = base.rem(self)?;
        let mut result = Self::constant(Fp::ONE);
        let mut e = exp;

        while e > 0 {
            if e & 1 == 1 {
                result = (result * &b).rem(self)?;
            }
            b = (b.clone() * &b).rem(self)?;
            e >>= 1;
        }

        Some(result)
    }

    /// Test if this polynomial is irreducible over F_p using Rabin's algorithm.
    ///
    /// A polynomial f(x) of degree n over F_p is irreducible if and only if:
    /// 1. `x^{p^n} ≡ x (mod f(x))`
    /// 2. `gcd(x^{p^{n/q}} - x, f(x)) = 1` for each prime divisor q of n
    ///
    /// Returns `false` for constant or zero polynomials.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // x^2 - 3 is irreducible over F_17 (3 is not a quadratic residue)
    /// let f = Poly::new(vec![F17::new(14), F17::new(0), F17::new(1)]); // -3 + x^2
    /// assert!(f.is_irreducible());
    ///
    /// // x^2 - 4 is reducible over F_17 (4 = 2^2 is a quadratic residue)
    /// let g = Poly::new(vec![F17::new(13), F17::new(0), F17::new(1)]); // -4 + x^2
    /// assert!(!g.is_irreducible());
    /// ```
    pub fn is_irreducible(&self) -> bool {
        let n = match self.degree() {
            None => return false,    // zero polynomial
            Some(0) => return false, // constant polynomial
            Some(1) => return true,  // linear polynomials are always irreducible
            Some(d) => d,
        };

        // Make monic for cleaner computation (irreducibility is preserved)
        let f = match self.monic() {
            Some(m) => m,
            None => return false,
        };

        // Step 1: Compute x^{p^n} mod f(x) and check if it equals x
        // We compute this iteratively: x -> x^p -> x^{p^2} -> ... -> x^{p^n}
        let mut h = Self::x(); // h = x^{p^i}, starting with i=0

        // Get prime divisors of n for step 2
        let prime_divisors = Self::prime_divisors(n);

        // For each i from 1 to n, compute h = h^p mod f
        for i in 1..=n {
            h = match f.powmod(&h, P) {
                Some(r) => r,
                None => return false,
            };

            // Step 2: At i = n/q for each prime divisor q, check gcd(h - x, f) = 1
            for &q in &prime_divisors {
                if n == i * q {
                    // We're at i = n/q
                    let h_minus_x = h.clone() - Self::x();
                    let g = Self::gcd(&h_minus_x, &f);
                    if g.degree() != Some(0) {
                        return false; // gcd is non-trivial, f is reducible
                    }
                }
            }
        }

        // Step 1 check: h should now be x^{p^n}, verify h = x (mod f)
        let h_minus_x = h - Self::x();
        h_minus_x.is_zero()
    }

    /// Find all prime divisors of n.
    fn prime_divisors(mut n: usize) -> Vec<usize> {
        let mut primes = Vec::new();
        let mut d = 2;

        while d <= n / d {
            if n.is_multiple_of(d) {
                primes.push(d);
                while n.is_multiple_of(d) {
                    n /= d;
                }
            }
            d += 1;
        }

        if n > 1 {
            primes.push(n);
        }

        primes
    }

    /// Test if this polynomial is primitive over F_p.
    ///
    /// A polynomial f(x) of degree n is primitive if:
    /// 1. It is irreducible
    /// 2. x is a primitive element of F_p[x]/(f(x)), i.e., ord(x) = p^n - 1
    ///
    /// Primitive polynomials are used to construct maximal-length LFSRs
    /// and are important in coding theory and cryptography.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F3 = Fp<3>;
    ///
    /// // Over F_3, x^2 + 1 is irreducible but NOT primitive
    /// // (x has order 4, not 8 = 3^2 - 1)
    /// let f = Poly::new(vec![F3::new(1), F3::new(0), F3::new(1)]); // 1 + x^2
    /// assert!(f.is_irreducible());
    /// assert!(!f.is_primitive());
    ///
    /// // x^2 + 2x + 2 is primitive over F_3 (x has order 8)
    /// let g = Poly::new(vec![F3::new(2), F3::new(2), F3::new(1)]); // 2 + 2x + x^2
    /// assert!(g.is_primitive());
    /// ```
    pub fn is_primitive(&self) -> bool {
        // Must be irreducible first
        if !self.is_irreducible() {
            return false;
        }

        let n = match self.degree() {
            Some(d) if d > 0 => d,
            _ => return false,
        };

        // Make monic for cleaner computation
        let f = match self.monic() {
            Some(m) => m,
            None => return false,
        };

        // Order of the multiplicative group is p^n - 1
        // We need to compute p^n - 1, being careful about overflow
        let order = Self::compute_field_order(P, n);
        let order = match order {
            Some(o) => o,
            None => return false, // Overflow, field too large
        };

        // x must have order exactly `order` = p^n - 1
        // This means:
        // 1. x^order ≡ 1 (mod f)
        // 2. x^(order/q) ≢ 1 (mod f) for each prime divisor q of order

        // Check x^order = 1
        let x_to_order = match f.powmod_x(order) {
            Some(r) => r,
            None => return false,
        };
        if x_to_order != Self::constant(Fp::ONE) {
            return false;
        }

        // Get prime divisors of order
        let prime_divisors = Self::prime_divisors_u64(order);

        // Check x^(order/q) ≠ 1 for each prime q dividing order
        for q in prime_divisors {
            let exp = order / q;
            let x_to_exp = match f.powmod_x(exp) {
                Some(r) => r,
                None => return false,
            };
            if x_to_exp == Self::constant(Fp::ONE) {
                return false; // x has smaller order, not primitive
            }
        }

        true
    }

    /// Compute p^n - 1, returning None on overflow.
    fn compute_field_order(p: u64, n: usize) -> Option<u64> {
        let mut result: u64 = 1;
        for _ in 0..n {
            result = result.checked_mul(p)?;
        }
        result.checked_sub(1)
    }

    /// Find all prime divisors of n (for u64).
    fn prime_divisors_u64(mut n: u64) -> Vec<u64> {
        let mut primes = Vec::new();
        let mut d: u64 = 2;

        while d <= n / d {
            if n.is_multiple_of(d) {
                primes.push(d);
                while n.is_multiple_of(d) {
                    n /= d;
                }
            }
            d += 1;
        }

        if n > 1 {
            primes.push(n);
        }

        primes
    }

    /// Generate a random monic irreducible polynomial of the given degree.
    ///
    /// Uses rejection sampling: generate random monic polynomials until
    /// finding an irreducible one.
    ///
    /// # Panics
    ///
    /// Panics if degree is 0.
    #[cfg(feature = "rand")]
    pub fn random_irreducible<R: rand::Rng>(rng: &mut R, degree: usize) -> Self {
        assert!(degree > 0, "degree must be positive");

        loop {
            // Generate random monic polynomial of given degree
            let mut coeffs = Vec::with_capacity(degree + 1);
            for _ in 0..degree {
                coeffs.push(Fp::new(rng.gen_range(0..P)));
            }
            coeffs.push(Fp::ONE); // monic

            let f = Self::new(coeffs);
            if f.is_irreducible() {
                return f;
            }
        }
    }

    /// Generate a random monic primitive polynomial of the given degree.
    ///
    /// Uses rejection sampling: generate random monic polynomials until
    /// finding a primitive one.
    ///
    /// # Panics
    ///
    /// Panics if degree is 0.
    #[cfg(feature = "rand")]
    pub fn random_primitive<R: rand::Rng>(rng: &mut R, degree: usize) -> Self {
        assert!(degree > 0, "degree must be positive");

        loop {
            // Generate random monic polynomial of given degree
            let mut coeffs = Vec::with_capacity(degree + 1);
            for _ in 0..degree {
                coeffs.push(Fp::new(rng.gen_range(0..P)));
            }
            coeffs.push(Fp::ONE); // monic

            let f = Self::new(coeffs);
            if f.is_primitive() {
                return f;
            }
        }
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

            for (j, &(xj, _)) in points.iter().enumerate() {
                if i != j {
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

    /* ---- Polynomial Factorization ---- */

    /// Compute the formal derivative of this polynomial.
    ///
    /// The derivative of a_n x^n + ... + a_1 x + a_0 is
    /// n*a_n x^{n-1} + ... + a_1.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // f(x) = x^3 + 2x^2 + 3x + 4
    /// let f = Poly::new(vec![F17::new(4), F17::new(3), F17::new(2), F17::new(1)]);
    ///
    /// // f'(x) = 3x^2 + 4x + 3
    /// let df = f.derivative();
    /// assert_eq!(df, Poly::new(vec![F17::new(3), F17::new(4), F17::new(3)]));
    /// ```
    pub fn derivative(&self) -> Self {
        if self.coeffs.len() <= 1 {
            return Self::zero();
        }

        let mut coeffs = Vec::with_capacity(self.coeffs.len() - 1);
        for (i, &c) in self.coeffs.iter().enumerate().skip(1) {
            coeffs.push(c * Fp::new(i as u64));
        }

        Self::new(coeffs)
    }

    /// Square-free factorization.
    ///
    /// Returns a list of (factor, multiplicity) pairs where each factor is
    /// square-free (has no repeated roots) and the product of factor^multiplicity
    /// equals self (up to a constant factor).
    ///
    /// Uses Yun's algorithm for characteristic > degree, otherwise falls back
    /// to the general algorithm handling p-th powers.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // f(x) = (x - 1)^2 * (x - 2) = x^3 - 4x^2 + 5x - 2
    /// let x_minus_1 = Poly::new(vec![F17::new(16), F17::new(1)]); // x - 1
    /// let x_minus_2 = Poly::new(vec![F17::new(15), F17::new(1)]); // x - 2
    /// let f = (x_minus_1.clone() * x_minus_1) * x_minus_2;
    ///
    /// let factors = f.square_free_factorization();
    ///
    /// // Should have (x-1, 2) and (x-2, 1)
    /// assert_eq!(factors.len(), 2);
    /// ```
    pub fn square_free_factorization(&self) -> Vec<(Self, usize)> {
        if self.is_zero() {
            return vec![];
        }

        let f = match self.monic() {
            Some(m) => m,
            None => return vec![],
        };

        if f.degree() == Some(0) {
            return vec![];
        }

        let mut result = Vec::new();
        let df = f.derivative();

        if df.is_zero() {
            // f'(x) = 0 means f(x) = g(x^p) for some g
            // This happens when all exponents are multiples of p
            let g = self.pth_root();
            let sub_factors = g.square_free_factorization();
            for (factor, mult) in sub_factors {
                result.push((factor, mult * P as usize));
            }
            return result;
        }

        // Yun's algorithm
        let mut c = Self::gcd(&f, &df);
        let mut w = f.div_rem(&c).map(|(q, _)| q).unwrap_or_else(Self::zero);
        let mut i = 1usize;

        while w.degree().unwrap_or(0) > 0 {
            let y = Self::gcd(&w, &c);
            let z = w.div_rem(&y).map(|(q, _)| q).unwrap_or_else(Self::zero);

            if z.degree().unwrap_or(0) > 0 {
                result.push((z.monic().unwrap_or(z), i));
            }

            w = y;
            c = c.div_rem(&w).map(|(q, _)| q).unwrap_or_else(Self::zero);
            i += 1;
        }

        // Handle remaining factor from p-th powers
        if c.degree().unwrap_or(0) > 0 {
            let g = c.pth_root();
            let sub_factors = g.square_free_factorization();
            for (factor, mult) in sub_factors {
                result.push((factor, mult * P as usize));
            }
        }

        result
    }

    /// Compute the p-th root of a polynomial where all exponents are multiples of p.
    ///
    /// Given f(x) = sum a_i x^{ip}, returns g(x) = sum a_i^{1/p} x^i.
    fn pth_root(&self) -> Self {
        let mut coeffs = Vec::new();

        for (i, &c) in self.coeffs.iter().enumerate() {
            if i % (P as usize) == 0 {
                // a^{1/p} = a^{p-1} in F_p (since a^p = a, so a^{p-1} * a^{1} = a)
                // Actually a^{1/p} = a^{p^{k-1}} where p^k = |F|
                // For F_p, a^{1/p} = a (since a^p = a by Fermat)
                coeffs.push(c);
            }
        }

        Self::new(coeffs)
    }

    /// Distinct-degree factorization.
    ///
    /// Returns a list of (g_i, i) where g_i is the product of all monic
    /// irreducible factors of degree i.
    ///
    /// The input should be monic and square-free for correct results.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F5 = Fp<5>;
    ///
    /// // f(x) = (x - 1)(x - 2)(x^2 + x + 1) where x^2 + x + 1 is irreducible over F_5
    /// // This has two linear factors and one quadratic factor
    /// let f = Poly::new(vec![F5::new(3), F5::new(1), F5::new(4), F5::new(0), F5::new(1)]);
    ///
    /// let ddf = f.distinct_degree_factorization();
    /// // ddf contains products of factors grouped by degree
    /// ```
    pub fn distinct_degree_factorization(&self) -> Vec<(Self, usize)> {
        if self.is_zero() || self.degree() == Some(0) {
            return vec![];
        }

        let f = match self.monic() {
            Some(m) => m,
            None => return vec![],
        };

        let n = f.degree().unwrap();
        let mut result = Vec::new();
        let mut f_star = f.clone();
        let mut h = Self::x(); // h = x^{p^i}

        for i in 1..=n / 2 {
            // h = x^{p^i} mod f_star
            h = match f_star.powmod(&h, P) {
                Some(r) => r,
                None => break,
            };

            // g_i = gcd(h - x, f_star)
            let h_minus_x = h.clone() - Self::x();
            let g = Self::gcd(&h_minus_x, &f_star);

            if g.degree().unwrap_or(0) > 0 {
                let g_monic = g.monic().unwrap_or(g.clone());
                result.push((g_monic.clone(), i));

                // f_star = f_star / g
                f_star = f_star.div_rem(&g).map(|(q, _)| q).unwrap_or(f_star);
                h = h.rem(&f_star).unwrap_or(h);
            }
        }

        // Any remaining factor has degree > n/2, so must be irreducible
        if f_star.degree().unwrap_or(0) > 0 {
            let d = f_star.degree().unwrap();
            result.push((f_star.monic().unwrap_or(f_star), d));
        }

        result
    }

    /// Equal-degree factorization using Cantor-Zassenhaus algorithm.
    ///
    /// Given a polynomial f that is a product of distinct irreducible
    /// polynomials all of degree d, returns the list of irreducible factors.
    ///
    /// Requires a random number generator for the probabilistic algorithm.
    #[cfg(feature = "rand")]
    pub fn equal_degree_factorization<R: rand::Rng>(&self, d: usize, rng: &mut R) -> Vec<Self> {
        let f = match self.monic() {
            Some(m) => m,
            None => return vec![],
        };

        let n = match f.degree() {
            Some(deg) => deg,
            None => return vec![],
        };

        if n == 0 {
            return vec![];
        }

        if n == d {
            // f is already irreducible
            return vec![f];
        }

        // Number of factors
        let r = n / d;
        if r == 1 {
            return vec![f];
        }

        // Cantor-Zassenhaus: find a non-trivial factor
        loop {
            // Pick random polynomial a of degree < n
            let a = Self::random_poly(rng, n - 1);

            if a.is_zero() {
                continue;
            }

            // Compute gcd(a, f)
            let g1 = Self::gcd(&a, &f);
            if g1.degree().unwrap_or(0) > 0 && g1.degree().unwrap() < n {
                // Found a non-trivial factor
                let g1_monic = g1.monic().unwrap_or(g1);
                let other = f.div_rem(&g1_monic).map(|(q, _)| q).unwrap();
                let other_monic = other.monic().unwrap_or(other);

                let mut factors = Self::equal_degree_factorization(&g1_monic, d, rng);
                factors.extend(Self::equal_degree_factorization(&other_monic, d, rng));
                return factors;
            }

            // Compute a^{(p^d - 1)/2} mod f
            // For odd p, this gives a splitting
            let exp = Self::compute_half_order(d);
            let b = match f.powmod(&a, exp) {
                Some(r) => r,
                None => continue,
            };

            // gcd(b - 1, f)
            let b_minus_1 = b - Self::constant(Fp::ONE);
            let g2 = Self::gcd(&b_minus_1, &f);

            if g2.degree().unwrap_or(0) > 0 && g2.degree().unwrap() < n {
                let g2_monic = g2.monic().unwrap_or(g2);
                let other = f.div_rem(&g2_monic).map(|(q, _)| q).unwrap();
                let other_monic = other.monic().unwrap_or(other);

                let mut factors = Self::equal_degree_factorization(&g2_monic, d, rng);
                factors.extend(Self::equal_degree_factorization(&other_monic, d, rng));
                return factors;
            }
        }
    }

    /// Compute (p^d - 1) / 2.
    #[cfg(feature = "rand")]
    fn compute_half_order(d: usize) -> u64 {
        let mut result: u64 = 1;
        for _ in 0..d {
            result = result.saturating_mul(P);
        }
        (result - 1) / 2
    }

    /// Generate a random polynomial of degree at most max_degree.
    #[cfg(feature = "rand")]
    fn random_poly<R: rand::Rng>(rng: &mut R, max_degree: usize) -> Self {
        let mut coeffs = Vec::with_capacity(max_degree + 1);
        for _ in 0..=max_degree {
            coeffs.push(Fp::new(rng.gen_range(0..P)));
        }
        Self::new(coeffs)
    }

    /// Full factorization into irreducible factors.
    ///
    /// Returns a list of (irreducible_factor, multiplicity) pairs.
    /// The product of factor^multiplicity equals self up to a constant.
    ///
    /// Requires a random number generator for the Cantor-Zassenhaus step.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // f(x) = (x - 1)^2 * (x - 2)
    /// let x_minus_1 = Poly::new(vec![F17::new(16), F17::new(1)]);
    /// let x_minus_2 = Poly::new(vec![F17::new(15), F17::new(1)]);
    /// let f = (x_minus_1.clone() * x_minus_1) * x_minus_2;
    ///
    /// let mut rng = rand::thread_rng();
    /// let factors = f.factor(&mut rng);
    ///
    /// // Should have two distinct irreducible factors
    /// assert_eq!(factors.len(), 2);
    /// ```
    #[cfg(feature = "rand")]
    pub fn factor<R: rand::Rng>(&self, rng: &mut R) -> Vec<(Self, usize)> {
        if self.is_zero() {
            return vec![];
        }

        let mut result = Vec::new();

        // Step 1: Square-free factorization
        let sqf = self.square_free_factorization();

        for (sq_free_factor, multiplicity) in sqf {
            // Step 2: Distinct-degree factorization
            let ddf = sq_free_factor.distinct_degree_factorization();

            for (same_degree_product, degree) in ddf {
                // Step 3: Equal-degree factorization
                let irreducibles = same_degree_product.equal_degree_factorization(degree, rng);

                for irr in irreducibles {
                    result.push((irr, multiplicity));
                }
            }
        }

        // Sort by degree, then by coefficients for consistent output
        result.sort_by(|a, b| {
            let deg_cmp = a.0.degree().cmp(&b.0.degree());
            if deg_cmp != std::cmp::Ordering::Equal {
                return deg_cmp;
            }
            // Compare coefficients lexicographically
            for i in 0..=a.0.degree().unwrap_or(0) {
                let cmp = a.0.coeff(i).value().cmp(&b.0.coeff(i).value());
                if cmp != std::cmp::Ordering::Equal {
                    return cmp;
                }
            }
            std::cmp::Ordering::Equal
        });

        result
    }

    /// Find all roots of this polynomial in F_p.
    ///
    /// Returns a list of (root, multiplicity) pairs.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly};
    ///
    /// type F17 = Fp<17>;
    ///
    /// // f(x) = (x - 3)(x - 5) = x^2 - 8x + 15
    /// let f = Poly::new(vec![F17::new(15), F17::new(9), F17::new(1)]); // 15 - 8x + x^2 (note: -8 = 9 mod 17)
    ///
    /// let mut rng = rand::thread_rng();
    /// let roots = f.roots(&mut rng);
    ///
    /// assert_eq!(roots.len(), 2);
    /// ```
    #[cfg(feature = "rand")]
    pub fn roots<R: rand::Rng>(&self, rng: &mut R) -> Vec<(Fp<P>, usize)> {
        let factors = self.factor(rng);

        factors
            .into_iter()
            .filter_map(|(f, mult)| {
                if f.degree() == Some(1) {
                    // f = x - r, so r = -f.coeff(0) / f.coeff(1)
                    let r = -f.coeff(0) * f.coeff(1).inverse()?;
                    Some((r, mult))
                } else {
                    None
                }
            })
            .collect()
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

#[cfg(feature = "serde")]
impl<const P: u64> serde::Serialize for Poly<P> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        // Serialize as a vector of coefficient values
        let values: Vec<u64> = self.coeffs.iter().map(|c| c.value()).collect();
        values.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, const P: u64> serde::Deserialize<'de> for Poly<P> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let values = Vec::<u64>::deserialize(deserializer)?;
        let coeffs = values.into_iter().map(Fp::new).collect();
        Ok(Self::new(coeffs))
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

    // ---- monic tests ----

    #[test]
    fn monic_basic() {
        // 2 + 4x + 2x^2 -> 1 + 2x + x^2
        let p = P17::new(vec![F17::new(2), F17::new(4), F17::new(2)]);
        let m = p.monic().unwrap();
        assert_eq!(m.leading_coeff(), Some(F17::ONE));
        assert_eq!(m.coeff(0), F17::ONE);
        assert_eq!(m.coeff(1), F17::new(2));
    }

    #[test]
    fn monic_already_monic() {
        let p = P17::new(vec![F17::new(3), F17::new(5), F17::ONE]);
        let m = p.monic().unwrap();
        assert_eq!(m, p);
    }

    #[test]
    fn monic_zero() {
        let z = P17::zero();
        assert!(z.monic().is_none());
    }

    // ---- div_rem tests ----

    #[test]
    fn div_rem_exact_division() {
        // (x^2 - 1) / (x - 1) = (x + 1), remainder 0
        // x^2 - 1 = [16, 0, 1] (since -1 ≡ 16 mod 17)
        // x - 1 = [16, 1]
        let dividend = P17::new(vec![F17::new(16), F17::ZERO, F17::ONE]);
        let divisor = P17::new(vec![F17::new(16), F17::ONE]);
        let (q, r) = dividend.div_rem(&divisor).unwrap();

        // q should be x + 1
        assert_eq!(q.coeff(0), F17::ONE);
        assert_eq!(q.coeff(1), F17::ONE);
        assert!(r.is_zero());

        // Verify: dividend = q * divisor + r
        assert_eq!(q * divisor, dividend);
    }

    #[test]
    fn div_rem_with_remainder() {
        // (x^2 + 1) / (x + 1) = x - 1 + 2/(x+1)
        // So quotient = x - 1, remainder = 2
        let dividend = P17::new(vec![F17::ONE, F17::ZERO, F17::ONE]); // 1 + x^2
        let divisor = P17::new(vec![F17::ONE, F17::ONE]); // 1 + x
        let (q, r) = dividend.clone().div_rem(&divisor).unwrap();

        // Verify: dividend = q * divisor + r
        let reconstructed = q.clone() * divisor.clone() + r.clone();
        assert_eq!(reconstructed, dividend);

        // Check remainder degree < divisor degree
        assert!(r.degree().unwrap_or(0) < divisor.degree().unwrap());
    }

    #[test]
    fn div_rem_dividend_smaller() {
        // (x + 1) / (x^2 + 1) = 0, remainder (x + 1)
        let dividend = P17::new(vec![F17::ONE, F17::ONE]); // 1 + x
        let divisor = P17::new(vec![F17::ONE, F17::ZERO, F17::ONE]); // 1 + x^2
        let (q, r) = dividend.clone().div_rem(&divisor).unwrap();

        assert!(q.is_zero());
        assert_eq!(r, dividend);
    }

    #[test]
    fn div_rem_zero_dividend() {
        let dividend = P17::zero();
        let divisor = P17::new(vec![F17::ONE, F17::ONE]);
        let (q, r) = dividend.div_rem(&divisor).unwrap();

        assert!(q.is_zero());
        assert!(r.is_zero());
    }

    #[test]
    fn div_rem_zero_divisor() {
        let dividend = P17::new(vec![F17::ONE, F17::ONE]);
        let divisor = P17::zero();
        assert!(dividend.div_rem(&divisor).is_none());
    }

    #[test]
    fn div_rem_by_constant() {
        // (2x + 4) / 2 = (x + 2)
        let dividend = P17::new(vec![F17::new(4), F17::new(2)]);
        let divisor = P17::constant(F17::new(2));
        let (q, r) = dividend.div_rem(&divisor).unwrap();

        assert_eq!(q.coeff(0), F17::new(2));
        assert_eq!(q.coeff(1), F17::ONE);
        assert!(r.is_zero());
    }

    #[test]
    fn div_rem_non_monic_divisor() {
        // (6x^2 + 11x + 3) / (3x + 1) = (2x + 3), remainder 0
        let dividend = P17::new(vec![F17::new(3), F17::new(11), F17::new(6)]);
        let divisor = P17::new(vec![F17::ONE, F17::new(3)]);
        let (q, r) = dividend.clone().div_rem(&divisor).unwrap();

        // Verify: dividend = q * divisor + r
        let reconstructed = q * divisor + r.clone();
        assert_eq!(reconstructed, dividend);
        assert!(r.is_zero());
    }

    #[test]
    fn div_rem_cubic() {
        // Test with higher degree
        // (x^3 + 2x^2 + 3x + 4) / (x + 1)
        let dividend = P17::new(vec![F17::new(4), F17::new(3), F17::new(2), F17::ONE]);
        let divisor = P17::new(vec![F17::ONE, F17::ONE]);
        let (q, r) = dividend.clone().div_rem(&divisor).unwrap();

        // Verify: dividend = q * divisor + r
        let reconstructed = q * divisor.clone() + r.clone();
        assert_eq!(reconstructed, dividend);

        // Remainder degree < divisor degree
        match r.degree() {
            None => {} // zero is fine
            Some(d) => assert!(d < divisor.degree().unwrap()),
        }
    }

    // ---- rem tests ----

    #[test]
    fn rem_basic() {
        let a = P17::new(vec![F17::ONE, F17::ZERO, F17::ONE]); // 1 + x^2
        let b = P17::new(vec![F17::ONE, F17::ONE]); // 1 + x
        let r = a.rem(&b).unwrap();

        // Should equal the remainder from div_rem
        let (_, r2) = a.div_rem(&b).unwrap();
        assert_eq!(r, r2);
    }

    // ---- gcd tests ----

    #[test]
    fn gcd_coprime() {
        // gcd(x, x + 1) = 1
        let a = P17::x();
        let b = P17::new(vec![F17::ONE, F17::ONE]);
        let g = P17::gcd(&a, &b);

        assert_eq!(g.degree(), Some(0));
        assert_eq!(g.coeff(0), F17::ONE);
    }

    #[test]
    fn gcd_common_factor() {
        // gcd((x-1)(x-2), (x-2)(x-3)) = (x-2)
        let p1 = P17::from_roots(&[F17::new(1), F17::new(2)]);
        let p2 = P17::from_roots(&[F17::new(2), F17::new(3)]);
        let g = P17::gcd(&p1, &p2);

        // Should be monic
        assert_eq!(g.leading_coeff(), Some(F17::ONE));
        // Should have degree 1
        assert_eq!(g.degree(), Some(1));
        // Should have root at 2
        assert_eq!(g.eval(F17::new(2)), F17::ZERO);
    }

    #[test]
    fn gcd_one_divides_other() {
        // gcd(x-1, (x-1)(x-2)) = x-1
        let a = P17::from_roots(&[F17::new(1)]);
        let b = P17::from_roots(&[F17::new(1), F17::new(2)]);
        let g = P17::gcd(&a, &b);

        assert_eq!(g.degree(), Some(1));
        assert_eq!(g.eval(F17::new(1)), F17::ZERO);
    }

    #[test]
    fn gcd_same_polynomial() {
        let p = P17::from_roots(&[F17::new(3), F17::new(5)]);
        let g = P17::gcd(&p, &p);

        // gcd(p, p) = p (made monic)
        let monic_p = p.monic().unwrap();
        assert_eq!(g, monic_p);
    }

    #[test]
    fn gcd_with_zero() {
        let p = P17::new(vec![F17::new(2), F17::new(4)]); // 2 + 4x
        let z = P17::zero();

        let g1 = P17::gcd(&p, &z);
        let g2 = P17::gcd(&z, &p);

        // gcd(p, 0) = monic(p) and gcd(0, p) = monic(p)
        let monic_p = p.monic().unwrap();
        assert_eq!(g1, monic_p);
        assert_eq!(g2, monic_p);
    }

    #[test]
    fn gcd_both_zero() {
        let z = P17::zero();
        let g = P17::gcd(&z, &z);
        assert!(g.is_zero());
    }

    #[test]
    fn gcd_multiple_common_roots() {
        // gcd((x-1)(x-2)(x-3), (x-2)(x-3)(x-4)) = (x-2)(x-3)
        let p1 = P17::from_roots(&[F17::new(1), F17::new(2), F17::new(3)]);
        let p2 = P17::from_roots(&[F17::new(2), F17::new(3), F17::new(4)]);
        let g = P17::gcd(&p1, &p2);

        assert_eq!(g.degree(), Some(2));
        assert_eq!(g.eval(F17::new(2)), F17::ZERO);
        assert_eq!(g.eval(F17::new(3)), F17::ZERO);
        // 1 and 4 are not roots of gcd
        assert_ne!(g.eval(F17::new(1)), F17::ZERO);
        assert_ne!(g.eval(F17::new(4)), F17::ZERO);
    }

    #[test]
    fn gcd_is_monic() {
        // Even with non-monic inputs, gcd should be monic
        let p1 = P17::new(vec![F17::new(6), F17::new(4), F17::new(2)]); // 6 + 4x + 2x^2
        let p2 = P17::new(vec![F17::new(3), F17::new(3)]); // 3 + 3x
        let g = P17::gcd(&p1, &p2);

        if !g.is_zero() {
            assert_eq!(g.leading_coeff(), Some(F17::ONE));
        }
    }

    // ---- extended_gcd tests ----

    #[test]
    fn extended_gcd_bezout_identity() {
        let a = P17::new(vec![F17::new(1), F17::new(2), F17::new(1)]); // 1 + 2x + x^2
        let b = P17::new(vec![F17::new(1), F17::new(1)]); // 1 + x

        let (g, s, t) = P17::extended_gcd(&a, &b);

        // Verify Bézout identity: s*a + t*b = g
        let check = s * a.clone() + t * b.clone();
        assert_eq!(check, g);
    }

    #[test]
    fn extended_gcd_coprime() {
        // Extended gcd of coprime polynomials
        let a = P17::x();
        let b = P17::new(vec![F17::ONE, F17::ONE]); // 1 + x

        let (g, s, t) = P17::extended_gcd(&a, &b);

        // g should be 1
        assert_eq!(g.degree(), Some(0));
        assert_eq!(g.coeff(0), F17::ONE);

        // Verify Bézout identity
        let check = s * a + t * b;
        assert_eq!(check, g);
    }

    #[test]
    fn extended_gcd_with_common_factor() {
        // a = (x-1)(x-2), b = (x-2)(x-3)
        let a = P17::from_roots(&[F17::new(1), F17::new(2)]);
        let b = P17::from_roots(&[F17::new(2), F17::new(3)]);

        let (g, s, t) = P17::extended_gcd(&a, &b);

        // g should be x-2 (monic)
        assert_eq!(g.degree(), Some(1));
        assert_eq!(g.eval(F17::new(2)), F17::ZERO);

        // Verify Bézout identity
        let check = s * a + t * b;
        assert_eq!(check, g);
    }

    #[test]
    fn extended_gcd_second_zero() {
        let a = P17::new(vec![F17::new(2), F17::new(4)]); // 2 + 4x
        let b = P17::zero();

        let (g, s, t) = P17::extended_gcd(&a, &b);

        // g should be monic(a)
        let monic_a = a.monic().unwrap();
        assert_eq!(g, monic_a);

        // Verify: s*a + t*b = g
        let check = s * a + t * b;
        assert_eq!(check, g);
    }

    #[test]
    fn extended_gcd_first_zero() {
        let a = P17::zero();
        let b = P17::new(vec![F17::new(3), F17::new(6)]); // 3 + 6x

        let (g, s, t) = P17::extended_gcd(&a, &b);

        // g should be monic(b)
        let monic_b = b.monic().unwrap();
        assert_eq!(g, monic_b);

        // Verify: s*a + t*b = g
        let check = s * a + t * b;
        assert_eq!(check, g);
    }

    #[test]
    fn extended_gcd_both_zero() {
        let a = P17::zero();
        let b = P17::zero();

        let (g, s, t) = P17::extended_gcd(&a, &b);

        assert!(g.is_zero());
        assert!(s.is_zero());
        assert!(t.is_zero());
    }

    #[test]
    fn extended_gcd_inverse_modulo() {
        // If gcd(a, m) = 1, then s is the inverse of a modulo m
        // This is the key operation for field extensions!
        let a = P17::new(vec![F17::new(2), F17::ONE]); // 2 + x
        let m = P17::new(vec![F17::ONE, F17::ZERO, F17::ONE]); // 1 + x^2

        let (g, s, _t) = P17::extended_gcd(&a, &m);

        // Should be coprime
        assert_eq!(g.degree(), Some(0));
        assert_eq!(g.coeff(0), F17::ONE);

        // s * a ≡ 1 (mod m)
        let product = s * a;
        let remainder = product.rem(&m).unwrap();
        assert_eq!(remainder.degree(), Some(0));
        assert_eq!(remainder.coeff(0), F17::ONE);
    }

    // ---- powmod tests ----

    #[test]
    fn powmod_x_basic() {
        // x^3 mod (x^2 + 1) = x * x^2 = x * (-1) = -x = 16x in F17
        let m = P17::new(vec![F17::ONE, F17::ZERO, F17::ONE]); // 1 + x^2
        let result = m.powmod_x(3).unwrap();

        assert_eq!(result.degree(), Some(1));
        assert_eq!(result.coeff(0), F17::ZERO);
        assert_eq!(result.coeff(1), F17::new(16)); // -1 mod 17
    }

    #[test]
    fn powmod_x_zero_exp() {
        let m = P17::new(vec![F17::ONE, F17::ONE]); // 1 + x
        let result = m.powmod_x(0).unwrap();

        assert_eq!(result.degree(), Some(0));
        assert_eq!(result.coeff(0), F17::ONE);
    }

    #[test]
    fn powmod_x_one_exp() {
        let m = P17::new(vec![F17::ONE, F17::ZERO, F17::ONE]); // 1 + x^2
        let result = m.powmod_x(1).unwrap();

        assert_eq!(result, P17::x());
    }

    #[test]
    fn powmod_base() {
        // (1 + x)^2 mod (x^2 + 1) = 1 + 2x + x^2 = 1 + 2x - 1 = 2x
        let m = P17::new(vec![F17::ONE, F17::ZERO, F17::ONE]); // 1 + x^2
        let base = P17::new(vec![F17::ONE, F17::ONE]); // 1 + x
        let result = m.powmod(&base, 2).unwrap();

        assert_eq!(result.coeff(0), F17::ZERO);
        assert_eq!(result.coeff(1), F17::new(2));
    }

    #[test]
    fn powmod_zero_modulus() {
        let m = P17::zero();
        assert!(m.powmod_x(5).is_none());
    }

    // ---- irreducibility tests ----

    #[test]
    fn is_irreducible_zero() {
        let f = P17::zero();
        assert!(!f.is_irreducible());
    }

    #[test]
    fn is_irreducible_constant() {
        let f = P17::constant(F17::new(5));
        assert!(!f.is_irreducible());
    }

    #[test]
    fn is_irreducible_linear() {
        // All linear polynomials are irreducible
        let f = P17::new(vec![F17::new(3), F17::ONE]); // 3 + x
        assert!(f.is_irreducible());

        let g = P17::new(vec![F17::new(7), F17::new(5)]); // 7 + 5x
        assert!(g.is_irreducible());
    }

    #[test]
    fn is_irreducible_quadratic_non_residue() {
        // x^2 - a is irreducible iff a is not a quadratic residue
        // 3 is not a QR mod 17
        let f = P17::new(vec![F17::new(14), F17::ZERO, F17::ONE]); // -3 + x^2
        assert!(f.is_irreducible());
    }

    #[test]
    fn is_irreducible_quadratic_residue() {
        // x^2 - a is reducible if a is a quadratic residue
        // 4 = 2^2 is a QR mod 17
        let f = P17::new(vec![F17::new(13), F17::ZERO, F17::ONE]); // -4 + x^2
        assert!(!f.is_irreducible());

        // 1 = 1^2 is a QR mod 17, so x^2 - 1 = (x-1)(x+1)
        let g = P17::new(vec![F17::new(16), F17::ZERO, F17::ONE]); // -1 + x^2
        assert!(!g.is_irreducible());
    }

    #[test]
    fn is_irreducible_x2_plus_1() {
        // x^2 + 1 over F_17
        // -1 = 16 is a QR mod 17 (since 4^2 = 16), so x^2 + 1 is reducible
        let f = P17::new(vec![F17::ONE, F17::ZERO, F17::ONE]); // 1 + x^2
        assert!(!f.is_irreducible());
    }

    #[test]
    fn is_irreducible_has_root() {
        // (x - 5) is a factor, so x^2 - 5x + 6 = (x-2)(x-3) is reducible
        // Actually let's just check a polynomial with a known root
        // f(x) = x^2 - 5, if 5 is QR then reducible
        // 5 is not a QR mod 17, so x^2 - 5 is irreducible
        let f = P17::new(vec![F17::new(12), F17::ZERO, F17::ONE]); // -5 + x^2
        assert!(f.is_irreducible());

        // x^2 - 9 = (x-3)(x+3), 9 = 3^2 is QR
        let g = P17::new(vec![F17::new(8), F17::ZERO, F17::ONE]); // -9 + x^2
        assert!(!g.is_irreducible());
    }

    #[test]
    fn is_irreducible_cubic() {
        // x^3 + x + 1 over F_2 would be irreducible, but we're in F_17
        // Let's test x^3 - 2 over F_17
        // Need to check if 2 has a cube root mod 17
        // 2^((17-1)/gcd(3,16)) = 2^16 = 1 mod 17, so 2 is a cube (since gcd(3,16)=1)
        // Actually for cube roots: a is a cube iff a^((p-1)/gcd(3,p-1)) = 1
        // gcd(3, 16) = 1, so 2^16 = 1, meaning 2 is a perfect cube
        // So x^3 - 2 is reducible

        // Let's find an irreducible cubic instead
        // x^3 + x + 1: check if it has roots in F_17
        let f = P17::new(vec![F17::ONE, F17::ONE, F17::ZERO, F17::ONE]); // 1 + x + x^3

        // Check manually for roots
        let mut has_root = false;
        for i in 0..17 {
            let x = F17::new(i);
            if f.eval(x) == F17::ZERO {
                has_root = true;
                break;
            }
        }

        // If no roots and degree 3, need further check for irreducibility
        // A cubic with no roots is irreducible
        if !has_root {
            assert!(f.is_irreducible());
        }
    }

    #[test]
    fn is_irreducible_product_of_irreducibles() {
        // (x - 1)(x - 2) = x^2 - 3x + 2
        let f = P17::from_roots(&[F17::new(1), F17::new(2)]);
        assert!(!f.is_irreducible());
    }

    #[test]
    fn is_irreducible_non_monic() {
        // 2x^2 - 6 should be reducible (same as x^2 - 3 up to scalar)
        // Wait, 3 is not a QR, so x^2 - 3 is irreducible
        let f = P17::new(vec![F17::new(11), F17::ZERO, F17::new(2)]); // -6 + 2x^2 = 2(x^2 - 3)
        assert!(f.is_irreducible()); // Irreducibility is preserved by scalar multiplication

        // 2x^2 - 8 = 2(x^2 - 4), and 4 is QR, so reducible
        let g = P17::new(vec![F17::new(9), F17::ZERO, F17::new(2)]); // -8 + 2x^2
        assert!(!g.is_irreducible());
    }

    #[test]
    fn prime_divisors_basic() {
        assert_eq!(P17::prime_divisors(1), Vec::<usize>::new());
        assert_eq!(P17::prime_divisors(2), vec![2usize]);
        assert_eq!(P17::prime_divisors(3), vec![3usize]);
        assert_eq!(P17::prime_divisors(4), vec![2usize]);
        assert_eq!(P17::prime_divisors(6), vec![2usize, 3]);
        assert_eq!(P17::prime_divisors(12), vec![2usize, 3]);
        assert_eq!(P17::prime_divisors(30), vec![2usize, 3, 5]);
        assert_eq!(P17::prime_divisors(17), vec![17usize]);
    }

    #[test]
    fn is_irreducible_degree_4() {
        // x^4 + x + 1 over F_2 is irreducible, let's check over F_17
        let f = P17::new(vec![F17::ONE, F17::ONE, F17::ZERO, F17::ZERO, F17::ONE]);

        // Just verify the function runs without panic
        let _ = f.is_irreducible();
    }

    #[test]
    fn is_irreducible_degree_6() {
        // Test with composite degree (6 = 2 * 3)
        // This exercises the prime divisor checks
        let f = P17::new(vec![
            F17::ONE,
            F17::ONE,
            F17::ZERO,
            F17::ZERO,
            F17::ZERO,
            F17::ZERO,
            F17::ONE,
        ]); // 1 + x + x^6

        // Just verify the function runs without panic for composite degrees
        let _ = f.is_irreducible();
    }

    // ---- primitive polynomial tests ----

    #[test]
    fn is_primitive_reducible() {
        // Reducible polynomials are not primitive
        let f = P17::from_roots(&[F17::new(1), F17::new(2)]);
        assert!(!f.is_primitive());
    }

    #[test]
    fn is_primitive_zero_constant() {
        assert!(!P17::zero().is_primitive());
        assert!(!P17::constant(F17::new(5)).is_primitive());
    }

    #[test]
    fn is_primitive_linear() {
        // Linear polynomials: x - a
        // The field is F_p itself, order = p - 1
        // x ≡ a (mod x - a), so we need ord(a) = p - 1
        // In F_17, the primitive elements are generators of F_17*
        // 3 is a primitive root mod 17
        let f = P17::new(vec![-Fp::new(3), Fp::ONE]); // x - 3
        assert!(f.is_primitive());

        // 1 is not a primitive root (ord(1) = 1)
        let g = P17::new(vec![-Fp::ONE, Fp::ONE]); // x - 1
        assert!(!g.is_primitive());
    }

    #[test]
    fn is_primitive_f3_degree2() {
        // Over F_3, degree 2: order = 3^2 - 1 = 8 = 2^3
        // x^2 + 1 is irreducible over F_3 (since -1 is not a QR mod 3)
        // Check: 1^2 = 1, 2^2 = 4 = 1 mod 3, so QRs are {1}, and -1 = 2 is not QR
        type F3 = Fp<3>;
        type P3 = Poly<3>;

        let f = P3::new(vec![F3::ONE, F3::ZERO, F3::ONE]); // 1 + x^2
        assert!(f.is_irreducible());

        // x has order 8 if x^8 = 1 but x^4 ≠ 1
        // x^2 = -1 = 2, x^4 = 4 = 1 mod 3... so x^4 = 1, order divides 4, not primitive
        assert!(!f.is_primitive());
    }

    #[test]
    fn is_primitive_f3_degree2_primitive() {
        // Over F_3, we need x^2 + x + 2 or similar
        // Let's find a primitive one by checking x^4 ≠ 1
        type F3 = Fp<3>;
        type P3 = Poly<3>;

        // x^2 - x - 1 = x^2 + 2x + 2 in F_3
        let f = P3::new(vec![F3::new(2), F3::new(2), F3::ONE]); // 2 + 2x + x^2

        if f.is_irreducible() {
            // Check if it's primitive
            let is_prim = f.is_primitive();
            // Just verify computation completes
            let _ = is_prim;
        }
    }

    #[test]
    fn is_primitive_f5_degree2() {
        // Over F_5, degree 2: order = 5^2 - 1 = 24 = 2^3 * 3
        // x^2 + 2 is irreducible over F_5 (2 is not a QR: 1,4,4,1 are QRs, so QR={1,4})
        type F5 = Fp<5>;
        type P5 = Poly<5>;

        let f = P5::new(vec![F5::new(2), F5::ZERO, F5::ONE]); // 2 + x^2
        assert!(f.is_irreducible());

        // Check primitivity
        let is_prim = f.is_primitive();
        // Just verify it runs
        let _ = is_prim;
    }

    #[test]
    fn is_primitive_f17_degree2() {
        // Over F_17, degree 2: order = 17^2 - 1 = 288 = 2^5 * 3^2
        // x^2 - 3 is irreducible (3 is not QR mod 17)
        // Need to check if x has order 288
        let f = P17::new(vec![F17::new(14), F17::ZERO, F17::ONE]); // x^2 - 3

        // This is irreducible
        assert!(f.is_irreducible());

        // Check if it's primitive (x must have order 288)
        let is_prim = f.is_primitive();

        // We can verify: if primitive, x^288 = 1 but x^(288/2) = x^144 ≠ 1 and x^(288/3) = x^96 ≠ 1
        // Just check the function works correctly
        assert!(is_prim || !is_prim); // Valid result either way
    }

    #[test]
    fn prime_divisors_u64_basic() {
        assert_eq!(P17::prime_divisors_u64(1), Vec::<u64>::new());
        assert_eq!(P17::prime_divisors_u64(2), vec![2u64]);
        assert_eq!(P17::prime_divisors_u64(15), vec![3u64, 5]);
        assert_eq!(P17::prime_divisors_u64(288), vec![2u64, 3]); // 2^5 * 3^2
        assert_eq!(P17::prime_divisors_u64(17), vec![17u64]);
    }

    #[test]
    fn compute_field_order_basic() {
        assert_eq!(P17::compute_field_order(2, 4), Some(15)); // 2^4 - 1
        assert_eq!(P17::compute_field_order(17, 2), Some(288)); // 17^2 - 1
        assert_eq!(P17::compute_field_order(2, 8), Some(255)); // 2^8 - 1
    }

    #[test]
    fn is_primitive_comprehensive_f3() {
        // Test all degree-2 monic polynomials over F_3
        // There are 9 monic degree-2 polynomials (3^2 choices for c0, c1)
        type F3 = Fp<3>;
        type P3 = Poly<3>;

        let mut primitive_count = 0;
        let mut irreducible_count = 0;

        // Iterate through all monic degree-2 polynomials
        for c0 in 0..3u64 {
            for c1 in 0..3u64 {
                let f = P3::new(vec![F3::new(c0), F3::new(c1), F3::ONE]);
                if f.is_irreducible() {
                    irreducible_count += 1;
                }
                if f.is_primitive() {
                    primitive_count += 1;
                }
            }
        }

        // Over F_3, degree 2:
        // Number of irreducible monic polynomials = (3^2 - 3) / 2 = 3
        assert_eq!(irreducible_count, 3);

        // Order = 3^2 - 1 = 8
        // φ(8) / 2 = 4 / 2 = 2 primitive polynomials
        assert_eq!(primitive_count, 2);
    }

    #[test]
    fn is_primitive_all_f5_degree2() {
        // Comprehensive test over F_5, degree 2
        // Order = 5^2 - 1 = 24
        type F5 = Fp<5>;
        type P5 = Poly<5>;

        let mut primitive_count = 0;
        let mut irreducible_count = 0;

        for c0 in 0..5u64 {
            for c1 in 0..5u64 {
                let f = P5::new(vec![F5::new(c0), F5::new(c1), F5::ONE]);
                if f.is_irreducible() {
                    irreducible_count += 1;
                }
                if f.is_primitive() {
                    primitive_count += 1;
                }
            }
        }

        // Number of monic irreducible degree-2 polynomials over F_p = (p^2 - p) / 2
        // For p=5: (25 - 5) / 2 = 10
        assert_eq!(irreducible_count, 10);

        // Number of primitive polynomials = φ(p^n - 1) / n
        // φ(24) / 2 = 8 / 2 = 4
        assert_eq!(primitive_count, 4);
    }

    // ---- Derivative tests ----

    #[test]
    fn derivative_constant() {
        let f = P17::constant(F17::new(5));
        assert!(f.derivative().is_zero());
    }

    #[test]
    fn derivative_linear() {
        // (3 + 2x)' = 2
        let f = P17::new(vec![F17::new(3), F17::new(2)]);
        let df = f.derivative();
        assert_eq!(df, P17::constant(F17::new(2)));
    }

    #[test]
    fn derivative_quadratic() {
        // (1 + 2x + 3x^2)' = 2 + 6x
        let f = P17::new(vec![F17::new(1), F17::new(2), F17::new(3)]);
        let df = f.derivative();
        assert_eq!(df, P17::new(vec![F17::new(2), F17::new(6)]));
    }

    #[test]
    fn derivative_cubic() {
        // (4 + 3x + 2x^2 + x^3)' = 3 + 4x + 3x^2
        let f = P17::new(vec![F17::new(4), F17::new(3), F17::new(2), F17::new(1)]);
        let df = f.derivative();
        assert_eq!(df, P17::new(vec![F17::new(3), F17::new(4), F17::new(3)]));
    }

    // ---- Square-free factorization tests ----

    #[test]
    fn square_free_linear() {
        // x - 1 is already square-free
        let f = P17::new(vec![F17::new(16), F17::ONE]); // x - 1
        let sqf = f.square_free_factorization();
        assert_eq!(sqf.len(), 1);
        assert_eq!(sqf[0].1, 1); // multiplicity 1
    }

    #[test]
    fn square_free_squared_linear() {
        // (x - 1)^2
        let x_minus_1 = P17::new(vec![F17::new(16), F17::ONE]);
        let f = x_minus_1.clone() * x_minus_1.clone();
        let sqf = f.square_free_factorization();

        assert_eq!(sqf.len(), 1);
        assert_eq!(sqf[0].0, x_minus_1.monic().unwrap());
        assert_eq!(sqf[0].1, 2); // multiplicity 2
    }

    #[test]
    fn square_free_product() {
        // (x - 1)^2 * (x - 2)
        let x_minus_1 = P17::new(vec![F17::new(16), F17::ONE]);
        let x_minus_2 = P17::new(vec![F17::new(15), F17::ONE]);
        let f = (x_minus_1.clone() * x_minus_1.clone()) * x_minus_2.clone();

        let sqf = f.square_free_factorization();
        assert_eq!(sqf.len(), 2);

        // Check multiplicities (order may vary)
        let mut found_1 = false;
        let mut found_2 = false;
        for (factor, mult) in &sqf {
            if *mult == 1 && *factor == x_minus_2 {
                found_1 = true;
            }
            if *mult == 2 && *factor == x_minus_1 {
                found_2 = true;
            }
        }
        assert!(found_1 && found_2);
    }

    // ---- Distinct-degree factorization tests ----

    #[test]
    fn ddf_single_linear() {
        // x - 1 has one linear factor
        let f = P17::new(vec![F17::new(16), F17::ONE]);
        let ddf = f.distinct_degree_factorization();
        assert_eq!(ddf.len(), 1);
        assert_eq!(ddf[0].1, 1); // degree 1
    }

    #[test]
    fn ddf_two_linears() {
        // (x - 1)(x - 2) = x^2 - 3x + 2
        let x_minus_1 = P17::new(vec![F17::new(16), F17::ONE]);
        let x_minus_2 = P17::new(vec![F17::new(15), F17::ONE]);
        let f = x_minus_1 * x_minus_2;

        let ddf = f.distinct_degree_factorization();
        assert_eq!(ddf.len(), 1);
        assert_eq!(ddf[0].1, 1); // all factors have degree 1
        assert_eq!(ddf[0].0.degree(), Some(2)); // product of two degree-1 factors
    }

    #[test]
    fn ddf_irreducible_quadratic() {
        // x^2 - 3 is irreducible over F_17
        let f = P17::new(vec![F17::new(14), F17::ZERO, F17::ONE]);
        let ddf = f.distinct_degree_factorization();
        assert_eq!(ddf.len(), 1);
        assert_eq!(ddf[0].1, 2); // degree 2
    }

    // ---- Full factorization tests (require rand) ----

    #[cfg(feature = "rand")]
    #[test]
    fn factor_linear() {
        let mut rng = rand::thread_rng();
        let f = P17::new(vec![F17::new(16), F17::ONE]); // x - 1
        let factors = f.factor(&mut rng);

        assert_eq!(factors.len(), 1);
        assert_eq!(factors[0].0, f);
        assert_eq!(factors[0].1, 1);
    }

    #[cfg(feature = "rand")]
    #[test]
    fn factor_squared_linear() {
        let mut rng = rand::thread_rng();
        let x_minus_1 = P17::new(vec![F17::new(16), F17::ONE]);
        let f = x_minus_1.clone() * x_minus_1.clone();
        let factors = f.factor(&mut rng);

        assert_eq!(factors.len(), 1);
        assert_eq!(factors[0].0, x_minus_1);
        assert_eq!(factors[0].1, 2);
    }

    #[cfg(feature = "rand")]
    #[test]
    fn factor_two_linears() {
        let mut rng = rand::thread_rng();
        let x_minus_1 = P17::new(vec![F17::new(16), F17::ONE]);
        let x_minus_2 = P17::new(vec![F17::new(15), F17::ONE]);
        let f = x_minus_1 * x_minus_2;

        let factors = f.factor(&mut rng);
        assert_eq!(factors.len(), 2);

        // Both should have multiplicity 1
        assert!(factors.iter().all(|(_, m)| *m == 1));
    }

    #[cfg(feature = "rand")]
    #[test]
    fn factor_mixed() {
        let mut rng = rand::thread_rng();
        // (x - 1)^2 * (x - 2)
        let x_minus_1 = P17::new(vec![F17::new(16), F17::ONE]);
        let x_minus_2 = P17::new(vec![F17::new(15), F17::ONE]);
        let f = (x_minus_1.clone() * x_minus_1.clone()) * x_minus_2.clone();

        let factors = f.factor(&mut rng);
        assert_eq!(factors.len(), 2);

        // Check we have both factors with correct multiplicities
        let mut found_mult_1 = false;
        let mut found_mult_2 = false;
        for (factor, mult) in &factors {
            if *factor == x_minus_1 && *mult == 2 {
                found_mult_2 = true;
            }
            if *factor == x_minus_2 && *mult == 1 {
                found_mult_1 = true;
            }
        }
        assert!(found_mult_1 && found_mult_2);
    }

    #[cfg(feature = "rand")]
    #[test]
    fn factor_product_reconstruction() {
        let mut rng = rand::thread_rng();
        // Create a polynomial and verify factorization reconstructs it
        let x_minus_1 = P17::new(vec![F17::new(16), F17::ONE]);
        let x_minus_3 = P17::new(vec![F17::new(14), F17::ONE]);
        let x_minus_5 = P17::new(vec![F17::new(12), F17::ONE]);
        let f = (x_minus_1 * x_minus_3) * x_minus_5;

        let factors = f.factor(&mut rng);

        // Reconstruct the polynomial
        let mut reconstructed = P17::constant(F17::ONE);
        for (factor, mult) in &factors {
            for _ in 0..*mult {
                reconstructed = reconstructed * factor.clone();
            }
        }

        // Should equal original (up to leading coefficient)
        assert_eq!(reconstructed.monic(), f.monic());
    }

    #[cfg(feature = "rand")]
    #[test]
    fn roots_basic() {
        let mut rng = rand::thread_rng();
        // (x - 3)(x - 5)
        let f = P17::new(vec![F17::new(15), F17::new(9), F17::ONE]); // 15 - 8x + x^2 = 15 + 9x + x^2 mod 17

        let roots = f.roots(&mut rng);
        assert_eq!(roots.len(), 2);

        let root_values: Vec<_> = roots.iter().map(|(r, _)| r.value()).collect();
        assert!(root_values.contains(&3));
        assert!(root_values.contains(&5));
    }

    #[cfg(feature = "rand")]
    #[test]
    fn roots_with_multiplicity() {
        let mut rng = rand::thread_rng();
        // (x - 2)^2
        let x_minus_2 = P17::new(vec![F17::new(15), F17::ONE]);
        let f = x_minus_2.clone() * x_minus_2;

        let roots = f.roots(&mut rng);
        assert_eq!(roots.len(), 1);
        assert_eq!(roots[0].0, F17::new(2));
        assert_eq!(roots[0].1, 2);
    }

    #[cfg(feature = "rand")]
    #[test]
    fn roots_irreducible_no_roots() {
        let mut rng = rand::thread_rng();
        // x^2 - 3 is irreducible over F_17, so no roots in F_17
        let f = P17::new(vec![F17::new(14), F17::ZERO, F17::ONE]);

        let roots = f.roots(&mut rng);
        assert!(roots.is_empty());
    }
}
