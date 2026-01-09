#[cfg(feature = "alloc")]
use alloc::vec;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;
use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};

use crate::algebra::field::Field;
use crate::algebra::ring::Ring;
use crate::structures::fp::Fp;
#[cfg(feature = "alloc")]
use crate::structures::poly::Poly;

/// Extension field F_p[x]/(f(x)) where f(x) is an irreducible polynomial of degree D.
///
/// Elements are represented as polynomials of degree < D, stored as an array of D coefficients.
/// The modulus polynomial must be provided separately when performing field operations.
///
/// # Type Parameters
/// - `P`: The prime modulus for the base field F_p
/// - `D`: The degree of the extension (degree of the irreducible polynomial)
///
/// # Example
///
/// ```
/// use kreep::{Fp, Poly, ExtField};
///
/// type F17 = Fp<17>;
///
/// // F_17[x]/(x^2 + 1) is a degree-2 extension
/// // x^2 + 1 is irreducible over F_17 since -1 is not a quadratic residue mod 17
/// let modulus = Poly::new(vec![F17::new(1), F17::new(0), F17::new(1)]); // 1 + x^2
///
/// // Create element (2 + 3x)
/// let a = ExtField::<17, 2>::new([F17::new(2), F17::new(3)]);
///
/// // Create element (1 + x)
/// let b = ExtField::<17, 2>::new([F17::new(1), F17::new(1)]);
///
/// // Multiply: (2 + 3x)(1 + x) = 2 + 5x + 3x^2
/// // Since x^2 ≡ -1 (mod x^2 + 1), we get: 2 + 5x - 3 = -1 + 5x = 16 + 5x
/// let c = a.mul_mod(&b, &modulus);
/// assert_eq!(c.coeffs()[0], F17::new(16));
/// assert_eq!(c.coeffs()[1], F17::new(5));
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct ExtField<const P: u64, const D: usize> {
    /// Coefficients in ascending order: coeffs[i] is the coefficient of x^i
    coeffs: [Fp<P>; D],
}

impl<const P: u64, const D: usize> ExtField<P, D> {
    /// Create a new extension field element from coefficients.
    ///
    /// Coefficients are in ascending order: `coeffs[i]` is the coefficient of `x^i`.
    pub fn new(coeffs: [Fp<P>; D]) -> Self {
        Self { coeffs }
    }

    /// Create the zero element.
    pub fn zero() -> Self {
        Self {
            coeffs: [Fp::ZERO; D],
        }
    }

    /// Create the multiplicative identity (1).
    pub fn one() -> Self {
        let mut coeffs = [Fp::ZERO; D];
        if D > 0 {
            coeffs[0] = Fp::ONE;
        }
        Self { coeffs }
    }

    /// Create an element representing the variable x.
    ///
    /// Returns `None` if D < 2 (since x would not be a valid element).
    pub fn x() -> Option<Self> {
        if D < 2 {
            return None;
        }
        let mut coeffs = [Fp::ZERO; D];
        coeffs[1] = Fp::ONE;
        Some(Self { coeffs })
    }

    /// Create an element from a base field element.
    pub fn from_base(c: Fp<P>) -> Self {
        let mut coeffs = [Fp::ZERO; D];
        if D > 0 {
            coeffs[0] = c;
        }
        Self { coeffs }
    }

    /// Get the coefficients.
    pub fn coeffs(&self) -> &[Fp<P>; D] {
        &self.coeffs
    }

    /// Get a specific coefficient.
    pub fn coeff(&self, i: usize) -> Fp<P> {
        if i < D {
            self.coeffs[i]
        } else {
            Fp::ZERO
        }
    }

    /// Check if this is the zero element.
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&c| c == Fp::ZERO)
    }

    /// Convert to a polynomial representation.
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn to_poly(&self) -> Poly<P> {
        Poly::new(self.coeffs.to_vec())
    }

    /// Create from a polynomial, reducing if necessary.
    ///
    /// The polynomial is reduced modulo the given modulus.
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn from_poly(p: &Poly<P>, modulus: &Poly<P>) -> Self {
        let reduced = p.rem(modulus).unwrap_or_else(Poly::zero);
        let mut coeffs = [Fp::ZERO; D];
        for (i, &c) in reduced.coefficients().iter().enumerate() {
            if i < D {
                coeffs[i] = c;
            }
        }
        Self { coeffs }
    }

    /// Multiply two extension field elements modulo the irreducible polynomial.
    ///
    /// This is the core operation for extension field arithmetic.
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn mul_mod(&self, other: &Self, modulus: &Poly<P>) -> Self {
        let a_poly = self.to_poly();
        let b_poly = other.to_poly();
        let product = a_poly * b_poly;
        Self::from_poly(&product, modulus)
    }

    /// Compute the multiplicative inverse modulo the irreducible polynomial.
    ///
    /// Uses the extended Euclidean algorithm: if gcd(a, m) = 1,
    /// then s*a + t*m = 1, so s*a ≡ 1 (mod m), meaning s is the inverse.
    ///
    /// Returns `None` if the element is zero.
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn inverse_mod(&self, modulus: &Poly<P>) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let a_poly = self.to_poly();
        let (g, s, _t) = Poly::extended_gcd(&a_poly, modulus);

        // For a field extension with irreducible modulus, gcd should be 1
        // (a constant polynomial) for any non-zero element
        if g.degree() != Some(0) {
            return None;
        }

        // Normalize: if gcd is c (constant), we need s/c as inverse
        let c = g.coeff(0);
        let c_inv = c.inverse()?;
        let inv_poly = s * c_inv;

        Some(Self::from_poly(&inv_poly, modulus))
    }

    /// Divide by another element modulo the irreducible polynomial.
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn div_mod(&self, other: &Self, modulus: &Poly<P>) -> Option<Self> {
        let inv = other.inverse_mod(modulus)?;
        Some(self.mul_mod(&inv, modulus))
    }

    /// Compute self^exp modulo the irreducible polynomial.
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn pow_mod(&self, mut exp: u64, modulus: &Poly<P>) -> Self {
        if exp == 0 {
            return Self::one();
        }

        let mut base = *self;
        let mut result = Self::one();

        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul_mod(&base, modulus);
            }
            base = base.mul_mod(&base, modulus);
            exp >>= 1;
        }

        result
    }

    /// Evaluate the element as a polynomial at a given point in the base field.
    pub fn eval(&self, x: Fp<P>) -> Fp<P> {
        // Horner's method
        let mut result = Fp::ZERO;
        for &coeff in self.coeffs.iter().rev() {
            result = result * x + coeff;
        }
        result
    }

    /// Frobenius endomorphism: x -> x^p.
    ///
    /// In characteristic p, this is an automorphism of the extension field.
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn frobenius(&self, modulus: &Poly<P>) -> Self {
        self.pow_mod(P, modulus)
    }

    /// Compute the norm from the extension field to the base field.
    ///
    /// For F_{p^d}/F_p, the norm is the product of all conjugates:
    /// N(a) = a * a^p * a^{p^2} * ... * a^{p^{d-1}}
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn norm(&self, modulus: &Poly<P>) -> Fp<P> {
        let mut result = *self;
        let mut conjugate = self.frobenius(modulus);

        for _ in 1..D {
            result = result.mul_mod(&conjugate, modulus);
            conjugate = conjugate.frobenius(modulus);
        }

        // The result should be in the base field (degree 0)
        result.coeffs[0]
    }

    /// Compute the trace from the extension field to the base field.
    ///
    /// For F_{p^d}/F_p, the trace is the sum of all conjugates:
    /// Tr(a) = a + a^p + a^{p^2} + ... + a^{p^{d-1}}
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn trace(&self, modulus: &Poly<P>) -> Fp<P> {
        let mut result = *self;
        let mut conjugate = self.frobenius(modulus);

        for _ in 1..D {
            result = result + conjugate;
            conjugate = conjugate.frobenius(modulus);
        }

        // The result should be in the base field (degree 0)
        result.coeffs[0]
    }

    /// Compute all conjugates of this element under the Frobenius automorphism.
    ///
    /// For an element α in F_{p^d}, the conjugates are α, α^p, α^{p^2}, ..., α^{p^{d-1}}.
    /// These are the roots of the minimal polynomial of α over F_p.
    ///
    /// Requires the `alloc` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly, ExtField, Ring};
    ///
    /// type F17 = Fp<17>;
    /// type Ext2 = ExtField<17, 2>;
    ///
    /// // x^2 + 14 = x^2 - 3 is irreducible over F_17
    /// let modulus = Poly::new(vec![F17::new(14), F17::ZERO, F17::ONE]);
    ///
    /// // α = x (the generator)
    /// let alpha = Ext2::new([F17::ZERO, F17::ONE]);
    /// let conjugates = alpha.conjugates(&modulus);
    ///
    /// assert_eq!(conjugates.len(), 2);
    /// assert_eq!(conjugates[0], alpha);
    /// // α^17 is the other conjugate
    /// ```
    #[cfg(feature = "alloc")]
    pub fn conjugates(&self, modulus: &Poly<P>) -> Vec<Self> {
        let mut result = Vec::with_capacity(D);
        let mut current = *self;

        for _ in 0..D {
            result.push(current);
            current = current.frobenius(modulus);
        }

        result
    }

    /// Compute the minimal polynomial of this element over F_p.
    ///
    /// The minimal polynomial is the monic polynomial of smallest degree
    /// with coefficients in F_p that has this element as a root.
    ///
    /// For an element α in F_{p^d}, the minimal polynomial divides x^{p^d} - x
    /// and has degree equal to the smallest k such that α^{p^k} = α.
    ///
    /// Requires the `alloc` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly, ExtField, Ring};
    ///
    /// type F17 = Fp<17>;
    /// type Ext2 = ExtField<17, 2>;
    ///
    /// // x^2 + 14 = x^2 - 3 is irreducible over F_17
    /// let modulus = Poly::new(vec![F17::new(14), F17::ZERO, F17::ONE]);
    ///
    /// // For α = x (the generator), minimal poly should be x^2 - 3
    /// let alpha = Ext2::new([F17::ZERO, F17::ONE]);
    /// let min_poly = alpha.minimal_polynomial(&modulus);
    ///
    /// assert_eq!(min_poly.degree(), Some(2));
    /// // α should be a root of its minimal polynomial
    /// assert_eq!(alpha.eval_poly(&min_poly, &modulus), Ext2::zero());
    /// ```
    #[cfg(feature = "alloc")]
    pub fn minimal_polynomial(&self, modulus: &Poly<P>) -> Poly<P> {
        // The minimal polynomial is the product (x - α)(x - α^p)(x - α^{p^2})...
        // over all distinct conjugates.

        // First find the distinct conjugates (the orbit under Frobenius)
        let mut conjugates = Vec::new();
        let mut current = *self;

        loop {
            conjugates.push(current);
            current = current.frobenius(modulus);
            if current == *self {
                break;
            }
            // Safety check to avoid infinite loop
            if conjugates.len() > D {
                break;
            }
        }

        let k = conjugates.len();

        // Build product (x - α_1)(x - α_2)...(x - α_k) iteratively
        // using polynomial multiplication over the extension field,
        // then extract base field coefficients.
        //
        // We represent polynomials as Vec<Self> where index = degree.

        // Start with (x - α_0)
        let mut poly_coeffs: Vec<Self> = vec![Self::zero() - conjugates[0], Self::one()];

        for i in 1..k {
            // Multiply by (x - α_i)
            let mut new_coeffs = vec![Self::zero(); poly_coeffs.len() + 1];

            // Shift (multiply by x)
            for (j, c) in poly_coeffs.iter().enumerate() {
                new_coeffs[j + 1] = new_coeffs[j + 1] + *c;
            }

            // Subtract α_i times old poly
            for (j, c) in poly_coeffs.iter().enumerate() {
                let term = c.mul_mod(&conjugates[i], modulus);
                new_coeffs[j] = new_coeffs[j] - term;
            }

            poly_coeffs = new_coeffs;
        }

        // Extract base field coefficients
        // The result should have coefficients in F_p (verified by construction
        // since we're multiplying conjugates under Galois action)
        let coeffs: Vec<Fp<P>> = poly_coeffs.iter().map(|c| c.coeffs[0]).collect();

        Poly::new(coeffs)
    }

    /// Evaluate a polynomial at this extension field element.
    ///
    /// This computes p(self) where p is a polynomial over F_p.
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn eval_poly(&self, poly: &Poly<P>, modulus: &Poly<P>) -> Self {
        if poly.is_zero() {
            return Self::zero();
        }

        // Horner's method
        let coefficients = poly.coefficients();
        let mut result = Self::from_base(coefficients[coefficients.len() - 1]);

        for &coeff in coefficients.iter().rev().skip(1) {
            result = result.mul_mod(self, modulus);
            result = result + Self::from_base(coeff);
        }

        result
    }

    /// Compute the multiplicative order of this element in the field.
    ///
    /// The order is the smallest positive integer n such that self^n = 1.
    /// Returns `None` if the element is zero (which has no multiplicative order).
    ///
    /// For F_{p^d}, the multiplicative group has order p^d - 1, so the order
    /// of any element divides p^d - 1.
    ///
    /// Requires the `alloc` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly, ExtField, Ring};
    ///
    /// type F17 = Fp<17>;
    /// type Ext2 = ExtField<17, 2>;
    ///
    /// // x^2 + 14 = x^2 - 3 is irreducible over F_17
    /// let modulus = Poly::new(vec![F17::new(14), F17::ZERO, F17::ONE]);
    ///
    /// // The identity has order 1
    /// assert_eq!(Ext2::one().order(&modulus), Some(1));
    ///
    /// // -1 has order 2
    /// let neg_one = Ext2::from_base(F17::new(16));
    /// assert_eq!(neg_one.order(&modulus), Some(2));
    /// ```
    #[cfg(feature = "alloc")]
    pub fn order(&self, modulus: &Poly<P>) -> Option<u64> {
        if self.is_zero() {
            return None;
        }
        if *self == Self::one() {
            return Some(1);
        }

        // The group order is p^d - 1
        // We need to compute this carefully to avoid overflow
        let mut group_order: u64 = 1;
        for _ in 0..D {
            group_order = group_order.checked_mul(P)?;
        }
        group_order -= 1;

        // Factor the group order and find the smallest divisor n such that self^n = 1
        let factors = Self::factor_u64(group_order);

        // Start with group_order and divide out prime factors while self^n = 1
        let mut order = group_order;
        for (prime, mut exp) in factors {
            while exp > 0 {
                let candidate = order / prime;
                if self.pow_mod(candidate, modulus) == Self::one() {
                    order = candidate;
                    exp -= 1;
                } else {
                    break;
                }
            }
        }

        Some(order)
    }

    /// Check if this element is a primitive element (generator) of the field.
    ///
    /// An element is primitive if its multiplicative order equals p^d - 1,
    /// meaning it generates the entire multiplicative group of the field.
    ///
    /// Requires the `alloc` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Poly, ExtField, Ring};
    ///
    /// type F17 = Fp<17>;
    /// type Ext2 = ExtField<17, 2>;
    ///
    /// // x^2 + 14 = x^2 - 3 is irreducible over F_17
    /// let modulus = Poly::new(vec![F17::new(14), F17::ZERO, F17::ONE]);
    ///
    /// // 1 is not primitive (order 1, not 288)
    /// assert!(!Ext2::one().is_primitive(&modulus));
    ///
    /// // Check if x is primitive
    /// let x = Ext2::new([F17::ZERO, F17::ONE]);
    /// // x has some order dividing 288 = 17^2 - 1
    /// ```
    #[cfg(feature = "alloc")]
    pub fn is_primitive(&self, modulus: &Poly<P>) -> bool {
        if self.is_zero() {
            return false;
        }

        // Compute group order p^d - 1
        let mut group_order: u64 = 1;
        for _ in 0..D {
            if let Some(next) = group_order.checked_mul(P) {
                group_order = next;
            } else {
                return false; // Overflow
            }
        }
        group_order -= 1;

        // An element is primitive iff self^(group_order/q) != 1 for all prime q | group_order
        let factors = Self::factor_u64(group_order);

        for (prime, _) in factors {
            let exp = group_order / prime;
            if self.pow_mod(exp, modulus) == Self::one() {
                return false;
            }
        }

        true
    }

    /// Factor a u64 into prime factors with their multiplicities.
    #[cfg(feature = "alloc")]
    fn factor_u64(mut n: u64) -> Vec<(u64, u32)> {
        let mut factors = Vec::new();

        // Handle factor of 2
        if n % 2 == 0 {
            let mut exp = 0;
            while n % 2 == 0 {
                n /= 2;
                exp += 1;
            }
            factors.push((2, exp));
        }

        // Check odd factors (use d <= n / d to avoid overflow)
        let mut d = 3u64;
        while d <= n / d {
            if n % d == 0 {
                let mut exp = 0;
                while n % d == 0 {
                    n /= d;
                    exp += 1;
                }
                factors.push((d, exp));
            }
            d += 2;
        }

        if n > 1 {
            factors.push((n, 1));
        }

        factors
    }
}

/* ---- Arithmetic operators (without modulus - for additive operations only) ---- */

impl<const P: u64, const D: usize> Add for ExtField<P, D> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut coeffs = [Fp::ZERO; D];
        for (i, c) in coeffs.iter_mut().enumerate() {
            *c = self.coeffs[i] + rhs.coeffs[i];
        }
        Self { coeffs }
    }
}

impl<const P: u64, const D: usize> Sub for ExtField<P, D> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut coeffs = [Fp::ZERO; D];
        for (i, c) in coeffs.iter_mut().enumerate() {
            *c = self.coeffs[i] - rhs.coeffs[i];
        }
        Self { coeffs }
    }
}

impl<const P: u64, const D: usize> Neg for ExtField<P, D> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut coeffs = [Fp::ZERO; D];
        for (i, c) in coeffs.iter_mut().enumerate() {
            *c = -self.coeffs[i];
        }
        Self { coeffs }
    }
}

/// Scalar multiplication by base field element.
impl<const P: u64, const D: usize> Mul<Fp<P>> for ExtField<P, D> {
    type Output = Self;

    fn mul(self, rhs: Fp<P>) -> Self::Output {
        let mut coeffs = [Fp::ZERO; D];
        for (i, c) in coeffs.iter_mut().enumerate() {
            *c = self.coeffs[i] * rhs;
        }
        Self { coeffs }
    }
}

impl<const P: u64, const D: usize> fmt::Debug for ExtField<P, D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Format as array of coefficients
        write!(f, "ExtField{:?}", self.coeffs)
    }
}

impl<const P: u64, const D: usize> fmt::Display for ExtField<P, D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

impl<const P: u64, const D: usize> Default for ExtField<P, D> {
    fn default() -> Self {
        Self::zero()
    }
}

/// A wrapper that bundles an extension field element with its modulus,
/// allowing implementation of Ring and Field traits.
///
/// This is useful when you want to work with extension field elements
/// using the standard Ring/Field trait operations.
///
/// Requires the `alloc` feature.
#[cfg(feature = "alloc")]
#[derive(Clone)]
pub struct ExtFieldWithModulus<'a, const P: u64, const D: usize> {
    elem: ExtField<P, D>,
    modulus: &'a Poly<P>,
}

#[cfg(feature = "alloc")]
impl<'a, const P: u64, const D: usize> ExtFieldWithModulus<'a, P, D> {
    /// Create a new extension field element with its modulus.
    pub fn new(elem: ExtField<P, D>, modulus: &'a Poly<P>) -> Self {
        Self { elem, modulus }
    }

    /// Get the underlying element.
    pub fn elem(&self) -> &ExtField<P, D> {
        &self.elem
    }

    /// Get the modulus.
    pub fn modulus(&self) -> &Poly<P> {
        self.modulus
    }
}

/// Tower extension field F_{p^D1}[y]/(g(y)) where g(y) is irreducible of degree D2.
///
/// This represents a degree-D2 extension over F_{p^D1}, giving a total extension
/// of degree D1*D2 over F_p.
///
/// Elements are represented as D2 coefficients from the base extension F_{p^D1}.
///
/// # Type Parameters
/// - `P`: The prime modulus for the base field F_p
/// - `D1`: The degree of the first extension (F_p → F_{p^D1})
/// - `D2`: The degree of the second extension (F_{p^D1} → F_{p^{D1*D2}})
///
/// # Example
///
/// ```
/// use kreep::{Fp, Poly, ExtField, TowerField};
///
/// type F17 = Fp<17>;
/// type F17_2 = ExtField<17, 2>;  // F_{17^2} = F_289
///
/// // First extension: F_17 → F_{17^2} via x^2 - 3
/// let mod1 = Poly::new(vec![F17::new(14), F17::new(0), F17::new(1)]);
///
/// // Second extension: F_{17^2} → F_{17^4} via y^2 - x
/// // where x is the generator of F_{17^2}
/// // Coefficients of y^2 - x are: [-x, 0, 1] in F_{17^2}
/// let neg_x = F17_2::new([F17::new(0), F17::new(16)]); // -x = 16x in F_17
/// let zero = F17_2::zero();
/// let one = F17_2::one();
///
/// // Create element in the tower: (1 + x) + (2 + 3x)y
/// let a = TowerField::<17, 2, 2>::new([
///     F17_2::new([F17::new(1), F17::new(1)]),   // 1 + x
///     F17_2::new([F17::new(2), F17::new(3)]),   // 2 + 3x
/// ]);
/// ```
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct TowerField<const P: u64, const D1: usize, const D2: usize> {
    /// Coefficients in ascending order: coeffs[i] is the coefficient of y^i
    /// Each coefficient is an element of ExtField<P, D1>
    coeffs: [ExtField<P, D1>; D2],
}

impl<const P: u64, const D1: usize, const D2: usize> TowerField<P, D1, D2> {
    /// Create a new tower field element from coefficients.
    ///
    /// Coefficients are elements of F_{p^D1}, in ascending order of the
    /// outer variable y.
    pub fn new(coeffs: [ExtField<P, D1>; D2]) -> Self {
        Self { coeffs }
    }

    /// Create the zero element.
    pub fn zero() -> Self {
        Self {
            coeffs: [ExtField::zero(); D2],
        }
    }

    /// Create the multiplicative identity (1).
    pub fn one() -> Self {
        let mut coeffs = [ExtField::zero(); D2];
        if D2 > 0 {
            coeffs[0] = ExtField::one();
        }
        Self { coeffs }
    }

    /// Create an element representing the outer variable y.
    ///
    /// Returns `None` if D2 < 2.
    pub fn y() -> Option<Self> {
        if D2 < 2 {
            return None;
        }
        let mut coeffs = [ExtField::zero(); D2];
        coeffs[1] = ExtField::one();
        Some(Self { coeffs })
    }

    /// Create an element from the base extension field.
    pub fn from_base(c: ExtField<P, D1>) -> Self {
        let mut coeffs = [ExtField::zero(); D2];
        if D2 > 0 {
            coeffs[0] = c;
        }
        Self { coeffs }
    }

    /// Create an element from the prime field F_p.
    pub fn from_prime(c: Fp<P>) -> Self {
        Self::from_base(ExtField::from_base(c))
    }

    /// Get the coefficients.
    pub fn coeffs(&self) -> &[ExtField<P, D1>; D2] {
        &self.coeffs
    }

    /// Get a specific coefficient.
    pub fn coeff(&self, i: usize) -> ExtField<P, D1> {
        if i < D2 {
            self.coeffs[i]
        } else {
            ExtField::zero()
        }
    }

    /// Check if this is the zero element.
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|c| c.is_zero())
    }

    /// Multiply two tower field elements.
    ///
    /// Requires both the inner modulus (for F_{p^D1}) and outer modulus (for the tower).
    ///
    /// # Parameters
    /// - `other`: The element to multiply with
    /// - `inner_mod`: Irreducible polynomial for F_p → F_{p^D1}
    /// - `outer_mod`: Coefficients of the irreducible polynomial for F_{p^D1} → F_{p^{D1*D2}}
    ///   Given as D2 elements representing the non-leading coefficients
    ///   (the polynomial is monic, so the leading coefficient y^D2 is implicitly 1)
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn mul_mod(
        &self,
        other: &Self,
        inner_mod: &Poly<P>,
        outer_mod: &[ExtField<P, D1>],
    ) -> Self {
        // Use a Vec for intermediate product to avoid const generic expressions
        // Result before reduction has degree up to 2*(D2-1), so 2*D2-1 coefficients
        let mut product = vec![ExtField::<P, D1>::zero(); 2 * D2 - 1];

        for i in 0..D2 {
            for j in 0..D2 {
                let term = self.coeffs[i].mul_mod(&other.coeffs[j], inner_mod);
                product[i + j] = product[i + j] + term;
            }
        }

        // Reduce modulo the outer polynomial
        // outer_mod represents g(y) = outer_mod[0] + outer_mod[1]*y + ... + y^D2
        // So y^D2 ≡ -(outer_mod[0] + ... + outer_mod[D2-1]*y^{D2-1})
        Self::reduce_tower(&product, inner_mod, outer_mod)
    }

    /// Reduce a polynomial of degree < 2*D2-1 modulo the outer polynomial.
    #[cfg(feature = "alloc")]
    fn reduce_tower(
        poly: &[ExtField<P, D1>],
        inner_mod: &Poly<P>,
        outer_mod: &[ExtField<P, D1>],
    ) -> Self {
        let mut result = [ExtField::<P, D1>::zero(); D2];

        // Copy lower coefficients
        let copy_len = D2.min(poly.len());
        result[..copy_len].clone_from_slice(&poly[..copy_len]);

        // Reduce higher coefficients from highest to lowest
        // y^D2 ≡ -(outer_mod[0] + ... + outer_mod[D2-1]*y^{D2-1})
        for i in (D2..poly.len()).rev() {
            // Get the coefficient we need to reduce (may have been modified by earlier iterations)
            let high_coeff = if i < D2 { result[i] } else { poly[i] };
            if !high_coeff.is_zero() {
                // y^i = y^{i-D2} * y^D2 ≡ -y^{i-D2} * (outer_mod[0] + ... + outer_mod[D2-1]*y^{D2-1})
                let shift = i - D2;
                for (j, om) in outer_mod.iter().enumerate().take(D2) {
                    let term = high_coeff.mul_mod(om, inner_mod);
                    let neg_term = ExtField::zero() - term;
                    let target = shift + j;
                    if target < D2 {
                        result[target] = result[target] + neg_term;
                    }
                    // If target >= D2, we'd need another pass, but for 2*D2-1 input
                    // and reducing from top down, this won't happen on the first pass
                }
            }
        }

        Self { coeffs: result }
    }

    /// Compute the multiplicative inverse.
    ///
    /// For a degree-2 tower extension, we use the formula:
    /// (a + by)^{-1} = (a - by) / (a^2 - b^2 * non_residue)
    /// where non_residue is the constant term of the outer modulus (negated).
    ///
    /// Returns `None` if the element is zero or if D2 != 2 (not yet implemented).
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn inverse_mod(&self, inner_mod: &Poly<P>, outer_mod: &[ExtField<P, D1>]) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        // For D2 = 2: (a + by)^{-1} where y^2 = γ (the non-residue)
        // = (a - by) / (a^2 - b^2 * γ)
        if D2 != 2 {
            // For general D2, we'd need to use the extended Euclidean algorithm
            // on polynomials over F_{p^D1}, which is more complex
            // For now, only D2=2 is supported
            return None;
        }

        let a = self.coeffs[0];
        let b = self.coeffs[1];

        // γ = -outer_mod[0] (since y^2 + outer_mod[0] = 0 means y^2 = -outer_mod[0])
        let gamma = ExtField::zero() - outer_mod[0];

        // norm = a^2 - b^2 * γ
        let a_sq = a.mul_mod(&a, inner_mod);
        let b_sq = b.mul_mod(&b, inner_mod);
        let b_sq_gamma = b_sq.mul_mod(&gamma, inner_mod);
        let norm = a_sq - b_sq_gamma;

        // Invert the norm in the base extension
        let norm_inv = norm.inverse_mod(inner_mod)?;

        // Result = (a - by) * norm_inv
        let neg_b = ExtField::zero() - b;
        let c0 = a.mul_mod(&norm_inv, inner_mod);
        let c1 = neg_b.mul_mod(&norm_inv, inner_mod);

        // Build result array - we know D2 == 2 here
        let mut coeffs = [ExtField::zero(); D2];
        coeffs[0] = c0;
        coeffs[1] = c1;
        Some(Self::new(coeffs))
    }

    /// Compute self^exp.
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn pow_mod(
        &self,
        mut exp: u64,
        inner_mod: &Poly<P>,
        outer_mod: &[ExtField<P, D1>],
    ) -> Self {
        if exp == 0 {
            return Self::one();
        }

        let mut base = *self;
        let mut result = Self::one();

        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul_mod(&base, inner_mod, outer_mod);
            }
            base = base.mul_mod(&base, inner_mod, outer_mod);
            exp >>= 1;
        }

        result
    }

    /// Frobenius endomorphism: x -> x^p.
    ///
    /// Requires the `alloc` feature.
    #[cfg(feature = "alloc")]
    pub fn frobenius(&self, inner_mod: &Poly<P>, outer_mod: &[ExtField<P, D1>]) -> Self {
        self.pow_mod(P, inner_mod, outer_mod)
    }
}

/* ---- Arithmetic operators for TowerField ---- */

impl<const P: u64, const D1: usize, const D2: usize> Add for TowerField<P, D1, D2> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut coeffs = [ExtField::zero(); D2];
        for (i, c) in coeffs.iter_mut().enumerate() {
            *c = self.coeffs[i] + rhs.coeffs[i];
        }
        Self { coeffs }
    }
}

impl<const P: u64, const D1: usize, const D2: usize> Sub for TowerField<P, D1, D2> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut coeffs = [ExtField::zero(); D2];
        for (i, c) in coeffs.iter_mut().enumerate() {
            *c = self.coeffs[i] - rhs.coeffs[i];
        }
        Self { coeffs }
    }
}

impl<const P: u64, const D1: usize, const D2: usize> Neg for TowerField<P, D1, D2> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut coeffs = [ExtField::zero(); D2];
        for (i, c) in coeffs.iter_mut().enumerate() {
            *c = -self.coeffs[i];
        }
        Self { coeffs }
    }
}

/// Scalar multiplication by base extension element.
impl<const P: u64, const D1: usize, const D2: usize> Mul<ExtField<P, D1>>
    for TowerField<P, D1, D2>
{
    type Output = Self;

    fn mul(self, rhs: ExtField<P, D1>) -> Self::Output {
        // Note: This is coefficient-wise, not field multiplication
        // For proper field multiplication, use mul_mod
        let mut coeffs = [ExtField::zero(); D2];
        for (i, c) in coeffs.iter_mut().enumerate() {
            // We can't properly multiply without the inner modulus
            // This is just for scaling by constants where modular reduction isn't needed
            *c = ExtField::new({
                let mut arr = [Fp::ZERO; D1];
                for (j, a) in arr.iter_mut().enumerate() {
                    *a = self.coeffs[i].coeff(j) * rhs.coeff(0);
                }
                arr
            });
        }
        Self { coeffs }
    }
}

impl<const P: u64, const D1: usize, const D2: usize> fmt::Debug for TowerField<P, D1, D2> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first = true;
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }

            if !first {
                write!(f, " + ")?;
            }
            first = false;

            match i {
                0 => write!(f, "({:?})", coeff)?,
                1 => write!(f, "({:?})*y", coeff)?,
                _ => write!(f, "({:?})*y^{}", coeff, i)?,
            }
        }

        if first {
            write!(f, "0")?;
        }

        Ok(())
    }
}

impl<const P: u64, const D1: usize, const D2: usize> fmt::Display for TowerField<P, D1, D2> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

impl<const P: u64, const D1: usize, const D2: usize> Default for TowerField<P, D1, D2> {
    fn default() -> Self {
        Self::zero()
    }
}

// ============================================================================
// Serde implementations
// ============================================================================

#[cfg(all(feature = "serde", feature = "alloc"))]
impl<const P: u64, const D: usize> serde::Serialize for ExtField<P, D> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        // Serialize as an array of coefficient values
        let values: Vec<u64> = self.coeffs.iter().map(|c| c.value()).collect();
        values.serialize(serializer)
    }
}

#[cfg(all(feature = "serde", feature = "alloc"))]
impl<'de, const P: u64, const D: usize> serde::Deserialize<'de> for ExtField<P, D> {
    fn deserialize<De>(deserializer: De) -> Result<Self, De::Error>
    where
        De: serde::Deserializer<'de>,
    {
        let values = Vec::<u64>::deserialize(deserializer)?;
        if values.len() != D {
            return Err(serde::de::Error::custom(alloc::format!(
                "expected {} coefficients, got {}",
                D,
                values.len()
            )));
        }
        let mut coeffs = [Fp::ZERO; D];
        for (i, v) in values.into_iter().enumerate() {
            coeffs[i] = Fp::new(v);
        }
        Ok(Self::new(coeffs))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::format;

    type F17 = Fp<17>;
    type Ext2 = ExtField<17, 2>;

    // x^2 - 3 is irreducible over F_17 (since 3 is not a QR mod 17)
    // In F_17: x^2 - 3 = x^2 + 14
    fn modulus_x2_minus_3() -> Poly<17> {
        Poly::new(vec![F17::new(14), F17::ZERO, F17::ONE]) // -3 + x^2 = 14 + x^2
    }

    #[test]
    fn new_and_coeffs() {
        let a = Ext2::new([F17::new(3), F17::new(5)]);
        assert_eq!(a.coeff(0), F17::new(3));
        assert_eq!(a.coeff(1), F17::new(5));
        assert_eq!(a.coeff(2), F17::ZERO); // out of range
    }

    #[test]
    fn zero_and_one() {
        let z = Ext2::zero();
        assert!(z.is_zero());
        assert_eq!(z.coeff(0), F17::ZERO);
        assert_eq!(z.coeff(1), F17::ZERO);

        let one = Ext2::one();
        assert!(!one.is_zero());
        assert_eq!(one.coeff(0), F17::ONE);
        assert_eq!(one.coeff(1), F17::ZERO);
    }

    #[test]
    fn x_element() {
        let x = Ext2::x().unwrap();
        assert_eq!(x.coeff(0), F17::ZERO);
        assert_eq!(x.coeff(1), F17::ONE);
    }

    #[test]
    fn from_base() {
        let a = Ext2::from_base(F17::new(7));
        assert_eq!(a.coeff(0), F17::new(7));
        assert_eq!(a.coeff(1), F17::ZERO);
    }

    #[test]
    fn add_sub() {
        let a = Ext2::new([F17::new(2), F17::new(3)]);
        let b = Ext2::new([F17::new(5), F17::new(7)]);

        let sum = a + b;
        assert_eq!(sum.coeff(0), F17::new(7));
        assert_eq!(sum.coeff(1), F17::new(10));

        let diff = a - b;
        assert_eq!(diff.coeff(0), F17::new(14)); // 2 - 5 = -3 = 14 mod 17
        assert_eq!(diff.coeff(1), F17::new(13)); // 3 - 7 = -4 = 13 mod 17
    }

    #[test]
    fn neg() {
        let a = Ext2::new([F17::new(2), F17::new(3)]);
        let neg_a = -a;

        assert_eq!(neg_a.coeff(0), F17::new(15)); // -2 = 15 mod 17
        assert_eq!(neg_a.coeff(1), F17::new(14)); // -3 = 14 mod 17
    }

    #[test]
    fn scalar_mul() {
        let a = Ext2::new([F17::new(2), F17::new(3)]);
        let c = F17::new(5);
        let result = a * c;

        assert_eq!(result.coeff(0), F17::new(10));
        assert_eq!(result.coeff(1), F17::new(15));
    }

    #[test]
    fn mul_mod_basic() {
        let m = modulus_x2_minus_3();

        // (1 + x) * (1 + x) = 1 + 2x + x^2
        // Since x^2 ≡ 3 (mod x^2 - 3): 1 + 2x + 3 = 4 + 2x
        let a = Ext2::new([F17::ONE, F17::ONE]);
        let result = a.mul_mod(&a, &m);

        assert_eq!(result.coeff(0), F17::new(4));
        assert_eq!(result.coeff(1), F17::new(2));
    }

    #[test]
    fn mul_mod_x_squared() {
        let m = modulus_x2_minus_3();

        // x * x = x^2 ≡ 3 (mod x^2 - 3)
        let x = Ext2::x().unwrap();
        let x_squared = x.mul_mod(&x, &m);

        assert_eq!(x_squared.coeff(0), F17::new(3));
        assert_eq!(x_squared.coeff(1), F17::ZERO);
    }

    #[test]
    fn mul_mod_associative() {
        let m = modulus_x2_minus_3();

        let a = Ext2::new([F17::new(2), F17::new(3)]);
        let b = Ext2::new([F17::new(5), F17::new(7)]);
        let c = Ext2::new([F17::new(11), F17::new(13)]);

        let ab_c = a.mul_mod(&b, &m).mul_mod(&c, &m);
        let a_bc = a.mul_mod(&b.mul_mod(&c, &m), &m);

        assert_eq!(ab_c, a_bc);
    }

    #[test]
    fn mul_mod_commutative() {
        let m = modulus_x2_minus_3();

        let a = Ext2::new([F17::new(2), F17::new(3)]);
        let b = Ext2::new([F17::new(5), F17::new(7)]);

        assert_eq!(a.mul_mod(&b, &m), b.mul_mod(&a, &m));
    }

    #[test]
    fn mul_mod_identity() {
        let m = modulus_x2_minus_3();

        let a = Ext2::new([F17::new(2), F17::new(3)]);
        let one = Ext2::one();

        assert_eq!(a.mul_mod(&one, &m), a);
        assert_eq!(one.mul_mod(&a, &m), a);
    }

    #[test]
    fn mul_mod_zero() {
        let m = modulus_x2_minus_3();

        let a = Ext2::new([F17::new(2), F17::new(3)]);
        let z = Ext2::zero();

        assert!(a.mul_mod(&z, &m).is_zero());
        assert!(z.mul_mod(&a, &m).is_zero());
    }

    #[test]
    fn inverse_mod_basic() {
        let m = modulus_x2_minus_3();

        let a = Ext2::new([F17::new(2), F17::new(3)]);
        let a_inv = a.inverse_mod(&m).unwrap();

        // a * a^{-1} = 1
        let product = a.mul_mod(&a_inv, &m);
        assert_eq!(product, Ext2::one());
    }

    #[test]
    fn inverse_mod_x() {
        let m = modulus_x2_minus_3();

        // x * x^{-1} = 1
        // Since x^2 = 3, we have x * (x/3) = x^2/3 = 3/3 = 1
        // So x^{-1} = x/3 = x * 3^{-1}
        // 3^{-1} mod 17 = 6 (since 3*6 = 18 = 1 mod 17)
        let x = Ext2::x().unwrap();
        let x_inv = x.inverse_mod(&m).unwrap();

        let product = x.mul_mod(&x_inv, &m);
        assert_eq!(product, Ext2::one());

        // Verify x^{-1} = 6x
        assert_eq!(x_inv.coeff(0), F17::ZERO);
        assert_eq!(x_inv.coeff(1), F17::new(6));
    }

    #[test]
    fn inverse_mod_one() {
        let m = modulus_x2_minus_3();

        let one = Ext2::one();
        let one_inv = one.inverse_mod(&m).unwrap();

        assert_eq!(one_inv, one);
    }

    #[test]
    fn inverse_mod_zero_fails() {
        let m = modulus_x2_minus_3();

        let z = Ext2::zero();
        assert!(z.inverse_mod(&m).is_none());
    }

    #[test]
    fn inverse_mod_all_nonzero() {
        let m = modulus_x2_minus_3();

        // Test several random elements
        let elements = [
            Ext2::new([F17::new(1), F17::new(1)]),
            Ext2::new([F17::new(2), F17::new(5)]),
            Ext2::new([F17::new(7), F17::new(11)]),
            Ext2::new([F17::new(0), F17::new(1)]), // x
            Ext2::new([F17::new(3), F17::new(0)]), // constant
        ];

        for a in &elements {
            let a_inv = a
                .inverse_mod(&m)
                .expect("non-zero element should have inverse");
            let product = a.mul_mod(&a_inv, &m);
            assert_eq!(product, Ext2::one(), "failed for {:?}", a);
        }
    }

    #[test]
    fn div_mod_basic() {
        let m = modulus_x2_minus_3();

        let a = Ext2::new([F17::new(6), F17::new(9)]);
        let b = Ext2::new([F17::new(2), F17::new(3)]);

        // a / b * b = a
        let quotient = a.div_mod(&b, &m).unwrap();
        let product = quotient.mul_mod(&b, &m);
        assert_eq!(product, a);
    }

    #[test]
    fn pow_mod_basic() {
        let m = modulus_x2_minus_3();

        let a = Ext2::new([F17::new(2), F17::new(3)]);

        // a^0 = 1
        assert_eq!(a.pow_mod(0, &m), Ext2::one());

        // a^1 = a
        assert_eq!(a.pow_mod(1, &m), a);

        // a^2 = a * a
        assert_eq!(a.pow_mod(2, &m), a.mul_mod(&a, &m));

        // a^3 = a * a * a
        let a2 = a.mul_mod(&a, &m);
        assert_eq!(a.pow_mod(3, &m), a2.mul_mod(&a, &m));
    }

    #[test]
    fn pow_mod_fermat() {
        let m = modulus_x2_minus_3();

        // In F_{p^2}, every non-zero element satisfies a^{p^2 - 1} = 1
        let a = Ext2::new([F17::new(5), F17::new(7)]);
        let order = 17u64 * 17 - 1; // p^2 - 1 = 288

        let result = a.pow_mod(order, &m);
        assert_eq!(result, Ext2::one());
    }

    #[test]
    fn frobenius_basic() {
        let m = modulus_x2_minus_3();

        // Frobenius: a -> a^p
        // For elements in base field, frobenius is identity
        let base = Ext2::from_base(F17::new(5));
        assert_eq!(base.frobenius(&m), base);

        // For x: x^p = x^17
        // In F_{17}[x]/(x^2-3), x^2 = 3
        // x^17 = x * (x^2)^8 = x * 3^8
        // 3^8 mod 17 = 6561 mod 17 = 16 = -1
        // So x^17 = -x
        let x = Ext2::x().unwrap();
        let x_frob = x.frobenius(&m);
        let neg_x = -x;
        assert_eq!(x_frob, neg_x);
    }

    #[test]
    fn to_poly_roundtrip() {
        let m = modulus_x2_minus_3();

        let a = Ext2::new([F17::new(3), F17::new(7)]);
        let p = a.to_poly();
        let b = Ext2::from_poly(&p, &m);

        assert_eq!(a, b);
    }

    #[test]
    fn from_poly_reduces() {
        let m = modulus_x2_minus_3();

        // Create a polynomial of degree >= 2
        // 5 + 3x + 2x^2 should reduce to 5 + 3x + 2*3 = 11 + 3x (since x^2 = 3)
        let p = Poly::new(vec![F17::new(5), F17::new(3), F17::new(2)]);
        let a = Ext2::from_poly(&p, &m);

        assert_eq!(a.coeff(0), F17::new(11)); // 5 + 6 = 11
        assert_eq!(a.coeff(1), F17::new(3));
    }

    #[test]
    fn distributive() {
        let m = modulus_x2_minus_3();

        let a = Ext2::new([F17::new(2), F17::new(3)]);
        let b = Ext2::new([F17::new(5), F17::new(7)]);
        let c = Ext2::new([F17::new(11), F17::new(13)]);

        // a * (b + c) = a*b + a*c
        let lhs = a.mul_mod(&(b + c), &m);
        let rhs = a.mul_mod(&b, &m) + a.mul_mod(&c, &m);
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn debug_format() {
        let a = Ext2::new([F17::new(3), F17::new(5)]);
        let s = format!("{:?}", a);
        assert_eq!(s, "ExtField[Fp<17>(3), Fp<17>(5)]");
    }

    #[test]
    fn debug_format_zero() {
        let z = Ext2::zero();
        let s = format!("{:?}", z);
        assert_eq!(s, "ExtField[Fp<17>(0), Fp<17>(0)]");
    }

    // Degree 3 extension tests
    type Ext3 = ExtField<17, 3>;

    // x^3 + x + 3 is irreducible over F_17
    // (verified: no roots in F_17)
    fn modulus_cubic() -> Poly<17> {
        Poly::new(vec![F17::new(3), F17::ONE, F17::ZERO, F17::ONE]) // 3 + x + x^3
    }

    #[test]
    fn cubic_extension_basic() {
        let m = modulus_cubic();

        let a = Ext3::new([F17::new(2), F17::new(3), F17::new(5)]);
        let b = Ext3::new([F17::new(7), F17::new(11), F17::new(13)]);

        // Basic operations
        let sum = a + b;
        assert_eq!(sum.coeff(0), F17::new(9));
        assert_eq!(sum.coeff(1), F17::new(14));
        assert_eq!(sum.coeff(2), F17::new(1)); // 18 mod 17

        // Multiplication
        let _product = a.mul_mod(&b, &m);
        // Verify a * a^{-1} = 1
        let a_inv = a.inverse_mod(&m).unwrap();
        assert_eq!(a.mul_mod(&a_inv, &m), Ext3::one());
    }

    #[test]
    fn cubic_extension_x_cubed() {
        let m = modulus_cubic();

        // x^3 ≡ -x - 3 (mod x^3 + x + 3)
        let x = Ext3::x().unwrap();
        let x2 = x.mul_mod(&x, &m);
        let x3 = x2.mul_mod(&x, &m);

        // x^3 = -x - 3 = 16x + 14 (mod 17)
        assert_eq!(x3.coeff(0), F17::new(14)); // -3 mod 17
        assert_eq!(x3.coeff(1), F17::new(16)); // -1 mod 17
        assert_eq!(x3.coeff(2), F17::ZERO);
    }

    // ---- Tower extension tests ----

    type Tower2x2 = TowerField<17, 2, 2>;

    // Inner modulus: x^2 - 3 (irreducible over F_17)
    // Outer modulus: y^2 - x (we need x to not be a square in F_{17^2})

    // For the outer modulus y^2 + c0 + c1*x = 0, we store [c0, c1]
    // y^2 = -c0 - c1*x means y^2 = γ where γ = -c0 - c1*x
    // For y^2 = x, we have c0 = 0, c1 = -1 = 16
    fn outer_modulus_y2_minus_x() -> [Ext2; 2] {
        // y^2 - x = 0, so y^2 + (-x) = 0
        // Coefficients: [0, -1] = [0, 16] in the monic form y^2 + 0 + 16x
        // Actually for our reduction: y^2 = -outer_mod[0] - outer_mod[1]*y
        // We want y^2 = x, so we need outer_mod = [-x, 0] but that's for y^2 + outer_mod[0] = 0
        // Let me reconsider: g(y) = y^2 - x is stored as coefficients [-x, 0, 1]
        // But we only store the non-leading coeffs: [-x, 0]
        [
            Ext2::new([F17::ZERO, F17::new(16)]), // -x = 16x in F_17
            Ext2::zero(),                         // coefficient of y is 0
        ]
    }

    #[test]
    fn tower_zero_one() {
        let z = Tower2x2::zero();
        assert!(z.is_zero());

        let one = Tower2x2::one();
        assert!(!one.is_zero());
        assert_eq!(one.coeff(0), Ext2::one());
        assert_eq!(one.coeff(1), Ext2::zero());
    }

    #[test]
    fn tower_from_base() {
        let base_elem = Ext2::new([F17::new(3), F17::new(5)]);
        let tower_elem = Tower2x2::from_base(base_elem);

        assert_eq!(tower_elem.coeff(0), base_elem);
        assert_eq!(tower_elem.coeff(1), Ext2::zero());
    }

    #[test]
    fn tower_from_prime() {
        let prime_elem = F17::new(7);
        let tower_elem = Tower2x2::from_prime(prime_elem);

        assert_eq!(tower_elem.coeff(0).coeff(0), prime_elem);
        assert_eq!(tower_elem.coeff(0).coeff(1), F17::ZERO);
        assert_eq!(tower_elem.coeff(1), Ext2::zero());
    }

    #[test]
    fn tower_y_element() {
        let y = Tower2x2::y().unwrap();

        assert_eq!(y.coeff(0), Ext2::zero());
        assert_eq!(y.coeff(1), Ext2::one());
    }

    #[test]
    fn tower_add_sub() {
        let a = Tower2x2::new([
            Ext2::new([F17::new(1), F17::new(2)]),
            Ext2::new([F17::new(3), F17::new(4)]),
        ]);
        let b = Tower2x2::new([
            Ext2::new([F17::new(5), F17::new(6)]),
            Ext2::new([F17::new(7), F17::new(8)]),
        ]);

        let sum = a + b;
        assert_eq!(sum.coeff(0), Ext2::new([F17::new(6), F17::new(8)]));
        assert_eq!(sum.coeff(1), Ext2::new([F17::new(10), F17::new(12)]));

        let diff = a - b;
        assert_eq!(diff.coeff(0), Ext2::new([F17::new(13), F17::new(13)])); // -4 = 13
        assert_eq!(diff.coeff(1), Ext2::new([F17::new(13), F17::new(13)]));
    }

    #[test]
    fn tower_neg() {
        let a = Tower2x2::new([
            Ext2::new([F17::new(1), F17::new(2)]),
            Ext2::new([F17::new(3), F17::new(4)]),
        ]);
        let neg_a = -a;

        assert_eq!(neg_a.coeff(0), Ext2::new([F17::new(16), F17::new(15)]));
        assert_eq!(neg_a.coeff(1), Ext2::new([F17::new(14), F17::new(13)]));
    }

    #[test]
    fn tower_mul_identity() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        let a = Tower2x2::new([
            Ext2::new([F17::new(2), F17::new(3)]),
            Ext2::new([F17::new(5), F17::new(7)]),
        ]);
        let one = Tower2x2::one();

        let result = a.mul_mod(&one, &inner_mod, &outer_mod);
        assert_eq!(result, a);

        let result2 = one.mul_mod(&a, &inner_mod, &outer_mod);
        assert_eq!(result2, a);
    }

    #[test]
    fn tower_mul_zero() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        let a = Tower2x2::new([
            Ext2::new([F17::new(2), F17::new(3)]),
            Ext2::new([F17::new(5), F17::new(7)]),
        ]);
        let zero = Tower2x2::zero();

        let result = a.mul_mod(&zero, &inner_mod, &outer_mod);
        assert!(result.is_zero());
    }

    #[test]
    fn tower_mul_y_squared() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        // y^2 should equal x (the generator of F_{17^2})
        let y = Tower2x2::y().unwrap();
        let y_squared = y.mul_mod(&y, &inner_mod, &outer_mod);

        // y^2 = x means coeff(0) = x = (0, 1) and coeff(1) = 0
        assert_eq!(y_squared.coeff(0), Ext2::new([F17::ZERO, F17::ONE]));
        assert_eq!(y_squared.coeff(1), Ext2::zero());
    }

    #[test]
    fn tower_mul_commutative() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        let a = Tower2x2::new([
            Ext2::new([F17::new(2), F17::new(3)]),
            Ext2::new([F17::new(5), F17::new(7)]),
        ]);
        let b = Tower2x2::new([
            Ext2::new([F17::new(11), F17::new(13)]),
            Ext2::new([F17::new(1), F17::new(2)]),
        ]);

        let ab = a.mul_mod(&b, &inner_mod, &outer_mod);
        let ba = b.mul_mod(&a, &inner_mod, &outer_mod);
        assert_eq!(ab, ba);
    }

    #[test]
    fn tower_mul_associative() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        let a = Tower2x2::new([
            Ext2::new([F17::new(2), F17::new(3)]),
            Ext2::new([F17::new(5), F17::new(7)]),
        ]);
        let b = Tower2x2::new([
            Ext2::new([F17::new(11), F17::new(13)]),
            Ext2::new([F17::new(1), F17::new(2)]),
        ]);
        let c = Tower2x2::new([
            Ext2::new([F17::new(4), F17::new(6)]),
            Ext2::new([F17::new(8), F17::new(9)]),
        ]);

        let ab_c = a
            .mul_mod(&b, &inner_mod, &outer_mod)
            .mul_mod(&c, &inner_mod, &outer_mod);
        let a_bc = a.mul_mod(
            &b.mul_mod(&c, &inner_mod, &outer_mod),
            &inner_mod,
            &outer_mod,
        );
        assert_eq!(ab_c, a_bc);
    }

    #[test]
    fn tower_inverse_one() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        let one = Tower2x2::one();
        let one_inv = one.inverse_mod(&inner_mod, &outer_mod).unwrap();

        assert_eq!(one_inv, one);
    }

    #[test]
    fn tower_inverse_zero_fails() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        let zero = Tower2x2::zero();
        assert!(zero.inverse_mod(&inner_mod, &outer_mod).is_none());
    }

    #[test]
    fn tower_inverse_basic() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        let a = Tower2x2::new([
            Ext2::new([F17::new(2), F17::new(3)]),
            Ext2::new([F17::new(5), F17::new(7)]),
        ]);

        let a_inv = a.inverse_mod(&inner_mod, &outer_mod).unwrap();
        let product = a.mul_mod(&a_inv, &inner_mod, &outer_mod);

        assert_eq!(product, Tower2x2::one());
    }

    #[test]
    fn tower_inverse_y() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        // y^{-1} where y^2 = x
        // y * y^{-1} = 1
        // y^{-1} = y / y^2 = y / x = y * x^{-1}
        let y = Tower2x2::y().unwrap();
        let y_inv = y.inverse_mod(&inner_mod, &outer_mod).unwrap();

        let product = y.mul_mod(&y_inv, &inner_mod, &outer_mod);
        assert_eq!(product, Tower2x2::one());
    }

    #[test]
    fn tower_pow_basic() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        let a = Tower2x2::new([
            Ext2::new([F17::new(2), F17::new(3)]),
            Ext2::new([F17::new(5), F17::new(7)]),
        ]);

        // a^0 = 1
        assert_eq!(a.pow_mod(0, &inner_mod, &outer_mod), Tower2x2::one());

        // a^1 = a
        assert_eq!(a.pow_mod(1, &inner_mod, &outer_mod), a);

        // a^2 = a * a
        let a_squared = a.mul_mod(&a, &inner_mod, &outer_mod);
        assert_eq!(a.pow_mod(2, &inner_mod, &outer_mod), a_squared);
    }

    #[test]
    fn tower_pow_fermat() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        // In F_{17^4}, every non-zero element satisfies a^{17^4 - 1} = 1
        let a = Tower2x2::new([
            Ext2::new([F17::new(3), F17::new(5)]),
            Ext2::new([F17::new(7), F17::new(11)]),
        ]);

        // 17^4 - 1 = 83520
        let order = 17u64 * 17 * 17 * 17 - 1;
        let result = a.pow_mod(order, &inner_mod, &outer_mod);

        assert_eq!(result, Tower2x2::one());
    }

    #[test]
    fn tower_distributive() {
        let inner_mod = modulus_x2_minus_3();
        let outer_mod = outer_modulus_y2_minus_x();

        let a = Tower2x2::new([
            Ext2::new([F17::new(2), F17::new(3)]),
            Ext2::new([F17::new(5), F17::new(7)]),
        ]);
        let b = Tower2x2::new([
            Ext2::new([F17::new(11), F17::new(13)]),
            Ext2::new([F17::new(1), F17::new(2)]),
        ]);
        let c = Tower2x2::new([
            Ext2::new([F17::new(4), F17::new(6)]),
            Ext2::new([F17::new(8), F17::new(9)]),
        ]);

        // a * (b + c) = a*b + a*c
        let lhs = a.mul_mod(&(b + c), &inner_mod, &outer_mod);
        let rhs = a.mul_mod(&b, &inner_mod, &outer_mod) + a.mul_mod(&c, &inner_mod, &outer_mod);
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn tower_debug_format() {
        let a = Tower2x2::new([
            Ext2::new([F17::new(1), F17::new(2)]),
            Ext2::new([F17::new(3), F17::new(4)]),
        ]);
        let s = format!("{:?}", a);
        // Should contain both coefficients
        assert!(s.contains("1"));
        assert!(s.contains("y"));
    }

    // ---- Extension field extras tests ----

    #[test]
    fn conjugates_degree2() {
        let modulus = modulus_x2_minus_3();
        let alpha = Ext2::new([F17::ZERO, F17::ONE]); // x

        let conjugates = alpha.conjugates(&modulus);
        assert_eq!(conjugates.len(), 2);
        assert_eq!(conjugates[0], alpha);
        // The second conjugate is α^17
        assert_eq!(conjugates[1], alpha.frobenius(&modulus));
    }

    #[test]
    fn conjugates_base_element() {
        let modulus = modulus_x2_minus_3();
        // An element in the base field has only one distinct conjugate (itself)
        let a = Ext2::from_base(F17::new(5));

        let conjugates = a.conjugates(&modulus);
        // All conjugates of a base field element are the same
        assert_eq!(conjugates[0], a);
        assert_eq!(conjugates[1], a); // a^p = a for a in F_p
    }

    #[test]
    fn minimal_poly_base_element() {
        let modulus = modulus_x2_minus_3();
        // For a ∈ F_p, the minimal polynomial is x - a
        let a = Ext2::from_base(F17::new(5));
        let min_poly = a.minimal_polynomial(&modulus);

        assert_eq!(min_poly.degree(), Some(1));
        // Check it's monic
        assert_eq!(min_poly.leading_coeff(), Some(F17::ONE));
        // Check a is a root
        assert_eq!(a.eval_poly(&min_poly, &modulus), Ext2::zero());
    }

    #[test]
    fn minimal_poly_generator() {
        let modulus = modulus_x2_minus_3();
        // For α = x, the minimal polynomial should be the modulus (x^2 - 3)
        let alpha = Ext2::new([F17::ZERO, F17::ONE]);
        let min_poly = alpha.minimal_polynomial(&modulus);

        assert_eq!(min_poly.degree(), Some(2));
        // Check α is a root
        assert_eq!(alpha.eval_poly(&min_poly, &modulus), Ext2::zero());
        // The minimal poly should be x^2 - 3 = x^2 + 14 (mod 17)
        assert_eq!(min_poly.coeff(0), F17::new(14)); // -3 mod 17
        assert_eq!(min_poly.coeff(1), F17::ZERO);
        assert_eq!(min_poly.coeff(2), F17::ONE);
    }

    #[test]
    fn eval_poly_basic() {
        let modulus = modulus_x2_minus_3();
        // p(x) = 1 + 2x + 3x^2
        let p = Poly::new(vec![F17::new(1), F17::new(2), F17::new(3)]);

        let alpha = Ext2::new([F17::new(1), F17::new(1)]); // 1 + x
        let result = alpha.eval_poly(&p, &modulus);

        // p(1 + x) = 1 + 2(1+x) + 3(1+x)^2
        // (1+x)^2 = 1 + 2x + x^2 = 1 + 2x + 3 = 4 + 2x (since x^2 = 3)
        // = 1 + 2 + 2x + 3(4 + 2x) = 3 + 2x + 12 + 6x = 15 + 8x
        assert_eq!(result.coeffs[0], F17::new(15));
        assert_eq!(result.coeffs[1], F17::new(8));
    }

    #[test]
    fn order_one() {
        let modulus = modulus_x2_minus_3();
        assert_eq!(Ext2::one().order(&modulus), Some(1));
    }

    #[test]
    fn order_minus_one() {
        let modulus = modulus_x2_minus_3();
        let neg_one = Ext2::from_base(F17::new(16)); // -1 mod 17
        assert_eq!(neg_one.order(&modulus), Some(2));
    }

    #[test]
    fn order_zero() {
        let modulus = modulus_x2_minus_3();
        assert_eq!(Ext2::zero().order(&modulus), None);
    }

    #[test]
    fn order_divides_group_order() {
        let modulus = modulus_x2_minus_3();
        let alpha = Ext2::new([F17::ZERO, F17::ONE]); // x
        let order = alpha.order(&modulus).unwrap();

        // Order must divide p^2 - 1 = 288
        assert_eq!(288 % order, 0);
        // And α^order = 1
        assert_eq!(alpha.pow_mod(order, &modulus), Ext2::one());
    }

    #[test]
    fn is_primitive_one() {
        let modulus = modulus_x2_minus_3();
        // 1 is not primitive
        assert!(!Ext2::one().is_primitive(&modulus));
    }

    #[test]
    fn is_primitive_zero() {
        let modulus = modulus_x2_minus_3();
        // 0 is not primitive
        assert!(!Ext2::zero().is_primitive(&modulus));
    }

    #[test]
    fn is_primitive_consistency() {
        let modulus = modulus_x2_minus_3();
        let alpha = Ext2::new([F17::ZERO, F17::ONE]); // x

        // An element is primitive iff its order equals p^d - 1
        let is_prim = alpha.is_primitive(&modulus);
        let order = alpha.order(&modulus).unwrap();

        // Group order is 17^2 - 1 = 288
        assert_eq!(is_prim, order == 288);
    }

    #[test]
    fn find_primitive_element() {
        let modulus = modulus_x2_minus_3();
        // Search for a primitive element
        // In F_{17^2}, the group order is 288 = 2^5 * 3^2
        // We need to find an element of order 288

        // Try x + 1
        let candidate = Ext2::new([F17::ONE, F17::ONE]);
        let order = candidate.order(&modulus).unwrap();

        // Verify the order is correct
        assert_eq!(candidate.pow_mod(order, &modulus), Ext2::one());
        // And no smaller power works
        if order > 1 {
            assert_ne!(candidate.pow_mod(order / 2, &modulus), Ext2::one());
        }
    }

    #[test]
    fn minimal_poly_degree3() {
        let modulus = modulus_cubic();
        type Ext3 = ExtField<17, 3>;

        // For the generator x in a degree-3 extension
        let alpha = Ext3::new([F17::ZERO, F17::ONE, F17::ZERO]);
        let min_poly = alpha.minimal_polynomial(&modulus);

        // Minimal poly should have degree 3 (same as extension degree for generator)
        assert_eq!(min_poly.degree(), Some(3));
        // α should be a root
        assert_eq!(alpha.eval_poly(&min_poly, &modulus), Ext3::zero());
    }
}

#[cfg(all(test, feature = "serde"))]
mod serde_tests {
    use super::*;

    type F17 = Fp<17>;
    type Ext2 = ExtField<17, 2>;
    type Ext3 = ExtField<17, 3>;

    #[test]
    fn serialize_extfield() {
        let a = Ext2::new([F17::new(2), F17::new(3)]);
        let json = serde_json::to_string(&a).unwrap();
        assert_eq!(json, "[2,3]");
    }

    #[test]
    fn deserialize_extfield() {
        let a: Ext2 = serde_json::from_str("[5,7]").unwrap();
        assert_eq!(a.coeff(0), F17::new(5));
        assert_eq!(a.coeff(1), F17::new(7));
    }

    #[test]
    fn roundtrip_extfield() {
        let a = Ext3::new([F17::new(1), F17::new(2), F17::new(3)]);
        let json = serde_json::to_string(&a).unwrap();
        let b: Ext3 = serde_json::from_str(&json).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn deserialize_wrong_length_fails() {
        // 3 elements for Ext2 which expects 2
        let result: Result<Ext2, _> = serde_json::from_str("[1,2,3]");
        assert!(result.is_err());
    }

    #[test]
    fn roundtrip_zero() {
        let a = Ext2::zero();
        let json = serde_json::to_string(&a).unwrap();
        assert_eq!(json, "[0,0]");
        let b: Ext2 = serde_json::from_str(&json).unwrap();
        assert_eq!(a, b);
    }
}
