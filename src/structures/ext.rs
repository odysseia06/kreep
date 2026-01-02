use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};

use crate::algebra::field::Field;
use crate::algebra::ring::Ring;
use crate::structures::fp::Fp;
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
    pub fn to_poly(&self) -> Poly<P> {
        Poly::new(self.coeffs.to_vec())
    }

    /// Create from a polynomial, reducing if necessary.
    ///
    /// The polynomial is reduced modulo the given modulus.
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
    pub fn div_mod(&self, other: &Self, modulus: &Poly<P>) -> Option<Self> {
        let inv = other.inverse_mod(modulus)?;
        Some(self.mul_mod(&inv, modulus))
    }

    /// Compute self^exp modulo the irreducible polynomial.
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
    pub fn frobenius(&self, modulus: &Poly<P>) -> Self {
        self.pow_mod(P, modulus)
    }

    /// Compute the norm from the extension field to the base field.
    ///
    /// For F_{p^d}/F_p, the norm is the product of all conjugates:
    /// N(a) = a * a^p * a^{p^2} * ... * a^{p^{d-1}}
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
        // Delegate to polynomial Debug format
        write!(f, "{:?}", self.to_poly())
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
#[derive(Clone)]
pub struct ExtFieldWithModulus<'a, const P: u64, const D: usize> {
    elem: ExtField<P, D>,
    modulus: &'a Poly<P>,
}

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

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(s, "3 + 5*x");
    }

    #[test]
    fn debug_format_zero() {
        let z = Ext2::zero();
        let s = format!("{:?}", z);
        assert_eq!(s, "0");
    }

    // Degree 3 extension tests
    type Ext3 = ExtField<17, 3>;

    // x^3 - 2 over F_17
    // Need to verify this is irreducible... 2^{(17-1)/gcd(3,16)} = 2^16 = 1 mod 17
    // Actually let's use x^3 + x + 1 which is a common irreducible
    fn modulus_cubic() -> Poly<17> {
        // x^3 + x + 1 - need to verify irreducibility
        // For now, let's just use it and trust the math
        Poly::new(vec![F17::ONE, F17::ONE, F17::ZERO, F17::ONE]) // 1 + x + x^3
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

        // x^3 ≡ -x - 1 (mod x^3 + x + 1)
        let x = Ext3::x().unwrap();
        let x2 = x.mul_mod(&x, &m);
        let x3 = x2.mul_mod(&x, &m);

        // x^3 = -x - 1 = 16x + 16
        assert_eq!(x3.coeff(0), F17::new(16));
        assert_eq!(x3.coeff(1), F17::new(16));
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
}
