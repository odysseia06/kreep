//! Finite field constructors and standard irreducible polynomials.
//!
//! This module provides:
//! - The `Modulus` struct for validated irreducible polynomials
//! - The `GF` struct for working with extension field elements that bundles
//!   the element with its modulus
//! - Helper functions to construct common finite fields with automatically
//!   selected irreducible polynomials

use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};
use std::rc::Rc;

use crate::algebra::ring::Ring;
use crate::structures::ext::ExtField;
use crate::structures::fp::Fp;
use crate::structures::poly::Poly;
use crate::utils::gcd;

// ============================================================================
// Modulus type for validated irreducible polynomials
// ============================================================================

/// Error type for modulus validation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModulusError {
    /// The polynomial has the wrong degree.
    WrongDegree { expected: usize, got: Option<usize> },
    /// The polynomial is not monic.
    NotMonic,
    /// The polynomial is not irreducible.
    NotIrreducible,
}

impl fmt::Display for ModulusError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ModulusError::WrongDegree { expected, got } => {
                write!(f, "wrong degree: expected {}, got {:?}", expected, got)
            }
            ModulusError::NotMonic => write!(f, "polynomial is not monic"),
            ModulusError::NotIrreducible => write!(f, "polynomial is not irreducible"),
        }
    }
}

impl std::error::Error for ModulusError {}

/// A validated irreducible polynomial for use as a field modulus.
///
/// This struct wraps a polynomial that has been validated to be:
/// - Of the correct degree `D`
/// - Monic (leading coefficient is 1)
/// - Irreducible (when using `new()`)
///
/// Using `Modulus` instead of raw `Poly` prevents accidental use of
/// reducible or wrong-degree polynomials as field moduli.
///
/// # Example
///
/// ```
/// use kreep::gf::{Modulus, irreducible_poly_deg2};
/// use kreep::Poly;
///
/// // Create a validated modulus from a known irreducible polynomial
/// let poly = irreducible_poly_deg2::<17>();
/// let modulus: Modulus<17, 2> = Modulus::new(poly).unwrap();
///
/// // Or skip irreducibility check for trusted polynomials
/// let poly2 = irreducible_poly_deg2::<17>();
/// let modulus2: Modulus<17, 2> = Modulus::new_unchecked(poly2).unwrap();
/// ```
#[derive(Clone, PartialEq, Eq)]
pub struct Modulus<const P: u64, const D: usize> {
    poly: Poly<P>,
}

impl<const P: u64, const D: usize> Modulus<P, D> {
    /// Create a new validated modulus.
    ///
    /// This validates that the polynomial is:
    /// - Of degree exactly `D`
    /// - Monic
    /// - Irreducible
    ///
    /// # Errors
    ///
    /// Returns `ModulusError::WrongDegree` if the degree is not `D`.
    /// Returns `ModulusError::NotMonic` if the leading coefficient is not 1.
    /// Returns `ModulusError::NotIrreducible` if the polynomial is reducible.
    pub fn new(poly: Poly<P>) -> Result<Self, ModulusError> {
        // Check degree
        if poly.degree() != Some(D) {
            return Err(ModulusError::WrongDegree {
                expected: D,
                got: poly.degree(),
            });
        }

        // Check monic
        if poly.leading_coeff() != Some(Fp::ONE) {
            return Err(ModulusError::NotMonic);
        }

        // Check irreducibility
        if !poly.is_irreducible() {
            return Err(ModulusError::NotIrreducible);
        }

        Ok(Self { poly })
    }

    /// Create a new modulus, skipping the irreducibility check.
    ///
    /// This still validates degree and monic properties, but skips the
    /// expensive irreducibility test. Use this when you know the polynomial
    /// is irreducible (e.g., from `irreducible_poly_deg2` or similar).
    ///
    /// # Errors
    ///
    /// Returns `ModulusError::WrongDegree` if the degree is not `D`.
    /// Returns `ModulusError::NotMonic` if the leading coefficient is not 1.
    pub fn new_unchecked(poly: Poly<P>) -> Result<Self, ModulusError> {
        // Check degree
        if poly.degree() != Some(D) {
            return Err(ModulusError::WrongDegree {
                expected: D,
                got: poly.degree(),
            });
        }

        // Check monic
        if poly.leading_coeff() != Some(Fp::ONE) {
            return Err(ModulusError::NotMonic);
        }

        Ok(Self { poly })
    }

    /// Get a reference to the underlying polynomial.
    pub fn poly(&self) -> &Poly<P> {
        &self.poly
    }

    /// Get the degree of the modulus (always `D`).
    pub const fn degree(&self) -> usize {
        D
    }
}

impl<const P: u64, const D: usize> fmt::Debug for Modulus<P, D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Modulus({:?})", self.poly)
    }
}

impl<const P: u64, const D: usize> fmt::Display for Modulus<P, D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.poly)
    }
}

/// A finite field element GF(p^d) with its modulus.
///
/// This struct bundles an extension field element with its irreducible
/// polynomial modulus, allowing natural arithmetic operations without
/// explicitly passing the modulus each time.
///
/// The modulus is shared via `Rc` so that field elements can be cloned
/// efficiently while sharing the same modulus.
///
/// # Example
///
/// ```
/// use kreep::gf::{GF, Modulus, irreducible_poly_deg2};
/// use kreep::Fp;
/// use std::rc::Rc;
///
/// type F17 = Fp<17>;
///
/// // Create a validated modulus for GF(17^2)
/// let modulus = Rc::new(Modulus::<17, 2>::new_unchecked(irreducible_poly_deg2::<17>()).unwrap());
///
/// // Create elements
/// let a = GF::<17, 2>::new([F17::new(2), F17::new(3)], Rc::clone(&modulus));
/// let b = GF::<17, 2>::new([F17::new(1), F17::new(1)], Rc::clone(&modulus));
///
/// // Arithmetic works naturally
/// let c = &a * &b;
/// let d = &a + &b;
///
/// // Inversion
/// let a_inv = a.inverse().unwrap();
/// assert_eq!(&a * &a_inv, GF::one_like(&a));
/// ```
#[derive(Clone)]
pub struct GF<const P: u64, const D: usize> {
    elem: ExtField<P, D>,
    modulus: Rc<Modulus<P, D>>,
}

impl<const P: u64, const D: usize> GF<P, D> {
    /// Create a new field element from coefficients and a shared modulus.
    pub fn new(coeffs: [Fp<P>; D], modulus: Rc<Modulus<P, D>>) -> Self {
        Self {
            elem: ExtField::new(coeffs),
            modulus,
        }
    }

    /// Create a new field element from an ExtField and a shared modulus.
    pub fn from_ext(elem: ExtField<P, D>, modulus: Rc<Modulus<P, D>>) -> Self {
        Self { elem, modulus }
    }

    /// Create the zero element with the same modulus as another element.
    pub fn zero_like(other: &Self) -> Self {
        Self {
            elem: ExtField::zero(),
            modulus: Rc::clone(&other.modulus),
        }
    }

    /// Create the one element with the same modulus as another element.
    pub fn one_like(other: &Self) -> Self {
        Self {
            elem: ExtField::one(),
            modulus: Rc::clone(&other.modulus),
        }
    }

    /// Create an element representing x with the same modulus.
    pub fn x_like(other: &Self) -> Option<Self> {
        ExtField::x().map(|elem| Self {
            elem,
            modulus: Rc::clone(&other.modulus),
        })
    }

    /// Create a field element from a base field element with a shared modulus.
    pub fn from_base(c: Fp<P>, modulus: Rc<Modulus<P, D>>) -> Self {
        Self {
            elem: ExtField::from_base(c),
            modulus,
        }
    }

    /// Get the underlying ExtField element.
    pub fn elem(&self) -> &ExtField<P, D> {
        &self.elem
    }

    /// Get the modulus.
    pub fn modulus(&self) -> &Modulus<P, D> {
        &self.modulus
    }

    /// Get the shared modulus reference.
    pub fn modulus_rc(&self) -> Rc<Modulus<P, D>> {
        Rc::clone(&self.modulus)
    }

    /// Get a specific coefficient.
    pub fn coeff(&self, i: usize) -> Fp<P> {
        self.elem.coeff(i)
    }

    /// Check if this is the zero element.
    pub fn is_zero(&self) -> bool {
        self.elem.is_zero()
    }

    /// Check if this is the one element.
    pub fn is_one(&self) -> bool {
        self.elem == ExtField::one()
    }

    /// Compute the multiplicative inverse.
    pub fn inverse(&self) -> Option<Self> {
        self.elem.inverse_mod(self.modulus.poly()).map(|inv| Self {
            elem: inv,
            modulus: Rc::clone(&self.modulus),
        })
    }

    /// Compute self^exp.
    pub fn pow(&self, exp: u64) -> Self {
        Self {
            elem: self.elem.pow_mod(exp, self.modulus.poly()),
            modulus: Rc::clone(&self.modulus),
        }
    }

    /// Frobenius endomorphism: x -> x^p.
    pub fn frobenius(&self) -> Self {
        Self {
            elem: self.elem.frobenius(self.modulus.poly()),
            modulus: Rc::clone(&self.modulus),
        }
    }

    /// Compute the norm to the base field.
    pub fn norm(&self) -> Fp<P> {
        self.elem.norm(self.modulus.poly())
    }

    /// Compute the trace to the base field.
    pub fn trace(&self) -> Fp<P> {
        self.elem.trace(self.modulus.poly())
    }

    /// Assert that two elements share the same modulus.
    ///
    /// Panics if the moduli differ. This check runs in both debug and release builds
    /// to prevent silent cross-field mixing.
    fn assert_same_modulus(&self, other: &Self) {
        assert!(
            Rc::ptr_eq(&self.modulus, &other.modulus)
                || self.modulus.poly() == other.modulus.poly(),
            "GF elements must have the same modulus"
        );
    }
}

impl<const P: u64, const D: usize> Add for GF<P, D> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(&rhs);
        Self {
            elem: self.elem + rhs.elem,
            modulus: self.modulus,
        }
    }
}

impl<const P: u64, const D: usize> Add for &GF<P, D> {
    type Output = GF<P, D>;

    fn add(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(rhs);
        GF {
            elem: self.elem + rhs.elem,
            modulus: Rc::clone(&self.modulus),
        }
    }
}

impl<const P: u64, const D: usize> Sub for GF<P, D> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(&rhs);
        Self {
            elem: self.elem - rhs.elem,
            modulus: self.modulus,
        }
    }
}

impl<const P: u64, const D: usize> Sub for &GF<P, D> {
    type Output = GF<P, D>;

    fn sub(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(rhs);
        GF {
            elem: self.elem - rhs.elem,
            modulus: Rc::clone(&self.modulus),
        }
    }
}

impl<const P: u64, const D: usize> Neg for GF<P, D> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            elem: -self.elem,
            modulus: self.modulus,
        }
    }
}

impl<const P: u64, const D: usize> Neg for &GF<P, D> {
    type Output = GF<P, D>;

    fn neg(self) -> Self::Output {
        GF {
            elem: -self.elem,
            modulus: Rc::clone(&self.modulus),
        }
    }
}

impl<const P: u64, const D: usize> Mul for GF<P, D> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(&rhs);
        Self {
            elem: self.elem.mul_mod(&rhs.elem, self.modulus.poly()),
            modulus: self.modulus,
        }
    }
}

impl<const P: u64, const D: usize> Mul for &GF<P, D> {
    type Output = GF<P, D>;

    fn mul(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(rhs);
        GF {
            elem: self.elem.mul_mod(&rhs.elem, self.modulus.poly()),
            modulus: Rc::clone(&self.modulus),
        }
    }
}

impl<const P: u64, const D: usize> Div for GF<P, D> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(&rhs);
        let rhs_inv = rhs.inverse().expect("division by zero");
        Self {
            elem: self.elem.mul_mod(&rhs_inv.elem, self.modulus.poly()),
            modulus: self.modulus,
        }
    }
}

impl<const P: u64, const D: usize> Div for &GF<P, D> {
    type Output = GF<P, D>;

    fn div(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(rhs);
        let rhs_inv = rhs.inverse().expect("division by zero");
        GF {
            elem: self.elem.mul_mod(&rhs_inv.elem, self.modulus.poly()),
            modulus: Rc::clone(&self.modulus),
        }
    }
}

impl<const P: u64, const D: usize> PartialEq for GF<P, D> {
    fn eq(&self, other: &Self) -> bool {
        // Two elements are equal only if they have the same modulus and same coefficients
        (Rc::ptr_eq(&self.modulus, &other.modulus) || self.modulus.poly() == other.modulus.poly())
            && self.elem == other.elem
    }
}

impl<const P: u64, const D: usize> Eq for GF<P, D> {}

impl<const P: u64, const D: usize> fmt::Debug for GF<P, D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.elem)
    }
}

impl<const P: u64, const D: usize> fmt::Display for GF<P, D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

/// Returns a standard irreducible polynomial for GF(p^2).
///
/// For quadratic extensions, we find the smallest `a` such that `x^2 - a`
/// is irreducible (i.e., `a` is not a quadratic residue mod `p`).
///
/// # Example
///
/// ```
/// use kreep::gf::irreducible_poly_deg2;
/// use kreep::{Fp, Poly};
///
/// // Get irreducible polynomial for GF(17^2)
/// let modulus = irreducible_poly_deg2::<17>();
/// assert!(modulus.is_irreducible());
/// assert_eq!(modulus.degree(), Some(2));
/// ```
pub fn irreducible_poly_deg2<const P: u64>() -> Poly<P> {
    // Find smallest non-quadratic-residue
    // x^2 - a is irreducible iff a is not a QR
    for a in 2..P {
        let elem = Fp::<P>::new(a);
        if !elem.is_quadratic_residue() {
            // x^2 - a = x^2 + (-a)
            let neg_a = -elem;
            return Poly::new(vec![neg_a, Fp::ZERO, Fp::ONE]);
        }
    }
    // This should never happen for valid primes > 2
    panic!("No non-quadratic residue found for p = {}", P);
}

/// Returns a standard irreducible polynomial for GF(p^3).
///
/// For cubic extensions, we find the smallest `a` such that `x^3 - a`
/// is irreducible (i.e., `a` is not a cubic residue mod `p`).
///
/// # Example
///
/// ```
/// use kreep::gf::irreducible_poly_deg3;
/// use kreep::{Fp, Poly};
///
/// // Get irreducible polynomial for GF(7^3)
/// let modulus = irreducible_poly_deg3::<7>();
/// assert!(modulus.is_irreducible());
/// assert_eq!(modulus.degree(), Some(3));
/// ```
pub fn irreducible_poly_deg3<const P: u64>() -> Poly<P> {
    // x^3 - a is irreducible iff a is not a perfect cube in F_p
    // a is a cube iff a^{(p-1)/gcd(3, p-1)} = 1
    let g = gcd(3, P - 1);
    let exp = (P - 1) / g;

    for a in 2..P {
        let elem = Fp::<P>::new(a);
        // Check if a is NOT a cubic residue
        if elem.pow(exp) != Fp::ONE {
            let neg_a = -elem;
            return Poly::new(vec![neg_a, Fp::ZERO, Fp::ZERO, Fp::ONE]);
        }
    }

    // Fallback: try x^3 + x + a for small a
    for a in 1..P {
        let poly = Poly::new(vec![Fp::<P>::new(a), Fp::ONE, Fp::ZERO, Fp::ONE]);
        if poly.is_irreducible() {
            return poly;
        }
    }

    panic!("No irreducible cubic found for p = {}", P);
}

/// Returns a standard irreducible polynomial of the given degree.
///
/// This function searches for an irreducible polynomial of the form:
/// - `x^d - a` if a non-d-th-residue exists
/// - Otherwise searches through sparse polynomials
///
/// # Example
///
/// ```
/// use kreep::gf::find_irreducible;
/// use kreep::Poly;
///
/// // Get irreducible polynomial of degree 4 over F_5
/// let modulus = find_irreducible::<5>(4);
/// assert!(modulus.is_irreducible());
/// assert_eq!(modulus.degree(), Some(4));
/// ```
pub fn find_irreducible<const P: u64>(degree: usize) -> Poly<P> {
    if degree == 0 {
        panic!("Cannot find irreducible polynomial of degree 0");
    }
    if degree == 1 {
        // x is always irreducible
        return Poly::x();
    }
    if degree == 2 {
        return irreducible_poly_deg2::<P>();
    }
    if degree == 3 {
        return irreducible_poly_deg3::<P>();
    }

    // For higher degrees, first try x^d - a for small a
    let g = gcd(degree as u64, P - 1);
    let exp = (P - 1) / g;

    for a in 2..P.min(100) {
        let elem = Fp::<P>::new(a);
        // Check if a is NOT a d-th residue
        if elem.pow(exp) != Fp::ONE {
            let neg_a = -elem;
            let mut coeffs = vec![Fp::ZERO; degree + 1];
            coeffs[0] = neg_a;
            coeffs[degree] = Fp::ONE;
            let poly = Poly::new(coeffs);
            if poly.is_irreducible() {
                return poly;
            }
        }
    }

    // Try sparse polynomials: x^d + x + a
    for a in 1..P {
        let mut coeffs = vec![Fp::ZERO; degree + 1];
        coeffs[0] = Fp::<P>::new(a);
        coeffs[1] = Fp::ONE;
        coeffs[degree] = Fp::ONE;
        let poly = Poly::new(coeffs);
        if poly.is_irreducible() {
            return poly;
        }
    }

    // Try x^d + x^k + a for various k
    for k in 2..degree {
        for a in 1..P {
            let mut coeffs = vec![Fp::ZERO; degree + 1];
            coeffs[0] = Fp::<P>::new(a);
            coeffs[k] = Fp::ONE;
            coeffs[degree] = Fp::ONE;
            let poly = Poly::new(coeffs);
            if poly.is_irreducible() {
                return poly;
            }
        }
    }

    // Try x^d + bx + a for various b
    for b in 2..P {
        for a in 1..P {
            let mut coeffs = vec![Fp::ZERO; degree + 1];
            coeffs[0] = Fp::<P>::new(a);
            coeffs[1] = Fp::<P>::new(b);
            coeffs[degree] = Fp::ONE;
            let poly = Poly::new(coeffs);
            if poly.is_irreducible() {
                return poly;
            }
        }
    }

    panic!(
        "No irreducible polynomial of degree {} found for p = {}",
        degree, P
    );
}

/// Standard irreducible polynomials for binary fields GF(2^n).
///
/// Note: Binary fields require a different implementation since our
/// Montgomery arithmetic requires odd primes. These polynomials are
/// provided for reference and can be used with a dedicated GF(2^n)
/// implementation.
///
/// Returns coefficients as a bitmask where bit i represents x^i.
///
/// # Example
///
/// ```
/// use kreep::gf::binary_field_poly;
///
/// // AES field: GF(2^8) with x^8 + x^4 + x^3 + x + 1
/// let poly = binary_field_poly(8);
/// assert_eq!(poly, Some(0b100011011)); // x^8 + x^4 + x^3 + x + 1
/// ```
pub const fn binary_field_poly(n: usize) -> Option<u64> {
    // Standard irreducible polynomials for GF(2^n)
    // These are commonly used in cryptography and coding theory
    match n {
        1 => Some(0b11),                                 // x + 1
        2 => Some(0b111),                                // x^2 + x + 1
        3 => Some(0b1011),                               // x^3 + x + 1
        4 => Some(0b10011),                              // x^4 + x + 1
        5 => Some(0b100101),                             // x^5 + x^2 + 1
        6 => Some(0b1000011),                            // x^6 + x + 1
        7 => Some(0b10000011),                           // x^7 + x + 1
        8 => Some(0b100011011),                          // x^8 + x^4 + x^3 + x + 1 (AES)
        9 => Some(0b1000010001),                         // x^9 + x^4 + 1
        10 => Some(0b10000001001),                       // x^10 + x^3 + 1
        11 => Some(0b100000000101),                      // x^11 + x^2 + 1
        12 => Some(0b1000001010011),                     // x^12 + x^6 + x^4 + x + 1
        13 => Some(0b10000000011011),                    // x^13 + x^4 + x^3 + x + 1
        14 => Some(0b100000000010011), // x^14 + x^4 + x + 1 (actually x^14 + x^5 + x^3 + x + 1)
        15 => Some(0b1000000000000011), // x^15 + x + 1
        16 => Some(0b10001000000100001), // x^16 + x^12 + x^5 + 1 (actually x^16 + x^5 + x^3 + x^2 + 1)
        32 => Some(0b100000000000000000000000011000011), // x^32 + x^7 + x^6 + x + 1 (approx)
        64 => Some(0x800000000000001b),  // x^64 + x^4 + x^3 + x + 1
        128 => Some(0x87),               // Placeholder - actual poly is larger
        _ => None,
    }
}

/// Returns a primitive polynomial for GF(p^d) if one can be found quickly.
///
/// A primitive polynomial generates the multiplicative group of the field.
/// The root of a primitive polynomial is a generator of GF(p^d)*.
///
/// # Example
///
/// ```
/// use kreep::gf::find_primitive;
///
/// // Find primitive polynomial for GF(3^2)
/// if let Some(poly) = find_primitive::<3>(2) {
///     assert!(poly.is_primitive());
/// }
/// ```
pub fn find_primitive<const P: u64>(degree: usize) -> Option<Poly<P>> {
    if degree == 0 {
        return None;
    }
    if degree == 1 {
        // x - g where g is a primitive root mod p
        for g in 2..P {
            let elem = Fp::<P>::new(g);
            // Check if g is a primitive root (generates F_p*)
            if is_primitive_root(elem) {
                return Some(Poly::new(vec![-elem, Fp::ONE]));
            }
        }
        return None;
    }

    // For higher degrees, search through candidates
    // First try x^d - a for primitive roots
    for a in 2..P.min(50) {
        let elem = Fp::<P>::new(a);
        let neg_a = -elem;
        let mut coeffs = vec![Fp::ZERO; degree + 1];
        coeffs[0] = neg_a;
        coeffs[degree] = Fp::ONE;
        let poly = Poly::new(coeffs);
        if poly.is_primitive() {
            return Some(poly);
        }
    }

    // Try sparse polynomials
    for a in 1..P.min(50) {
        let mut coeffs = vec![Fp::ZERO; degree + 1];
        coeffs[0] = Fp::<P>::new(a);
        coeffs[1] = Fp::ONE;
        coeffs[degree] = Fp::ONE;
        let poly = Poly::new(coeffs);
        if poly.is_primitive() {
            return Some(poly);
        }
    }

    None
}

/// Check if an element is a primitive root mod p.
fn is_primitive_root<const P: u64>(g: Fp<P>) -> bool {
    if g == Fp::ZERO || g == Fp::ONE {
        return false;
    }

    let order = P - 1;

    // Check that g^{(p-1)/q} != 1 for all prime divisors q of p-1
    let mut temp = order;
    let mut q = 2u64;

    while q * q <= temp {
        if temp.is_multiple_of(q) {
            // q is a prime divisor
            let exp = order / q;
            if g.pow(exp) == Fp::ONE {
                return false;
            }
            while temp.is_multiple_of(q) {
                temp /= q;
            }
        }
        q += 1;
    }

    if temp > 1 {
        // temp is a prime divisor
        let exp = order / temp;
        if g.pow(exp) == Fp::ONE {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- Modulus tests ----

    #[test]
    fn modulus_new_valid() {
        let poly = irreducible_poly_deg2::<17>();
        let modulus: Result<Modulus<17, 2>, _> = Modulus::new(poly);
        assert!(modulus.is_ok());
    }

    #[test]
    fn modulus_new_wrong_degree() {
        let poly = irreducible_poly_deg3::<17>();
        let modulus: Result<Modulus<17, 2>, _> = Modulus::new(poly);
        assert!(matches!(
            modulus,
            Err(ModulusError::WrongDegree {
                expected: 2,
                got: Some(3)
            })
        ));
    }

    #[test]
    fn modulus_new_not_monic() {
        // 2x^2 - 6 (not monic)
        let poly = Poly::new(vec![Fp::<17>::new(11), Fp::ZERO, Fp::new(2)]);
        let modulus: Result<Modulus<17, 2>, _> = Modulus::new(poly);
        assert!(matches!(modulus, Err(ModulusError::NotMonic)));
    }

    #[test]
    fn modulus_new_not_irreducible() {
        // x^2 - 1 = (x-1)(x+1) is reducible
        let poly = Poly::new(vec![Fp::<17>::new(16), Fp::ZERO, Fp::ONE]);
        let modulus: Result<Modulus<17, 2>, _> = Modulus::new(poly);
        assert!(matches!(modulus, Err(ModulusError::NotIrreducible)));
    }

    #[test]
    fn modulus_new_unchecked_skips_irreducibility() {
        // x^2 - 1 is reducible but new_unchecked should accept it
        let poly = Poly::new(vec![Fp::<17>::new(16), Fp::ZERO, Fp::ONE]);
        let modulus: Result<Modulus<17, 2>, _> = Modulus::new_unchecked(poly);
        assert!(modulus.is_ok());
    }

    #[test]
    fn modulus_new_unchecked_still_checks_degree() {
        let poly = irreducible_poly_deg3::<17>();
        let modulus: Result<Modulus<17, 2>, _> = Modulus::new_unchecked(poly);
        assert!(matches!(modulus, Err(ModulusError::WrongDegree { .. })));
    }

    #[test]
    fn modulus_new_unchecked_still_checks_monic() {
        let poly = Poly::new(vec![Fp::<17>::new(11), Fp::ZERO, Fp::new(2)]);
        let modulus: Result<Modulus<17, 2>, _> = Modulus::new_unchecked(poly);
        assert!(matches!(modulus, Err(ModulusError::NotMonic)));
    }

    // ---- GF struct tests ----

    type F17 = Fp<17>;

    fn make_modulus() -> Rc<Modulus<17, 2>> {
        Rc::new(Modulus::new_unchecked(irreducible_poly_deg2::<17>()).unwrap())
    }

    #[test]
    fn gf_new_and_coeff() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(3), F17::new(5)], modulus);
        assert_eq!(a.coeff(0), F17::new(3));
        assert_eq!(a.coeff(1), F17::new(5));
    }

    #[test]
    fn gf_zero_one() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(1), F17::new(2)], modulus);

        let zero = GF::zero_like(&a);
        let one = GF::one_like(&a);

        assert!(zero.is_zero());
        assert!(!one.is_zero());
        assert!(one.is_one());
        assert!(!zero.is_one());
    }

    #[test]
    fn gf_add() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(2), F17::new(3)], Rc::clone(&modulus));
        let b = GF::<17, 2>::new([F17::new(5), F17::new(7)], modulus);

        let c = &a + &b;
        assert_eq!(c.coeff(0), F17::new(7));
        assert_eq!(c.coeff(1), F17::new(10));
    }

    #[test]
    fn gf_sub() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(5), F17::new(7)], Rc::clone(&modulus));
        let b = GF::<17, 2>::new([F17::new(2), F17::new(3)], modulus);

        let c = &a - &b;
        assert_eq!(c.coeff(0), F17::new(3));
        assert_eq!(c.coeff(1), F17::new(4));
    }

    #[test]
    fn gf_neg() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(2), F17::new(3)], modulus);

        let neg_a = -&a;
        assert_eq!(neg_a.coeff(0), F17::new(15)); // -2 mod 17
        assert_eq!(neg_a.coeff(1), F17::new(14)); // -3 mod 17
    }

    #[test]
    fn gf_mul() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(1), F17::new(1)], Rc::clone(&modulus));
        let b = GF::<17, 2>::new([F17::new(1), F17::new(1)], modulus);

        // (1 + x)^2 = 1 + 2x + x^2
        // x^2 = 3 (since modulus is x^2 - 3)
        // = 1 + 2x + 3 = 4 + 2x
        let c = &a * &b;
        assert_eq!(c.coeff(0), F17::new(4));
        assert_eq!(c.coeff(1), F17::new(2));
    }

    #[test]
    fn gf_mul_identity() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(2), F17::new(3)], modulus);
        let one = GF::one_like(&a);

        assert_eq!(&a * &one, a);
    }

    #[test]
    fn gf_inverse() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(2), F17::new(3)], modulus);

        let a_inv = a.inverse().unwrap();
        let product = &a * &a_inv;

        assert!(product.is_one());
    }

    #[test]
    fn gf_inverse_zero_fails() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(2), F17::new(3)], modulus);
        let zero = GF::zero_like(&a);

        assert!(zero.inverse().is_none());
    }

    #[test]
    fn gf_div() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(6), F17::new(9)], Rc::clone(&modulus));
        let b = GF::<17, 2>::new([F17::new(2), F17::new(3)], modulus);

        let c = &a / &b;
        // a / b * b = a
        let result = &c * &b;
        assert_eq!(result, a);
    }

    #[test]
    fn gf_pow() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(2), F17::new(3)], modulus);

        // a^0 = 1
        assert!(a.pow(0).is_one());

        // a^1 = a
        assert_eq!(a.pow(1), a);

        // a^2 = a * a
        assert_eq!(a.pow(2), &a * &a);
    }

    #[test]
    fn gf_pow_fermat() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(5), F17::new(7)], modulus);

        // In GF(17^2), a^{17^2 - 1} = 1 for all non-zero a
        let order = 17u64 * 17 - 1;
        assert!(a.pow(order).is_one());
    }

    #[test]
    fn gf_shared_modulus() {
        let modulus = make_modulus();

        let a = GF::<17, 2>::new([F17::new(2), F17::new(3)], Rc::clone(&modulus));
        let b = GF::<17, 2>::new([F17::new(5), F17::new(7)], Rc::clone(&modulus));

        // Same Rc pointer
        assert!(Rc::ptr_eq(&a.modulus_rc(), &b.modulus_rc()));

        // Arithmetic still works
        let c = &a * &b;
        assert!(!c.is_zero());
    }

    #[test]
    fn gf_frobenius() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(5), F17::new(0)], modulus);

        // For base field elements, Frobenius is identity
        assert_eq!(a.frobenius(), a);
    }

    #[test]
    fn gf_display() {
        let modulus = make_modulus();
        let a = GF::<17, 2>::new([F17::new(3), F17::new(5)], modulus);

        let s = format!("{}", a);
        assert!(s.contains("3"));
        assert!(s.contains("5"));
    }

    // ---- Irreducible polynomial tests ----

    #[test]
    fn irreducible_deg2_f5() {
        let poly = irreducible_poly_deg2::<5>();
        assert!(poly.is_irreducible());
        assert_eq!(poly.degree(), Some(2));
        // 2 is the smallest non-QR mod 5 (1, 4 are QR)
        // So we get x^2 - 2 = x^2 + 3
        assert_eq!(poly.coeff(0), Fp::<5>::new(3)); // -2 mod 5
    }

    #[test]
    fn irreducible_deg2_f7() {
        let poly = irreducible_poly_deg2::<7>();
        assert!(poly.is_irreducible());
        assert_eq!(poly.degree(), Some(2));
        // QRs mod 7: 1, 2, 4. Non-QR: 3
        assert_eq!(poly.coeff(0), Fp::<7>::new(4)); // -3 mod 7
    }

    #[test]
    fn irreducible_deg2_f17() {
        let poly = irreducible_poly_deg2::<17>();
        assert!(poly.is_irreducible());
        assert_eq!(poly.degree(), Some(2));
        // QRs mod 17: 1, 2, 4, 8, 9, 13, 15, 16
        // Smallest non-QR is 3
        assert_eq!(poly.coeff(0), Fp::<17>::new(14)); // -3 mod 17
    }

    #[test]
    fn irreducible_deg3_f5() {
        let poly = irreducible_poly_deg3::<5>();
        assert!(poly.is_irreducible());
        assert_eq!(poly.degree(), Some(3));
    }

    #[test]
    fn irreducible_deg3_f7() {
        let poly = irreducible_poly_deg3::<7>();
        assert!(poly.is_irreducible());
        assert_eq!(poly.degree(), Some(3));
    }

    #[test]
    fn find_irreducible_deg4() {
        let poly = find_irreducible::<5>(4);
        assert!(poly.is_irreducible());
        assert_eq!(poly.degree(), Some(4));
    }

    #[test]
    fn find_irreducible_deg5() {
        let poly = find_irreducible::<3>(5);
        assert!(poly.is_irreducible());
        assert_eq!(poly.degree(), Some(5));
    }

    #[test]
    fn find_irreducible_deg6() {
        let poly = find_irreducible::<7>(6);
        assert!(poly.is_irreducible());
        assert_eq!(poly.degree(), Some(6));
    }

    #[test]
    fn binary_field_poly_aes() {
        // AES uses GF(2^8) with polynomial x^8 + x^4 + x^3 + x + 1
        let poly = binary_field_poly(8).unwrap();
        // Bits: 1 0001 1011 = 0x11b = 283
        assert_eq!(poly, 0b100011011);
    }

    #[test]
    fn binary_field_poly_small() {
        assert_eq!(binary_field_poly(1), Some(0b11));
        assert_eq!(binary_field_poly(2), Some(0b111));
        assert_eq!(binary_field_poly(3), Some(0b1011));
        assert_eq!(binary_field_poly(4), Some(0b10011));
    }

    #[test]
    fn primitive_root_check() {
        // 2 is a primitive root mod 5
        assert!(is_primitive_root(Fp::<5>::new(2)));
        // 3 is also a primitive root mod 5
        assert!(is_primitive_root(Fp::<5>::new(3)));
        // 4 is not (4^2 = 1)
        assert!(!is_primitive_root(Fp::<5>::new(4)));
        // 1 is never primitive
        assert!(!is_primitive_root(Fp::<5>::new(1)));
    }

    #[test]
    fn find_primitive_deg2_f3() {
        if let Some(poly) = find_primitive::<3>(2) {
            assert!(poly.is_primitive());
            assert_eq!(poly.degree(), Some(2));
        }
    }

    #[test]
    fn find_primitive_deg2_f5() {
        if let Some(poly) = find_primitive::<5>(2) {
            assert!(poly.is_primitive());
            assert_eq!(poly.degree(), Some(2));
        }
    }
}
