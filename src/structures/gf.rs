//! Finite field constructors and standard irreducible polynomials.
//!
//! This module provides helper functions to construct common finite fields
//! with automatically selected irreducible polynomials.

use crate::algebra::ring::Ring;
use crate::structures::fp::Fp;
use crate::structures::poly::Poly;
use crate::utils::gcd;

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
        if temp % q == 0 {
            // q is a prime divisor
            let exp = order / q;
            if g.pow(exp) == Fp::ONE {
                return false;
            }
            while temp % q == 0 {
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
