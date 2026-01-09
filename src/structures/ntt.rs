//! Number Theoretic Transform (NTT) for fast polynomial multiplication.
//!
//! NTT is the finite field analog of FFT, enabling O(n log n) polynomial
//! multiplication instead of O(n²).
//!
//! # Requirements
//!
//! For NTT to work efficiently, the prime `p` must be "NTT-friendly":
//! - p = k · 2^n + 1 for some k, n
//! - This ensures primitive 2^m-th roots of unity exist for m ≤ n
//!
//! # Common NTT-friendly primes
//!
//! - 998244353 = 119 · 2²³ + 1 (max NTT size: 2²³)
//! - 2013265921 = 15 · 2²⁷ + 1 (max NTT size: 2²⁷)
//! - 2281701377 = 17 · 2²⁷ + 1 (max NTT size: 2²⁷)
//! - 3221225473 = 3 · 2³⁰ + 1 (max NTT size: 2³⁰)

use crate::algebra::field::Field;
use crate::algebra::ring::Ring;
use crate::structures::fp::Fp;
use crate::structures::poly::Poly;

/// Information about NTT support for a prime.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NttInfo {
    /// The largest power of 2 that divides p-1.
    /// This is the maximum NTT size supported.
    pub max_log2: u32,
    /// A primitive root of unity of order 2^max_log2.
    pub primitive_root: u64,
}

/// Check if a prime p supports NTT and return relevant information.
///
/// Returns `Some(NttInfo)` if p is NTT-friendly (i.e., p-1 has a large power of 2).
/// Returns `None` if p-1 is odd (no NTT support).
///
/// # Example
///
/// ```
/// use kreep::ntt::ntt_info;
///
/// // 998244353 = 119 * 2^23 + 1
/// let info = ntt_info::<998244353>().unwrap();
/// assert_eq!(info.max_log2, 23);
/// ```
pub fn ntt_info<const P: u64>() -> Option<NttInfo> {
    let mut p_minus_1 = P - 1;
    let mut max_log2 = 0u32;

    // Factor out powers of 2
    while p_minus_1.is_multiple_of(2) {
        p_minus_1 /= 2;
        max_log2 += 1;
    }

    if max_log2 == 0 {
        return None;
    }

    // Find a primitive root modulo p
    // g is a primitive root if g^((p-1)/q) != 1 for all prime q dividing p-1
    let primitive_root = find_primitive_root::<P>()?;

    // Compute the 2^max_log2-th root of unity
    // ω = g^((p-1) / 2^max_log2)
    let exp = (P - 1) >> max_log2;
    let omega = Fp::<P>::new(primitive_root).pow(exp);

    Some(NttInfo {
        max_log2,
        primitive_root: omega.value(),
    })
}

/// Find a primitive root modulo p.
fn find_primitive_root<const P: u64>() -> Option<u64> {
    let p_minus_1 = P - 1;

    // Find prime factors of p-1
    let factors = prime_factors(p_minus_1);

    // Test candidates starting from 2
    'outer: for g in 2..P {
        let g_fp = Fp::<P>::new(g);

        // Check g^((p-1)/q) != 1 for all prime factors q
        for &q in &factors {
            let exp = p_minus_1 / q;
            if g_fp.pow(exp) == Fp::ONE {
                continue 'outer;
            }
        }

        return Some(g);
    }

    None
}

/// Get prime factors of n (without multiplicity).
fn prime_factors(mut n: u64) -> Vec<u64> {
    let mut factors = Vec::new();

    let mut d = 2;
    while d <= n / d {
        if n.is_multiple_of(d) {
            factors.push(d);
            while n.is_multiple_of(d) {
                n /= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        factors.push(n);
    }

    factors
}

/// Compute the NTT of a sequence in-place.
///
/// # Arguments
/// * `a` - The input/output array (must have length that is a power of 2)
/// * `omega` - A primitive n-th root of unity where n = a.len()
///
/// After this function, `a[i]` contains the evaluation of the polynomial
/// at ω^(bit_reverse(i)).
pub fn ntt<const P: u64>(a: &mut [Fp<P>], omega: Fp<P>) {
    let n = a.len();
    debug_assert!(n.is_power_of_two(), "NTT length must be a power of 2");

    if n <= 1 {
        return;
    }

    // Bit-reversal permutation
    bit_reverse_permutation(a);

    // Cooley-Tukey iterative NTT
    let mut len = 2;
    while len <= n {
        // ω_len = ω^(n/len) is a primitive len-th root of unity
        let omega_len = omega.pow((n / len) as u64);

        for start in (0..n).step_by(len) {
            let mut w = Fp::ONE;
            for j in 0..len / 2 {
                let u = a[start + j];
                let v = a[start + j + len / 2] * w;
                a[start + j] = u + v;
                a[start + j + len / 2] = u - v;
                w = w * omega_len;
            }
        }
        len *= 2;
    }
}

/// Compute the inverse NTT of a sequence in-place.
///
/// # Arguments
/// * `a` - The input/output array (must have length that is a power of 2)
/// * `omega` - A primitive n-th root of unity where n = a.len()
///
/// This recovers the original coefficients from NTT values.
pub fn intt<const P: u64>(a: &mut [Fp<P>], omega: Fp<P>) {
    let n = a.len();
    debug_assert!(n.is_power_of_two(), "NTT length must be a power of 2");

    if n <= 1 {
        return;
    }

    // Inverse NTT uses ω^(-1) instead of ω
    let omega_inv = omega.inverse().expect("omega must be invertible");
    ntt(a, omega_inv);

    // Scale by 1/n
    let n_inv = Fp::<P>::new(n as u64)
        .inverse()
        .expect("n must be invertible mod p");
    for x in a.iter_mut() {
        *x = *x * n_inv;
    }
}

/// Bit-reversal permutation.
fn bit_reverse_permutation<T: Copy>(a: &mut [T]) {
    let n = a.len();
    let log_n = n.trailing_zeros();

    for i in 0..n {
        let j = bit_reverse(i, log_n);
        if i < j {
            a.swap(i, j);
        }
    }
}

/// Reverse the lower `bits` bits of `x`.
fn bit_reverse(x: usize, bits: u32) -> usize {
    x.reverse_bits() >> (usize::BITS - bits)
}

/// Multiply two polynomials using NTT.
///
/// This is O(n log n) where n is the sum of the degrees.
///
/// # Panics
///
/// Panics if the prime doesn't support NTT for the required size.
///
/// # Example
///
/// ```
/// use kreep::{Fp, Poly};
/// use kreep::ntt::mul_ntt;
///
/// // Use NTT-friendly prime 998244353
/// type F = Fp<998244353>;
///
/// let a = Poly::new(vec![F::new(1), F::new(2), F::new(3)]); // 1 + 2x + 3x²
/// let b = Poly::new(vec![F::new(4), F::new(5)]);            // 4 + 5x
///
/// let c = mul_ntt(&a, &b);
///
/// // (1 + 2x + 3x²)(4 + 5x) = 4 + 13x + 22x² + 15x³
/// assert_eq!(c, a.clone() * b.clone());
/// ```
pub fn mul_ntt<const P: u64>(a: &Poly<P>, b: &Poly<P>) -> Poly<P> {
    if a.is_zero() || b.is_zero() {
        return Poly::zero();
    }

    let result_len = a.degree().unwrap() + b.degree().unwrap() + 1;

    // Find the smallest power of 2 >= result_len
    let n = result_len.next_power_of_two();

    // Get NTT info for this prime
    let info = ntt_info::<P>().expect("Prime does not support NTT");
    let log_n = n.trailing_zeros();

    assert!(
        log_n <= info.max_log2,
        "Polynomial too large for NTT with this prime. Need 2^{} but max is 2^{}",
        log_n,
        info.max_log2
    );

    // Compute the n-th root of unity
    // ω_n = ω_{2^max_log2}^{2^{max_log2 - log_n}}
    let omega_max = Fp::<P>::new(info.primitive_root);
    let omega = omega_max.pow(1u64 << (info.max_log2 - log_n));

    // Pad polynomials to length n
    let mut a_coeffs: Vec<Fp<P>> = a.coefficients().to_vec();
    let mut b_coeffs: Vec<Fp<P>> = b.coefficients().to_vec();
    a_coeffs.resize(n, Fp::ZERO);
    b_coeffs.resize(n, Fp::ZERO);

    // Forward NTT
    ntt(&mut a_coeffs, omega);
    ntt(&mut b_coeffs, omega);

    // Pointwise multiplication
    for i in 0..n {
        a_coeffs[i] = a_coeffs[i] * b_coeffs[i];
    }

    // Inverse NTT
    intt(&mut a_coeffs, omega);

    // Trim to actual result length
    a_coeffs.truncate(result_len);

    Poly::new(a_coeffs)
}

/// Get the primitive n-th root of unity for NTT.
///
/// Returns `Some(ω)` where ω^n = 1 and ω^k ≠ 1 for 0 < k < n.
/// Returns `None` if the prime doesn't support this NTT size.
///
/// # Arguments
/// * `n` - Must be a power of 2
pub fn get_root_of_unity<const P: u64>(n: usize) -> Option<Fp<P>> {
    if !n.is_power_of_two() || n == 0 {
        return None;
    }

    let info = ntt_info::<P>()?;
    let log_n = n.trailing_zeros();

    if log_n > info.max_log2 {
        return None;
    }

    let omega_max = Fp::<P>::new(info.primitive_root);
    Some(omega_max.pow(1u64 << (info.max_log2 - log_n)))
}

#[cfg(test)]
mod tests {
    use super::*;

    // NTT-friendly prime: 998244353 = 119 * 2^23 + 1
    const NTT_PRIME: u64 = 998244353;
    type F = Fp<NTT_PRIME>;
    type P = Poly<NTT_PRIME>;

    #[test]
    fn ntt_info_998244353() {
        let info = ntt_info::<NTT_PRIME>().unwrap();
        assert_eq!(info.max_log2, 23);

        // Verify the primitive root
        let omega = Fp::<NTT_PRIME>::new(info.primitive_root);

        // ω^(2^23) should equal 1
        let order = 1u64 << 23;
        assert_eq!(omega.pow(order), Fp::ONE);

        // ω^(2^22) should not equal 1
        assert_ne!(omega.pow(order / 2), Fp::ONE);
    }

    #[test]
    fn ntt_info_small_prime() {
        // 17 = 1 * 2^4 + 1, so max_log2 = 4
        let info = ntt_info::<17>().unwrap();
        assert_eq!(info.max_log2, 4);
    }

    #[test]
    fn ntt_info_another_prime() {
        // 2013265921 = 15 * 2^27 + 1
        let info = ntt_info::<2013265921>().unwrap();
        assert_eq!(info.max_log2, 27);
    }

    #[test]
    fn bit_reverse_test() {
        assert_eq!(bit_reverse(0b000, 3), 0b000);
        assert_eq!(bit_reverse(0b001, 3), 0b100);
        assert_eq!(bit_reverse(0b010, 3), 0b010);
        assert_eq!(bit_reverse(0b011, 3), 0b110);
        assert_eq!(bit_reverse(0b100, 3), 0b001);
        assert_eq!(bit_reverse(0b101, 3), 0b101);
        assert_eq!(bit_reverse(0b110, 3), 0b011);
        assert_eq!(bit_reverse(0b111, 3), 0b111);
    }

    #[test]
    fn ntt_roundtrip() {
        let omega = get_root_of_unity::<NTT_PRIME>(8).unwrap();

        let original: Vec<F> = (0..8).map(|i| F::new(i + 1)).collect();
        let mut a = original.clone();

        ntt(&mut a, omega);
        intt(&mut a, omega);

        assert_eq!(a, original);
    }

    #[test]
    fn ntt_roundtrip_large() {
        let n = 1024;
        let omega = get_root_of_unity::<NTT_PRIME>(n).unwrap();

        let original: Vec<F> = (0..n).map(|i| F::new(i as u64)).collect();
        let mut a = original.clone();

        ntt(&mut a, omega);
        intt(&mut a, omega);

        assert_eq!(a, original);
    }

    #[test]
    fn mul_ntt_simple() {
        // (1 + x) * (1 + x) = 1 + 2x + x²
        let a = P::new(vec![F::new(1), F::new(1)]);
        let b = P::new(vec![F::new(1), F::new(1)]);

        let c = mul_ntt(&a, &b);
        let expected = P::new(vec![F::new(1), F::new(2), F::new(1)]);

        assert_eq!(c, expected);
    }

    #[test]
    fn mul_ntt_different_sizes() {
        // (1 + 2x + 3x²) * (4 + 5x) = 4 + 13x + 22x² + 15x³
        let a = P::new(vec![F::new(1), F::new(2), F::new(3)]);
        let b = P::new(vec![F::new(4), F::new(5)]);

        let c = mul_ntt(&a, &b);
        let expected = P::new(vec![F::new(4), F::new(13), F::new(22), F::new(15)]);

        assert_eq!(c, expected);
    }

    #[test]
    fn mul_ntt_matches_naive() {
        let a = P::new(vec![F::new(1), F::new(2), F::new(3), F::new(4), F::new(5)]);
        let b = P::new(vec![F::new(6), F::new(7), F::new(8), F::new(9)]);

        let ntt_result = mul_ntt(&a, &b);
        let naive_result = a * b;

        assert_eq!(ntt_result, naive_result);
    }

    #[test]
    fn mul_ntt_with_zero() {
        let a = P::new(vec![F::new(1), F::new(2)]);
        let zero = P::zero();

        assert_eq!(mul_ntt(&a, &zero), P::zero());
        assert_eq!(mul_ntt(&zero, &a), P::zero());
    }

    #[test]
    fn mul_ntt_constant() {
        let a = P::new(vec![F::new(1), F::new(2), F::new(3)]);
        let c = P::constant(F::new(5));

        let result = mul_ntt(&a, &c);
        let expected = P::new(vec![F::new(5), F::new(10), F::new(15)]);

        assert_eq!(result, expected);
    }

    #[test]
    fn mul_ntt_large() {
        // Test with larger polynomials
        let n = 100;
        let a: P = Poly::new((1..=n).map(|i| F::new(i)).collect());
        let b: P = Poly::new((1..=n).map(|i| F::new(i * 2)).collect());

        let ntt_result = mul_ntt(&a, &b);
        let naive_result = a * b;

        assert_eq!(ntt_result, naive_result);
    }

    #[test]
    fn get_root_of_unity_powers() {
        let n = 8usize;
        let omega = get_root_of_unity::<NTT_PRIME>(n).unwrap();

        // ω^n = 1
        assert_eq!(omega.pow(n as u64), Fp::ONE);

        // ω^(n/2) = -1
        assert_eq!(omega.pow((n / 2) as u64), -Fp::ONE);

        // ω^k ≠ 1 for 0 < k < n
        for k in 1..n {
            assert_ne!(omega.pow(k as u64), Fp::ONE);
        }
    }

    #[test]
    fn get_root_of_unity_none_for_non_power_of_2() {
        assert!(get_root_of_unity::<NTT_PRIME>(3).is_none());
        assert!(get_root_of_unity::<NTT_PRIME>(5).is_none());
        assert!(get_root_of_unity::<NTT_PRIME>(6).is_none());
    }

    #[test]
    fn prime_factors_test() {
        assert_eq!(prime_factors(12), vec![2u64, 3]);
        assert_eq!(prime_factors(100), vec![2u64, 5]);
        assert_eq!(prime_factors(17), vec![17u64]);
        assert_eq!(prime_factors(1), Vec::<u64>::new());
        assert_eq!(prime_factors(2), vec![2u64]);
    }
}
