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

/// A precomputed NTT plan for a specific size.
///
/// This caches twiddle factors for efficient repeated NTT operations of the
/// same size. Use this when performing many NTT operations with the same size
/// to avoid recomputing roots of unity each time.
///
/// # Example
///
/// ```
/// use kreep::ntt::NttPlan;
/// use kreep::{Fp, Poly};
///
/// type F = Fp<998244353>;
///
/// // Create a plan for size 1024
/// let plan = NttPlan::<998244353>::new(1024).unwrap();
///
/// // Use the plan for multiple multiplications
/// let a = Poly::new(vec![F::new(1), F::new(2), F::new(3)]);
/// let b = Poly::new(vec![F::new(4), F::new(5)]);
/// let c = plan.mul(&a, &b);
///
/// assert_eq!(c, a * b);
/// ```
#[derive(Clone)]
pub struct NttPlan<const P: u64> {
    /// The size of the NTT (power of 2)
    n: usize,
    /// Precomputed twiddle factors for forward NTT
    /// twiddles[k] contains the twiddle factors for stage k (length 2^(k+1))
    twiddles: Vec<Vec<Fp<P>>>,
    /// Precomputed twiddle factors for inverse NTT
    inv_twiddles: Vec<Vec<Fp<P>>>,
    /// Precomputed 1/n for inverse NTT scaling
    n_inv: Fp<P>,
}

impl<const P: u64> NttPlan<P> {
    /// Create a new NTT plan for the given size.
    ///
    /// # Arguments
    /// * `n` - The NTT size (must be a power of 2)
    ///
    /// # Returns
    /// * `Some(plan)` if the prime supports NTT for this size
    /// * `None` if `n` is not a power of 2 or the prime doesn't support this size
    pub fn new(n: usize) -> Option<Self> {
        if n == 0 || !n.is_power_of_two() {
            return None;
        }

        let log_n = n.trailing_zeros();
        let info = ntt_info::<P>()?;

        if log_n > info.max_log2 {
            return None;
        }

        // Get the n-th primitive root of unity
        let omega_max = Fp::<P>::new(info.primitive_root);
        let omega = omega_max.pow(1u64 << (info.max_log2 - log_n));
        let omega_inv = omega.inverse()?;

        // Precompute twiddle factors for each stage
        let mut twiddles = Vec::with_capacity(log_n as usize);
        let mut inv_twiddles = Vec::with_capacity(log_n as usize);

        for k in 0..log_n {
            let len = 1 << (k + 1);
            let half_len = len / 2;

            // ω_len = ω^(n/len) is a primitive len-th root of unity
            let omega_len = omega.pow((n / len) as u64);
            let omega_len_inv = omega_inv.pow((n / len) as u64);

            let mut stage_twiddles = Vec::with_capacity(half_len);
            let mut stage_inv_twiddles = Vec::with_capacity(half_len);

            let mut w = Fp::ONE;
            let mut w_inv = Fp::ONE;

            for _ in 0..half_len {
                stage_twiddles.push(w);
                stage_inv_twiddles.push(w_inv);
                w = w * omega_len;
                w_inv = w_inv * omega_len_inv;
            }

            twiddles.push(stage_twiddles);
            inv_twiddles.push(stage_inv_twiddles);
        }

        let n_inv = Fp::<P>::new(n as u64).inverse()?;

        Some(Self {
            n,
            twiddles,
            inv_twiddles,
            n_inv,
        })
    }

    /// Get the NTT size.
    pub fn size(&self) -> usize {
        self.n
    }

    /// Perform forward NTT in-place using precomputed twiddles.
    ///
    /// # Panics
    /// Panics if `a.len() != self.size()`.
    pub fn ntt(&self, a: &mut [Fp<P>]) {
        assert_eq!(a.len(), self.n, "Input length must match plan size");

        if self.n <= 1 {
            return;
        }

        // Bit-reversal permutation
        bit_reverse_permutation(a);

        // Cooley-Tukey iterative NTT with precomputed twiddles
        for (k, stage_twiddles) in self.twiddles.iter().enumerate() {
            let len = 1 << (k + 1);
            let half_len = len / 2;

            for start in (0..self.n).step_by(len) {
                for (j, &w) in stage_twiddles.iter().enumerate() {
                    let u = a[start + j];
                    let v = a[start + j + half_len] * w;
                    a[start + j] = u + v;
                    a[start + j + half_len] = u - v;
                }
            }
        }
    }

    /// Perform inverse NTT in-place using precomputed twiddles.
    ///
    /// # Panics
    /// Panics if `a.len() != self.size()`.
    pub fn intt(&self, a: &mut [Fp<P>]) {
        assert_eq!(a.len(), self.n, "Input length must match plan size");

        if self.n <= 1 {
            return;
        }

        // Bit-reversal permutation
        bit_reverse_permutation(a);

        // Cooley-Tukey iterative NTT with inverse twiddles
        for (k, stage_twiddles) in self.inv_twiddles.iter().enumerate() {
            let len = 1 << (k + 1);
            let half_len = len / 2;

            for start in (0..self.n).step_by(len) {
                for (j, &w) in stage_twiddles.iter().enumerate() {
                    let u = a[start + j];
                    let v = a[start + j + half_len] * w;
                    a[start + j] = u + v;
                    a[start + j + half_len] = u - v;
                }
            }
        }

        // Scale by 1/n
        for x in a.iter_mut() {
            *x = *x * self.n_inv;
        }
    }

    /// Multiply two polynomials using this NTT plan.
    ///
    /// The result polynomial degree is `deg(a) + deg(b)`, so the required
    /// NTT size is at least `deg(a) + deg(b) + 1`. If the polynomials are
    /// too large for this plan, use a larger plan.
    ///
    /// # Panics
    /// Panics if the result would require a larger NTT size than this plan supports.
    pub fn mul(&self, a: &Poly<P>, b: &Poly<P>) -> Poly<P> {
        if a.is_zero() || b.is_zero() {
            return Poly::zero();
        }

        let result_len = a.degree().unwrap() + b.degree().unwrap() + 1;
        assert!(
            result_len <= self.n,
            "Polynomials too large for this plan. Need {} but plan size is {}",
            result_len,
            self.n
        );

        // Pad polynomials to plan size
        let mut a_coeffs: Vec<Fp<P>> = a.coefficients().to_vec();
        let mut b_coeffs: Vec<Fp<P>> = b.coefficients().to_vec();
        a_coeffs.resize(self.n, Fp::ZERO);
        b_coeffs.resize(self.n, Fp::ZERO);

        // Forward NTT
        self.ntt(&mut a_coeffs);
        self.ntt(&mut b_coeffs);

        // Pointwise multiplication
        for i in 0..self.n {
            a_coeffs[i] = a_coeffs[i] * b_coeffs[i];
        }

        // Inverse NTT
        self.intt(&mut a_coeffs);

        // Trim to actual result length
        a_coeffs.truncate(result_len);

        Poly::new(a_coeffs)
    }

    /// Accumulate a product into an existing NTT-domain buffer.
    ///
    /// Computes `acc += a * b` where all inputs are in NTT domain.
    /// This is useful for computing sums of products efficiently:
    /// `sum_i (a_i * b_i)` can be computed as pointwise operations in NTT domain.
    ///
    /// # Arguments
    /// * `acc` - Accumulator in NTT domain (modified in-place)
    /// * `a` - First operand in NTT domain
    /// * `b` - Second operand in NTT domain
    ///
    /// # Panics
    /// Panics if any slice length doesn't match the plan size.
    pub fn mul_accumulate(&self, acc: &mut [Fp<P>], a: &[Fp<P>], b: &[Fp<P>]) {
        assert_eq!(acc.len(), self.n, "Accumulator length must match plan size");
        assert_eq!(a.len(), self.n, "First operand length must match plan size");
        assert_eq!(
            b.len(),
            self.n,
            "Second operand length must match plan size"
        );

        for i in 0..self.n {
            acc[i] = acc[i] + a[i] * b[i];
        }
    }

    /// Transform a polynomial to NTT domain.
    ///
    /// Returns a vector of NTT values. The input polynomial is padded to the
    /// plan size with zeros.
    ///
    /// # Panics
    /// Panics if the polynomial degree is >= plan size.
    pub fn forward(&self, p: &Poly<P>) -> Vec<Fp<P>> {
        if let Some(deg) = p.degree() {
            assert!(
                deg < self.n,
                "Polynomial degree {} too large for plan size {}",
                deg,
                self.n
            );
        }

        let mut coeffs: Vec<Fp<P>> = p.coefficients().to_vec();
        coeffs.resize(self.n, Fp::ZERO);
        self.ntt(&mut coeffs);
        coeffs
    }

    /// Transform NTT values back to a polynomial.
    ///
    /// # Panics
    /// Panics if `values.len() != self.size()`.
    pub fn backward(&self, values: &[Fp<P>]) -> Poly<P> {
        assert_eq!(values.len(), self.n, "Values length must match plan size");

        let mut coeffs = values.to_vec();
        self.intt(&mut coeffs);
        Poly::new(coeffs)
    }
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

    // ---- NttPlan tests ----

    #[test]
    fn ntt_plan_new() {
        let plan = NttPlan::<NTT_PRIME>::new(1024);
        assert!(plan.is_some());
        assert_eq!(plan.unwrap().size(), 1024);
    }

    #[test]
    fn ntt_plan_invalid_size() {
        // Not a power of 2
        assert!(NttPlan::<NTT_PRIME>::new(100).is_none());
        assert!(NttPlan::<NTT_PRIME>::new(0).is_none());
    }

    #[test]
    fn ntt_plan_too_large() {
        // 2^24 > 2^23, so this should fail for 998244353
        let plan = NttPlan::<NTT_PRIME>::new(1 << 24);
        assert!(plan.is_none());
    }

    #[test]
    fn ntt_plan_roundtrip() {
        let plan = NttPlan::<NTT_PRIME>::new(8).unwrap();

        let original: Vec<F> = (0..8).map(|i| F::new(i + 1)).collect();
        let mut a = original.clone();

        plan.ntt(&mut a);
        plan.intt(&mut a);

        assert_eq!(a, original);
    }

    #[test]
    fn ntt_plan_roundtrip_large() {
        let plan = NttPlan::<NTT_PRIME>::new(1024).unwrap();

        let original: Vec<F> = (0..1024).map(|i| F::new(i as u64)).collect();
        let mut a = original.clone();

        plan.ntt(&mut a);
        plan.intt(&mut a);

        assert_eq!(a, original);
    }

    #[test]
    fn ntt_plan_mul() {
        let plan = NttPlan::<NTT_PRIME>::new(8).unwrap();

        let a = P::new(vec![F::new(1), F::new(2), F::new(3)]);
        let b = P::new(vec![F::new(4), F::new(5)]);

        let c = plan.mul(&a, &b);
        let expected = a.clone() * b.clone();

        assert_eq!(c, expected);
    }

    #[test]
    fn ntt_plan_mul_matches_naive() {
        let plan = NttPlan::<NTT_PRIME>::new(256).unwrap();

        let a = P::new((1..=50).map(F::new).collect());
        let b = P::new((1..=40).map(|i| F::new(i * 2)).collect());

        let plan_result = plan.mul(&a, &b);
        let naive_result = a * b;

        assert_eq!(plan_result, naive_result);
    }

    #[test]
    fn ntt_plan_mul_with_zero() {
        let plan = NttPlan::<NTT_PRIME>::new(8).unwrap();

        let a = P::new(vec![F::new(1), F::new(2)]);
        let zero = P::zero();

        assert_eq!(plan.mul(&a, &zero), P::zero());
        assert_eq!(plan.mul(&zero, &a), P::zero());
    }

    #[test]
    fn ntt_plan_forward_backward() {
        let plan = NttPlan::<NTT_PRIME>::new(16).unwrap();

        let p = P::new(vec![F::new(1), F::new(2), F::new(3), F::new(4)]);
        let ntt_values = plan.forward(&p);
        let recovered = plan.backward(&ntt_values);

        assert_eq!(p, recovered);
    }

    #[test]
    fn ntt_plan_mul_accumulate() {
        let plan = NttPlan::<NTT_PRIME>::new(16).unwrap();

        let a = P::new(vec![F::new(1), F::new(2)]);
        let b = P::new(vec![F::new(3), F::new(4)]);
        let c = P::new(vec![F::new(5), F::new(6)]);
        let d = P::new(vec![F::new(7), F::new(8)]);

        // Compute a*b + c*d using NTT domain accumulation
        let a_ntt = plan.forward(&a);
        let b_ntt = plan.forward(&b);
        let c_ntt = plan.forward(&c);
        let d_ntt = plan.forward(&d);

        let mut acc = vec![F::ZERO; plan.size()];
        plan.mul_accumulate(&mut acc, &a_ntt, &b_ntt);
        plan.mul_accumulate(&mut acc, &c_ntt, &d_ntt);

        let result = plan.backward(&acc);

        // Compare with naive computation
        let expected = a * b + c * d;
        assert_eq!(result, expected);
    }

    #[test]
    fn ntt_plan_matches_regular_ntt() {
        let plan = NttPlan::<NTT_PRIME>::new(8).unwrap();
        let omega = get_root_of_unity::<NTT_PRIME>(8).unwrap();

        let original: Vec<F> = (0..8).map(|i| F::new(i + 1)).collect();

        // Using plan
        let mut a_plan = original.clone();
        plan.ntt(&mut a_plan);

        // Using regular ntt
        let mut a_regular = original.clone();
        ntt(&mut a_regular, omega);

        assert_eq!(a_plan, a_regular);
    }
}
