#[cfg(feature = "alloc")]
use alloc::vec;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;
use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::algebra::field::Field;
use crate::algebra::ring::Ring;
use crate::utils::is_prime;

/// Compute -P^(-1) mod 2^64 using the Newton-Raphson method.
///
/// This is used for Montgomery reduction.
const fn compute_p_inv(p: u64) -> u64 {
    // We want x such that p * x ≡ -1 (mod 2^64)
    // Start with x = 1, then iterate: x = x * (2 - p * x)
    // Each iteration doubles the number of correct bits.
    let mut x: u64 = 1;
    let mut i = 0;
    while i < 6 {
        // 6 iterations: 1 -> 2 -> 4 -> 8 -> 16 -> 32 -> 64 bits
        x = x.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(x)));
        i += 1;
    }
    x.wrapping_neg() // return -x = -p^(-1) mod 2^64
}

/// Compute R mod P where R = 2^64.
const fn compute_r_mod_p(p: u64) -> u64 {
    // R mod P = 2^64 mod P
    // Since 2^64 = (2^64 - P) + P, we have 2^64 mod P = (2^64 - P) mod P
    // But we can't represent 2^64 directly, so we compute it as:
    // 2^64 mod P = (2^63 mod P) * 2 mod P, done carefully
    let r = (1u128 << 64) % (p as u128);
    r as u64
}

/// Compute R^2 mod P where R = 2^64.
const fn compute_r2_mod_p(p: u64) -> u64 {
    // R^2 mod P = 2^128 mod P
    let r2 = (1u128 << 64) % (p as u128);
    let r2 = (r2 * r2) % (p as u128);
    r2 as u64
}

/// Montgomery constants for a given prime P.
struct MontgomeryParams<const P: u64>;

impl<const P: u64> MontgomeryParams<P> {
    /// -P^(-1) mod 2^64
    const P_INV: u64 = compute_p_inv(P);

    /// R mod P where R = 2^64
    const R: u64 = compute_r_mod_p(P);

    /// R^2 mod P
    const R2: u64 = compute_r2_mod_p(P);
}

/// Montgomery reduction: given T < P * R, compute T * R^(-1) mod P.
#[inline]
const fn montgomery_reduce<const P: u64>(t: u128) -> u64 {
    // m = (T * P_INV) mod R
    let m = (t as u64).wrapping_mul(MontgomeryParams::<P>::P_INV);
    // t = (T + m * P) / R
    let t = (t + (m as u128) * (P as u128)) >> 64;
    // Conditional subtraction
    let t = t as u64;
    if t >= P {
        t - P
    } else {
        t
    }
}

/// Convert a value to Montgomery form: a -> aR mod P
#[inline]
const fn to_montgomery<const P: u64>(a: u64) -> u64 {
    // aR mod P = (a * R^2) * R^(-1) mod P = montgomery_reduce(a * R^2)
    montgomery_reduce::<P>((a as u128) * (MontgomeryParams::<P>::R2 as u128))
}

/// Convert from Montgomery form: aR -> a mod P
#[inline]
const fn from_montgomery<const P: u64>(a_mont: u64) -> u64 {
    montgomery_reduce::<P>(a_mont as u128)
}

/// Prime field GF(p) where `p` is a `u64`-sized modulus.
///
/// Internally uses Montgomery representation for efficient multiplication.
/// For correct field behavior, `P` must be an odd prime. Use [`Fp::validate_prime`]
/// at startup to verify, or rely on `debug_assert!` checks during development.
#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct Fp<const P: u64> {
    /// Value stored in Montgomery form: value = aR mod P
    mont: u64,
}

#[cfg(feature = "rand")]
impl<const P: u64> rand::distributions::Distribution<Fp<P>> for rand::distributions::Standard {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Fp<P> {
        Fp::new(rng.gen_range(0..P))
    }
}

#[cfg(feature = "subtle")]
impl<const P: u64> subtle::ConstantTimeEq for Fp<P> {
    fn ct_eq(&self, other: &Self) -> subtle::Choice {
        self.mont.ct_eq(&other.mont)
    }
}

#[cfg(feature = "subtle")]
impl<const P: u64> subtle::ConditionallySelectable for Fp<P> {
    fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
        Self {
            mont: u64::conditional_select(&a.mont, &b.mont, choice),
        }
    }
}

#[cfg(feature = "serde")]
impl<const P: u64> serde::Serialize for Fp<P> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.value().serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, const P: u64> serde::Deserialize<'de> for Fp<P> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let value = u64::deserialize(deserializer)?;
        Ok(Self::new(value))
    }
}

impl<const P: u64> Fp<P> {
    /// Create a new field element from a standard integer.
    ///
    /// In debug builds, this asserts that `P` is an odd prime.
    pub fn new(value: u64) -> Self {
        debug_assert!(is_prime(P), "Fp modulus P={} is not prime", P);
        debug_assert!(P % 2 == 1, "Fp modulus P={} must be odd for Montgomery", P);
        let reduced = value % P;
        Self {
            mont: to_montgomery::<P>(reduced),
        }
    }

    /// Create from a value already in Montgomery form.
    #[inline]
    const fn from_mont(mont: u64) -> Self {
        Self { mont }
    }

    /// Get the representative in `[0, P-1]`.
    pub const fn value(self) -> u64 {
        from_montgomery::<P>(self.mont)
    }

    /// The modulus `p`.
    pub const fn modulus() -> u64 {
        P
    }

    /// Validate that the modulus `P` is a valid odd prime.
    ///
    /// Returns `Ok(())` if `P` is an odd prime, or an error message otherwise.
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
        if P == 2 {
            return Err("modulus P=2 is not supported (must be odd for Montgomery)");
        }
        if !is_prime(P) {
            return Err("modulus P is not prime");
        }
        Ok(())
    }

    /// Compute `self^exp` using square-and-multiply.
    ///
    /// Time complexity: O(log exp) multiplications.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Ring};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let a = F17::new(3);
    /// assert_eq!(a.pow(0), F17::ONE);
    /// assert_eq!(a.pow(1), a);
    /// assert_eq!(a.pow(2), a * a);
    /// assert_eq!(a.pow(16), F17::ONE);  // Fermat: a^(p-1) = 1
    /// ```
    #[inline]
    pub fn pow(self, exp: u64) -> Self {
        if exp == 0 {
            return Self::ONE;
        }

        let mut base = self;
        let mut result = Self::ONE;
        let mut e = exp;

        while e > 0 {
            if e & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            e >>= 1;
        }
        result
    }

    /// Compute `self^exp` where `exp` can be negative.
    ///
    /// For negative exponents, computes `(self^(-1))^|exp|`.
    /// Returns `None` if `self` is zero and `exp` is negative.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Ring, Field};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let a = F17::new(3);
    /// assert_eq!(a.pow_signed(-1), a.inverse());
    /// assert_eq!(a.pow_signed(-2).unwrap(), a.inverse().unwrap() * a.inverse().unwrap());
    /// ```
    pub fn pow_signed(self, exp: i64) -> Option<Self> {
        if exp >= 0 {
            Some(self.pow(exp as u64))
        } else {
            let inv = self.inverse()?;
            Some(inv.pow((-exp) as u64))
        }
    }

    /// Constant-time exponentiation using square-and-multiply.
    ///
    /// Always performs exactly 64 iterations regardless of the exponent value.
    /// Uses conditional selection to avoid timing side channels.
    ///
    /// Requires the `subtle` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Ring};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let a = F17::new(3);
    /// assert_eq!(a.pow_ct(2), a * a);
    /// assert_eq!(a.pow_ct(16), F17::ONE);  // Fermat: a^(p-1) = 1
    /// ```
    #[cfg(feature = "subtle")]
    pub fn pow_ct(self, exp: u64) -> Self {
        use subtle::{Choice, ConditionallySelectable};

        let mut base = self;
        let mut result = Self::ONE;

        for i in 0..64 {
            let bit = Choice::from(((exp >> i) & 1) as u8);
            result = Self::conditional_select(&result, &(result * base), bit);
            base = base * base;
        }
        result
    }

    /// Constant-time multiplicative inverse using Fermat's little theorem.
    ///
    /// Computes `a^(p-2) mod p` which equals `a^(-1)` for non-zero `a`.
    /// Returns `None` for zero (checked in constant time).
    ///
    /// This is slower than the extended Euclidean algorithm but runs in
    /// constant time, making it suitable for cryptographic applications.
    ///
    /// Requires the `subtle` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Ring};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let a = F17::new(3);
    /// let a_inv = a.inverse_ct().unwrap();
    /// assert_eq!(a * a_inv, F17::ONE);
    ///
    /// assert!(F17::ZERO.inverse_ct().is_none());
    /// ```
    #[cfg(feature = "subtle")]
    pub fn inverse_ct(self) -> Option<Self> {
        use subtle::ConstantTimeEq;

        // Check for zero in constant time
        let is_zero = self.mont.ct_eq(&0u64);

        // Compute a^(p-2) mod p using constant-time exponentiation
        let result = self.pow_ct(P - 2);

        // Return None if input was zero
        if bool::from(is_zero) {
            None
        } else {
            Some(result)
        }
    }

    /// Compute the Legendre symbol (a/p).
    ///
    /// Returns:
    /// - `1` if `a` is a quadratic residue (has a square root) and `a != 0`
    /// - `-1` if `a` is a non-residue (no square root)
    /// - `0` if `a == 0`
    ///
    /// Uses Euler's criterion: `a^((p-1)/2) = (a/p) mod p`
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::Fp;
    ///
    /// type F17 = Fp<17>;
    ///
    /// // 2 is a quadratic residue mod 17 (6^2 = 36 = 2 mod 17)
    /// assert_eq!(F17::new(2).legendre(), 1);
    ///
    /// // 3 is not a quadratic residue mod 17
    /// assert_eq!(F17::new(3).legendre(), -1);
    ///
    /// assert_eq!(F17::new(0).legendre(), 0);
    /// ```
    pub fn legendre(self) -> i8 {
        if self.mont == 0 {
            return 0;
        }

        // Euler's criterion: a^((p-1)/2) mod p
        let result = self.pow((P - 1) / 2);

        if result == Self::ONE {
            1
        } else {
            // result == P - 1 (i.e., -1 mod P)
            -1
        }
    }

    /// Check if this element is a quadratic residue (has a square root).
    ///
    /// Returns `true` if there exists some `x` such that `x^2 = self`.
    /// Zero is considered a quadratic residue (with sqrt = 0).
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::Fp;
    ///
    /// type F17 = Fp<17>;
    ///
    /// assert!(F17::new(0).is_quadratic_residue());
    /// assert!(F17::new(1).is_quadratic_residue());
    /// assert!(F17::new(2).is_quadratic_residue());
    /// assert!(!F17::new(3).is_quadratic_residue());
    /// ```
    pub fn is_quadratic_residue(self) -> bool {
        self.legendre() >= 0
    }

    /// Compute a square root using the Tonelli-Shanks algorithm.
    ///
    /// Returns `Some(r)` where `r^2 = self`, or `None` if no square root exists.
    /// When two roots exist (±r), returns the smaller one.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::Fp;
    ///
    /// type F17 = Fp<17>;
    ///
    /// let a = F17::new(2);
    /// let r = a.sqrt().unwrap();
    /// assert_eq!(r * r, a);
    ///
    /// // 3 has no square root mod 17
    /// assert!(F17::new(3).sqrt().is_none());
    /// ```
    pub fn sqrt(self) -> Option<Self> {
        // Handle zero
        if self.mont == 0 {
            return Some(Self::ZERO);
        }

        // Check if quadratic residue
        if self.legendre() != 1 {
            return None;
        }

        // Special case: p ≡ 3 (mod 4)
        // sqrt(a) = a^((p+1)/4)
        if P % 4 == 3 {
            let r = self.pow((P + 1) / 4);
            return Some(self.normalize_sqrt(r));
        }

        // General case: Tonelli-Shanks algorithm
        // Write p - 1 = Q * 2^S where Q is odd
        let mut q = P - 1;
        let mut s = 0u32;
        while q.is_multiple_of(2) {
            q /= 2;
            s += 1;
        }

        // Find a quadratic non-residue z
        let mut z = Self::new(2);
        while z.legendre() != -1 {
            z = z + Self::ONE;
        }

        let mut m = s;
        let mut c = z.pow(q);
        let mut t = self.pow(q);
        let mut r = self.pow(q.div_ceil(2));

        loop {
            if t == Self::ONE {
                return Some(self.normalize_sqrt(r));
            }

            // Find the least i such that t^(2^i) = 1
            let mut i = 1u32;
            let mut t_pow = t * t;
            while t_pow != Self::ONE {
                t_pow = t_pow * t_pow;
                i += 1;
            }

            // Update values
            let b = c.pow(1u64 << (m - i - 1));
            m = i;
            c = b * b;
            t = t * c;
            r = r * b;
        }
    }

    /// Return the canonical (smaller) square root.
    #[inline]
    fn normalize_sqrt(self, r: Self) -> Self {
        let neg_r = -r;
        if r.value() <= neg_r.value() {
            r
        } else {
            neg_r
        }
    }

    /// Batch inversion using Montgomery's trick.
    ///
    /// Computes the multiplicative inverse of each element in the slice
    /// using only one field inversion plus 3(n-1) multiplications.
    ///
    /// Returns `None` if any element is zero.
    ///
    /// Requires the `alloc` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Field, Ring};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let elements = vec![F17::new(2), F17::new(3), F17::new(5)];
    /// let inverses = F17::batch_inverse(&elements).unwrap();
    ///
    /// for (a, a_inv) in elements.iter().zip(inverses.iter()) {
    ///     assert_eq!(*a * *a_inv, F17::ONE);
    /// }
    /// ```
    #[cfg(feature = "alloc")]
    pub fn batch_inverse(elements: &[Self]) -> Option<Vec<Self>> {
        let n = elements.len();
        if n == 0 {
            return Some(Vec::new());
        }

        // Check for zeros and compute partial products
        let mut partials = Vec::with_capacity(n);
        let mut acc = Self::ONE;

        for &elem in elements {
            if elem.mont == 0 {
                return None;
            }
            partials.push(acc);
            acc = acc * elem;
        }

        // Invert the accumulated product
        let mut acc_inv = acc.inverse()?;

        // Work backwards to compute individual inverses
        let mut result = vec![Self::ZERO; n];
        for i in (0..n).rev() {
            result[i] = acc_inv * partials[i];
            acc_inv = acc_inv * elements[i];
        }

        Some(result)
    }

    /// In-place batch inversion using Montgomery's trick.
    ///
    /// Replaces each element with its multiplicative inverse using only
    /// one field inversion plus 3(n-1) multiplications.
    ///
    /// Returns `false` if any element is zero (slice unchanged).
    ///
    /// Requires the `alloc` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::{Fp, Field, Ring};
    ///
    /// type F17 = Fp<17>;
    ///
    /// let mut elements = [F17::new(2), F17::new(3), F17::new(5)];
    /// let originals = elements.clone();
    ///
    /// assert!(F17::batch_inverse_in_place(&mut elements));
    ///
    /// for (a, a_inv) in originals.iter().zip(elements.iter()) {
    ///     assert_eq!(*a * *a_inv, F17::ONE);
    /// }
    /// ```
    #[cfg(feature = "alloc")]
    pub fn batch_inverse_in_place(elements: &mut [Self]) -> bool {
        let n = elements.len();
        if n == 0 {
            return true;
        }

        // Check for zeros first
        for elem in elements.iter() {
            if elem.mont == 0 {
                return false;
            }
        }

        // Compute partial products in-place
        // partials[i] = elements[0] * ... * elements[i-1]
        let mut partials = Vec::with_capacity(n);
        let mut acc = Self::ONE;

        for &elem in elements.iter() {
            partials.push(acc);
            acc = acc * elem;
        }

        // Invert the accumulated product
        let mut acc_inv = match acc.inverse() {
            Some(inv) => inv,
            None => return false,
        };

        // Work backwards to compute individual inverses
        for i in (0..n).rev() {
            let orig = elements[i];
            elements[i] = acc_inv * partials[i];
            acc_inv = acc_inv * orig;
        }

        true
    }

    /// Compute the discrete logarithm of `self` to base `base`.
    ///
    /// Finds the smallest non-negative integer `x` such that `base^x = self`,
    /// if it exists. Returns `None` if no such `x` exists (i.e., `self` is not
    /// in the subgroup generated by `base`).
    ///
    /// Uses baby-step giant-step algorithm with O(√n) time and space,
    /// where n is the group order.
    ///
    /// Requires the `std` feature (uses `HashMap` internally).
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::Fp;
    ///
    /// type F17 = Fp<17>;
    ///
    /// let base = F17::new(3); // 3 is a primitive root mod 17
    /// let target = F17::new(15);
    /// let x = target.discrete_log(base).unwrap();
    /// assert_eq!(base.pow(x), target);
    /// ```
    #[cfg(feature = "std")]
    pub fn discrete_log(self, base: Self) -> Option<u64> {
        self.discrete_log_with_order(base, P - 1)
    }

    /// Compute discrete logarithm when the group order is known.
    ///
    /// This is more efficient when working in a proper subgroup of F_p^*.
    ///
    /// Requires the `std` feature (uses `HashMap` internally).
    ///
    /// # Arguments
    /// * `base` - The base of the logarithm
    /// * `order` - The order of `base` in the multiplicative group
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::Fp;
    ///
    /// type F17 = Fp<17>;
    ///
    /// // 2^4 = 16 ≡ -1, so 2 has order 8 (not 16)
    /// let base = F17::new(2);
    /// let target = F17::new(8); // 2^3 = 8
    /// let x = target.discrete_log_with_order(base, 8).unwrap();
    /// assert_eq!(x, 3);
    /// ```
    #[cfg(feature = "std")]
    pub fn discrete_log_with_order(self, base: Self, order: u64) -> Option<u64> {
        // Handle edge cases
        if base == Self::ZERO {
            return if self == Self::ZERO { Some(1) } else { None };
        }
        if self == Self::ONE {
            return Some(0);
        }
        if self == base {
            return Some(1);
        }

        // Baby-step giant-step algorithm
        // We want to find x such that base^x = self
        // Write x = i*m + j where m = ceil(sqrt(order))
        // Then base^(i*m + j) = self
        // base^j = self * (base^(-m))^i
        //
        // Baby step: compute base^j for j = 0, 1, ..., m-1
        // Giant step: compute self * (base^(-m))^i for i = 0, 1, ..., m-1
        // Find collision

        let m = crate::utils::ceil_sqrt_u64(order);

        // Baby step: build table of (base^j, j) for j = 0..m
        // Use or_insert to keep the smallest j on collisions (when base has small order)
        let mut table = std::collections::HashMap::with_capacity(m as usize);
        let mut power = Self::ONE;
        for j in 0..m {
            table.entry(power.value()).or_insert(j);
            power = power * base;
        }

        // Compute base^(-m)
        let base_inv = base.inverse()?;
        let factor = base_inv.pow(m);

        // Giant step: find i such that self * factor^i is in the table
        let mut gamma = self;
        for i in 0..m {
            if let Some(&j) = table.get(&gamma.value()) {
                let x = i * m + j;
                if x < order {
                    return Some(x);
                }
            }
            gamma = gamma * factor;
        }

        None
    }

    /// Compute the multiplicative order of `self` in F_p^*.
    ///
    /// Returns the smallest positive integer `n` such that `self^n = 1`.
    /// Returns `None` if `self` is zero.
    ///
    /// Requires the `alloc` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::Fp;
    ///
    /// type F17 = Fp<17>;
    ///
    /// assert_eq!(F17::new(2).multiplicative_order(), Some(8));  // 2^8 = 256 ≡ 1
    /// assert_eq!(F17::new(3).multiplicative_order(), Some(16)); // 3 is primitive root
    /// assert_eq!(F17::new(1).multiplicative_order(), Some(1));
    /// assert_eq!(F17::new(16).multiplicative_order(), Some(2)); // 16 ≡ -1
    /// ```
    #[cfg(feature = "alloc")]
    pub fn multiplicative_order(self) -> Option<u64> {
        if self == Self::ZERO {
            return None;
        }
        if self == Self::ONE {
            return Some(1);
        }

        let group_order = P - 1;

        // The order must divide p-1, so we only need to check divisors
        // Factor p-1 and check divisors in increasing order
        let factors = Self::factor_u64(group_order);

        // Start with group_order and divide out prime factors while maintaining self^n = 1
        let mut order = group_order;
        for (prime, mut exp) in factors {
            while exp > 0 {
                let candidate = order / prime;
                if self.pow(candidate) == Self::ONE {
                    order = candidate;
                    exp -= 1;
                } else {
                    break;
                }
            }
        }

        Some(order)
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

        // Check odd factors
        let mut d = 3u64;
        while d * d <= n {
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

    /// Check if `self` is a primitive root modulo P.
    ///
    /// A primitive root has multiplicative order p-1, generating the
    /// entire multiplicative group F_p^*.
    ///
    /// Requires the `alloc` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::Fp;
    ///
    /// type F17 = Fp<17>;
    ///
    /// assert!(F17::new(3).is_primitive_root());
    /// assert!(!F17::new(2).is_primitive_root()); // order 8, not 16
    /// ```
    #[cfg(feature = "alloc")]
    pub fn is_primitive_root(self) -> bool {
        if self == Self::ZERO {
            return false;
        }

        let group_order = P - 1;

        // Check that self^(group_order/q) != 1 for each prime q dividing group_order
        let factors = Self::factor_u64(group_order);

        for (prime, _) in factors {
            let exp = group_order / prime;
            if self.pow(exp) == Self::ONE {
                return false;
            }
        }

        true
    }

    /// Find the smallest primitive root modulo P.
    ///
    /// Requires the `alloc` feature.
    ///
    /// # Example
    ///
    /// ```
    /// use kreep::Fp;
    ///
    /// type F17 = Fp<17>;
    ///
    /// let g = F17::primitive_root().unwrap();
    /// assert!(g.is_primitive_root());
    /// ```
    #[cfg(feature = "alloc")]
    pub fn primitive_root() -> Option<Self> {
        for a in 2..P {
            let elem = Self::new(a);
            if elem.is_primitive_root() {
                return Some(elem);
            }
        }
        None
    }
}

impl<const P: u64> fmt::Debug for Fp<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Fp<{}>({})", P, self.value())
    }
}

impl<const P: u64> fmt::Display for Fp<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value())
    }
}

/* ---- From/Into conversions ---- */

impl<const P: u64> From<u64> for Fp<P> {
    fn from(value: u64) -> Self {
        Self::new(value)
    }
}

impl<const P: u64> From<u32> for Fp<P> {
    fn from(value: u32) -> Self {
        Self::new(value as u64)
    }
}

impl<const P: u64> From<u16> for Fp<P> {
    fn from(value: u16) -> Self {
        Self::new(value as u64)
    }
}

impl<const P: u64> From<u8> for Fp<P> {
    fn from(value: u8) -> Self {
        Self::new(value as u64)
    }
}

impl<const P: u64> From<Fp<P>> for u64 {
    fn from(fp: Fp<P>) -> Self {
        fp.value()
    }
}

/* ---- standard arithmetic operators ---- */

impl<const P: u64> Add for Fp<P> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        let mut sum = self.mont + rhs.mont;
        if sum >= P {
            sum -= P;
        }
        Self::from_mont(sum)
    }
}

impl<const P: u64> Sub for Fp<P> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        if self.mont >= rhs.mont {
            Self::from_mont(self.mont - rhs.mont)
        } else {
            Self::from_mont(self.mont + P - rhs.mont)
        }
    }
}

impl<const P: u64> Mul for Fp<P> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        // (aR) * (bR) = abR^2, then reduce to get abR
        let prod = (self.mont as u128) * (rhs.mont as u128);
        Self::from_mont(montgomery_reduce::<P>(prod))
    }
}

impl<const P: u64> Neg for Fp<P> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        if self.mont == 0 {
            self
        } else {
            Self::from_mont(P - self.mont)
        }
    }
}

/// Division implemented via multiplicative inverse.
impl<const P: u64> Div for Fp<P> {
    type Output = Self;

    #[inline]
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inverse().expect("division by zero in Fp")
    }
}

/* ---- implement Ring ---- */

impl<const P: u64> Ring for Fp<P> {
    const ZERO: Self = Self {
        mont: 0, // 0 in Montgomery form is still 0
    };
    const ONE: Self = Self {
        mont: MontgomeryParams::<P>::R, // 1 in Montgomery form is R mod P
    };
}

/* ---- implement Field ---- */

impl<const P: u64> Field for Fp<P> {
    fn inverse(self) -> Option<Self> {
        if self.mont == 0 {
            return None;
        }

        // Convert out of Montgomery form, compute inverse, convert back
        let a = self.value();
        let m = P as i128;

        let (g, x, _) = egcd(a as i128, m);
        if g != 1 {
            return None;
        }

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
    fn montgomery_constants() {
        // Verify Montgomery constants for P=17
        // R = 2^64 mod 17
        let r = (1u128 << 64) % 17;
        assert_eq!(MontgomeryParams::<17>::R, r as u64);

        // R^2 mod 17
        let r2 = (r * r) % 17;
        assert_eq!(MontgomeryParams::<17>::R2, r2 as u64);

        // Verify P_INV: P * P_INV ≡ -1 (mod 2^64)
        let check = 17u64.wrapping_mul(MontgomeryParams::<17>::P_INV);
        assert_eq!(check, u64::MAX); // -1 mod 2^64
    }

    #[test]
    fn to_from_montgomery() {
        for x in 0u64..17 {
            let fp = F17::new(x);
            assert_eq!(fp.value(), x);
        }
    }

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
    fn zero_and_one() {
        assert_eq!(F17::ZERO.value(), 0);
        assert_eq!(F17::ONE.value(), 1);
    }

    #[test]
    fn validate_prime_ok() {
        assert!(Fp::<3>::validate_prime().is_ok());
        assert!(Fp::<17>::validate_prime().is_ok());
        assert!(Fp::<101>::validate_prime().is_ok());
    }

    #[test]
    fn validate_prime_err() {
        assert!(Fp::<1>::validate_prime().is_err());
        assert!(Fp::<2>::validate_prime().is_err()); // P=2 not supported
        assert!(Fp::<4>::validate_prime().is_err());
        assert!(Fp::<15>::validate_prime().is_err());
        assert!(Fp::<100>::validate_prime().is_err());
    }

    #[test]
    fn larger_prime() {
        type F101 = Fp<101>;
        let a = F101::new(50);
        let b = F101::new(60);
        assert_eq!((a + b).value(), 9); // 50 + 60 = 110 ≡ 9 (mod 101)
        assert_eq!((a * b).value(), (3000 % 101) as u64);
    }

    #[test]
    fn batch_inverse_basic() {
        let elements: Vec<F17> = (1u64..17).map(F17::new).collect();
        let inverses = F17::batch_inverse(&elements).unwrap();

        assert_eq!(elements.len(), inverses.len());
        for (a, a_inv) in elements.iter().zip(inverses.iter()) {
            assert_eq!(*a * *a_inv, F17::ONE);
        }
    }

    #[test]
    fn batch_inverse_empty() {
        let elements: Vec<F17> = vec![];
        let inverses = F17::batch_inverse(&elements).unwrap();
        assert!(inverses.is_empty());
    }

    #[test]
    fn batch_inverse_single() {
        let elements = vec![F17::new(5)];
        let inverses = F17::batch_inverse(&elements).unwrap();
        assert_eq!(inverses.len(), 1);
        assert_eq!(elements[0] * inverses[0], F17::ONE);
    }

    #[test]
    fn batch_inverse_with_zero_fails() {
        let elements = vec![F17::new(3), F17::ZERO, F17::new(5)];
        assert!(F17::batch_inverse(&elements).is_none());
    }

    #[test]
    fn batch_inverse_in_place_basic() {
        let originals: Vec<F17> = (1u64..17).map(F17::new).collect();
        let mut elements = originals.clone();

        assert!(F17::batch_inverse_in_place(&mut elements));

        for (a, a_inv) in originals.iter().zip(elements.iter()) {
            assert_eq!(*a * *a_inv, F17::ONE);
        }
    }

    #[test]
    fn batch_inverse_in_place_empty() {
        let mut elements: [F17; 0] = [];
        assert!(F17::batch_inverse_in_place(&mut elements));
    }

    #[test]
    fn batch_inverse_in_place_with_zero_fails() {
        let mut elements = [F17::new(3), F17::ZERO, F17::new(5)];
        let original = elements.clone();

        assert!(!F17::batch_inverse_in_place(&mut elements));
        // Elements should be unchanged
        assert_eq!(elements, original);
    }

    #[test]
    fn batch_inverse_matches_individual() {
        let elements: Vec<F17> = vec![F17::new(2), F17::new(7), F17::new(11)];
        let batch_inv = F17::batch_inverse(&elements).unwrap();
        let individual_inv: Vec<F17> = elements.iter().map(|x| x.inverse().unwrap()).collect();

        assert_eq!(batch_inv, individual_inv);
    }

    #[test]
    fn pow_basic() {
        let a = F17::new(3);

        assert_eq!(a.pow(0), F17::ONE);
        assert_eq!(a.pow(1), a);
        assert_eq!(a.pow(2), a * a);
        assert_eq!(a.pow(3), a * a * a);
    }

    #[test]
    fn pow_zero_base() {
        assert_eq!(F17::ZERO.pow(0), F17::ONE); // 0^0 = 1 by convention
        assert_eq!(F17::ZERO.pow(1), F17::ZERO);
        assert_eq!(F17::ZERO.pow(100), F17::ZERO);
    }

    #[test]
    fn pow_one_base() {
        assert_eq!(F17::ONE.pow(0), F17::ONE);
        assert_eq!(F17::ONE.pow(1), F17::ONE);
        assert_eq!(F17::ONE.pow(1000), F17::ONE);
    }

    #[test]
    fn pow_fermat_little_theorem() {
        // a^(p-1) = 1 for all nonzero a in Fp
        for x in 1u64..17 {
            let a = F17::new(x);
            assert_eq!(a.pow(16), F17::ONE);
        }
    }

    #[test]
    fn pow_large_exponent() {
        let a = F17::new(5);
        // 5^100 mod 17
        // Using Fermat: 5^16 = 1, so 5^100 = 5^(16*6 + 4) = 5^4
        let expected = a.pow(4);
        assert_eq!(a.pow(100), expected);
    }

    #[test]
    fn pow_signed_positive() {
        let a = F17::new(3);
        assert_eq!(a.pow_signed(0), Some(F17::ONE));
        assert_eq!(a.pow_signed(1), Some(a));
        assert_eq!(a.pow_signed(2), Some(a * a));
    }

    #[test]
    fn pow_signed_negative() {
        let a = F17::new(3);
        let a_inv = a.inverse().unwrap();

        assert_eq!(a.pow_signed(-1), Some(a_inv));
        assert_eq!(a.pow_signed(-2), Some(a_inv * a_inv));
        assert_eq!(a.pow_signed(-3), Some(a_inv * a_inv * a_inv));
    }

    #[test]
    fn pow_signed_zero_base() {
        assert_eq!(F17::ZERO.pow_signed(0), Some(F17::ONE));
        assert_eq!(F17::ZERO.pow_signed(5), Some(F17::ZERO));
        assert_eq!(F17::ZERO.pow_signed(-1), None); // 0^(-1) undefined
    }

    #[test]
    fn pow_inverse_relationship() {
        // a^(-n) * a^n = 1
        let a = F17::new(7);
        for n in 1i64..10 {
            let pos = a.pow_signed(n).unwrap();
            let neg = a.pow_signed(-n).unwrap();
            assert_eq!(pos * neg, F17::ONE);
        }
    }

    #[test]
    fn legendre_zero() {
        assert_eq!(F17::ZERO.legendre(), 0);
    }

    #[test]
    fn legendre_residues() {
        // Quadratic residues mod 17: 1, 2, 4, 8, 9, 13, 15, 16
        let residues = [1, 2, 4, 8, 9, 13, 15, 16];
        for &r in &residues {
            assert_eq!(F17::new(r).legendre(), 1, "{} should be a residue", r);
        }
    }

    #[test]
    fn legendre_non_residues() {
        // Non-residues mod 17: 3, 5, 6, 7, 10, 11, 12, 14
        let non_residues = [3, 5, 6, 7, 10, 11, 12, 14];
        for &n in &non_residues {
            assert_eq!(F17::new(n).legendre(), -1, "{} should be a non-residue", n);
        }
    }

    #[test]
    fn is_quadratic_residue_basic() {
        assert!(F17::ZERO.is_quadratic_residue());
        assert!(F17::ONE.is_quadratic_residue());
        assert!(F17::new(4).is_quadratic_residue()); // 2^2 = 4
        assert!(!F17::new(3).is_quadratic_residue());
    }

    #[test]
    fn sqrt_zero() {
        assert_eq!(F17::ZERO.sqrt(), Some(F17::ZERO));
    }

    #[test]
    fn sqrt_one() {
        assert_eq!(F17::ONE.sqrt(), Some(F17::ONE));
    }

    #[test]
    fn sqrt_perfect_squares() {
        // Test all perfect squares in F17
        for x in 0u64..17 {
            let a = F17::new(x);
            let a_sq = a * a;
            let r = a_sq.sqrt().expect("perfect square should have sqrt");
            assert_eq!(
                r * r,
                a_sq,
                "sqrt({})^2 should equal {}",
                a_sq.value(),
                a_sq.value()
            );
        }
    }

    #[test]
    fn sqrt_non_residues() {
        // Non-residues should return None
        let non_residues = [3, 5, 6, 7, 10, 11, 12, 14];
        for &n in &non_residues {
            assert!(F17::new(n).sqrt().is_none(), "{} should have no sqrt", n);
        }
    }

    #[test]
    fn sqrt_canonical() {
        // sqrt should return the smaller of the two roots
        for x in 1u64..17 {
            let a = F17::new(x);
            if let Some(r) = a.sqrt() {
                let neg_r = -r;
                assert!(
                    r.value() <= neg_r.value(),
                    "sqrt should return smaller root"
                );
            }
        }
    }

    #[test]
    fn sqrt_p_mod_4_eq_1() {
        // Test with a prime where p ≡ 1 (mod 4), which uses full Tonelli-Shanks
        // 41 ≡ 1 (mod 4)
        type F41 = Fp<41>;

        for x in 0u64..41 {
            let a = F41::new(x);
            let a_sq = a * a;
            let r = a_sq.sqrt().expect("perfect square should have sqrt");
            assert_eq!(r * r, a_sq);
        }
    }

    #[test]
    fn sqrt_p_mod_4_eq_3() {
        // Test with a prime where p ≡ 3 (mod 4), which uses fast path
        // 23 ≡ 3 (mod 4)
        type F23 = Fp<23>;

        for x in 0u64..23 {
            let a = F23::new(x);
            let a_sq = a * a;
            let r = a_sq.sqrt().expect("perfect square should have sqrt");
            assert_eq!(r * r, a_sq);
        }
    }

    #[test]
    fn discrete_log_basic() {
        // 3 is a primitive root mod 17
        let base = F17::new(3);

        // Test all elements in F17*
        for x in 0u64..16 {
            let target = base.pow(x);
            let result = target.discrete_log(base);
            assert_eq!(result, Some(x), "3^{} = {}", x, target.value());
        }
    }

    #[test]
    fn discrete_log_with_known_order() {
        // 2 has order 8 in F17 (not 16)
        let base = F17::new(2);

        for x in 0u64..8 {
            let target = base.pow(x);
            let result = target.discrete_log_with_order(base, 8);
            assert_eq!(result, Some(x));
        }
    }

    #[test]
    fn discrete_log_not_in_subgroup() {
        // 2 has order 8, so 3 is not in <2>
        let base = F17::new(2);
        let target = F17::new(3);

        // 3 is not a power of 2 mod 17
        let result = target.discrete_log_with_order(base, 8);
        assert!(result.is_none());
    }

    #[test]
    fn discrete_log_one() {
        let base = F17::new(3);
        let one = F17::ONE;
        assert_eq!(one.discrete_log(base), Some(0));
    }

    #[test]
    fn discrete_log_base() {
        let base = F17::new(3);
        assert_eq!(base.discrete_log(base), Some(1));
    }

    #[test]
    fn multiplicative_order_primitive_root() {
        // 3 is a primitive root mod 17, so order is 16
        assert_eq!(F17::new(3).multiplicative_order(), Some(16));
    }

    #[test]
    fn multiplicative_order_non_primitive() {
        // 2^8 = 256 = 15*17 + 1 ≡ 1 mod 17
        assert_eq!(F17::new(2).multiplicative_order(), Some(8));
    }

    #[test]
    fn multiplicative_order_minus_one() {
        // -1 has order 2
        assert_eq!(F17::new(16).multiplicative_order(), Some(2));
    }

    #[test]
    fn multiplicative_order_one() {
        assert_eq!(F17::ONE.multiplicative_order(), Some(1));
    }

    #[test]
    fn multiplicative_order_zero() {
        assert_eq!(F17::ZERO.multiplicative_order(), None);
    }

    #[test]
    fn is_primitive_root_true() {
        // 3 is a primitive root mod 17
        assert!(F17::new(3).is_primitive_root());
    }

    #[test]
    fn is_primitive_root_false() {
        // 2 has order 8, not primitive
        assert!(!F17::new(2).is_primitive_root());
        // 1 has order 1
        assert!(!F17::new(1).is_primitive_root());
        // -1 has order 2
        assert!(!F17::new(16).is_primitive_root());
    }

    #[test]
    fn is_primitive_root_zero() {
        assert!(!F17::ZERO.is_primitive_root());
    }

    #[test]
    fn primitive_root_exists() {
        let g = F17::primitive_root().unwrap();
        assert!(g.is_primitive_root());
        assert_eq!(g.multiplicative_order(), Some(16));
    }

    #[test]
    fn primitive_root_is_smallest() {
        // The smallest primitive root mod 17 is 3
        let g = F17::primitive_root().unwrap();
        assert_eq!(g.value(), 3);
    }

    #[test]
    fn discrete_log_larger_prime() {
        type F101 = Fp<101>;
        // Find a primitive root
        let g = F101::primitive_root().unwrap();

        // Test a few values
        for x in [0, 1, 10, 50, 99] {
            let target = g.pow(x);
            let result = target.discrete_log(g);
            assert_eq!(result, Some(x));
        }
    }

    #[test]
    fn discrete_log_non_primitive_base_collision() {
        // Regression test: BSGS with non-primitive base that causes baby-step collisions.
        // Base with small order means base^j cycles, causing hash table overwrites.
        // The fix (using or_insert) ensures we keep the smallest j.

        // 2 has order 8 in F17 (not 16), so base^8 = 1 and values repeat
        let base = F17::new(2);
        assert_eq!(base.multiplicative_order(), Some(8));

        // Test all elements in <2> = {1, 2, 4, 8, 16, 15, 13, 9}
        for x in 0u64..8 {
            let target = base.pow(x);
            let result = target.discrete_log_with_order(base, 8);
            assert_eq!(
                result,
                Some(x),
                "discrete_log of {}^{} should be {}",
                base.value(),
                x,
                x
            );
        }

        // Specifically test that we get the SMALLEST exponent
        // base^0 = 1, and with order 8, base^8 = 1 too
        // We must get 0, not 8
        let one = F17::ONE;
        assert_eq!(one.discrete_log_with_order(base, 8), Some(0));

        // base^1 = 2, base^9 = 2 (since 9 mod 8 = 1)
        // We must get 1
        assert_eq!(base.discrete_log_with_order(base, 8), Some(1));
    }

    #[test]
    fn discrete_log_large_order() {
        // Test with a larger prime to exercise ceil_sqrt_u64 with bigger values
        type F65537 = Fp<65537>; // Fermat prime F4

        let g = F65537::primitive_root().unwrap();

        // Test several exponents including large ones
        for x in [0, 1, 100, 1000, 32768, 65535] {
            let target = g.pow(x);
            let result = target.discrete_log(g);
            assert_eq!(result, Some(x), "discrete_log failed for exponent {}", x);
        }
    }

    #[test]
    fn discrete_log_very_small_order() {
        // Edge case: base = -1 has order 2
        let base = F17::new(16); // -1 mod 17
        assert_eq!(base.multiplicative_order(), Some(2));

        // Only two elements in <-1>: {1, -1}
        assert_eq!(F17::ONE.discrete_log_with_order(base, 2), Some(0));
        assert_eq!(base.discrete_log_with_order(base, 2), Some(1));
    }
}

#[cfg(all(test, feature = "serde"))]
mod serde_tests {
    use super::*;

    type F17 = Fp<17>;

    #[test]
    fn serialize_json() {
        let a = F17::new(5);
        let json = serde_json::to_string(&a).unwrap();
        assert_eq!(json, "5");
    }

    #[test]
    fn deserialize_json() {
        let a: F17 = serde_json::from_str("5").unwrap();
        assert_eq!(a.value(), 5);
    }

    #[test]
    fn roundtrip() {
        for x in 0u64..17 {
            let a = F17::new(x);
            let json = serde_json::to_string(&a).unwrap();
            let b: F17 = serde_json::from_str(&json).unwrap();
            assert_eq!(a, b);
        }
    }

    #[test]
    fn deserialize_reduces_mod_p() {
        // Values >= P should be reduced
        let a: F17 = serde_json::from_str("20").unwrap();
        assert_eq!(a.value(), 3); // 20 mod 17 = 3
    }

    #[test]
    fn serialize_vec() {
        let elements: Vec<F17> = vec![F17::new(1), F17::new(5), F17::new(16)];
        let json = serde_json::to_string(&elements).unwrap();
        assert_eq!(json, "[1,5,16]");
    }

    #[test]
    fn deserialize_vec() {
        let elements: Vec<F17> = serde_json::from_str("[1,5,16]").unwrap();
        assert_eq!(elements.len(), 3);
        assert_eq!(elements[0].value(), 1);
        assert_eq!(elements[1].value(), 5);
        assert_eq!(elements[2].value(), 16);
    }
}

#[cfg(all(test, feature = "rand"))]
mod rand_tests {
    use super::*;
    use rand::Rng;

    type F17 = Fp<17>;

    #[test]
    fn random_in_range() {
        let mut rng = rand::thread_rng();
        for _ in 0..100 {
            let a: F17 = rng.gen();
            assert!(a.value() < 17);
        }
    }

    #[test]
    fn random_distribution() {
        // Generate many random elements and check we get variety
        let mut rng = rand::thread_rng();
        let mut seen = std::collections::HashSet::new();
        for _ in 0..1000 {
            let a: F17 = rng.gen();
            seen.insert(a.value());
        }
        // With 1000 samples from 17 elements, we should see most of them
        assert!(seen.len() >= 15, "should see most field elements");
    }

    #[test]
    fn random_nonzero() {
        // Helper to generate non-zero elements
        let mut rng = rand::thread_rng();
        for _ in 0..100 {
            let mut a: F17 = rng.gen();
            while a == F17::ZERO {
                a = rng.gen();
            }
            assert_ne!(a, F17::ZERO);
        }
    }
}

#[cfg(all(test, feature = "subtle"))]
mod subtle_tests {
    use super::*;
    use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

    type F17 = Fp<17>;

    #[test]
    fn ct_eq_equal() {
        let a = F17::new(5);
        let b = F17::new(5);
        assert!(bool::from(a.ct_eq(&b)));
    }

    #[test]
    fn ct_eq_not_equal() {
        let a = F17::new(5);
        let b = F17::new(7);
        assert!(!bool::from(a.ct_eq(&b)));
    }

    #[test]
    fn conditional_select_false() {
        let a = F17::new(5);
        let b = F17::new(7);
        let result = F17::conditional_select(&a, &b, Choice::from(0));
        assert_eq!(result, a);
    }

    #[test]
    fn conditional_select_true() {
        let a = F17::new(5);
        let b = F17::new(7);
        let result = F17::conditional_select(&a, &b, Choice::from(1));
        assert_eq!(result, b);
    }

    #[test]
    fn pow_ct_basic() {
        let a = F17::new(3);
        assert_eq!(a.pow_ct(0), F17::ONE);
        assert_eq!(a.pow_ct(1), a);
        assert_eq!(a.pow_ct(2), a * a);
        assert_eq!(a.pow_ct(3), a * a * a);
    }

    #[test]
    fn pow_ct_matches_pow() {
        for base in 0u64..17 {
            let a = F17::new(base);
            for exp in 0u64..20 {
                assert_eq!(
                    a.pow_ct(exp),
                    a.pow(exp),
                    "pow_ct({}, {}) mismatch",
                    base,
                    exp
                );
            }
        }
    }

    #[test]
    fn pow_ct_fermat() {
        // a^(p-1) = 1 for all nonzero a
        for x in 1u64..17 {
            let a = F17::new(x);
            assert_eq!(a.pow_ct(16), F17::ONE);
        }
    }

    #[test]
    fn inverse_ct_basic() {
        for x in 1u64..17 {
            let a = F17::new(x);
            let inv = a.inverse_ct().unwrap();
            assert_eq!(a * inv, F17::ONE, "inverse_ct({}) failed", x);
        }
    }

    #[test]
    fn inverse_ct_zero() {
        assert!(F17::ZERO.inverse_ct().is_none());
    }

    #[test]
    fn inverse_ct_matches_inverse() {
        use crate::algebra::field::Field;
        for x in 1u64..17 {
            let a = F17::new(x);
            assert_eq!(a.inverse_ct(), a.inverse());
        }
    }
}
