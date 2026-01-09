/// Compute ceil(sqrt(n)) using integer arithmetic only.
///
/// Uses Newton's method (Heron's method) for integer square root,
/// then adjusts to get the ceiling. This avoids f64 precision loss
/// for large values of n.
pub const fn ceil_sqrt_u64(n: u64) -> u64 {
    if n <= 1 {
        return n;
    }

    // Newton's method for floor(sqrt(n))
    // Initial guess: use bit manipulation to get close to sqrt
    // For n with k bits, sqrt(n) has about k/2 bits
    let mut x = 1u64 << ((64 - n.leading_zeros() + 1) / 2);

    loop {
        let y = (x + n / x) / 2;
        if y >= x {
            break;
        }
        x = y;
    }
    // x is now floor(sqrt(n))

    // Return ceil: check if x*x == n without overflow
    // x*x == n iff n / x == x && n % x == 0
    if n / x == x && n % x == 0 {
        x
    } else {
        x + 1
    }
}

/// Compute the greatest common divisor of two numbers.
///
/// Uses the Euclidean algorithm.
pub const fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// Check if `n` is a prime number.
///
/// Uses trial division up to sqrt(n). Suitable for validating
/// moduli at startup, not for high-performance primality testing.
pub const fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 {
        return true;
    }
    if n.is_multiple_of(2) {
        return false;
    }

    let mut i = 3;
    while i <= n / i {
        if n.is_multiple_of(i) {
            return false;
        }
        i += 2;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_primes() {
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        assert!(is_prime(2));
        assert!(is_prime(3));
        assert!(!is_prime(4));
        assert!(is_prime(5));
        assert!(!is_prime(6));
        assert!(is_prime(7));
        assert!(!is_prime(8));
        assert!(!is_prime(9));
        assert!(!is_prime(10));
        assert!(is_prime(11));
        assert!(is_prime(13));
        assert!(is_prime(17));
        assert!(is_prime(19));
        assert!(is_prime(23));
    }

    #[test]
    fn composites() {
        assert!(!is_prime(15));
        assert!(!is_prime(21));
        assert!(!is_prime(25));
        assert!(!is_prime(100));
        assert!(!is_prime(1000));
    }

    #[test]
    fn larger_primes() {
        assert!(is_prime(101));
        assert!(is_prime(1009));
        assert!(is_prime(10007));
        assert!(is_prime(104729)); // 10000th prime
    }

    #[test]
    fn gcd_basic() {
        assert_eq!(gcd(12, 8), 4);
        assert_eq!(gcd(17, 5), 1);
        assert_eq!(gcd(100, 25), 25);
        assert_eq!(gcd(7, 0), 7);
        assert_eq!(gcd(0, 5), 5);
        assert_eq!(gcd(48, 18), 6);
    }

    #[test]
    fn ceil_sqrt_small() {
        assert_eq!(ceil_sqrt_u64(0), 0);
        assert_eq!(ceil_sqrt_u64(1), 1);
        assert_eq!(ceil_sqrt_u64(2), 2);
        assert_eq!(ceil_sqrt_u64(3), 2);
        assert_eq!(ceil_sqrt_u64(4), 2);
        assert_eq!(ceil_sqrt_u64(5), 3);
        assert_eq!(ceil_sqrt_u64(8), 3);
        assert_eq!(ceil_sqrt_u64(9), 3);
        assert_eq!(ceil_sqrt_u64(10), 4);
    }

    #[test]
    fn ceil_sqrt_perfect_squares() {
        assert_eq!(ceil_sqrt_u64(16), 4);
        assert_eq!(ceil_sqrt_u64(25), 5);
        assert_eq!(ceil_sqrt_u64(100), 10);
        assert_eq!(ceil_sqrt_u64(10000), 100);
        assert_eq!(ceil_sqrt_u64(1_000_000), 1000);
    }

    #[test]
    fn ceil_sqrt_non_perfect() {
        assert_eq!(ceil_sqrt_u64(17), 5);
        assert_eq!(ceil_sqrt_u64(99), 10);
        assert_eq!(ceil_sqrt_u64(101), 11);
    }

    #[test]
    fn ceil_sqrt_large() {
        // Test values near u64::MAX where f64 loses precision
        let large = (1u64 << 62) - 1;
        let result = ceil_sqrt_u64(large);
        // result should satisfy: (result-1)^2 < large <= result^2
        assert!(result > 0);
        assert!((result - 1).saturating_mul(result - 1) < large);
        // Can't check result^2 directly due to overflow, but verify it's close
        assert_eq!(result, 2147483648); // 2^31

        // Test u64::MAX
        let max_sqrt = ceil_sqrt_u64(u64::MAX);
        assert_eq!(max_sqrt, 1u64 << 32); // ceil(sqrt(2^64 - 1)) = 2^32
    }
}
