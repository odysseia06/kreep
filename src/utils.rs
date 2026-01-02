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
    if n % 2 == 0 {
        return false;
    }

    let mut i = 3;
    while i * i <= n {
        if n % i == 0 {
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
}
