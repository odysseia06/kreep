//! NTT-Based Polynomial Multiplication
//!
//! This example demonstrates:
//! - NTT-friendly primes and their properties
//! - Fast O(n log n) polynomial multiplication using NTT
//! - Comparison with naive O(n^2) multiplication
//!
//! Run with: cargo run --example ntt_mul

use kreep::{ntt, Fp, Poly, Ring};

// NTT-friendly prime: 998244353 = 119 * 2^23 + 1
// This supports NTT sizes up to 2^23 = 8,388,608
const NTT_PRIME: u64 = 998244353;

type F = Fp<NTT_PRIME>;
type P = Poly<NTT_PRIME>;

fn main() {
    println!("=== NTT-Based Polynomial Multiplication ===\n");

    ntt_friendly_primes();
    basic_ntt_multiplication();
    compare_with_naive();
    roots_of_unity();
}

/// Explain NTT-friendly primes
fn ntt_friendly_primes() {
    println!("--- NTT-Friendly Primes ---\n");

    println!("For NTT to work, the prime p must have p - 1 divisible by a large power of 2.");
    println!("This ensures primitive 2^k-th roots of unity exist.\n");

    // Check our prime
    let info = ntt::ntt_info::<NTT_PRIME>().unwrap();
    println!("Prime: {} = 119 * 2^23 + 1", NTT_PRIME);
    println!(
        "  max_log2: {} (supports NTT up to size 2^{} = {})",
        info.max_log2,
        info.max_log2,
        1u64 << info.max_log2
    );
    println!(
        "  primitive 2^{}-th root of unity: {}",
        info.max_log2, info.primitive_root
    );
    println!();

    // Other common NTT primes
    println!("Other common NTT-friendly primes:");

    if let Some(info) = ntt::ntt_info::<2013265921>() {
        println!("  2013265921 = 15 * 2^27 + 1, max size 2^{}", info.max_log2);
    }

    if let Some(info) = ntt::ntt_info::<2281701377>() {
        println!("  2281701377 = 17 * 2^27 + 1, max size 2^{}", info.max_log2);
    }

    // Primes with limited NTT support
    println!("\nPrimes with limited NTT support:");
    if let Some(info) = ntt::ntt_info::<17>() {
        println!(
            "  17: 17-1 = 16 = 2^4, max_log2 = {} (max NTT size 16)",
            info.max_log2
        );
    }
    if let Some(info) = ntt::ntt_info::<19>() {
        println!(
            "  19: 19-1 = 18 = 2*9, max_log2 = {} (max NTT size 2)",
            info.max_log2
        );
    }
    println!();
}

/// Demonstrate basic NTT multiplication
fn basic_ntt_multiplication() {
    println!("--- Basic NTT Multiplication ---\n");

    // Create two polynomials
    // a(x) = 1 + 2x + 3x^2
    let a = P::new(vec![F::new(1), F::new(2), F::new(3)]);
    // b(x) = 4 + 5x
    let b = P::new(vec![F::new(4), F::new(5)]);

    println!("a(x) = {:?}", a);
    println!("b(x) = {:?}", b);
    println!();

    // Multiply using NTT (O(n log n))
    let c_ntt = ntt::mul_ntt(&a, &b);
    println!("a * b (NTT) = {:?}", c_ntt);

    // Expected: (1 + 2x + 3x^2)(4 + 5x) = 4 + 13x + 22x^2 + 15x^3
    println!("Expected:    4 + 13x + 22x^2 + 15x^3");
    println!();

    // Verify coefficients
    println!("Coefficient verification:");
    println!("  [0]: {} (expected 4)", c_ntt.coeff(0));
    println!("  [1]: {} (expected 13)", c_ntt.coeff(1));
    println!("  [2]: {} (expected 22)", c_ntt.coeff(2));
    println!("  [3]: {} (expected 15)", c_ntt.coeff(3));
    println!();
}

/// Compare NTT with naive multiplication
fn compare_with_naive() {
    println!("--- Comparison: NTT vs Naive ---\n");

    // Create larger polynomials
    let degree = 100;
    let a: P = Poly::new((1..=degree).map(|i| F::new(i as u64)).collect());
    let b: P = Poly::new((1..=degree).map(|i| F::new((i * 2) as u64)).collect());

    println!("Multiplying two degree-{} polynomials...", degree - 1);
    println!();

    // NTT multiplication
    let c_ntt = ntt::mul_ntt(&a, &b);

    // Naive multiplication (using Poly's * operator)
    let c_naive = a.clone() * b.clone();

    // Verify they match
    println!("Results match: {}", c_ntt == c_naive);
    println!("Result degree: {:?}", c_ntt.degree());
    println!();

    // Show a few coefficients
    println!("Sample coefficients:");
    for i in [0, 1, 50, 100, 150, 197] {
        if let Some(d) = c_ntt.degree() {
            if i <= d {
                println!("  [{}]: {}", i, c_ntt.coeff(i));
            }
        }
    }
    println!();

    println!("Complexity comparison:");
    println!("  Naive: O(n^2) = O({}) operations", degree * degree);
    println!(
        "  NTT:   O(n log n) ≈ O({}) operations",
        degree * (degree as f64).log2() as u64
    );
    println!();
}

/// Demonstrate roots of unity
fn roots_of_unity() {
    println!("--- Roots of Unity ---\n");

    // Get 8th root of unity
    let n = 8;
    let omega = ntt::get_root_of_unity::<NTT_PRIME>(n).unwrap();

    println!("8th primitive root of unity: ω = {}", omega);
    println!();

    // Show powers of omega
    println!("Powers of ω:");
    let mut power = F::ONE;
    for i in 0..=n {
        println!("  ω^{} = {}", i, power);
        power = power * omega;
    }
    println!();

    // Verify key properties
    println!("Verification:");
    println!("  ω^8 = {} (should be 1)", omega.pow(8));
    println!(
        "  ω^4 = {} (should be -1 = {})",
        omega.pow(4),
        F::new(NTT_PRIME - 1)
    );

    // ω^k ≠ 1 for 0 < k < 8
    let mut primitive = true;
    for k in 1..8 {
        if omega.pow(k) == F::ONE {
            primitive = false;
            println!("  ω^{} = 1 (not primitive!)", k);
        }
    }
    if primitive {
        println!("  ω is primitive (ω^k ≠ 1 for 0 < k < 8)");
    }
}
