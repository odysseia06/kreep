//! Extension Fields
//!
//! This example demonstrates:
//! - Low-level `ExtField` with explicit modulus passing
//! - High-level `GF` wrapper with natural arithmetic syntax
//! - Field operations: multiplication, inversion, exponentiation
//! - Frobenius endomorphism, norm, and trace
//!
//! Run with: cargo run --example extension_fields

use kreep::{gf, ExtField, Fp, Poly, Ring};

type F17 = Fp<17>;

fn main() {
    println!("=== Extension Fields ===\n");

    extfield_basics();
    gf_wrapper();
    frobenius_norm_trace();
}

/// Demonstrate low-level ExtField with explicit modulus
fn extfield_basics() {
    println!("--- ExtField (Low-Level API) ---\n");

    // Create the irreducible modulus: x^2 - 3 over F_17
    // Since 3 is not a QR mod 17, this is irreducible
    let modulus = Poly::new(vec![F17::new(14), F17::ZERO, F17::ONE]); // -3 + x^2
    println!("Modulus: {:?}  (x^2 - 3)", modulus);
    println!("This defines F_17[x]/(x^2 - 3) ≅ F_289\n");

    // Create elements: coeffs[0] + coeffs[1]*x
    type Ext2 = ExtField<17, 2>;

    let a = Ext2::new([F17::new(2), F17::new(3)]); // 2 + 3x
    let b = Ext2::new([F17::new(1), F17::new(1)]); // 1 + x

    println!("a = {:?}", a);
    println!("b = {:?}", b);
    println!();

    // Addition (no modulus needed)
    println!("a + b = {:?}", a + b);
    println!("a - b = {:?}", a - b);
    println!();

    // Multiplication requires the modulus
    let ab = a.mul_mod(&b, &modulus);
    println!("a * b = {:?}", ab);

    // Let's verify: (2 + 3x)(1 + x) = 2 + 5x + 3x^2
    // Since x^2 = 3: = 2 + 5x + 9 = 11 + 5x
    println!("  (verify: (2+3x)(1+x) = 2 + 5x + 3x^2 = 2 + 5x + 9 = 11 + 5x)");
    println!();

    // Inversion
    let a_inv = a.inverse_mod(&modulus).unwrap();
    println!("a^(-1) = {:?}", a_inv);

    // Verify: a * a^(-1) = 1
    let product = a.mul_mod(&a_inv, &modulus);
    println!("a * a^(-1) = {:?} (should be 1)", product);
    println!();

    // Exponentiation
    let a_cubed = a.pow_mod(3, &modulus);
    println!("a^3 = {:?}", a_cubed);

    // Verify manually
    let a_sq = a.mul_mod(&a, &modulus);
    let a_cubed_manual = a_sq.mul_mod(&a, &modulus);
    println!("a * a * a = {:?}", a_cubed_manual);
    println!();
}

/// Demonstrate high-level GF wrapper with natural syntax
fn gf_wrapper() {
    println!("--- GF Wrapper (High-Level API) ---\n");

    // Get standard irreducible polynomial for F_17^2
    let modulus = gf::irreducible_poly_deg2::<17>();
    println!("Modulus: {:?}", modulus);
    println!();

    // Create GF elements - modulus is bundled with the element
    type GF17_2 = gf::GF<17, 2>;

    let a = GF17_2::new([F17::new(2), F17::new(3)], modulus.clone());
    let b = GF17_2::new([F17::new(1), F17::new(1)], modulus.clone());

    println!("a = {}", a);
    println!("b = {}", b);
    println!();

    // Natural arithmetic syntax - no need to pass modulus!
    println!("a + b = {}", &a + &b);
    println!("a - b = {}", &a - &b);
    println!("a * b = {}", &a * &b);
    println!("a / b = {}", &a / &b);
    println!();

    // Verify division: (a / b) * b = a
    let quotient = &a / &b;
    let product = &quotient * &b;
    println!("Verify: (a / b) * b = {}", product);
    println!("Equals a? {}", product == a);
    println!();

    // Exponentiation
    println!("a^5 = {}", a.pow(5));
    println!();

    // Fermat's little theorem in extension field: a^(q-1) = 1 where q = 17^2 = 289
    let order = 17u64 * 17 - 1; // 288
    let a_to_order = a.pow(order);
    println!("a^288 = {} (should be 1, by Fermat)", a_to_order);
    println!();
}

/// Demonstrate Frobenius endomorphism, norm, and trace
fn frobenius_norm_trace() {
    println!("--- Frobenius, Norm, and Trace ---\n");

    let modulus = gf::irreducible_poly_deg2::<17>();
    type GF17_2 = gf::GF<17, 2>;

    // Base field element
    let base = GF17_2::from_base(F17::new(5), modulus.clone());
    println!("Base field element: c = {}", base);
    println!(
        "  Frobenius(c) = {} (identity for base field)",
        base.frobenius()
    );
    println!("  Norm(c) = {} (= c^2 = 25 = 8 mod 17)", base.norm());
    println!("  Trace(c) = {} (= 2c = 10)", base.trace());
    println!();

    // Extension field element
    let a = GF17_2::new([F17::new(2), F17::new(3)], modulus.clone());
    println!("Extension element: a = {}", a);

    // Frobenius: x -> x^p
    let a_frob = a.frobenius();
    println!("  Frobenius(a) = a^17 = {}", a_frob);

    // For degree-2 extension, Frobenius^2 = identity
    let a_frob2 = a_frob.frobenius();
    println!("  Frobenius^2(a) = {} (should equal a)", a_frob2);
    println!("  Equals a? {}", a_frob2 == a);
    println!();

    // Norm: product of conjugates = a * Frobenius(a)
    let norm = a.norm();
    println!("  Norm(a) = {} (in base field F_17)", norm);

    // Verify: norm = a * a^p, and the result should be in the base field
    let norm_check = &a * &a_frob;
    let norm_as_gf = GF17_2::from_base(norm, modulus.clone());
    println!("  Verify: a * Frobenius(a) = {}", norm_check);
    println!("  Matches Norm(a)? {}", norm_check == norm_as_gf);
    println!();

    // Trace: sum of conjugates = a + Frobenius(a)
    let trace = a.trace();
    println!("  Trace(a) = {} (in base field F_17)", trace);

    // Verify: trace = a + a^p, and the result should be in the base field
    let trace_check = &a + &a_frob;
    let trace_as_gf = GF17_2::from_base(trace, modulus.clone());
    println!("  Verify: a + Frobenius(a) = {}", trace_check);
    println!("  Matches Trace(a)? {}", trace_check == trace_as_gf);
}
