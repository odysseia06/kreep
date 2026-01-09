//! Irreducible and Primitive Polynomials
//!
//! This example demonstrates:
//! - Testing polynomials for irreducibility
//! - Testing polynomials for primitivity
//! - Finding irreducible polynomials of a given degree
//! - Standard irreducible polynomial generators
//!
//! Run with: cargo run --example irreducible_polys

use kreep::{gf, Fp, Poly, Ring};

type F17 = Fp<17>;
type P17 = Poly<17>;

fn main() {
    println!("=== Irreducible and Primitive Polynomials ===\n");

    irreducibility_tests();
    primitivity_tests();
    finding_irreducibles();
    standard_generators();
}

/// Test specific polynomials for irreducibility
fn irreducibility_tests() {
    println!("--- Irreducibility Tests ---\n");

    // x^2 - 3 over F_17
    // 3 is NOT a quadratic residue mod 17, so x^2 - 3 is irreducible
    let irred = P17::new(vec![F17::new(14), F17::ZERO, F17::ONE]); // -3 + x^2
    println!("f(x) = x^2 - 3 = {:?}", irred);
    println!("  is_irreducible: {}", irred.is_irreducible());
    println!("  (3 is not a quadratic residue mod 17)\n");

    // x^2 - 4 over F_17
    // 4 = 2^2 is a quadratic residue, so x^2 - 4 = (x-2)(x+2) is reducible
    let reducible = P17::new(vec![F17::new(13), F17::ZERO, F17::ONE]); // -4 + x^2
    println!("g(x) = x^2 - 4 = {:?}", reducible);
    println!("  is_irreducible: {}", reducible.is_irreducible());
    println!("  (factors as (x-2)(x+2))\n");

    // Linear polynomials are always irreducible
    let linear = P17::new(vec![F17::new(5), F17::ONE]); // x + 5
    println!("h(x) = x + 5 = {:?}", linear);
    println!(
        "  is_irreducible: {} (linear polynomials are always irreducible)\n",
        linear.is_irreducible()
    );

    // A degree 3 example: x^3 + x + 1
    let cubic = P17::new(vec![F17::ONE, F17::ONE, F17::ZERO, F17::ONE]); // 1 + x + x^3
    println!("p(x) = x^3 + x + 1 = {:?}", cubic);
    println!("  is_irreducible: {}\n", cubic.is_irreducible());
}

/// Test polynomials for primitivity
fn primitivity_tests() {
    println!("--- Primitivity Tests ---\n");

    println!("A primitive polynomial is irreducible AND its root generates");
    println!("the multiplicative group of the extension field.\n");

    // Over F_3
    type F3 = Fp<3>;
    type P3 = Poly<3>;

    // x^2 + 1 over F_3: irreducible but NOT primitive
    // The root has order 4, not 8 = 3^2 - 1
    let irred_not_prim = P3::new(vec![F3::new(1), F3::ZERO, F3::new(1)]); // 1 + x^2
    println!("Over F_3:");
    println!("  f(x) = x^2 + 1 = {:?}", irred_not_prim);
    println!("  is_irreducible: {}", irred_not_prim.is_irreducible());
    println!("  is_primitive: {}", irred_not_prim.is_primitive());
    println!("  (root has order 4, not 8 = 3^2 - 1)\n");

    // x^2 + 2x + 2 over F_3: primitive
    let primitive = P3::new(vec![F3::new(2), F3::new(2), F3::new(1)]); // 2 + 2x + x^2
    println!("  g(x) = x^2 + 2x + 2 = {:?}", primitive);
    println!("  is_irreducible: {}", primitive.is_irreducible());
    println!("  is_primitive: {}", primitive.is_primitive());
    println!("  (root has order 8 = 3^2 - 1)\n");
}

/// Find irreducible polynomials of various degrees
fn finding_irreducibles() {
    println!("--- Finding Irreducible Polynomials ---\n");

    // Find irreducible polynomials over F_5 for various degrees
    println!("Irreducible polynomials over F_5:");
    for degree in 2..=5 {
        let irred = gf::find_irreducible::<5>(degree);
        println!("  degree {}: {:?}", degree, irred);
        assert!(irred.is_irreducible());
    }
    println!();

    // Find irreducible polynomials over F_7
    println!("Irreducible polynomials over F_7:");
    for degree in 2..=4 {
        let irred = gf::find_irreducible::<7>(degree);
        println!("  degree {}: {:?}", degree, irred);
        assert!(irred.is_irreducible());
    }
    println!();
}

/// Demonstrate standard irreducible polynomial generators
fn standard_generators() {
    println!("--- Standard Generators ---\n");

    // Degree 2 irreducible polynomials (x^2 - a where a is non-QR)
    println!("Quadratic irreducibles (x^2 - a for smallest non-QR a):");

    let mod17 = gf::irreducible_poly_deg2::<17>();
    println!("  F_17: {:?}", mod17);
    println!("        (3 is smallest non-QR mod 17)");

    let mod7 = gf::irreducible_poly_deg2::<7>();
    println!("  F_7:  {:?}", mod7);
    println!("        (3 is smallest non-QR mod 7)");

    let mod5 = gf::irreducible_poly_deg2::<5>();
    println!("  F_5:  {:?}", mod5);
    println!("        (2 is smallest non-QR mod 5)");
    println!();

    // Degree 3 irreducible polynomials
    println!("Cubic irreducibles:");
    let cubic7 = gf::irreducible_poly_deg3::<7>();
    println!("  F_7: {:?}", cubic7);

    let cubic5 = gf::irreducible_poly_deg3::<5>();
    println!("  F_5: {:?}", cubic5);
    println!();

    // Binary field polynomials (for reference)
    println!("Binary field polynomials (for reference):");
    for n in [4, 8, 16] {
        if let Some(poly) = gf::binary_field_poly(n) {
            println!("  GF(2^{}): 0x{:x} = 0b{:b}", n, poly, poly);
        }
    }
    println!("  (Note: GF(2^8) with 0x11b is used in AES)");
}
