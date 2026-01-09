//! Finite Field Basics
//!
//! This example demonstrates core `Fp<P>` operations including:
//! - Basic arithmetic (add, mul, pow, inverse)
//! - Quadratic residues and square roots
//! - Discrete logarithm and multiplicative order
//!
//! Run with: cargo run --example field_basics

use kreep::{Field, Fp, Ring};

// Define a prime field with modulus 17
type F17 = Fp<17>;

fn main() {
    println!("=== Finite Field Basics over F_17 ===\n");

    basic_arithmetic();
    quadratic_residues();
    discrete_logarithm();
    field_identities();
}

/// Demonstrate basic field arithmetic
fn basic_arithmetic() {
    println!("--- Basic Arithmetic ---\n");

    let a = F17::new(5);
    let b = F17::new(9);

    println!("a = {}", a);
    println!("b = {}", b);
    println!();

    // Addition and subtraction
    println!("a + b = {}", a + b);
    println!("a - b = {}", a - b);

    // Multiplication
    println!("a * b = {}", a * b);

    // Exponentiation
    println!("a^3 = {}", a.pow(3));

    // Multiplicative inverse
    let a_inv = a.inverse().unwrap();
    println!("a^(-1) = {}", a_inv);
    println!("a * a^(-1) = {} (should be 1)", a * a_inv);

    // Division
    println!("a / b = {}", a / b);

    // Fermat's little theorem: a^(p-1) = 1
    println!("a^16 = {} (Fermat's little theorem)", a.pow(16));

    println!();
}

/// Demonstrate quadratic residues and square roots
fn quadratic_residues() {
    println!("--- Quadratic Residues and Square Roots ---\n");

    // Find all quadratic residues mod 17
    print!("Quadratic residues mod 17: ");
    for x in 0..17 {
        let elem = F17::new(x);
        if elem.is_quadratic_residue() {
            print!("{} ", x);
        }
    }
    println!("\n");

    // Compute Legendre symbols
    println!("Legendre symbols (a/17):");
    for a in 1..=5 {
        let elem = F17::new(a);
        let legendre = elem.legendre();
        let status = match legendre {
            1 => "quadratic residue",
            -1 => "non-residue",
            _ => "zero",
        };
        println!("  ({}/17) = {} ({})", a, legendre, status);
    }
    println!();

    // Compute square roots
    println!("Square roots:");
    let residue = F17::new(2); // 2 is a QR mod 17
    if let Some(sqrt) = residue.sqrt() {
        println!(
            "  sqrt(2) = {} (verify: {}^2 = {})",
            sqrt,
            sqrt,
            sqrt * sqrt
        );
    }

    let non_residue = F17::new(3); // 3 is not a QR mod 17
    match non_residue.sqrt() {
        Some(r) => println!("  sqrt(3) = {}", r),
        None => println!("  sqrt(3) = None (3 is not a quadratic residue)"),
    }

    println!();
}

/// Demonstrate discrete logarithm and multiplicative order
fn discrete_logarithm() {
    println!("--- Discrete Logarithm ---\n");

    // 3 is a primitive root mod 17 (generates all of F_17*)
    let g = F17::new(3);
    let order = g.multiplicative_order().unwrap();
    println!("g = {} has multiplicative order {}", g, order);
    println!("(3 is a primitive root since order = p-1 = 16)\n");

    // Compute discrete logs for various elements
    println!("Discrete logarithms base 3:");
    for target_val in [1, 2, 5, 10, 16] {
        let target = F17::new(target_val);
        if let Some(x) = target.discrete_log(g) {
            println!(
                "  log_3({}) = {} (verify: 3^{} = {})",
                target_val,
                x,
                x,
                g.pow(x)
            );
        }
    }
    println!();

    // Demonstrate with a non-primitive base
    let base = F17::new(2); // 2 has order 8 in F_17
    let base_order = base.multiplicative_order().unwrap();
    println!("base = {} has order {} (not primitive)", base, base_order);

    let target = base.pow(5);
    let x = target.discrete_log_with_order(base, base_order).unwrap();
    println!(
        "  log_2({}) = {} in subgroup of order {}",
        target, x, base_order
    );

    println!();
}

/// Verify basic field identities
fn field_identities() {
    println!("--- Field Identities Check ---\n");

    let mut all_passed = true;

    // Check for all non-zero elements
    for x in 1..17 {
        let a = F17::new(x);

        // a * a^(-1) = 1
        let a_inv = a.inverse().unwrap();
        if a * a_inv != F17::ONE {
            println!("FAIL: {} * {}^(-1) != 1", a, a);
            all_passed = false;
        }

        // Distributivity: a * (b + c) = a*b + a*c
        for y in 0..17 {
            for z in 0..17 {
                let b = F17::new(y);
                let c = F17::new(z);
                if a * (b + c) != a * b + a * c {
                    println!("FAIL: distributivity for a={}, b={}, c={}", a, b, c);
                    all_passed = false;
                }
            }
        }
    }

    if all_passed {
        println!("All field identities verified for F_17!");
    }
}
