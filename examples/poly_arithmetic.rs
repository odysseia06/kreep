//! Polynomial Arithmetic
//!
//! This example demonstrates `Poly<P>` operations including:
//! - Basic arithmetic (add, mul, scalar mul)
//! - Division with remainder
//! - GCD and extended GCD
//! - Derivatives and evaluation
//! - Interpolation
//!
//! Run with: cargo run --example poly_arithmetic

use kreep::{Fp, Poly, Ring};

type F17 = Fp<17>;
type P17 = Poly<17>;

fn main() {
    println!("=== Polynomial Arithmetic over F_17[x] ===\n");

    basic_operations();
    division_and_gcd();
    derivatives_and_eval();
    interpolation();
}

/// Demonstrate basic polynomial operations
fn basic_operations() {
    println!("--- Basic Operations ---\n");

    // Create polynomials: p(x) = 1 + 2x + x^2, q(x) = 1 + x
    let p = P17::new(vec![F17::new(1), F17::new(2), F17::new(1)]);
    let q = P17::new(vec![F17::new(1), F17::new(1)]);

    println!("p(x) = {:?}", p);
    println!("q(x) = {:?}", q);
    println!("degree(p) = {:?}", p.degree());
    println!("degree(q) = {:?}", q.degree());
    println!();

    // Addition
    println!("p + q = {:?}", p.clone() + q.clone());

    // Subtraction
    println!("p - q = {:?}", p.clone() - q.clone());

    // Multiplication
    let product = p.clone() * q.clone();
    println!("p * q = {:?}", product);

    // Scalar multiplication
    let scaled = p.clone() * F17::new(3);
    println!("3 * p = {:?}", scaled);

    // Special constructors
    println!("\nSpecial polynomials:");
    println!("  zero = {:?}", P17::zero());
    println!("  x = {:?}", P17::x());
    println!("  constant(5) = {:?}", P17::constant(F17::new(5)));
    println!("  monomial(3, 2) = {:?}", P17::monomial(F17::new(3), 2)); // 3x^2

    println!();
}

/// Demonstrate division and GCD
fn division_and_gcd() {
    println!("--- Division and GCD ---\n");

    // p(x) = x^3 - 1, q(x) = x - 1
    // x^3 - 1 = (x-1)(x^2 + x + 1)
    let p = P17::new(vec![F17::new(16), F17::ZERO, F17::ZERO, F17::ONE]); // -1 + x^3
    let q = P17::new(vec![F17::new(16), F17::ONE]); // -1 + x = x - 1

    println!("p(x) = {:?}  (x^3 - 1)", p);
    println!("q(x) = {:?}  (x - 1)", q);
    println!();

    // Division with remainder
    let (quotient, remainder) = p.div_rem(&q).unwrap();
    println!("p / q:");
    println!("  quotient  = {:?}", quotient);
    println!("  remainder = {:?}", remainder);
    println!();

    // Verify: p = q * quotient + remainder
    let reconstructed = q.clone() * quotient.clone() + remainder.clone();
    println!("Verify: q * quotient + remainder = {:?}", reconstructed);
    println!("Equals p? {}", reconstructed == p);
    println!();

    // GCD example
    // gcd((x-1)(x-2), (x-2)(x-3)) = (x-2)
    let poly1 = P17::from_roots(&[F17::new(1), F17::new(2)]); // (x-1)(x-2)
    let poly2 = P17::from_roots(&[F17::new(2), F17::new(3)]); // (x-2)(x-3)

    println!("poly1 = {:?}  [(x-1)(x-2)]", poly1);
    println!("poly2 = {:?}  [(x-2)(x-3)]", poly2);

    let gcd = P17::gcd(&poly1, &poly2);
    println!("gcd(poly1, poly2) = {:?}  [should be (x-2)]", gcd);
    println!();

    // Extended GCD: find s, t such that s*a + t*b = gcd(a, b)
    let a = P17::new(vec![F17::new(1), F17::new(2), F17::new(1)]); // 1 + 2x + x^2
    let b = P17::new(vec![F17::new(1), F17::new(1)]); // 1 + x

    let (g, s, t) = P17::extended_gcd(&a, &b);
    println!("Extended GCD:");
    println!("  a = {:?}", a);
    println!("  b = {:?}", b);
    println!("  gcd = {:?}", g);
    println!("  s = {:?}", s);
    println!("  t = {:?}", t);

    // Verify: s*a + t*b = g
    let check = s.clone() * a.clone() + t.clone() * b.clone();
    println!("  s*a + t*b = {:?}", check);
    println!("  Equals gcd? {}", check == g);

    println!();
}

/// Demonstrate derivatives and evaluation
fn derivatives_and_eval() {
    println!("--- Derivatives and Evaluation ---\n");

    // p(x) = 1 + 2x + 3x^2 + 4x^3
    let p = P17::new(vec![F17::new(1), F17::new(2), F17::new(3), F17::new(4)]);

    println!("p(x) = {:?}", p);
    println!();

    // Derivative: p'(x) = 2 + 6x + 12x^2
    let dp = p.derivative();
    println!("p'(x) = {:?}", dp);

    // Second derivative
    let ddp = dp.derivative();
    println!("p''(x) = {:?}", ddp);
    println!();

    // Evaluation at various points
    println!("Evaluation:");
    for x_val in [0, 1, 2, 5] {
        let x = F17::new(x_val);
        println!("  p({}) = {}", x_val, p.eval(x));
    }

    println!();
}

/// Demonstrate Lagrange interpolation
fn interpolation() {
    println!("--- Lagrange Interpolation ---\n");

    // Find polynomial passing through (0, 1), (1, 3), (2, 7)
    let points = [
        (F17::new(0), F17::new(1)),
        (F17::new(1), F17::new(3)),
        (F17::new(2), F17::new(7)),
    ];

    println!("Points: (0, 1), (1, 3), (2, 7)");

    let p = P17::interpolate(&points).unwrap();
    println!("Interpolating polynomial: {:?}", p);
    println!();

    // Verify it passes through all points
    println!("Verification:");
    for (x, y) in &points {
        let eval = p.eval(*x);
        println!("  p({}) = {} (expected {})", x, eval, y);
    }

    // Build polynomial from roots
    println!("\nPolynomial from roots:");
    let roots = [F17::new(1), F17::new(5), F17::new(10)];
    let from_roots = P17::from_roots(&roots);
    println!("  roots: 1, 5, 10");
    println!("  polynomial: {:?}", from_roots);

    // Verify roots
    for &r in &roots {
        println!("  p({}) = {}", r, from_roots.eval(r));
    }
}
