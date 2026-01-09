//! Polynomial Factorization
//!
//! This example demonstrates:
//! - Square-free factorization
//! - Distinct-degree factorization
//! - Full factorization into irreducible factors
//! - Root finding
//!
//! Run with: cargo run --example factorization --features rand

#[cfg(not(feature = "rand"))]
fn main() {
    eprintln!("This example requires the `rand` feature.");
    eprintln!("Run with: cargo run --example factorization --features rand");
    std::process::exit(1);
}

#[cfg(feature = "rand")]
fn main() {
    use kreep::{Fp, Poly, Ring};

    type F17 = Fp<17>;
    type P17 = Poly<17>;

    println!("=== Polynomial Factorization over F_17 ===\n");

    // --- Square-Free Factorization ---
    println!("--- Square-Free Factorization ---\n");

    println!("Square-free factorization separates repeated factors.");
    println!("f = f1 * f2^2 * f3^3 * ... where each fi is square-free.\n");

    // Create f(x) = (x - 1)^2 * (x - 2) * (x - 3)^3
    let x_minus_1 = P17::new(vec![F17::new(16), F17::ONE]); // x - 1
    let x_minus_2 = P17::new(vec![F17::new(15), F17::ONE]); // x - 2
    let x_minus_3 = P17::new(vec![F17::new(14), F17::ONE]); // x - 3

    let f = x_minus_1.clone()
        * x_minus_1.clone()
        * x_minus_2.clone()
        * x_minus_3.clone()
        * x_minus_3.clone()
        * x_minus_3.clone();

    println!("f(x) = (x-1)^2 * (x-2) * (x-3)^3");
    println!("f(x) = {:?}", f);
    println!();

    let sqf = f.square_free_factorization();
    println!("Square-free factorization:");
    for (factor, multiplicity) in &sqf {
        println!("  {:?} with multiplicity {}", factor, multiplicity);
    }
    println!();

    // --- Distinct-Degree Factorization ---
    println!("--- Distinct-Degree Factorization ---\n");

    println!("Groups factors by their degree.");
    println!("Returns products of all irreducible factors of each degree.\n");

    // Create a polynomial with factors of different degrees
    // Linear factors: (x - 1)(x - 2)
    // Quadratic factor: x^2 - 3 (irreducible since 3 is non-QR mod 17)
    let linear1 = P17::new(vec![F17::new(16), F17::ONE]); // x - 1
    let linear2 = P17::new(vec![F17::new(15), F17::ONE]); // x - 2
    let quadratic = P17::new(vec![F17::new(14), F17::ZERO, F17::ONE]); // x^2 - 3

    let g = linear1 * linear2 * quadratic;

    println!("g(x) = (x-1)(x-2)(x^2-3)");
    println!("g(x) = {:?}", g);
    println!();

    let ddf = g.distinct_degree_factorization();
    println!("Distinct-degree factorization:");
    for (product, degree) in &ddf {
        println!("  Degree {}: {:?}", degree, product);
    }
    println!();

    // --- Full Factorization ---
    println!("--- Full Factorization ---\n");

    let mut rng = rand::thread_rng();

    // Create a polynomial with known factors
    // h(x) = (x - 1)^2 * (x - 5) * (x^2 - 3)
    let x_minus_1 = P17::new(vec![F17::new(16), F17::ONE]);
    let x_minus_5 = P17::new(vec![F17::new(12), F17::ONE]);
    let x2_minus_3 = P17::new(vec![F17::new(14), F17::ZERO, F17::ONE]);

    let h = x_minus_1.clone() * x_minus_1 * x_minus_5 * x2_minus_3;

    println!("h(x) = (x-1)^2 * (x-5) * (x^2-3)");
    println!("h(x) = {:?}", h);
    println!();

    let factors = h.factor(&mut rng);
    println!("Complete factorization:");
    for (factor, multiplicity) in &factors {
        println!("  {:?}  (multiplicity {})", factor, multiplicity);
    }
    println!();

    // Verify by multiplying back
    let mut reconstructed = P17::constant(F17::ONE);
    for (factor, mult) in &factors {
        for _ in 0..*mult {
            reconstructed = reconstructed * factor.clone();
        }
    }

    // Make both monic for comparison
    let h_monic = h.monic().unwrap();
    let r_monic = reconstructed.monic().unwrap();
    println!("Reconstruction matches: {}", h_monic == r_monic);
    println!();

    // --- Root Finding ---
    println!("--- Root Finding ---\n");

    // Create polynomial with known roots
    // p(x) = (x - 3)^2 * (x - 7) * (x - 11)
    let roots_expected: [(u64, u32); 3] = [(3, 2), (7, 1), (11, 1)];

    let x_minus_3 = P17::new(vec![F17::new(14), F17::ONE]); // x - 3 = x + 14
    let x_minus_7 = P17::new(vec![F17::new(10), F17::ONE]); // x - 7 = x + 10
    let x_minus_11 = P17::new(vec![F17::new(6), F17::ONE]); // x - 11 = x + 6

    let p = x_minus_3.clone() * x_minus_3 * x_minus_7 * x_minus_11;

    println!("p(x) = (x-3)^2 * (x-7) * (x-11)");
    println!("p(x) = {:?}", p);
    println!();

    // Verify roots evaluate to zero
    println!("Verification that these are roots:");
    for (r, _) in &roots_expected {
        let val = p.eval(F17::new(*r));
        println!("  p({}) = {}", r, val);
    }
    println!();

    // Find roots using factorization
    let roots = p.roots(&mut rng);
    println!("Found roots:");
    for (root, multiplicity) in &roots {
        println!("  {} with multiplicity {}", root, multiplicity);
    }
    println!();

    // Polynomial with no roots in F_17
    let no_roots = P17::new(vec![F17::new(14), F17::ZERO, F17::ONE]); // x^2 - 3
    println!("q(x) = x^2 - 3 (irreducible, no roots in F_17)");
    let roots = no_roots.roots(&mut rng);
    println!("Roots of q: {:?}", roots);
}
