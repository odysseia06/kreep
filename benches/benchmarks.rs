//! Benchmarks for kreep finite field operations.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use kreep::{ExtField, Field, Fp, Poly, Ring};

// Use the common NTT-friendly prime
type F = Fp<998244353>;

fn bench_fp_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("Fp Operations");

    let a = F::new(123456789);
    let b = F::new(987654321);

    group.bench_function("add", |bencher| {
        bencher.iter(|| black_box(a) + black_box(b))
    });

    group.bench_function("mul", |bencher| {
        bencher.iter(|| black_box(a) * black_box(b))
    });

    group.bench_function("inverse", |bencher| bencher.iter(|| black_box(a).inverse()));

    group.bench_function("pow_small", |bencher| {
        bencher.iter(|| black_box(a).pow(1000))
    });

    group.bench_function("pow_large", |bencher| {
        bencher.iter(|| black_box(a).pow(998244352))
    });

    group.finish();
}

fn bench_poly_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("Poly Operations");

    // Create polynomials of various sizes
    let coeffs_small: Vec<F> = (0..16).map(|i| F::new(i)).collect();
    let coeffs_medium: Vec<F> = (0..64).map(|i| F::new(i)).collect();
    let coeffs_large: Vec<F> = (0..256).map(|i| F::new(i)).collect();

    let p_small = Poly::new(coeffs_small.clone());
    let p_medium = Poly::new(coeffs_medium.clone());
    let p_large = Poly::new(coeffs_large.clone());

    // Polynomial multiplication (naive)
    group.bench_function("mul_naive_16x16", |bencher| {
        bencher.iter(|| black_box(p_small.clone()) * black_box(p_small.clone()))
    });

    group.bench_function("mul_naive_64x64", |bencher| {
        bencher.iter(|| black_box(p_medium.clone()) * black_box(p_medium.clone()))
    });

    // Polynomial evaluation
    let x = F::new(42);
    group.bench_function("eval_16", |bencher| {
        bencher.iter(|| black_box(&p_small).eval(black_box(x)))
    });

    group.bench_function("eval_64", |bencher| {
        bencher.iter(|| black_box(&p_medium).eval(black_box(x)))
    });

    group.bench_function("eval_256", |bencher| {
        bencher.iter(|| black_box(&p_large).eval(black_box(x)))
    });

    group.finish();
}

fn bench_ntt_operations(c: &mut Criterion) {
    use kreep::ntt::{mul_ntt, NttPlan};

    let mut group = c.benchmark_group("NTT Operations");

    // Create polynomials for NTT multiplication
    for size in [16, 64, 256, 1024] {
        let coeffs: Vec<F> = (0..size).map(|i| F::new(i as u64)).collect();
        let p = Poly::new(coeffs);

        group.bench_with_input(BenchmarkId::new("mul_ntt", size), &p, |bencher, p| {
            bencher.iter(|| mul_ntt(black_box(p), black_box(p)))
        });
    }

    // NttPlan benchmarks
    let plan = NttPlan::<998244353>::new(1024).unwrap();
    let coeffs: Vec<F> = (0..512).map(|i| F::new(i as u64)).collect();
    let p = Poly::new(coeffs);

    group.bench_function("plan_mul_512", |bencher| {
        bencher.iter(|| plan.mul(black_box(&p), black_box(&p)))
    });

    // Compare naive vs NTT for various sizes
    group.finish();

    let mut group = c.benchmark_group("Naive vs NTT");
    for size in [16, 32, 64, 128] {
        let coeffs: Vec<F> = (0..size).map(|i| F::new(i as u64)).collect();
        let p = Poly::new(coeffs);

        group.bench_with_input(BenchmarkId::new("naive", size), &p, |bencher, p| {
            bencher.iter(|| black_box(p.clone()) * black_box(p.clone()))
        });

        group.bench_with_input(BenchmarkId::new("ntt", size), &p, |bencher, p| {
            bencher.iter(|| mul_ntt(black_box(p), black_box(p)))
        });
    }
    group.finish();
}

fn bench_poly_algorithms(c: &mut Criterion) {
    let mut group = c.benchmark_group("Poly Algorithms");

    // GCD benchmark
    let p1 = Poly::new(vec![F::new(1), F::new(2), F::new(1), F::new(1)]); // x^3 + x^2 + 2x + 1
    let p2 = Poly::new(vec![F::new(1), F::new(1), F::new(1)]); // x^2 + x + 1

    group.bench_function("gcd", |bencher| {
        bencher.iter(|| Poly::gcd(black_box(&p1), black_box(&p2)))
    });

    group.bench_function("extended_gcd", |bencher| {
        bencher.iter(|| Poly::extended_gcd(black_box(&p1), black_box(&p2)))
    });

    // Irreducibility check
    // x^4 + x + 1 over F_2 is irreducible, but we use a larger prime
    let irred = Poly::new(vec![F::new(1), F::new(1), F::ZERO, F::ZERO, F::ONE]);
    group.bench_function("is_irreducible_deg4", |bencher| {
        bencher.iter(|| black_box(&irred).is_irreducible())
    });

    // Interpolation
    let points: Vec<(F, F)> = (0..16).map(|i| (F::new(i), F::new(i * i))).collect();
    group.bench_function("interpolate_16", |bencher| {
        bencher.iter(|| Poly::interpolate(black_box(&points)))
    });

    group.finish();
}

fn bench_ext_field(c: &mut Criterion) {
    let mut group = c.benchmark_group("ExtField Operations");

    type Ext2 = ExtField<998244353, 2>;

    // x^2 - 3 (3 is not a QR mod 998244353)
    let modulus = Poly::new(vec![F::new(998244353 - 3), F::ZERO, F::ONE]);

    let a = Ext2::new([F::new(123), F::new(456)]);
    let b = Ext2::new([F::new(789), F::new(101112)]);

    group.bench_function("add", |bencher| {
        bencher.iter(|| black_box(a) + black_box(b))
    });

    group.bench_function("mul_mod", |bencher| {
        bencher.iter(|| black_box(a).mul_mod(black_box(&b), black_box(&modulus)))
    });

    group.bench_function("inverse_mod", |bencher| {
        bencher.iter(|| black_box(a).inverse_mod(black_box(&modulus)))
    });

    group.bench_function("pow_mod_1000", |bencher| {
        bencher.iter(|| black_box(a).pow_mod(1000, black_box(&modulus)))
    });

    group.bench_function("frobenius", |bencher| {
        bencher.iter(|| black_box(a).frobenius(black_box(&modulus)))
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_fp_operations,
    bench_poly_operations,
    bench_ntt_operations,
    bench_poly_algorithms,
    bench_ext_field,
);
criterion_main!(benches);
