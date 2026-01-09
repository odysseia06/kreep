# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test Commands

```bash
cargo build                    # Build the library
cargo test                     # Run all tests
cargo test <test_name>         # Run a specific test
cargo test --features serde    # Test with serde feature
cargo test --features rand     # Test with rand feature
cargo test --features subtle   # Test with subtle feature (constant-time ops)
cargo test --all-features      # Test with all features enabled
cargo clippy                   # Run linter
cargo run --example f17        # Run example
```

## Architecture

This is a Rust library for abstract algebra and finite field arithmetic, organized into two main modules:

### `algebra/` - Abstract Trait Definitions
- `Ring` trait: Defines `ZERO`, `ONE`, and arithmetic operators (`+`, `-`, `*`, `neg`)
- `Field` trait: Extends `Ring` with `inverse()` for multiplicative inverses
- `Group` trait: Abstract group with `identity()`, `op()`, `inverse()`, and `pow()`

### `structures/` - Concrete Implementations

**`Fp<P>`** - Prime field elements (core type)
- Uses Montgomery representation internally for efficient modular multiplication
- `P` is a const generic `u64` prime modulus
- Supports `pow()`, `sqrt()` (Tonelli-Shanks), `legendre()`, `discrete_log()`, `multiplicative_order()`
- Batch operations: `batch_inverse()`, `batch_inverse_in_place()`
- Constant-time operations with `subtle` feature: `pow_ct()`, `inverse_ct()`

**`Poly<P>`** - Polynomials over Fp
- Coefficients in ascending order: `coeffs[i]` is coefficient of x^i
- `div_rem()`, `gcd()`, `extended_gcd()` for polynomial arithmetic
- `is_irreducible()` (Rabin's algorithm), `is_primitive()`
- `factor()` - Full factorization via square-free + distinct-degree + Cantor-Zassenhaus
- `interpolate()` - Lagrange interpolation

**`ExtField<P, D>`** - Extension field F_p[x]/(f(x))
- Fixed-size array of D coefficients from Fp
- Requires explicit modulus polynomial for `mul_mod()`, `inverse_mod()`, `pow_mod()`
- `frobenius()`, `norm()`, `trace()` operations

**`TowerField<P, D1, D2>`** - Tower extensions F_{p^D1}[y]/(g(y))
- Two-level extension: degree D2 over F_{p^D1}
- Currently only D2=2 supported for inversion

**`gf::GF<P, D>`** - Extension field with bundled modulus
- Wraps `ExtField` with `Rc<Poly>` modulus for natural arithmetic syntax
- Use `gf::irreducible_poly_deg2()`, `gf::find_irreducible()` for standard moduli

**`ntt`** - Number Theoretic Transform
- O(n log n) polynomial multiplication via `mul_ntt()`
- Requires NTT-friendly primes (e.g., 998244353 = 119 * 2^23 + 1)
- `ntt_info::<P>()` to check prime compatibility

## Optional Features

- `serde` - Serialization for `Fp` (serializes as raw u64 value)
- `rand` - Random element generation, `random_irreducible()`, `random_primitive()`
- `subtle` - Constant-time operations (`pow_ct`, `inverse_ct`, `ConstantTimeEq`, `ConditionallySelectable`)

## Key Design Patterns

- Zero polynomial represented as empty coefficient vector
- Montgomery form: `Fp::new()` converts to Montgomery, `Fp::value()` converts back
- Extension field operations require passing modulus explicitly (except `GF` wrapper)
- Polynomials auto-normalize (strip trailing zeros) on construction
