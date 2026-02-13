# kreep

`kreep` is a Rust 2021 library focused on finite-field and algebraic tooling, built with `no_std`-first design and feature-gated ergonomics. The implemented scope today is prime fields and related finite-field math (polynomials, extension/tower fields, and NTT utilities), with a staged direction toward a broader cryptography/math library.

## Current Features (Verified)

- Core algebra traits in `src/algebra/`: `Ring`, `Field`, `Group`, `AbelianGroup`.
- Prime field arithmetic via `Fp<P>` (`src/structures/fp.rs`):
  - Montgomery-backed arithmetic over `u64` prime moduli.
  - `+`, `-`, `*`, `/`, negation, exponentiation, inverse.
  - Legendre symbol, quadratic residue checks, square roots.
  - Batch inverse, discrete log, multiplicative order, primitive-root helpers.
- Polynomial arithmetic via `Poly<P>` (`src/structures/poly.rs`, requires `alloc`):
  - Construction, evaluation/interpolation, add/sub/mul.
  - Division/remainder, GCD/extended GCD.
  - Irreducibility/primitivity checks.
  - Square-free, distinct-degree, equal-degree, and full factorization (randomized parts behind `rand`).
- Extension field support (`src/structures/ext.rs`):
  - `ExtField<P, D>` arithmetic, modular inverse/pow.
  - Frobenius, norm, trace, conjugates, minimal polynomial.
  - `TowerField<P, D1, D2>` support.
- Validated wrappers (`src/structures/gf.rs`, requires `alloc`):
  - `Modulus`, `GF`, `TowerModulus`, `TowerGF`.
  - Irreducible/primitive polynomial finders.
  - `GFWithModulus` for serde round-trips (`serde` feature).
- NTT utilities (`src/structures/ntt.rs`, requires `alloc`):
  - Root-of-unity helpers, forward/inverse NTT, `mul_ntt`, `NttPlan`.
- Numeric helpers (`src/utils.rs`): `is_prime`, `gcd`, `ceil_sqrt_u64`.

## Non-Goals / Safety Notes

- This crate does **not** claim formal security proofs or production hardening across the full API surface.
- Constant-time behavior is **not guaranteed globally**. The `subtle` feature provides constant-time helpers for `Fp` (for example `ct_eq`, `pow_ct`, `inverse_ct`), but other operations are not presented as constant-time audited.
- `Fp<P>` semantics require `P` to be an odd prime. Use `Fp::<P>::validate_prime()` for explicit runtime validation.
- Some operator overloads can panic on invalid division (for example division by zero).
- NTT APIs require NTT-friendly prime/size constraints; some invalid inputs panic.
- `TowerField`/`TowerGF` inverse support is currently limited for tower degree `D2 = 2`.

## Quick Start

Add the crate and run a small prime-field program:

```rust
use kreep::{Field, Fp};

type F17 = Fp<17>;

fn main() {
    let a = F17::new(5);
    let b = F17::new(9);

    let sum = a + b;
    let product = a * b;
    let inv_a = a.inverse().unwrap();

    println!("sum={sum}, product={product}, inv_a={inv_a}");
}
```

## Examples

Only examples that currently exist in `examples/`:

- `cargo run --example field_basics`
- `cargo run --example poly_arithmetic`
- `cargo run --example irreducible_polys`
- `cargo run --example extension_fields`
- `cargo run --example ntt_mul`
- `cargo run --example factorization --features rand`

## Testing and Feature Flags

### Common Commands

- `cargo build`
- `cargo test`
- `cargo test --all-features`
- `cargo test --features serde --test serde_tests`
- `cargo clippy --all-targets --all-features`
- `cargo fmt --check`
- `cargo bench`

Note: there are currently no CI workflows under `.github/workflows/`; validation is run locally via the commands above.

### Feature Flags (`Cargo.toml`)

- `default = ["std"]`
- `std = ["alloc"]`
- `alloc = []`
- `serde = ["dep:serde", "alloc"]`
- `rand = ["dep:rand"]`
- `subtle = ["dep:subtle"]`

Practical gating:

- `Poly`, `gf`, and `ntt` modules require `alloc`.
- `serde` support requires `alloc`.
- Randomized factorization/root-finding helpers require `rand`.
- Constant-time `Fp` helper APIs are behind `subtle`.

## Roadmap (Planned)

The current implemented scope is finite fields and supporting algebra. Planned expansion is staged and test-backed:

- **Planned:** larger-prime backends using multi-limb/big-integer representations.
- **Planned:** stronger `no_std` coverage and tighter alloc/no-alloc boundaries.
- **Planned:** broader constant-time hardening and audit coverage beyond current `Fp` helpers.
- **Planned:** more ergonomic extension-field and tower-field constructors/builders.
- **Planned:** ECC layer and group abstractions on top of field primitives.
- **Planned:** expanded benchmark coverage and performance optimizations.

## License

MIT License (see `LICENSE`).

## Maintainer Checklist

When updating this README, verify each claim against code before merging:

- [ ] Feature list matches `src/` modules and public exports in `src/lib.rs`.
- [ ] Example commands match files in `examples/` and required feature flags.
- [ ] Test/lint/build commands match current `Cargo.toml` and repository layout.
- [ ] Feature-flag section matches `[features]` in `Cargo.toml` exactly.
- [ ] Safety notes still match current panic behavior and implementation limits.
- [ ] Roadmap items are clearly labeled as **Planned** and not implied as shipped.
- [ ] License section matches `LICENSE`.
