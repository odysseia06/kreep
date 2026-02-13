# AGENTS Guide

## What this repo is
`kreep` is a `no_std`-first Rust cryptography-math library centered on finite-field arithmetic over prime fields (`Fp<P>`) and higher-level finite-field tooling (polynomials, extension fields, tower fields, and NTT-based multiplication). Today the implemented scope is finite fields and related algebra/number-theory utilities; the intended direction is to grow from this base into broader crypto utilities in clearly staged, test-backed increments.

## Current capabilities (verified from code/tests/examples)
- Core algebra traits in `src/algebra/`: `Ring`, `Field`, `Group`, and `AbelianGroup`.
- Prime field `Fp<P>` in `src/structures/fp.rs` (Montgomery representation over `u64` primes):
  `+`, `-`, `*`, `/`, negation, exponentiation, inverse, Legendre symbol, quadratic residue/sqrt, batch inverse, discrete log, multiplicative order, primitive-root helpers.
- Polynomial arithmetic `Poly<P>` in `src/structures/poly.rs` (`alloc` required):
  constructors, eval/interpolation, add/sub/mul, Euclidean div/rem, GCD/extended GCD, irreducibility and primitivity checks, derivative, square-free/distinct-degree/equal-degree/full factorization, root extraction (`rand`-gated pieces where applicable).
- Extension fields in `src/structures/ext.rs`:
  `ExtField<P, D>` arithmetic, modular multiplication/inversion/pow, Frobenius/norm/trace/conjugates/minimal polynomial (`alloc`-gated methods where needed), plus `TowerField<P, D1, D2>` support.
- Validated field wrappers and constructors in `src/structures/gf.rs` (`alloc`):
  `Modulus`, `GF`, `TowerModulus`, `TowerGF`, irreducible/primitive polynomial finders, and `GFWithModulus` (`serde`) serialization helper.
- NTT support in `src/structures/ntt.rs` (`alloc`):
  roots-of-unity helpers, forward/inverse NTT, `mul_ntt`, and reusable `NttPlan`.
- Utilities in `src/utils.rs`: `is_prime`, `gcd`, `ceil_sqrt_u64`.
- Tests/examples present and active in:
  `tests/field_properties.rs`, `tests/serde_tests.rs`, module unit tests, and examples under `examples/`.

## Workspace and module layout
- `Cargo.toml`: crate metadata, features, deps, bench target.
- `src/lib.rs`: crate entrypoint, feature-gated exports, `no_std` setup.
- `src/algebra/`: abstract algebra traits.
- `src/structures/fp.rs`: prime-field implementation.
- `src/structures/poly.rs`: polynomial algorithms over `Fp`.
- `src/structures/ext.rs`: extension and tower field element types.
- `src/structures/gf.rs`: validated modulus types, ergonomic wrappers, constructors.
- `src/structures/ntt.rs`: NTT primitives and planning APIs.
- `src/utils.rs`: numeric helpers.
- `tests/`: integration tests (`field_properties`, `serde_tests`).
- `examples/`: runnable demos:
  `field_basics`, `poly_arithmetic`, `irreducible_polys`, `extension_fields`, `ntt_mul`, `factorization`.
- `benches/benchmarks.rs`: Criterion benchmarks.
- `.github/ISSUE_TEMPLATE/`: issue templates only (no CI workflow files currently).

## Build, test, examples (exact commands)
- Build:
  `cargo build`
- Full test suite:
  `cargo test`
- Targeted test name:
  `cargo test <test_name>`
- Feature-specific tests:
  `cargo test --features serde`
  `cargo test --features rand`
  `cargo test --features subtle`
- Serde integration test explicitly:
  `cargo test --features serde --test serde_tests`
- All optional features together:
  `cargo test --all-features`
- Lint:
  `cargo clippy --all-targets --all-features`
- Format check:
  `cargo fmt --check`
- Run examples:
  `cargo run --example field_basics`
  `cargo run --example poly_arithmetic`
  `cargo run --example irreducible_polys`
  `cargo run --example extension_fields`
  `cargo run --example ntt_mul`
  `cargo run --example factorization --features rand`
- Benchmarks:
  `cargo bench`

## Feature flags (from `Cargo.toml`)
- `default = ["std"]`
- `std = ["alloc"]`
- `alloc = []`
- `serde = ["dep:serde", "alloc"]`
- `rand = ["dep:rand"]`
- `subtle = ["dep:subtle"]`

Practical impact:
- `poly`, `gf`, and `ntt` modules are `alloc`-gated.
- `serde` support needs `alloc` (and integration serde tests are `std`-gated).
- `rand` gates random sampling/factorization helpers.
- `subtle` gates constant-time equality/select and constant-time `Fp` helpers.

## Coding conventions
- Rust 2021 + standard `rustfmt` style.
- Naming conventions in code:
  `snake_case` (functions/modules), `UpperCamelCase` (types), `UPPER_SNAKE_CASE` (constants).
- Keep trait APIs in `src/algebra/*` and concrete arithmetic in `src/structures/*`.
- Prefer small arithmetic helpers and explicit invariants.
- Run `cargo fmt --check` and `cargo clippy --all-targets --all-features` before landing changes.

## Error handling style
- Fallible arithmetic APIs generally return `Option<T>` (for example `inverse`, `try_div`, many modular operations).
- Validation APIs use typed `Result` errors (for example `Modulus::new(...) -> Result<_, ModulusError>`).
- Operator overloads (`Div` on `Fp`, `GF`, `TowerGF`) may panic on invalid division (for example division by zero); use explicit fallible APIs when you need non-panicking behavior.
- Some invariants are debug-only checks (for example `Fp::new` uses `debug_assert!` for prime/odd modulus checks); explicit validation helpers exist (`Fp::validate_prime`).

## How to add a new module
1. Pick the right layer:
   `src/algebra/` for traits, `src/structures/` for concrete implementations.
2. Add the file and wire module exports:
   update `src/algebra/mod.rs` or `src/structures/mod.rs`.
3. Re-export in `src/lib.rs` (feature-gated where needed).
4. Feature-gate correctly:
   `alloc` for heap-backed types, `serde`/`rand`/`subtle` for optional behavior.
5. Follow existing error-handling patterns:
   `Option`/`Result` for fallible operations, panic only where operator semantics already do so.
6. Add tests:
   unit tests in-module; integration tests in `tests/` when API-level behavior matters.
7. Add or update an example under `examples/` for user-facing behavior.
8. Verify with:
   `cargo fmt --check`
   `cargo clippy --all-targets --all-features`
   `cargo test --all-features`

## Roadmap (Planned, not yet implemented)
- Larger-prime backends using limb-based/big-integer representations (beyond current `u64` modulus model).
- Stronger `no_std` coverage and tighter alloc/no-alloc surface guarantees across modules.
- Broader constant-time hardening/audit across non-`Fp` operations.
- Elliptic-curve structures and group operations built on top of current field primitives.
- More ergonomic constructors/builders for extension/tower fields and modulus selection.
- Continued evolution from finite-field core into broader reusable cryptographic utilities.

## Known pitfalls / gotchas
- `Fp<P>` requires `P` to be an odd prime for correct field behavior (Montgomery math); `P = 2` is rejected by `validate_prime`.
- `Fp::new` checks prime/oddness with `debug_assert!`; call `Fp::<P>::validate_prime()` at startup for explicit runtime validation.
- Division operators panic on invalid division:
  `Fp`/`GF` divide-by-zero, and `TowerGF` division can also fail when inverse is unsupported (`D2 != 2` currently).
- `alloc` feature is required for `Poly`, `gf`, and `ntt` modules.
- `factorization`/randomized polynomial routines need the `rand` feature.
- `GF` serialization is serialize-only; use `GFWithModulus` for round-trip serialization.
- `GFWithModulus::to_gf` validates degree/monic, but intentionally uses unchecked irreducibility validation (`Modulus::new_unchecked`).
- NTT APIs require NTT-friendly primes and power-of-two sizes; some calls panic if constraints are violated.
- `TowerField`/`TowerGF` inverse support is currently limited for tower degree `D2 = 2`.
- No CI workflows are currently configured under `.github/workflows/`.

## What changed
- Replaced outdated/partial guidance with a code-verified map of the current crate surface.
- Corrected commands and examples to match actual files/features in this repository.
- Added explicit feature-flag behavior from `Cargo.toml`.
- Added concrete module-add workflow aligned with existing `algebra`/`structures` patterns.
- Added clearly labeled planned roadmap items for expansion beyond current finite-field scope.
- Added gotchas section covering panic behavior, feature gating, and current implementation limits.

## Unknown/TODOs
- Decide and document preferred CI matrix once `.github/workflows/` is introduced.
- Confirm roadmap ordering/priorities for limbs vs constant-time hardening vs ECC work.
- Decide whether additional non-panicking division APIs should be standardized across wrappers.
