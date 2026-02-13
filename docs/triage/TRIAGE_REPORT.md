# Triage Report

## Scope and Method
This report was generated from repository inspection and local command verification on 2026-02-13.

Primary checks run:
- `cargo test -q`
- `cargo test -q --all-features`
- `cargo build -q --no-default-features`
- `cargo build -q --no-default-features --features alloc`
- `cargo build -q --no-default-features --features "alloc serde"`
- `cargo test -q --no-default-features` (fails; see backlog)
- `cargo test -q --doc --no-default-features` (fails; see backlog)
- `cargo clippy -q --all-targets --all-features` (fails; see backlog)

## Repo Overview
- Crate: `kreep`
- Version: `0.1.0`
- Edition: Rust `2021` (`Cargo.toml`)
- MSRV: not explicitly specified (no `rust-version` key in `Cargo.toml`, no `rust-toolchain*` file)
- CI/workflows: none currently checked in (no `.github/workflows/*`)

### Feature Flags
From `Cargo.toml`:
- `default = ["std"]`
- `std = ["alloc"]`
- `alloc = []`
- `serde = ["dep:serde", "alloc"]`
- `rand = ["dep:rand"]`
- `subtle = ["dep:subtle"]`

### Crate Layout
- `src/lib.rs`: crate docs/re-exports; `#![no_std]` with `std`/`alloc` feature gates
- `src/algebra/`: abstract traits (`Ring`, `Field`, `Group`, `AbelianGroup`)
- `src/structures/fp.rs`: prime field `Fp<P>` (Montgomery arithmetic)
- `src/structures/poly.rs`: `Poly<P>` + interpolation, irreducibility, primitivity, factorization, roots
- `src/structures/ntt.rs`: NTT info/plans and polynomial multiplication via NTT
- `src/structures/ext.rs`: `ExtField<P,D>`, `TowerField<P,D1,D2>`
- `src/structures/gf.rs`: modulus-validated wrappers (`GF`, `TowerGF`), constructors/helpers
- `src/utils.rs`: `ceil_sqrt_u64`, `gcd`, `is_prime`

### Examples
- `examples/field_basics.rs`
- `examples/poly_arithmetic.rs`
- `examples/irreducible_polys.rs`
- `examples/extension_fields.rs`
- `examples/ntt_mul.rs`
- `examples/factorization.rs` (feature-gated at runtime for `rand`)

### Tests and Benchmarks
- Integration tests: `tests/field_properties.rs`, `tests/serde_tests.rs`
- Large unit-test suites embedded in `src/structures/*`
- Benchmarks: `benches/benchmarks.rs` (criterion)

Observed test/quality status:
- `cargo test -q`: pass
- `cargo test -q --all-features`: pass
- `cargo test -q --no-default-features`: fails (feature-gating/test-target mismatches)
- `cargo test -q --doc --no-default-features`: fails (alloc-dependent docs not feature-gated)
- `cargo clippy -q --all-targets --all-features`: fails (tautological assert in `poly.rs` tests)

## Module Map (Implemented Surface)

### Algebraic Traits (`src/algebra/*`)
- `Ring`: additive/multiplicative identities and ring ops
- `Field`: `inverse`, `try_div`
- `Group`, `AbelianGroup`

### Prime Field (`src/structures/fp.rs`)
Implemented:
- `Fp<P>` with Montgomery representation and arithmetic operators
- Inverse, pow, signed pow, sqrt, Legendre symbol
- `validate_prime` helper
- Optional features:
  - `serde` serialization
  - `rand` distribution
  - `subtle` ct-eq/select + `pow_ct`/`inverse_ct`
  - `std`-gated discrete log (`discrete_log`, BSGS)
  - `alloc`-gated multiplicative order / primitive root / batch inverse

### Polynomials (`src/structures/poly.rs`)
Implemented:
- Core arithmetic/eval/division/GCD/extended GCD
- Interpolation and derivatives
- Irreducibility and primitivity checks
- Factorization pipeline (square-free, distinct-degree, equal-degree/Cantor-Zassenhaus)
- Root finding (feature `rand` for randomized factorization step)
- `serde` support under feature

### NTT/FFT (`src/structures/ntt.rs`)
Implemented:
- NTT capability detection (`ntt_info`)
- Root-of-unity lookup
- In-place NTT/INTT
- `mul_ntt`
- Precomputed `NttPlan` with forward/backward/mul APIs

### Extension Fields (`src/structures/ext.rs` + `src/structures/gf.rs`)
Implemented:
- `ExtField<P,D>` low-level type with explicit modulus-passing ops
- `TowerField<P,D1,D2>` low-level tower type
- `GF<P,D>` and `TowerGF<P,D1,D2>` wrappers bundling element + modulus
- `Modulus` validation type
- Constructors/helpers for irreducible/primitive polynomials (`gf::*`)
- Binary polynomial reference table for `GF(2^n)` (reference-only, not a binary-field arithmetic backend)

## Roadmap Seed Verification (from maintainer notes)

Status legend:
- `Implemented`: verified in code
- `Partial`: present but incomplete/limited
- `Planned`: not implemented (kept as roadmap issue)

- More tests/property-based: `Partial` (strong Fp proptests; less proptest coverage elsewhere)
- Benchmarks: `Partial` (criterion exists; missing some core algorithms)
- More `#[inline]`: `Partial`
- `no_std` support: `Partial` (library builds, but no-default test/doc paths currently fail)
- Constant-time operations: `Partial` (Fp has subtle-gated pieces; broader API not constant-time)
- Serialization (`serde`), randomness (`rand`): `Implemented` (feature-gated)
- Primality validation helpers / safer division API: `Partial` (`validate_prime`, `Field::try_div`; panic-y `/` still common)
- Multi-limb (`[u64; N]`) large-prime fields: `Planned`
- Tower extensions / extension constructors: `Partial` (tower/ext exist; ergonomics/completeness gaps)
- ECC suites + scalar multiplication: `Planned`
- FFT/NTT improvements: `Partial` (good base; fallible APIs/perf work pending)
- Polynomial factorization (Cantor-Zassenhaus): `Implemented`
- Discrete logarithm algorithms: `Partial` (Fp BSGS implemented)
- Specialized quadratic extension helper `QuadExt<P,N>`: `Planned`

## Proposed Label Taxonomy
Minimal taxonomy applied in scripts/templates.

### Type (exactly one per issue)
- `bug`: incorrect behavior or broken workflows
- `feature`: net-new capability
- `perf`: performance-specific improvements
- `tech-debt`: maintainability/quality/tooling cleanup
- `docs`: documentation/examples clarity and correctness
- `security`: side-channel, misuse-hazard, and safety hardening

### Priority (exactly one per issue)
- `P0-correctness/blocker`: correctness/blocked development path
- `P1-important`: near-term quality and developer UX
- `P2-later`: planned expansion and longer-horizon work

### Area (1-3 per issue)
- `area:traits-algebra`
- `area:field-fp`
- `area:ext-fields`
- `area:polynomials`
- `area:ntt-fft`
- `area:factorization`
- `area:discrete-log`
- `area:serialization`
- `area:randomness`
- `area:constant-time`
- `area:no-std`
- `area:benchmarks`
- `area:docs-examples`
- `area:infra-ci`

### Status (optional)
- `status:needs-investigation`
- `status:blocked`

## Proposed Milestones
- `v0.1-hardening`: unblock no_std/test matrix, fix clippy/CI, panic-safety improvements
- `v0.2-math-quality`: extension/tower correctness hardening, perf baseline expansion, API cleanup
- `v0.3-crypto-roadmap`: large-field backends, advanced dlog, ECC scaffolding, specialized extension helpers

## Triage Policy (to avoid issue bloat)
Use this intake rule:
1. Require concrete evidence before opening as `bug` (failing test, reproducible command, or precise file/line mismatch).
2. If evidence is weak, open `status:needs-investigation` with a narrow hypothesis and explicit exit criteria.
3. Reject “mega issues”; split by one outcome and one subsystem.
4. For roadmap ideas, open as `feature + P2-later` unless immediate blockers exist.
5. Any issue without a crisp Definition of Done should be refined before backlog inclusion.
6. Close stale exploratory issues that have no reproducer and no owner after one milestone.

## Automation Notes
Artifacts added for automation:
- `scripts/triage_create_labels.sh`
- `scripts/triage_create_issues.sh`
- `.github/ISSUE_TEMPLATE/*.yml`
- `docs/triage/BACKLOG.md` as seed source for issue creation

Execution run (authenticated `gh`) completed:
- Labels synced: 25 labels present (idempotent re-run verified).
- Milestones created: `v0.1-hardening`, `v0.2-math-quality`, `v0.3-crypto-roadmap`.
- Issues created from backlog: 23 (`#1` through `#23` on `odysseia06/kreep`), with idempotent re-run showing `created=0, skipped=23`.

If `gh` is unavailable or unauthenticated:
1. Run `scripts/triage_create_labels.sh` later in an authenticated shell.
2. Run `scripts/triage_create_issues.sh` after labels exist.
3. Or manually create labels/issues from `docs/triage/BACKLOG.md` using the same titles and labels.
